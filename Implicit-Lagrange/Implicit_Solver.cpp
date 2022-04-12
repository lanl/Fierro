
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>  // fmin, fmax, abs note: fminl is long
#include <sys/stat.h>
#include <mpi.h>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_Parallel_Reduce.hpp>
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"
#include "Tpetra_Details_FixedHashTable.hpp"
#include "Tpetra_Import.hpp"
#include <set>

#include "elements.h"
#include "swage.h"
#include "matar.h"
#include "utilities.h"
#include "node_combination.h"
#include "Simulation_Parameter_Headers.h"
#include "FEA_Module_Headers.h"
#include "Implicit_Solver.h"


//Optimization Package
#include "ROL_Algorithm.hpp"
#include "ROL_Solver.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Types.hpp"
#include "ROL_Elementwise_Reduce.hpp"
#include "ROL_Stream.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_StdBoundConstraint.hpp"
#include "ROL_ParameterList.hpp"
#include <ROL_TpetraMultiVector.hpp>

//Objective Functions and Constraint Functions
#include "Topology_Optimization_Function_Headers.h"

//debug and performance includes
#include <sys/time.h>
#include <sys/resource.h>
#include <ctime>

#define BUFFER_LINES 1000
#define MAX_WORD 30
#define MAX_ELEM_NODES 32
#define STRAIN_EPSILON 0.000000001
#define DENSITY_EPSILON 0.0001
#define BC_EPSILON 1.0e-8

using namespace utils;

/*

Swage is a reference to a swage block used in blacksmithing.  
Its a large metal block that has multiple shaps carved into 
each surface to use for hammering metal into to form it. 

*/

Implicit_Solver::Implicit_Solver() : Solver(){
  //create parameter objects
  simparam = new Simulation_Parameters_Topology_Optimization(this);
  // ---- Read input file, define state and boundary conditions ---- //
  simparam->Simulation_Parameters::input();
  simparam->input();
  //create ref element object
  ref_elem = new elements::ref_element();
  //create mesh objects
  //init_mesh = new swage::mesh_t(simparam);
  //mesh = new swage::mesh_t(simparam);

  element_select = new elements::element_selector();
  num_nodes = 0;

  //boundary condition data
  current_bdy_id = 0;

  //Trilinos output stream
  std::ostream &out = std::cout;
  fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  (*fos).setOutputToRootOnly(0);

  //FEA module data init
  nfea_modules = 0;
  displacement_module = -1;
}

Implicit_Solver::~Implicit_Solver(){
   delete simparam;
   delete ref_elem;
   delete element_select;
   if(myrank==0)
   delete in;
}

//==============================================================================
//    Primary simulation runtime routine
//==============================================================================


void Implicit_Solver::run(int argc, char *argv[]){
    
    //MPI info
    world = MPI_COMM_WORLD; //used for convenience to represent all the ranks in the job
    MPI_Comm_rank(world,&myrank);
    MPI_Comm_size(world,&nranks);
    
    if(myrank == 0){
      std::cout << "Running TO Solver" << std::endl;
       // check to see of a mesh was supplied when running the code
      if (argc == 1) {
        std::cout << "\n\n**********************************\n\n";
        std::cout << " ERROR:\n";
        std::cout << " Please supply a mesh file as the second command line argument \n";
        std::cout << "**********************************\n\n" << std::endl;
        return;
      }
    }
    //initialize Trilinos communicator class
    comm = Tpetra::getDefaultComm();

    //error handle for file input name
    //if(argc < 2)

    // ---- Read intial mesh, refine, and build connectivity ---- //
    if(simparam->tecplot_input)
      read_mesh_tecplot(argv[1]);
    else if(simparam->ansys_dat_input)
      read_mesh_ansys_dat(argv[1]);
    else
      read_mesh_ensight(argv[1]);
    
    //construct list of FEA modules requested
    simparam->FEA_module_setup();

    //process process list of requested FEA modules to construct list of objects
    FEA_module_setup();

    //debug
    //return;
    init_maps();
    
    std::cout << "Num elements on process " << myrank << " = " << rnum_elem << std::endl;
    
    //initialize timing
    if(simparam->report_runtime_flag)
      init_clock();

    //initialize runtime counters and timers
    int hessvec_count = 0;
    int update_count = 0;
    file_index = 0;
    real_t linear_solve_time = 0;
    real_t hessvec_time = 0;
    real_t hessvec_linear_time = 0;
    
    // ---- Find Boundaries on mesh ---- //
    init_boundaries();

    //set boundary conditions
    generate_tcs();

    //initialize TO design variable storage
    init_design();

    //Have modules read in boundary/loading conditions if file format provides it
    for(int imodule = 0; imodule < nfea_modules; imodule++){
      if(fea_module_must_read[imodule]){
        fea_modules[imodule]->read_conditions_ansys_dat(in, before_condition_header);
      }
      else{
        fea_modules[imodule]->init_boundaries();

        //set boundary conditions for FEA modules
        fea_modules[imodule]->generate_bcs();

        //set applied loading conditions for FEA modules
        fea_modules[imodule]->generate_applied_loads();
      }
    }

    //std::cout << "FEA MODULES " << nfea_modules << " " << simparam->nfea_modules << std::endl;
    //call boundary routines on fea modules

    if(myrank == 0)
    std::cout << "Starting init assembly" << std::endl <<std::flush;
    //allocate and fill sparse structures needed for global solution in each FEA module
    for(int imodule = 0; imodule < nfea_modules; imodule++)
      fea_modules[imodule]->init_assembly();
    
    //assemble the global solution (stiffness matrix etc. and nodal forces)
    for(int imodule = 0; imodule < nfea_modules; imodule++)
      fea_modules[imodule]->assemble_matrix();

    if(myrank == 0)
    std::cout << "Finished matrix assembly" << std::endl <<std::flush;
    
    for(int imodule = 0; imodule < nfea_modules; imodule++)
      fea_modules[imodule]->assemble_vector();

    for(int imodule = 0; imodule < nfea_modules; imodule++)
      fea_modules[imodule]->linear_solver_parameters();
    
    if(myrank == 0)
    std::cout << "Starting First Solve" << std::endl <<std::flush;
    
    for(int imodule = 0; imodule < nfea_modules; imodule++){
      int solver_exit = fea_modules[imodule]->solve();
      if(solver_exit == EXIT_SUCCESS){
        std::cout << "Linear Solver Error" << std::endl <<std::flush;
        return;
      }
    }
    /*
    //debug print
    
    Teuchos::RCP<MV> design_gradients_distributed = Teuchos::rcp(new MV(map, 1));
    const_host_vec_array node_densities = design_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    host_vec_array design_gradients = design_gradients_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    compute_adjoint_gradients(node_densities, design_gradients);
    std::ostream &out = std::cout;
    Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    if(myrank==0)
    *fos << "Strain Energy gradients :" << std::endl;
    design_gradients_distributed->describe(*fos,Teuchos::VERB_EXTREME);
    *fos << std::endl;
    std::fflush(stdout);
    */
    //return;
    setup_optimization_problem();
    
    //solver_exit = solve();
    //if(solver_exit == EXIT_SUCCESS){
      //std::cout << "Linear Solver Error" << std::endl <<std::flush;
      //return;
    //}
    
    //CPU time
    double current_cpu = CPU_Time();
    for(int imodule = 0; imodule < nfea_modules; imodule++){
      linear_solve_time += fea_modules[imodule]->linear_solve_time;
      hessvec_linear_time += fea_modules[imodule]->hessvec_linear_time;
    }

    std::cout << " RUNTIME OF CODE ON TASK " << myrank << " is "<< current_cpu-initial_CPU_time << " update solve time "
              << linear_solve_time << " hess solve time " << hessvec_linear_time <<std::endl;

    // Data writers
    tecplot_writer();
    // vtk_writer();
    if(myrank==0){
      std::cout << "Total number of solves and assembly " << fea_modules[0]->update_count <<std::endl;
      std::cout << "Total number of hessvec counts " << fea_modules[0]->hessvec_count <<std::endl;
      std::cout << "End of Optimization" << std::endl;
    }
}

/* ----------------------------------------------------------------------
   Read Ensight format mesh file
------------------------------------------------------------------------- */
void Implicit_Solver::read_mesh_ensight(char *MESH){

  char ch;
  int num_dim = simparam->num_dim;
  int p_order = simparam->p_order;
  real_t unit_scaling = simparam->unit_scaling;
  int local_node_index, current_column_index;
  size_t strain_count;
  std::string skip_line, read_line, substring;
  std::stringstream line_parse;
  CArrayKokkos<char, array_layout, HostSpace, memory_traits> read_buffer;
  int buffer_loop, buffer_iteration, scan_loop;
  size_t read_index_start, node_rid, elem_gid;
  GO node_gid;
  real_t dof_value;
  //Nodes_Per_Element_Type =  elements::elem_types::Nodes_Per_Element_Type;

  //read the mesh
  //PLACEHOLDER: ensight_format(MESH);
  // abaqus_format(MESH);
  // vtk_format(MESH)

  //task 0 reads file
  if(myrank==0){
  in = new std::ifstream();
  in->open(MESH);  
    
  
  if(simparam->tecplot_input){
    //skip 2 lines
    for (int j = 1; j <= 2; j++) {
      getline(*in, skip_line);
      std::cout << skip_line << std::endl;
    } //for
  }
  else{
    //skip 8 lines
    for (int j = 1; j <= 8; j++) {
      getline(*in, skip_line);
      std::cout << skip_line << std::endl;
    } //for
  }
  }

  // --- Read the number of nodes in the mesh --- //
  if(myrank==0){
    getline(*in, read_line);
    line_parse.str(read_line);
    line_parse >> num_nodes;
    std::cout << "declared node count: " << num_nodes << std::endl;
  }
  
  //broadcast number of nodes
  MPI_Bcast(&num_nodes,1,MPI_LONG_LONG_INT,0,world);
  
  //construct contiguous parallel row map now that we know the number of nodes
  map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(num_nodes,0,comm));

  // set the vertices in the mesh read in
  global_size_t local_nrows = map->getNodeNumElements();
  nlocal_nodes = local_nrows;
  //populate local row offset data from global data
  global_size_t min_gid = map->getMinGlobalIndex();
  global_size_t max_gid = map->getMaxGlobalIndex();
  global_size_t index_base = map->getIndexBase();
  //debug print
  //std::cout << "local node count on task: " << " " << nlocal_nodes << std::endl;
  //construct dof map that follows from the node map (used for distributed matrix and vector objects later)
  CArrayKokkos<GO, array_layout, device_type, memory_traits> local_dof_indices(nlocal_nodes*num_dim, "local_dof_indices");
  for(int i = 0; i < nlocal_nodes; i++){
    for(int j = 0; j < num_dim; j++)
    local_dof_indices(i*num_dim + j) = map->getGlobalElement(i)*num_dim + j;
  }
  
  local_dof_map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(num_nodes*num_dim,local_dof_indices.get_kokkos_view(),0,comm) );

  //allocate node storage with dual view
  dual_node_coords = dual_vec_array("dual_node_coords", nlocal_nodes,num_dim);

  //local variable for host view in the dual view
  host_vec_array node_coords = dual_node_coords.view_host();
  //notify that the host view is going to be modified in the file readin
  dual_node_coords.modify_host();

  //old swage method
  //mesh->init_nodes(local_nrows); // add 1 for index starting at 1
    
  std::cout << "Num nodes assigned to task " << myrank << " = " << nlocal_nodes << std::endl;

  // read the initial mesh coordinates
  // x-coords
  /*only task 0 reads in nodes and elements from the input file
  stores node data in a buffer and communicates once the buffer cap is reached
  or the data ends*/

  words_per_line = simparam->words_per_line;
  elem_words_per_line = simparam->elem_words_per_line;

  //allocate read buffer
  read_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(BUFFER_LINES,words_per_line,MAX_WORD);

  int dof_limit = num_nodes;
  int buffer_iterations = dof_limit/BUFFER_LINES;
  if(dof_limit%BUFFER_LINES!=0) buffer_iterations++;
  
  //x-coords
  read_index_start = 0;
  for(buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
    //pack buffer on rank 0
    if(myrank==0&&buffer_iteration<buffer_iterations-1){
      for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        
        for(int iword = 0; iword < words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //debug print
        //std::cout<<" "<< substring <<std::endl;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
      }
    }
    else if(myrank==0){
      buffer_loop=0;
      while(buffer_iteration*BUFFER_LINES+buffer_loop < num_nodes) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
        buffer_loop++;
      }
      
    }

    //broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
    MPI_Bcast(read_buffer.pointer(),BUFFER_LINES*words_per_line*MAX_WORD,MPI_CHAR,0,world);
    //broadcast how many nodes were read into this buffer iteration
    MPI_Bcast(&buffer_loop,1,MPI_INT,0,world);

    //debug_print
    //std::cout << "NODE BUFFER LOOP IS: " << buffer_loop << std::endl;
    //for(int iprint=0; iprint < buffer_loop; iprint++)
      //std::cout<<"buffer packing: " << std::string(&read_buffer(iprint,0,0)) << std::endl;
    //return;

    //determine which data to store in the swage mesh members (the local node data)
    //loop through read buffer
    for(scan_loop = 0; scan_loop < buffer_loop; scan_loop++){
      //set global node id (ensight specific order)
      node_gid = read_index_start + scan_loop;
      //let map decide if this node id belongs locally; if yes store data
      if(map->isNodeGlobalElement(node_gid)){
        //set local node index in this mpi rank
        node_rid = map->getLocalElement(node_gid);
        //extract nodal position from the read buffer
        //for ensight format this is just one coordinate per line
        dof_value = atof(&read_buffer(scan_loop,0,0));
        node_coords(node_rid, 0) = dof_value * unit_scaling;
      }
    }
    read_index_start+=BUFFER_LINES;
  }

  
  // y-coords
  read_index_start = 0;
  for(buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
    
    //pack buffer on rank 0
    if(myrank==0&&buffer_iteration<buffer_iterations-1){
      for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
      }
    }
    else if(myrank==0){
      buffer_loop=0;
      while(buffer_iteration*BUFFER_LINES+buffer_loop < num_nodes) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
        buffer_loop++;
        //std::cout<<" "<< node_coords(node_gid, 0)<<std::endl;
      }
    }

    //broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
    MPI_Bcast(read_buffer.pointer(),BUFFER_LINES*words_per_line*MAX_WORD,MPI_CHAR,0,world);
    //broadcast how many nodes were read into this buffer iteration
    MPI_Bcast(&buffer_loop,1,MPI_INT,0,world);

    //determine which data to store in the swage mesh members (the local node data)
    //loop through read buffer
    for(scan_loop = 0; scan_loop < buffer_loop; scan_loop++){
      //set global node id (ensight specific order)
      node_gid = read_index_start + scan_loop;
      //let map decide if this node id belongs locally; if yes store data
      if(map->isNodeGlobalElement(node_gid)){
        //set local node index in this mpi rank
        node_rid = map->getLocalElement(node_gid);
        //extract nodal position from the read buffer
        //for ensight format this is just one coordinate per line
        dof_value = atof(&read_buffer(scan_loop,0,0));
        node_coords(node_rid, 1) = dof_value * unit_scaling;
      }
    }
    read_index_start+=BUFFER_LINES;
  }

  // z-coords
  read_index_start = 0;
  if(num_dim==3)
  for(buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
    
    //pack buffer on rank 0
    if(myrank==0&&buffer_iteration<buffer_iterations-1){
      for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
      }
    }
    else if(myrank==0){
      buffer_loop=0;
      while(buffer_iteration*BUFFER_LINES+buffer_loop < num_nodes) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
        buffer_loop++;
        //std::cout<<" "<< node_coords(node_gid, 0)<<std::endl;
      }
    }

    //broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
    MPI_Bcast(read_buffer.pointer(),BUFFER_LINES*words_per_line*MAX_WORD,MPI_CHAR,0,world);
    //broadcast how many nodes were read into this buffer iteration
    MPI_Bcast(&buffer_loop,1,MPI_INT,0,world);

    //determine which data to store in the swage mesh members (the local node data)
    //loop through read buffer
    for(scan_loop = 0; scan_loop < buffer_loop; scan_loop++){
      //set global node id (ensight specific order)
      node_gid = read_index_start + scan_loop;
      //let map decide if this node id belongs locally; if yes store data
      if(map->isNodeGlobalElement(node_gid)){
        //set local node index in this mpi rank
        node_rid = map->getLocalElement(node_gid);
        //extract nodal position from the read buffer
        //for ensight format this is just one coordinate per line
        dof_value = atof(&read_buffer(scan_loop,0,0));
        node_coords(node_rid, 2) = dof_value * unit_scaling;
      }
    }
    read_index_start+=BUFFER_LINES;
  }
  
  
  //debug print of nodal data
  
  //debug print nodal positions and indices
  /*
  std::cout << " ------------NODAL POSITIONS ON TASK " << myrank << " --------------"<<std::endl;
  for (int inode = 0; inode < local_nrows; inode++){
      std::cout << "node: " << map->getGlobalElement(inode) + 1 << " { ";
    for (int istride = 0; istride < num_dim; istride++){
        std::cout << node_coords(inode,istride) << " , ";
    }
    std::cout << " }"<< std::endl;
  }
  */

  //check that local assignments match global total

  
  //read in element info (ensight file format is organized in element type sections)
  //loop over this later for several element type sections

  num_elem = 0;
  rnum_elem = 0;
  CArrayKokkos<int, array_layout, HostSpace, memory_traits> node_store(elem_words_per_line);

  if(myrank==0){
  //skip element type name line
    getline(*in, skip_line);
    std::cout << skip_line << std::endl;
  }
    
  // --- read the number of cells in the mesh ---
  // --- Read the number of vertices in the mesh --- //
  if(myrank==0){
    getline(*in, read_line);
    line_parse.clear();
    line_parse.str(read_line);
    line_parse >> num_elem;
    std::cout << "declared element count: " << num_elem << std::endl;
    if(num_elem <= 0) std::cout << "ERROR, NO ELEMENTS IN MESH" << std::endl;
  }
  
  //broadcast number of elements
  MPI_Bcast(&num_elem,1,MPI_LONG_LONG_INT,0,world);
  
  if(myrank == 0)
  std::cout<<"before mesh initialization"<<std::endl;
  
  //read in element connectivity
  //we're gonna reallocate for the words per line expected for the element connectivity
  read_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(BUFFER_LINES,elem_words_per_line,MAX_WORD);

  //calculate buffer iterations to read number of lines
  buffer_iterations = num_elem/BUFFER_LINES;
  int assign_flag;

  //dynamic buffer used to store elements before we know how many this rank needs
  std::vector<size_t> element_temp(BUFFER_LINES*elem_words_per_line);
  std::vector<size_t> global_indices_temp(BUFFER_LINES);
  size_t buffer_max = BUFFER_LINES*elem_words_per_line;
  size_t indices_buffer_max = BUFFER_LINES;

  if(num_elem%BUFFER_LINES!=0) buffer_iterations++;
  read_index_start = 0;
  //std::cout << "ELEMENT BUFFER ITERATIONS: " << buffer_iterations << std::endl;
  rnum_elem = 0;
  for(buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
    //pack buffer on rank 0
    if(myrank==0&&buffer_iteration<buffer_iterations-1){
      for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < elem_words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
      }
    }
    else if(myrank==0){
      buffer_loop=0;
      while(buffer_iteration*BUFFER_LINES+buffer_loop < num_elem) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < elem_words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
        buffer_loop++;
        //std::cout<<" "<< node_coords(node_gid, 0)<<std::endl;
      }
    }

    //broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
    MPI_Bcast(read_buffer.pointer(),BUFFER_LINES*elem_words_per_line*MAX_WORD,MPI_CHAR,0,world);
    //broadcast how many nodes were read into this buffer iteration
    MPI_Bcast(&buffer_loop,1,MPI_INT,0,world);
    
    //store element connectivity that belongs to this rank
    //loop through read buffer
    for(scan_loop = 0; scan_loop < buffer_loop; scan_loop++){
      //set global node id (ensight specific order)
      elem_gid = read_index_start + scan_loop;
      //add this element to the local list if any of its nodes belong to this rank according to the map
      //get list of nodes for each element line and check if they belong to the map
      assign_flag = 0;
      for(int inode = 0; inode < elem_words_per_line; inode++){
        //as we loop through the nodes belonging to this element we store them
        //if any of these nodes belongs to this rank this list is used to store the element locally
        node_gid = atoi(&read_buffer(scan_loop,inode,0));
        node_store(inode) = node_gid - 1; //subtract 1 since file index start is 1 but code expects 0
        //first we add the elements to a dynamically allocated list
        if(map->isNodeGlobalElement(node_gid-1)&&!assign_flag){
          assign_flag = 1;
          rnum_elem++;
        }
      }

      if(assign_flag){
        for(int inode = 0; inode < elem_words_per_line; inode++){
          if((rnum_elem-1)*elem_words_per_line + inode>=buffer_max){ 
            element_temp.resize((rnum_elem-1)*elem_words_per_line + inode + BUFFER_LINES*elem_words_per_line);
            buffer_max = (rnum_elem-1)*elem_words_per_line + inode + BUFFER_LINES*elem_words_per_line;
          }
          element_temp[(rnum_elem-1)*elem_words_per_line + inode] = node_store(inode); 
          //std::cout << "VECTOR STORAGE FOR ELEM " << rnum_elem << " ON TASK " << myrank << " NODE " << inode+1 << " IS " << node_store(inode) + 1 << std::endl;
        }
        //assign global element id to temporary list
        if(rnum_elem-1>=indices_buffer_max){ 
          global_indices_temp.resize(rnum_elem-1 + BUFFER_LINES);
          indices_buffer_max = rnum_elem-1 + BUFFER_LINES;
        }
        global_indices_temp[rnum_elem-1] = elem_gid;
      }
    }
    read_index_start+=BUFFER_LINES;
  }

  // Close mesh input file
  if(myrank==0)
  in->close();
  
  //std::cout << "RNUM ELEMENTS IS: " << rnum_elem << std::endl;
  
  Element_Types = CArrayKokkos<elements::elem_types::elem_type, array_layout, HostSpace, memory_traits>(rnum_elem);
  
  elements::elem_types::elem_type mesh_element_type;
  if(simparam->element_type == "Hex8"){
    mesh_element_type = elements::elem_types::Hex8;
  }
  else if(simparam->element_type == "Hex20"){
    mesh_element_type = elements::elem_types::Hex20;
  }
  else if(simparam->element_type == "Hex32"){
    mesh_element_type = elements::elem_types::Hex32;
  }

  //set element object pointer
  if(simparam->num_dim==2){
    element_select->choose_2Delem_type(mesh_element_type, elem2D);
     max_nodes_per_element = elem2D->num_nodes();
  }
  else if(simparam->num_dim==3){
    element_select->choose_3Delem_type(mesh_element_type, elem);
     max_nodes_per_element = elem->num_nodes();
  }

  //1 type per mesh for now
  for(int ielem = 0; ielem < rnum_elem; ielem++)
    Element_Types(ielem) = mesh_element_type;
  
  //copy temporary element storage to multivector storage
  dual_nodes_in_elem = dual_elem_conn_array("dual_nodes_in_elem", rnum_elem, max_nodes_per_element);
  nodes_in_elem = dual_nodes_in_elem.view_host();
  dual_nodes_in_elem.modify_host();

  for(int ielem = 0; ielem < rnum_elem; ielem++)
    for(int inode = 0; inode < elem_words_per_line; inode++)
      nodes_in_elem(ielem, inode) = element_temp[ielem*elem_words_per_line + inode];

  //view storage for all local elements connected to local nodes on this rank
  CArrayKokkos<GO, array_layout, device_type, memory_traits> All_Element_Global_Indices(rnum_elem);

  //copy temporary global indices storage to view storage
  for(int ielem = 0; ielem < rnum_elem; ielem++)
    All_Element_Global_Indices(ielem) = global_indices_temp[ielem];

  //delete temporary element connectivity and index storage
  std::vector<size_t>().swap(element_temp);
  std::vector<size_t>().swap(global_indices_temp);
  
  //construct overlapping element map (since different ranks can own the same elements due to the local node map)
  all_element_map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),All_Element_Global_Indices.get_kokkos_view(),0,comm));

  //element type selection (subject to change)
  // ---- Set Element Type ---- //
  // allocate element type memory
  //elements::elem_type_t* elem_choice;

  int NE = 1; // number of element types in problem
    

  // Convert ijk index system to the finite element numbering convention
  // for vertices in cell
  CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> convert_ensight_to_ijk(max_nodes_per_element);
  CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> tmp_ijk_indx(max_nodes_per_element);
  convert_ensight_to_ijk(0) = 0;
  convert_ensight_to_ijk(1) = 1;
  convert_ensight_to_ijk(2) = 3;
  convert_ensight_to_ijk(3) = 2;
  convert_ensight_to_ijk(4) = 4;
  convert_ensight_to_ijk(5) = 5;
  convert_ensight_to_ijk(6) = 7;
  convert_ensight_to_ijk(7) = 6;
    
  int nodes_per_element;
  
  if(num_dim==2)
  for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++) {
    //set nodes per element
    element_select->choose_2Delem_type(Element_Types(cell_rid), elem2D);
    nodes_per_element = elem2D->num_nodes();
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      tmp_ijk_indx(node_lid) = nodes_in_elem(cell_rid, convert_ensight_to_ijk(node_lid));
    }   
        
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      nodes_in_elem(cell_rid, node_lid) = tmp_ijk_indx(node_lid);
    }
  }

  if(num_dim==3)
  for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++) {
    //set nodes per element
    element_select->choose_3Delem_type(Element_Types(cell_rid), elem);
    nodes_per_element = elem->num_nodes();
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      tmp_ijk_indx(node_lid) = nodes_in_elem(cell_rid, convert_ensight_to_ijk(node_lid));
    }   
        
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      nodes_in_elem(cell_rid, node_lid) = tmp_ijk_indx(node_lid);
    }
  }
  
  //synchronize device data
  dual_node_coords.sync_device();
  dual_node_coords.modify_device();
  //debug print element edof
  /*
  std::cout << " ------------ELEMENT EDOF ON TASK " << myrank << " --------------"<<std::endl;

  for (int ielem = 0; ielem < rnum_elem; ielem++){
    std::cout << "elem:  " << ielem+1 << std::endl;
    for (int lnode = 0; lnode < 8; lnode++){
        std::cout << "{ ";
          std::cout << lnode+1 << " = " << nodes_in_elem(ielem,lnode) + 1 << " ";
        
        std::cout << " }"<< std::endl;
    }
    std::cout << std::endl;
  }
  */
 
} // end read_mesh

/* ----------------------------------------------------------------------
   Read Tecplot format mesh file
------------------------------------------------------------------------- */
void Implicit_Solver::read_mesh_tecplot(char *MESH){

  char ch;
  int num_dim = simparam->num_dim;
  int p_order = simparam->p_order;
  real_t unit_scaling = simparam->unit_scaling;
  bool restart_file = simparam->restart_file;
  int local_node_index, current_column_index;
  size_t strain_count;
  std::string skip_line, read_line, substring;
  std::stringstream line_parse;
  CArrayKokkos<char, array_layout, HostSpace, memory_traits> read_buffer;
  int buffer_loop, buffer_iteration, scan_loop;
  size_t read_index_start, node_rid, elem_gid;
  GO node_gid;
  real_t dof_value;
  host_vec_array node_densities;
  //Nodes_Per_Element_Type =  elements::elem_types::Nodes_Per_Element_Type;

  //read the mesh
  //PLACEHOLDER: ensight_format(MESH);
  // abaqus_format(MESH);
  // vtk_format(MESH)

  //task 0 reads file
  if(myrank==0){
  in = new std::ifstream();
  in->open(MESH);  
    //skip 2 lines
    for (int j = 1; j <= 2; j++) {
      getline(*in, skip_line);
      std::cout << skip_line << std::endl;
    } //for
  }
  

  // --- Read the number of nodes in the mesh --- //
  if(myrank==0){
    getline(*in, read_line);
    line_parse.str(read_line);
    //stop when the NODES= string is reached
    while (!line_parse.eof()){
      line_parse >> substring;
      if(!substring.compare("NODES=")){
        line_parse >> num_nodes;
      }
      if(!substring.compare("ELEMENTS=")){
        line_parse >> num_elem;
      }
    } //while
    std::cout << "declared node count: " << num_nodes << std::endl;
    std::cout << "declared element count: " << num_elem << std::endl;
    if(num_elem <= 0) std::cout << "ERROR, NO ELEMENTS IN MESH!!!!" << std::endl;
  }
  
  //broadcast number of nodes
  MPI_Bcast(&num_nodes,1,MPI_LONG_LONG_INT,0,world);
  
  //construct contiguous parallel row map now that we know the number of nodes
  map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(num_nodes,0,comm));

  // set the vertices in the mesh read in
  global_size_t local_nrows = map->getNodeNumElements();
  nlocal_nodes = local_nrows;
  //populate local row offset data from global data
  global_size_t min_gid = map->getMinGlobalIndex();
  global_size_t max_gid = map->getMaxGlobalIndex();
  global_size_t index_base = map->getIndexBase();
  //debug print
  //std::cout << "local node count on task: " << " " << nlocal_nodes << std::endl;
  //construct dof map that follows from the node map (used for distributed matrix and vector objects later)
  CArrayKokkos<GO, array_layout, device_type, memory_traits> local_dof_indices(nlocal_nodes*num_dim, "local_dof_indices");
  for(int i = 0; i < nlocal_nodes; i++){
    for(int j = 0; j < num_dim; j++)
    local_dof_indices(i*num_dim + j) = map->getGlobalElement(i)*num_dim + j;
  }
  
  local_dof_map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(num_nodes*num_dim,local_dof_indices.get_kokkos_view(),0,comm) );

  //allocate node storage with dual view
  dual_node_coords = dual_vec_array("dual_node_coords", nlocal_nodes,num_dim);
  if(restart_file)
    dual_node_densities = dual_vec_array("dual_node_densities", nlocal_nodes,1);

  //local variable for host view in the dual view
  host_vec_array node_coords = dual_node_coords.view_host();
  if(restart_file)
    node_densities = dual_node_densities.view_host();
  //notify that the host view is going to be modified in the file readin
  dual_node_coords.modify_host();
  if(restart_file)
    dual_node_densities.modify_host();

  //old swage method
  //mesh->init_nodes(local_nrows); // add 1 for index starting at 1
    
  std::cout << "Num nodes assigned to task " << myrank << " = " << nlocal_nodes << std::endl;

  // read the initial mesh coordinates
  // x-coords
  /*only task 0 reads in nodes and elements from the input file
  stores node data in a buffer and communicates once the buffer cap is reached
  or the data ends*/

  words_per_line = simparam->tecplot_words_per_line;
  if(restart_file) words_per_line++;
  elem_words_per_line = simparam->elem_words_per_line;

  //allocate read buffer
  read_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(BUFFER_LINES,words_per_line,MAX_WORD);

  int dof_limit = num_nodes;
  int buffer_iterations = dof_limit/BUFFER_LINES;
  if(dof_limit%BUFFER_LINES!=0) buffer_iterations++;
  
  //read coords, also density if restarting
  read_index_start = 0;
  for(buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
    //pack buffer on rank 0
    if(myrank==0&&buffer_iteration<buffer_iterations-1){
      for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        
        for(int iword = 0; iword < words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //debug print
        //std::cout<<" "<< substring <<std::endl;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
      }
    }
    else if(myrank==0){
      buffer_loop=0;
      while(buffer_iteration*BUFFER_LINES+buffer_loop < num_nodes) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
        buffer_loop++;
      }
      
    }

    //broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
    MPI_Bcast(read_buffer.pointer(),BUFFER_LINES*words_per_line*MAX_WORD,MPI_CHAR,0,world);
    //broadcast how many nodes were read into this buffer iteration
    MPI_Bcast(&buffer_loop,1,MPI_INT,0,world);

    //debug_print
    //std::cout << "NODE BUFFER LOOP IS: " << buffer_loop << std::endl;
    //for(int iprint=0; iprint < buffer_loop; iprint++)
      //std::cout<<"buffer packing: " << std::string(&read_buffer(iprint,0,0)) << std::endl;
    //return;

    //determine which data to store in the swage mesh members (the local node data)
    //loop through read buffer
    for(scan_loop = 0; scan_loop < buffer_loop; scan_loop++){
      //set global node id (ensight specific order)
      node_gid = read_index_start + scan_loop;
      //let map decide if this node id belongs locally; if yes store data
      if(map->isNodeGlobalElement(node_gid)){
        //set local node index in this mpi rank
        node_rid = map->getLocalElement(node_gid);
        //extract nodal position from the read buffer
        //for tecplot format this is the three coords in the same line
        dof_value = atof(&read_buffer(scan_loop,0,0));
        node_coords(node_rid, 0) = dof_value * unit_scaling;
        dof_value = atof(&read_buffer(scan_loop,1,0));
        node_coords(node_rid, 1) = dof_value * unit_scaling;
        dof_value = atof(&read_buffer(scan_loop,2,0));
        node_coords(node_rid, 2) = dof_value * unit_scaling;
        if(restart_file){
          dof_value = atof(&read_buffer(scan_loop,3,0));
          node_densities(node_rid, 0) = dof_value;
        }
        //extract density if restarting
      }
    }
    read_index_start+=BUFFER_LINES;
  }

  
  //debug print of nodal data
  
  //debug print nodal positions and indices
  
  //std::cout << " ------------NODAL POSITIONS ON TASK " << myrank << " --------------"<<std::endl;
  //for (int inode = 0; inode < local_nrows; inode++){
      //std::cout << "node: " << map->getGlobalElement(inode) + 1 << " { ";
    //for (int istride = 0; istride < num_dim; istride++){
       //std::cout << node_coords(inode,istride) << " , ";
    //}
    //std::cout << node_densities(inode,0);
    //std::cout << " }"<< std::endl;
  //}
  

  //check that local assignments match global total

  
  //read in element info (supported tecplot format currently assumes one type)

  CArrayKokkos<int, array_layout, HostSpace, memory_traits> node_store(elem_words_per_line);
  
  //broadcast number of elements
  MPI_Bcast(&num_elem,1,MPI_LONG_LONG_INT,0,world);
  //std::cout<<"before initial mesh initialization"<<std::endl;
  
  //read in element connectivity
  //we're gonna reallocate for the words per line expected for the element connectivity
  read_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(BUFFER_LINES,elem_words_per_line,MAX_WORD);

  //calculate buffer iterations to read number of lines
  buffer_iterations = num_elem/BUFFER_LINES;
  int assign_flag;

  //dynamic buffer used to store elements before we know how many this rank needs
  std::vector<size_t> element_temp(BUFFER_LINES*elem_words_per_line);
  std::vector<size_t> global_indices_temp(BUFFER_LINES);
  size_t buffer_max = BUFFER_LINES*elem_words_per_line;
  size_t indices_buffer_max = BUFFER_LINES;

  if(num_elem%BUFFER_LINES!=0) buffer_iterations++;
  read_index_start = 0;
  //std::cout << "ELEMENT BUFFER ITERATIONS: " << buffer_iterations << std::endl;
  rnum_elem = 0;
  //std::cout << "BUFFER ITERATIONS IS: " << buffer_iterations << std::endl;
  for(buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
    //pack buffer on rank 0
    if(myrank==0&&buffer_iteration<buffer_iterations-1){
      for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < elem_words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
      }
    }
    else if(myrank==0){
      buffer_loop=0;
      while(buffer_iteration*BUFFER_LINES+buffer_loop < num_elem) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < elem_words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
        buffer_loop++;
        //std::cout<<" "<< node_coords(node_gid, 0)<<std::endl;
      }
    }

    //broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
    MPI_Bcast(read_buffer.pointer(),BUFFER_LINES*elem_words_per_line*MAX_WORD,MPI_CHAR,0,world);
    //broadcast how many nodes were read into this buffer iteration
    MPI_Bcast(&buffer_loop,1,MPI_INT,0,world);
    
    //store element connectivity that belongs to this rank
    //loop through read buffer
    for(scan_loop = 0; scan_loop < buffer_loop; scan_loop++){
      //set global node id (ensight specific order)
      elem_gid = read_index_start + scan_loop;
      //add this element to the local list if any of its nodes belong to this rank according to the map
      //get list of nodes for each element line and check if they belong to the map
      assign_flag = 0;
      for(int inode = 0; inode < elem_words_per_line; inode++){
        //as we loop through the nodes belonging to this element we store them
        //if any of these nodes belongs to this rank this list is used to store the element locally
        node_gid = atoi(&read_buffer(scan_loop,inode,0));
        node_store(inode) = node_gid - 1; //subtract 1 since file index start is 1 but code expects 0
        //first we add the elements to a dynamically allocated list
        if(map->isNodeGlobalElement(node_gid-1)&&!assign_flag){
          assign_flag = 1;
          rnum_elem++;
        }
      }

      if(assign_flag){
        for(int inode = 0; inode < elem_words_per_line; inode++){
          if((rnum_elem-1)*elem_words_per_line + inode>=buffer_max){ 
            element_temp.resize((rnum_elem-1)*elem_words_per_line + inode + BUFFER_LINES*elem_words_per_line);
            buffer_max = (rnum_elem-1)*elem_words_per_line + inode + BUFFER_LINES*elem_words_per_line;
          }
          element_temp[(rnum_elem-1)*elem_words_per_line + inode] = node_store(inode); 
          //std::cout << "VECTOR STORAGE FOR ELEM " << rnum_elem << " ON TASK " << myrank << " NODE " << inode+1 << " IS " << node_store(inode) + 1 << std::endl;
        }
        //assign global element id to temporary list
        if(rnum_elem-1>=indices_buffer_max){ 
          global_indices_temp.resize(rnum_elem-1 + BUFFER_LINES);
          indices_buffer_max = rnum_elem-1 + BUFFER_LINES;
        }
        global_indices_temp[rnum_elem-1] = elem_gid;
      }
    }
    read_index_start+=BUFFER_LINES;
  }

  // Close mesh input file
  if(myrank==0)
  in->close();
  
  std::cout << "RNUM ELEMENTS IS: " << rnum_elem << std::endl;
  //copy temporary element storage to multivector storage
  Element_Types = CArrayKokkos<elements::elem_types::elem_type, array_layout, HostSpace, memory_traits>(rnum_elem);
  
  elements::elem_types::elem_type mesh_element_type;
  if(simparam->element_type == "Hex8"){
    mesh_element_type = elements::elem_types::Hex8;
  }
  else if(simparam->element_type == "Hex20"){
    mesh_element_type = elements::elem_types::Hex20;
  }
  else if(simparam->element_type == "Hex32"){
    mesh_element_type = elements::elem_types::Hex32;
  }

  //1 type per mesh for now
  for(int ielem = 0; ielem < rnum_elem; ielem++)
    Element_Types(ielem) = mesh_element_type;

  //set element object pointer
  if(simparam->num_dim==2){
    element_select->choose_2Delem_type(mesh_element_type, elem2D);
     max_nodes_per_element = elem2D->num_nodes();
  }
  else if(simparam->num_dim==3){
    element_select->choose_3Delem_type(mesh_element_type, elem);
     max_nodes_per_element = elem->num_nodes();
  }

  dual_nodes_in_elem = dual_elem_conn_array("dual_nodes_in_elem", rnum_elem, max_nodes_per_element);
  nodes_in_elem = dual_nodes_in_elem.view_host();
  dual_nodes_in_elem.modify_host();

  for(int ielem = 0; ielem < rnum_elem; ielem++)
    for(int inode = 0; inode < elem_words_per_line; inode++)
      nodes_in_elem(ielem, inode) = element_temp[ielem*elem_words_per_line + inode];

  //view storage for all local elements connected to local nodes on this rank
  CArrayKokkos<GO, array_layout, device_type, memory_traits> All_Element_Global_Indices(rnum_elem);

  //copy temporary global indices storage to view storage
  for(int ielem = 0; ielem < rnum_elem; ielem++)
    All_Element_Global_Indices(ielem) = global_indices_temp[ielem];

  //delete temporary element connectivity and index storage
  std::vector<size_t>().swap(element_temp);
  std::vector<size_t>().swap(global_indices_temp);
  
  //construct overlapping element map (since different ranks can own the same elements due to the local node map)
  all_element_map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),All_Element_Global_Indices.get_kokkos_view(),0,comm));


  //element type selection (subject to change)
  // ---- Set Element Type ---- //
  // allocate element type memory
  //elements::elem_type_t* elem_choice;

  int NE = 1; // number of element types in problem

  // Convert ijk index system to the finite element numbering convention
  // for vertices in cell
  CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> convert_ensight_to_ijk(max_nodes_per_element);
  CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> tmp_ijk_indx(max_nodes_per_element);
  convert_ensight_to_ijk(0) = 0;
  convert_ensight_to_ijk(1) = 1;
  convert_ensight_to_ijk(2) = 3;
  convert_ensight_to_ijk(3) = 2;
  convert_ensight_to_ijk(4) = 4;
  convert_ensight_to_ijk(5) = 5;
  convert_ensight_to_ijk(6) = 7;
  convert_ensight_to_ijk(7) = 6;
    
  int nodes_per_element;
  
  if(num_dim==2)
  for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++) {
    //set nodes per element
    element_select->choose_2Delem_type(Element_Types(cell_rid), elem2D);
    nodes_per_element = elem2D->num_nodes();
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      tmp_ijk_indx(node_lid) = nodes_in_elem(cell_rid, convert_ensight_to_ijk(node_lid));
    }   
        
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      nodes_in_elem(cell_rid, node_lid) = tmp_ijk_indx(node_lid);
    }
  }

  if(num_dim==3)
  for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++) {
    //set nodes per element
    element_select->choose_3Delem_type(Element_Types(cell_rid), elem);
    nodes_per_element = elem->num_nodes();
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      tmp_ijk_indx(node_lid) = nodes_in_elem(cell_rid, convert_ensight_to_ijk(node_lid));
    }   
        
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      nodes_in_elem(cell_rid, node_lid) = tmp_ijk_indx(node_lid);
    }
  }

  //synchronize device data
  dual_node_coords.sync_device();
  dual_node_coords.modify_device();
  if(restart_file){
    dual_node_densities.sync_device();
    dual_node_densities.modify_device();
  }

  //debug print element edof
  
  //std::cout << " ------------ELEMENT EDOF ON TASK " << myrank << " --------------"<<std::endl;

  //for (int ielem = 0; ielem < rnum_elem; ielem++){
    //std::cout << "elem:  " << ielem+1 << std::endl;
    //for (int lnode = 0; lnode < 8; lnode++){
        //std::cout << "{ ";
          //std::cout << lnode+1 << " = " << nodes_in_elem(ielem,lnode) + 1 << " ";
        
        //std::cout << " }"<< std::endl;
    //}
    //std::cout << std::endl;
  //}
  
 
} // end read_mesh

/* ----------------------------------------------------------------------
   Read ANSYS dat format mesh file
------------------------------------------------------------------------- */
void Implicit_Solver::read_mesh_ansys_dat(char *MESH){

  char ch;
  int num_dim = simparam->num_dim;
  int p_order = simparam->p_order;
  real_t unit_scaling = simparam->unit_scaling;
  bool restart_file = simparam->restart_file;
  int local_node_index, current_column_index;
  size_t strain_count;
  std::string skip_line, read_line, substring, token;
  std::stringstream line_parse;
  CArrayKokkos<char, array_layout, HostSpace, memory_traits> read_buffer;
  int buffer_loop, buffer_iteration, scan_loop, nodes_per_element;
  size_t read_index_start, node_rid, elem_gid;
  GO node_gid;
  real_t dof_value;
  host_vec_array node_densities;
  //Nodes_Per_Element_Type =  elements::elem_types::Nodes_Per_Element_Type;

  //read the mesh
  //PLACEHOLDER: ensight_format(MESH);
  // abaqus_format(MESH);
  // vtk_format(MESH)

  //task 0 reads file
  if(myrank==0){
    in = new std::ifstream();
    in->open(MESH);
  }

  //ANSYS dat file doesn't specify total number of nodes, which is needed for the node map.
  //First pass reads in node section to determine the maximum number of nodes, second pass distributes node data
  //The elements section header does specify element count
  num_nodes = 0;
  if(myrank==0){
    bool searching_for_nodes = true;
    //skip lines at the top with nonessential info; stop skipping when "Nodes for the whole assembly" string is reached
    while (searching_for_nodes&&in->good()) {
      getline(*in, skip_line);
      //std::cout << skip_line << std::endl;
      line_parse.clear();
      line_parse.str(skip_line);
      //stop when the NODES= string is reached
      while (!line_parse.eof()){
        line_parse >> substring;
        //std::cout << substring << std::endl;
        if(!substring.compare("Nodes")){
          searching_for_nodes = false;
          break;
        }
      } //while
      
    }
    if(searching_for_nodes){
      std::cout << "FILE FORMAT ERROR" << std::endl;
    }
    //skip 2 lines
    for (int j = 0; j < 2; j++) {
      getline(*in, skip_line);
      std::cout << skip_line << std::endl;
    } //for

    //tally node count (bug in dat files seems to print wrong node count so read is done in two passes)
    //stop when apparent "-1" zone delimiter is reacher
    searching_for_nodes = true;
    int node_tally = 0;
    while (searching_for_nodes) {
      getline(*in, read_line);
      //std::cout << read_line << std::endl;
      line_parse.clear();
      line_parse.str(read_line);
      line_parse >> substring;
      
      //std::cout << substring << std::endl;
      if(substring == "-1"){
        searching_for_nodes = false;
          break;
      }
      else{
        node_tally++;
      }
      
    }
    num_nodes = node_tally;
    std::cout << "declared node count: " << num_nodes << std::endl;
  }

  //broadcast number of nodes
  MPI_Bcast(&num_nodes,1,MPI_LONG_LONG_INT,0,world);
  
  //construct contiguous parallel row map now that we know the number of nodes
  map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(num_nodes,0,comm));
  
  //close and reopen file for second pass now that global node count is known
  if(myrank==0){
    in->close();
    //in = new std::ifstream();
    in->open(MESH);
  }

  // set the vertices in the mesh read in
  global_size_t local_nrows = map->getNodeNumElements();
  nlocal_nodes = local_nrows;
  //populate local row offset data from global data
  global_size_t min_gid = map->getMinGlobalIndex();
  global_size_t max_gid = map->getMaxGlobalIndex();
  global_size_t index_base = map->getIndexBase();
  //debug print
  //std::cout << "local node count on task: " << " " << nlocal_nodes << std::endl;
  //construct dof map that follows from the node map (used for distributed matrix and vector objects later)
  CArrayKokkos<GO, array_layout, device_type, memory_traits> local_dof_indices(nlocal_nodes*num_dim, "local_dof_indices");
  for(int i = 0; i < nlocal_nodes; i++){
    for(int j = 0; j < num_dim; j++)
    local_dof_indices(i*num_dim + j) = map->getGlobalElement(i)*num_dim + j;
  }
  
  local_dof_map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(num_nodes*num_dim,local_dof_indices.get_kokkos_view(),0,comm) );

  //allocate node storage with dual view
  dual_node_coords = dual_vec_array("dual_node_coords", nlocal_nodes,num_dim);
  if(restart_file)
    dual_node_densities = dual_vec_array("dual_node_densities", nlocal_nodes,1);

  //local variable for host view in the dual view
  host_vec_array node_coords = dual_node_coords.view_host();
  if(restart_file)
    node_densities = dual_node_densities.view_host();
  //notify that the host view is going to be modified in the file readin
  dual_node_coords.modify_host();
  if(restart_file)
    dual_node_densities.modify_host();

  //old swage method
  //mesh->init_nodes(local_nrows); // add 1 for index starting at 1
    
  std::cout << "Num nodes assigned to task " << myrank << " = " << nlocal_nodes << std::endl;

  // read the initial mesh coordinates
  // x-coords
  /*only task 0 reads in nodes and elements from the input file
  stores node data in a buffer and communicates once the buffer cap is reached
  or the data ends*/

  words_per_line = simparam->ansys_dat_node_words_per_line;
  //if(restart_file) words_per_line++;
  elem_words_per_line = simparam->ansys_dat_elem_words_per_line;

  //allocate read buffer
  read_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(BUFFER_LINES,words_per_line,MAX_WORD);

  int dof_limit = num_nodes;
  int buffer_iterations = dof_limit/BUFFER_LINES;
  if(dof_limit%BUFFER_LINES!=0) buffer_iterations++;
  
  //second pass to now read node coords with global node map defines
  if(myrank==0){
    bool searching_for_nodes = true;
    //skip lines at the top with nonessential info; stop skipping when "Nodes for the whole assembly" string is reached
    while (searching_for_nodes&&in->good()) {
      getline(*in, skip_line);
      //std::cout << skip_line << std::endl;
      line_parse.clear();
      line_parse.str(skip_line);
      //stop when the NODES= string is reached
      while (!line_parse.eof()){
        line_parse >> substring;
        //std::cout << substring << std::endl;
        if(!substring.compare("Nodes")){
          searching_for_nodes = false;
          break;
        }
      } //while
      
    }
    if(searching_for_nodes){
      std::cout << "FILE FORMAT ERROR" << std::endl;
    }
    //skip 2 lines
    for (int j = 0; j < 2; j++) {
      getline(*in, skip_line);
      std::cout << skip_line << std::endl;
    } //for
  }
  
  //read coords, also density if restarting
  read_index_start = 0;
  for(buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
    //pack buffer on rank 0
    if(myrank==0&&buffer_iteration<buffer_iterations-1){
      for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        
        for(int iword = 0; iword < words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //debug print
        //std::cout<<" "<< substring <<std::endl;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
      }
    }
    else if(myrank==0){
      buffer_loop=0;
      while(buffer_iteration*BUFFER_LINES+buffer_loop < num_nodes) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
        buffer_loop++;
      }
      
    }

    //broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
    MPI_Bcast(read_buffer.pointer(),BUFFER_LINES*words_per_line*MAX_WORD,MPI_CHAR,0,world);
    //broadcast how many nodes were read into this buffer iteration
    MPI_Bcast(&buffer_loop,1,MPI_INT,0,world);

    //debug_print
    //std::cout << "NODE BUFFER LOOP IS: " << buffer_loop << std::endl;
    //for(int iprint=0; iprint < buffer_loop; iprint++)
      //std::cout<<"buffer packing: " << std::string(&read_buffer(iprint,0,0)) << std::endl;
    //return;

    //determine which data to store in the swage mesh members (the local node data)
    //loop through read buffer
    for(scan_loop = 0; scan_loop < buffer_loop; scan_loop++){
      //set global node id (ensight specific order)
      node_gid = read_index_start + scan_loop;
      //let map decide if this node id belongs locally; if yes store data
      if(map->isNodeGlobalElement(node_gid)){
        //set local node index in this mpi rank
        node_rid = map->getLocalElement(node_gid);
        //extract nodal position from the read buffer
        //for tecplot format this is the three coords in the same line
        dof_value = atof(&read_buffer(scan_loop,1,0));
        node_coords(node_rid, 0) = dof_value * unit_scaling;
        dof_value = atof(&read_buffer(scan_loop,2,0));
        node_coords(node_rid, 1) = dof_value * unit_scaling;
        dof_value = atof(&read_buffer(scan_loop,3,0));
        node_coords(node_rid, 2) = dof_value * unit_scaling;
        //extract density if restarting
      }
    }
    read_index_start+=BUFFER_LINES;
  }

  
  //debug print of nodal data
  
  //debug print nodal positions and indices
  /*
  std::cout << " ------------NODAL POSITIONS ON TASK " << myrank << " --------------"<<std::endl;
  for (int inode = 0; inode < local_nrows; inode++){
      std::cout << "node: " << map->getGlobalElement(inode) + 1 << " { ";
    for (int istride = 0; istride < num_dim; istride++){
       std::cout << node_coords(inode,istride) << " , ";
    }
    //std::cout << node_densities(inode,0);
    std::cout << " }"<< std::endl;
  }
  */

  //check that local assignments match global total

  
  //read in element info
  //seek element connectivity zone
  int etype_index = 0;
  if(myrank==0){
    bool searching_for_elements = true;
    //skip lines at the top with nonessential info; stop skipping when "Nodes for the whole assembly" string is reached
    while (searching_for_elements&&in->good()) {
      getline(*in, skip_line);
      //std::cout << skip_line << std::endl;
      line_parse.clear();
      line_parse.str(skip_line);
      //stop when the NODES= string is reached
      while (!line_parse.eof()){
        line_parse >> substring;
        //std::cout << substring << std::endl;
        if(!substring.compare("Elements")){
          searching_for_elements = false;
          break;
        }
      } //while
      
    }
    if(searching_for_elements){
      std::cout << "FILE FORMAT ERROR" << std::endl;
    }

    //read in element type from following line
    getline(*in, read_line);
    std::cout << read_line << std::endl;
    line_parse.clear();
    line_parse.str(read_line);
    line_parse >> substring;
    //std::cout << substring << std::endl;
    if(!substring.compare("et,1,185")){
      //Hex8 type
      etype_index = 1;

    }
    else if(!substring.compare("et,1,186")){
      //Hex20 type
      etype_index = 2;
    }
    else{
      etype_index = 0;
    }
     //for
    
    //seek element count line
    //read in element count from the following line
    getline(*in, read_line);
    std::cout << read_line << std::endl;
    line_parse.clear();
    line_parse.str(read_line);
    line_parse >> substring;
    //parse element line out of jumble of comma delimited entries
    line_parse.clear();
    line_parse.str(substring);
    while(line_parse.good()){
      getline(line_parse, token, ',');
    }
    //element count should be the last token read in
    num_elem = std::stoi(token);

    //skip line
    for (int j = 0; j < 1; j++) {
      getline(*in, skip_line);
      std::cout << skip_line << std::endl;
    }
  }

  //broadcast element type
  MPI_Bcast(&etype_index,1,MPI_INT,0,world);

  elements::elem_types::elem_type mesh_element_type;
  int elem_words_per_line_no_nodes = elem_words_per_line;
  if(etype_index==1){
    mesh_element_type = elements::elem_types::Hex8;
    nodes_per_element = 8;
    elem_words_per_line += 8;
  }
  else if(etype_index==2){
    mesh_element_type = elements::elem_types::Hex20;
    nodes_per_element = 20;
    elem_words_per_line += 20;
  }
  else if(etype_index==3){
    mesh_element_type = elements::elem_types::Hex32;
    nodes_per_element = 32;
    elem_words_per_line += 32;
  }
  else{
    *fos << "ERROR: ANSYS ELEMENT TYPE NOT FOUND OR RECOGNIZED" << std::endl;
    exit_solver(0);
  }
  
  //broadcast number of elements
  MPI_Bcast(&num_elem,1,MPI_LONG_LONG_INT,0,world);
  
  *fos << "declared element count: " << num_elem << std::endl;
  //std::cout<<"before initial mesh initialization"<<std::endl;
  
  //read in element connectivity
  //we're gonna reallocate for the words per line expected for the element connectivity
  read_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(BUFFER_LINES,elem_words_per_line,MAX_WORD); 
  CArrayKokkos<int, array_layout, HostSpace, memory_traits> node_store(nodes_per_element);

  //calculate buffer iterations to read number of lines
  buffer_iterations = num_elem/BUFFER_LINES;
  int assign_flag;

  //dynamic buffer used to store elements before we know how many this rank needs
  std::vector<size_t> element_temp(BUFFER_LINES*elem_words_per_line);
  std::vector<size_t> global_indices_temp(BUFFER_LINES);
  size_t buffer_max = BUFFER_LINES*elem_words_per_line;
  size_t indices_buffer_max = BUFFER_LINES;

  if(num_elem%BUFFER_LINES!=0) buffer_iterations++;
  read_index_start = 0;
  //std::cout << "ELEMENT BUFFER ITERATIONS: " << buffer_iterations << std::endl;
  rnum_elem = 0;
  //std::cout << "BUFFER ITERATIONS IS: " << buffer_iterations << std::endl;
  for(buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
    //pack buffer on rank 0
    if(myrank==0&&buffer_iteration<buffer_iterations-1){
      for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < elem_words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
      }
    }
    else if(myrank==0){
      buffer_loop=0;
      while(buffer_iteration*BUFFER_LINES+buffer_loop < num_elem) {
        getline(*in,read_line);
        line_parse.clear();
        line_parse.str(read_line);
        for(int iword = 0; iword < elem_words_per_line; iword++){
        //read portions of the line into the substring variable
        line_parse >> substring;
        //assign the substring variable as a word of the read buffer
        strcpy(&read_buffer(buffer_loop,iword,0),substring.c_str());
        }
        buffer_loop++;
        //std::cout<<" "<< node_coords(node_gid, 0)<<std::endl;
      }
    }

    //broadcast buffer to all ranks; each rank will determine which nodes in the buffer belong
    MPI_Bcast(read_buffer.pointer(),BUFFER_LINES*elem_words_per_line*MAX_WORD,MPI_CHAR,0,world);
    //broadcast how many nodes were read into this buffer iteration
    MPI_Bcast(&buffer_loop,1,MPI_INT,0,world);
    
    //store element connectivity that belongs to this rank
    //loop through read buffer
    for(scan_loop = 0; scan_loop < buffer_loop; scan_loop++){
      //set global node id (ensight specific order)
      elem_gid = read_index_start + scan_loop;
      //add this element to the local list if any of its nodes belong to this rank according to the map
      //get list of nodes for each element line and check if they belong to the map
      assign_flag = 0;
      for(int inode = elem_words_per_line_no_nodes; inode < elem_words_per_line; inode++){
        //as we loop through the nodes belonging to this element we store them
        //if any of these nodes belongs to this rank this list is used to store the element locally
        node_gid = atoi(&read_buffer(scan_loop,inode,0));
        node_store(inode-elem_words_per_line_no_nodes) = node_gid - 1; //subtract 1 since file index start is 1 but code expects 0
        //first we add the elements to a dynamically allocated list
        if(map->isNodeGlobalElement(node_gid-1)&&!assign_flag){
          assign_flag = 1;
          rnum_elem++;
        }
      }

      if(assign_flag){
        for(int inode = 0; inode < nodes_per_element; inode++){
          if((rnum_elem-1)*nodes_per_element + inode>=buffer_max){ 
            element_temp.resize((rnum_elem-1)*nodes_per_element + inode + BUFFER_LINES*nodes_per_element);
            buffer_max = (rnum_elem-1)*nodes_per_element + inode + BUFFER_LINES*nodes_per_element;
          }
          element_temp[(rnum_elem-1)*nodes_per_element + inode] = node_store(inode); 
          //std::cout << "VECTOR STORAGE FOR ELEM " << rnum_elem << " ON TASK " << myrank << " NODE " << inode+1 << " IS " << node_store(inode) + 1 << std::endl;
        }
        //assign global element id to temporary list
        if(rnum_elem-1>=indices_buffer_max){ 
          global_indices_temp.resize(rnum_elem-1 + BUFFER_LINES);
          indices_buffer_max = rnum_elem-1 + BUFFER_LINES;
        }
        global_indices_temp[rnum_elem-1] = elem_gid;
      }
    }
    read_index_start+=BUFFER_LINES;
  }
  
  //check if ANSYS file has boundary and loading condition zones
  bool No_Conditions = true;
  if(myrank==0){
    if(in->good())
      before_condition_header = in->tellg();
    bool searching_for_conditions = true;
    //skip lines at the top with nonessential info; stop skipping when "Fixed Supports or Pressure" string is reached
    while (searching_for_conditions&&in->good()) {
      getline(*in, skip_line);
      //std::cout << skip_line << std::endl;
      line_parse.clear();
      line_parse.str(skip_line);
      //stop when the NODES= string is reached
      while (!line_parse.eof()){
        line_parse >> substring;
        //std::cout << substring << std::endl;
        if(!substring.compare("Supports")||!substring.compare("Pressure")){
          No_Conditions = searching_for_conditions = false;
          break;
        }
      } //while
      
    } //while
  }
  
  //broadcast search condition
  MPI_Bcast(&No_Conditions,1,MPI_CXX_BOOL,0,world);

  //flag elasticity fea module for boundary/loading conditions readin that remains
  if(!No_Conditions){
    //look for elasticity module in Simulation Parameters data; if not declared add the module
    int nfea_modules = simparam->nfea_modules;
    bool elasticity_found = false;
    std::vector<std::string> FEA_Module_List = simparam->FEA_Module_List;
    for(int imodule = 0; imodule < nfea_modules; imodule++){
      if(FEA_Module_List[imodule]=="Elasticity"){ 
        elasticity_found = true;
        simparam->fea_module_must_read[imodule] = true;
      }
    }
    
    //add Elasticity module to requested modules in the Simulation Parameters data
    if(!elasticity_found){
      if(nfea_modules==simparam->FEA_Module_List.capacity()){
        simparam->FEA_Module_List.push_back("Elasticity");
        simparam->fea_module_must_read.push_back(true);
        simparam->nfea_modules++;
      }
      else{
        simparam->FEA_Module_List[nfea_modules] = "Elasticity";
        simparam->fea_module_must_read[nfea_modules]= true;
        simparam->nfea_modules++;
      }
    }

    
  }

  // Close mesh input file if no further readin is done by FEA modules for conditions
  if(myrank==0&&No_Conditions){
    in->close();
  }

  std::cout << "RNUM ELEMENTS IS: " << rnum_elem << std::endl;
  //copy temporary element storage to multivector storage
  Element_Types = CArrayKokkos<elements::elem_types::elem_type, array_layout, HostSpace, memory_traits>(rnum_elem);

  //set element object pointer
  if(simparam->num_dim==2){
    element_select->choose_2Delem_type(mesh_element_type, elem2D);
     max_nodes_per_element = elem2D->num_nodes();
  }
  else if(simparam->num_dim==3){
    element_select->choose_3Delem_type(mesh_element_type, elem);
     max_nodes_per_element = elem->num_nodes();
  }

  //1 type per mesh for now
  for(int ielem = 0; ielem < rnum_elem; ielem++)
    Element_Types(ielem) = mesh_element_type;

  dual_nodes_in_elem = dual_elem_conn_array("dual_nodes_in_elem", rnum_elem, max_nodes_per_element);
  nodes_in_elem = dual_nodes_in_elem.view_host();
  dual_nodes_in_elem.modify_host();

  for(int ielem = 0; ielem < rnum_elem; ielem++)
    for(int inode = 0; inode < nodes_per_element; inode++)
      nodes_in_elem(ielem, inode) = element_temp[ielem*nodes_per_element + inode];

  //view storage for all local elements connected to local nodes on this rank
  CArrayKokkos<GO, array_layout, device_type, memory_traits> All_Element_Global_Indices(rnum_elem);

  //copy temporary global indices storage to view storage
  for(int ielem = 0; ielem < rnum_elem; ielem++)
    All_Element_Global_Indices(ielem) = global_indices_temp[ielem];

  //debug print element edof
  /*
  std::cout << " ------------ELEMENT EDOF ON TASK " << myrank << " --------------"<<std::endl;
  
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    std::cout << "elem:  " << All_Element_Global_Indices(ielem)+1 << std::endl;
    for (int lnode = 0; lnode < 8; lnode++){
        std::cout << "{ ";
          std::cout << lnode+1 << " = " << nodes_in_elem(ielem,lnode) + 1 << " ";
        
        std::cout << " }"<< std::endl;
    }
    std::cout << std::endl;
  }
  */

  //delete temporary element connectivity and index storage
  std::vector<size_t>().swap(element_temp);
  std::vector<size_t>().swap(global_indices_temp);
  
  //construct overlapping element map (since different ranks can own the same elements due to the local node map)
  all_element_map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),All_Element_Global_Indices.get_kokkos_view(),0,comm));


  //element type selection (subject to change)
  // ---- Set Element Type ---- //
  // allocate element type memory
  //elements::elem_type_t* elem_choice;

  int NE = 1; // number of element types in problem

  // Convert ijk index system to the finite element numbering convention
  // for vertices in cell
  CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> convert_ensight_to_ijk(max_nodes_per_element);
  CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> tmp_ijk_indx(max_nodes_per_element);
  convert_ensight_to_ijk(0) = 0;
  convert_ensight_to_ijk(1) = 1;
  convert_ensight_to_ijk(2) = 3;
  convert_ensight_to_ijk(3) = 2;
  convert_ensight_to_ijk(4) = 4;
  convert_ensight_to_ijk(5) = 5;
  convert_ensight_to_ijk(6) = 7;
  convert_ensight_to_ijk(7) = 6;
  
  if(num_dim==2)
  for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++) {
    //set nodes per element
    element_select->choose_2Delem_type(Element_Types(cell_rid), elem2D);
    nodes_per_element = elem2D->num_nodes();
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      tmp_ijk_indx(node_lid) = nodes_in_elem(cell_rid, convert_ensight_to_ijk(node_lid));
    }   
        
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      nodes_in_elem(cell_rid, node_lid) = tmp_ijk_indx(node_lid);
    }
  }

  if(num_dim==3)
  for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++) {
    //set nodes per element
    element_select->choose_3Delem_type(Element_Types(cell_rid), elem);
    nodes_per_element = elem->num_nodes();
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      tmp_ijk_indx(node_lid) = nodes_in_elem(cell_rid, convert_ensight_to_ijk(node_lid));
    }   
        
    for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
      nodes_in_elem(cell_rid, node_lid) = tmp_ijk_indx(node_lid);
    }
  }

  //synchronize device data
  dual_node_coords.sync_device();
  dual_node_coords.modify_device();
  if(restart_file){
    dual_node_densities.sync_device();
    dual_node_densities.modify_device();
  }
 
} // end read_mesh

/* ----------------------------------------------------------------------
   Initialize Ghost and Non-Overlapping Element Maps
------------------------------------------------------------------------- */

void Implicit_Solver::init_maps(){
  char ch;
  int num_dim = simparam->num_dim;
  int p_order = simparam->p_order;
  real_t unit_scaling = simparam->unit_scaling;
  int local_node_index, current_column_index;
  size_t strain_count;
  std::stringstream line_parse;
  CArrayKokkos<char, array_layout, HostSpace, memory_traits> read_buffer;
  int nodes_per_element;
  GO node_gid;
  
  if(rnum_elem >= 1) {

    //Construct set of ghost nodes; start with a buffer with upper limit
    size_t buffer_limit = 0;
    if(num_dim==2)
    for(int ielem = 0; ielem < rnum_elem; ielem++){
      element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
      buffer_limit += elem2D->num_nodes();
    }

    if(num_dim==3)
    for(int ielem = 0; ielem < rnum_elem; ielem++){
      element_select->choose_3Delem_type(Element_Types(ielem), elem);
      buffer_limit += elem->num_nodes();
    }

    CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> ghost_node_buffer(buffer_limit);
    std::set<GO> ghost_node_set;

    //search through local elements for global node indices not owned by this MPI rank
    if(num_dim==2)
    for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++) {
      //set nodes per element
      element_select->choose_2Delem_type(Element_Types(cell_rid), elem2D);
      nodes_per_element = elem2D->num_nodes();  
      for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
        node_gid = nodes_in_elem(cell_rid, node_lid);
        if(!map->isNodeGlobalElement(node_gid)) ghost_node_set.insert(node_gid);
      }
    }

    if(num_dim==3)
    for (int cell_rid = 0; cell_rid < rnum_elem; cell_rid++) {
      //set nodes per element
      element_select->choose_3Delem_type(Element_Types(cell_rid), elem);
      nodes_per_element = elem->num_nodes();  
      for (int node_lid = 0; node_lid < nodes_per_element; node_lid++){
        node_gid = nodes_in_elem(cell_rid, node_lid);
        if(!map->isNodeGlobalElement(node_gid)) ghost_node_set.insert(node_gid);
      }
    }

    //by now the set contains, with no repeats, all the global node indices that are ghosts for this rank
    //now pass the contents of the set over to a CArrayKokkos, then create a map to find local ghost indices from global ghost indices
    nghost_nodes = ghost_node_set.size();
    ghost_nodes = CArrayKokkos<GO, Kokkos::LayoutLeft, node_type::device_type>(nghost_nodes, "ghost_nodes");
    ghost_node_ranks = CArrayKokkos<int, array_layout, device_type, memory_traits>(nghost_nodes, "ghost_node_ranks");
    int ighost = 0;
    auto it = ghost_node_set.begin();
    while(it!=ghost_node_set.end()){
      ghost_nodes(ighost++) = *it;
      it++;
    }

    //debug print of ghost nodes
    //std::cout << " GHOST NODE SET ON TASK " << myrank << std::endl;
    //for(int i = 0; i < nghost_nodes; i++)
      //std::cout << "{" << i + 1 << "," << ghost_nodes(i) + 1 << "}" << std::endl;

    //find which mpi rank each ghost node belongs to and store the information in a CArrayKokkos
    //allocate Teuchos Views since they are the only input available at the moment in the map definitions
    Teuchos::ArrayView<const GO> ghost_nodes_pass(ghost_nodes.get_kokkos_view().data(), nghost_nodes);
    Teuchos::ArrayView<int> ghost_node_ranks_pass(ghost_node_ranks.get_kokkos_view().data(), nghost_nodes);
    map->getRemoteIndexList(ghost_nodes_pass, ghost_node_ranks_pass);
    
    //debug print of ghost nodes
    //std::cout << " GHOST NODE MAP ON TASK " << myrank << std::endl;
    //for(int i = 0; i < nghost_nodes; i++)
      //std::cout << "{" << i + 1 << "," << global2local_map.get(ghost_nodes(i)) + 1 << "}" << std::endl;

  }

  // create a Map for ghost node indices
  ghost_node_map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),ghost_nodes.get_kokkos_view(),0,comm));
    
  // Create reference element
  //ref_elem->init(p_order, num_dim, elem->num_basis());
  //std::cout<<"done with ref elem"<<std::endl;

  //communicate ghost node positions; construct multivector distributed object using local node data

  //construct array for all indices (ghost + local)
  nall_nodes = nlocal_nodes + nghost_nodes;
  CArrayKokkos<GO, array_layout, device_type, memory_traits> all_node_indices(nall_nodes, "all_node_indices");
  for(int i = 0; i < nall_nodes; i++){
    if(i<nlocal_nodes) all_node_indices(i) = map->getGlobalElement(i);
    else all_node_indices(i) = ghost_nodes(i-nlocal_nodes);
  }
  
  //debug print of node indices
  //for(int inode=0; inode < index_counter; inode++)
  //std::cout << " my_reduced_global_indices " << my_reduced_global_indices(inode) <<std::endl;
  
  // create a Map for all the node indices (ghost + local)
  all_node_map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),all_node_indices.get_kokkos_view(),0,comm));

  //remove elements from the local set so that each rank has a unique set of global ids
  
  //local elements belonging to the non-overlapping element distribution to each rank with buffer
  CArrayKokkos<GO, array_layout, device_type, memory_traits> Initial_Element_Global_Indices(rnum_elem);
  size_t nonoverlapping_count = 0;
  int my_element_flag;
  //loop through local element set
  if(num_dim == 2)
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
    nodes_per_element = elem2D->num_nodes();
    my_element_flag = 1;
    for (int lnode = 0; lnode < nodes_per_element; lnode++){
      node_gid = nodes_in_elem(ielem, lnode);
      if(ghost_node_map->isNodeGlobalElement(node_gid)){
        local_node_index = ghost_node_map->getLocalElement(node_gid);
        if(ghost_node_ranks(local_node_index) < myrank) my_element_flag = 0;
      }
    }
    if(my_element_flag){
      Initial_Element_Global_Indices(nonoverlapping_count++) = all_element_map->getGlobalElement(ielem);
    }
  }

  if(num_dim == 3)
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_3Delem_type(Element_Types(ielem), elem);
    nodes_per_element = elem->num_nodes();
    my_element_flag = 1;
    for (int lnode = 0; lnode < nodes_per_element; lnode++){
      node_gid = nodes_in_elem(ielem, lnode);
      if(ghost_node_map->isNodeGlobalElement(node_gid)){
        local_node_index = ghost_node_map->getLocalElement(node_gid);
        if(ghost_node_ranks(local_node_index) < myrank) my_element_flag = 0;
      }
    }
    if(my_element_flag){
      Initial_Element_Global_Indices(nonoverlapping_count++) = all_element_map->getGlobalElement(ielem);
    }
  }

  //copy over from buffer to compressed storage
  CArrayKokkos<GO, array_layout, device_type, memory_traits> Element_Global_Indices(nonoverlapping_count);
  for(int ibuffer = 0; ibuffer < nonoverlapping_count; ibuffer++)
  Element_Global_Indices(ibuffer) = Initial_Element_Global_Indices(ibuffer);

  //create nonoverlapping element map
  element_map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),Element_Global_Indices.get_kokkos_view(),0,comm));

  //construct dof map that follows from the all_node map (used for distributed matrix and vector objects later)
  CArrayKokkos<GO, array_layout, device_type, memory_traits> all_dof_indices(nall_nodes*num_dim, "all_dof_indices");
  for(int i = 0; i < nall_nodes; i++){
    for(int j = 0; j < num_dim; j++)
    all_dof_indices(i*num_dim + j) = all_node_map->getGlobalElement(i)*num_dim + j;
  }
  
  //pass invalid global count so the map reduces the global count automatically
  all_dof_map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),all_dof_indices.get_kokkos_view(),0,comm) );

  //debug print of map
  //debug print
  std::ostream &out = std::cout;
  Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  //if(myrank==0)
  //*fos << "Ghost Node Map :" << std::endl;
  //all_node_map->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  //std::fflush(stdout);

  //Count how many elements connect to each local node
  node_nconn_distributed = Teuchos::rcp(new MCONN(map, 1));
  host_elem_conn_array node_nconn = node_nconn_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  for(int inode = 0; inode < nlocal_nodes; inode++)
    node_nconn(inode,0) = 0;

  for(int ielem = 0; ielem < rnum_elem; ielem++){
    for(int inode = 0; inode < nodes_per_element; inode++){
      node_gid = nodes_in_elem(ielem, inode);
      if(map->isNodeGlobalElement(node_gid))
        node_nconn(map->getLocalElement(node_gid),0)++;
    }
  }

  //create distributed multivector of the local node data and all (local + ghost) node storage
  node_coords_distributed = Teuchos::rcp(new MV(map, dual_node_coords));
  all_node_coords_distributed = Teuchos::rcp(new MV(all_node_map, num_dim));
  
  //debug print
  //std::ostream &out = std::cout;
  //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  //if(myrank==0)
  //*fos << "Node Data :" << std::endl;
  //node_coords_distributed->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  //std::fflush(stdout);

  //create import object using local node indices map and all indices map
  Tpetra::Import<LO, GO> importer(map, all_node_map);

  //comms to get ghosts
  all_node_coords_distributed->doImport(*node_coords_distributed, importer, Tpetra::INSERT);
  //all_node_nconn_distributed->doImport(*node_nconn_distributed, importer, Tpetra::INSERT);
  
  dual_nodes_in_elem.sync_device();
  dual_nodes_in_elem.modify_device();
  //construct distributed element connectivity multivector
  nodes_in_elem_distributed = Teuchos::rcp(new MCONN(all_element_map, dual_nodes_in_elem));

  //debug print
  //std::ostream &out = std::cout;
  //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  //if(myrank==0)
  //*fos << "Element Connectivity :" << std::endl;
  //nodes_in_elem_distributed->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  //std::fflush(stdout);

  //debug print
  //std::ostream &out = std::cout;
  //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  //if(myrank==0)
  //*fos << "Number of Elements Connected to Each Node :" << std::endl;
  //all_node_nconn_distributed->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  //std::fflush(stdout);

  //debug print
  //std::ostream &out = std::cout;
  //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  //if(myrank==0)
  //*fos << "Node Data with Ghosts :" << std::endl;
  //all_node_coords_distributed->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  //std::fflush(stdout);

  //debug print of views node indices
  //std::cout << "Local View of All Nodes on Task " << myrank <<std::endl;
  //for(int inode=0; inode < all_node_map->getNodeNumElements(); inode++){
    //std::cout << "node "<<all_node_map->getGlobalElement(inode) << " } " ;
    //std::cout << dual_all_node_coords.view_host()(inode,0) << " " << dual_all_node_coords.view_host()(inode,1) << " " << dual_all_node_coords.view_host()(inode,2) << " " << std::endl;
  //}
     
  //std::cout << "number of patches = " << mesh->num_patches() << std::endl;
  if(myrank == 0)
  std::cout << "End of map setup " << std::endl;
}

/* ----------------------------------------------------------------------
   Construct list of objects for FEA modules 
------------------------------------------------------------------------- */

void Implicit_Solver::FEA_module_setup(){
  nfea_modules = simparam->nfea_modules;
  std::vector<std::string> FEA_Module_List = simparam->FEA_Module_List;
  fea_module_must_read = simparam->fea_module_must_read;
  //allocate lists to size
  fea_module_types = std::vector<std::string>(nfea_modules);
  fea_modules = std::vector<FEA_Module*>(nfea_modules);
  bool module_found = false;
  
  //list should not have repeats since that was checked by simulation parameters setups
  for(int imodule = 0; imodule < nfea_modules; imodule++){
    //decides which FEA module objects to setup based on string.
    //automate selection list later; use std::map maybe?
    if(FEA_Module_List[imodule] == "Elasticity"){
      fea_module_types[imodule] = "Elasticity";
      fea_modules[imodule] = new FEA_Module_Elasticity(this);
      module_found = true;
      displacement_module = imodule;
      //debug print
      *fos << " ELASTICITY MODULE ALLOCATED AS " <<imodule << std::endl;
      
    }
    else if(FEA_Module_List[imodule] == "Inertial"){
      fea_module_types[imodule] = "Inertial";
      fea_modules[imodule] = new FEA_Module_Inertial(this);
      module_found = true;
      displacement_module = imodule;
      //debug print
      *fos << " INERTIAL MODULE ALLOCATED AS " <<imodule << std::endl;
      
    }
    else if(FEA_Module_List[imodule] == "Heat_Conduction"){
      fea_module_types[imodule] = "Heat_Conduction";
      fea_modules[imodule] = new FEA_Module_Heat_Conduction(this);
      module_found = true; 
      //debug print
      *fos << " HEAT MODULE ALLOCATED AS " <<imodule << std::endl;
    }
    else{
      *fos << "PROGRAM IS ENDING DUE TO ERROR; UNDEFINED FEA MODULE REQUESTED WITH NAME \"" <<FEA_Module_List[imodule]<<"\"" << std::endl;
      exit_solver(0);
    }
  }
}

/* ----------------------------------------------------------------------
   Setup Optimization Problem Object, Relevant Objective, and Constraints
------------------------------------------------------------------------- */

void Implicit_Solver::setup_optimization_problem(){
  int num_dim = simparam->num_dim;
  bool nodal_density_flag = simparam->nodal_density_flag;
  int nTO_modules = simparam->nTO_modules;
  int nmulti_objective_modules = simparam->nmulti_objective_modules;
  std::vector<std::string> TO_Module_List = simparam->TO_Module_List;
  std::vector<std::string> FEA_Module_List = simparam->FEA_Module_List;
  std::vector<int> TO_Module_My_FEA_Module = simparam->TO_Module_My_FEA_Module;
  std::vector<int> Multi_Objective_Modules = simparam->Multi_Objective_Modules;
  std::vector<real_t> Multi_Objective_Weights = simparam->Multi_Objective_Weights;
  std::vector<std::vector<real_t>> Function_Arguments = simparam->Function_Arguments;
  std::vector<Simulation_Parameters_Topology_Optimization::function_type> TO_Function_Type = simparam->TO_Function_Type;
  std::vector<ROL::Ptr<ROL::Objective<real_t>>> Multi_Objective_Terms;

  std::string constraint_base, constraint_name;
  std::stringstream number_union;
  CArrayKokkos<GO, array_layout, device_type, memory_traits> Surface_Nodes;
  GO current_node_index;
  LO local_node_index;
  int num_bdy_patches_in_set;
  size_t node_id, patch_id, module_id;
  int num_boundary_sets;
  int current_element_index, local_surface_id;
  const_host_vec_array design_densities;
  typedef ROL::TpetraMultiVector<real_t,LO,GO,node_type> ROL_MV;
  
  // fill parameter list with desired algorithmic options or leave as default
  // Read optimization input parameter list.
  std::string filename = "input_ex01.xml";
  auto parlist = ROL::getParametersFromXmlFile( filename );
  //ROL::ParameterList parlist;

  //Design variables to optimize
  ROL::Ptr<ROL::Vector<real_t>> x;
  if(nodal_density_flag)
    x = ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO>>(design_node_densities_distributed);
  else
    x = ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO>>(Global_Element_Densities);
  
  //Instantiate (the one) objective function for the problem
  ROL::Ptr<ROL::Objective<real_t>> obj;
  bool objective_declared = false;
  for(int imodule = 0; imodule < nTO_modules; imodule++){
    if(TO_Function_Type[imodule] == Simulation_Parameters_Topology_Optimization::OBJECTIVE){
      //check if previous module already defined an objective, there must be one objective module
      if(objective_declared){
        *fos << "PROGRAM IS ENDING DUE TO ERROR; ANOTHER OBJECTIVE FUNCTION WITH NAME \"" <<TO_Module_List[imodule] <<"\" ATTEMPTED TO REPLACE A PREVIOUS OBJECTIVE; THERE MUST BE ONE OBJECTIVE." << std::endl;
          exit_solver(0);
      }
      if(TO_Module_List[imodule] == "Strain_Energy_Minimize"){
        //debug print
        *fos << " STRAIN ENERGY OBJECTIVE EXPECTS FEA MODULE INDEX " <<TO_Module_My_FEA_Module[imodule] << std::endl;
        obj = ROL::makePtr<StrainEnergyMinimize_TopOpt>(fea_modules[TO_Module_My_FEA_Module[imodule]], nodal_density_flag);
      }
      else if(TO_Module_List[imodule] == "Heat_Capacity_Potential_Minimize"){
        //debug print
        *fos << " HEAT CAPACITY POTENTIAL OBJECTIVE EXPECTS FEA MODULE INDEX " <<TO_Module_My_FEA_Module[imodule] << std::endl;
        obj = ROL::makePtr<HeatCapacityPotentialMinimize_TopOpt>(fea_modules[TO_Module_My_FEA_Module[imodule]], nodal_density_flag);
      }
      //Multi-Objective case
      else if(TO_Module_List[imodule] == "Multi_Objective"){
        //allocate vector of Objective Functions to pass
        Multi_Objective_Terms = std::vector<ROL::Ptr<ROL::Objective<real_t>>>(nmulti_objective_modules);
        for(int imulti = 0; imulti < nmulti_objective_modules; imulti++){
          //get module index for objective term
          module_id = Multi_Objective_Modules[imulti];
          if(TO_Module_List[module_id] == "Strain_Energy_Minimize"){
            //debug print
            *fos << " STRAIN ENERGY OBJECTIVE EXPECTS FEA MODULE INDEX " <<TO_Module_My_FEA_Module[module_id] << std::endl;
            Multi_Objective_Terms[imulti] = ROL::makePtr<StrainEnergyMinimize_TopOpt>(fea_modules[TO_Module_My_FEA_Module[module_id]], nodal_density_flag);
          }
          else if(TO_Module_List[module_id] == "Heat_Capacity_Potential_Minimize"){
            //debug print
            *fos << " HEAT CAPACITY POTENTIAL OBJECTIVE EXPECTS FEA MODULE INDEX " <<TO_Module_My_FEA_Module[module_id] << std::endl;
            Multi_Objective_Terms[imulti] = ROL::makePtr<HeatCapacityPotentialMinimize_TopOpt>(fea_modules[TO_Module_My_FEA_Module[module_id]], nodal_density_flag);
          }
        }
        //allocate multi objective function
        obj = ROL::makePtr<MultiObjective_TopOpt>(Multi_Objective_Terms, Multi_Objective_Weights);
      }
      else{
        *fos << "PROGRAM IS ENDING DUE TO ERROR; UNDEFINED OBJECTIVE FUNCTION REQUESTED WITH NAME \"" <<TO_Module_List[imodule] <<"\"" << std::endl;
        exit_solver(0);
      }
      objective_declared = true;
    }
  }
  
  //optimization problem interface that can have constraints added to it before passing to solver object
   ROL::Ptr<ROL::Problem<real_t>> problem = ROL::makePtr<ROL::Problem<real_t>>(obj,x);
  
  //ROL::Ptr<ROL::Constraint<double>>     lin_icon = ROL::makePtr<MyLinearInequalityConstraint<double>>();
  //ROL::Ptr<ROL::Vector<double>>         lin_imul = ROL::makePtr<MyLinearInequalityConstraintMultiplier<double>>();
  //ROL::Ptr<ROL:BoundConstraint<double>> lin_ibnd = ROL::makePtr<MyLinearInequalityConstraintBound<double>>();
  //problem.addLinearConstraint("Linear Inequality Constraint",lin_icon,lin_imul,lin_ibnd);
    
  // TypeG (generally constrained) specification
  //ROL::Ptr<ROL::Constraint<double>> econ = ROL::makePtr<MyEqualityConstraint<double>>();
  //ROL::Ptr<ROL::Vector<double>>     emul = ROL::makePtr<MyEqualityConstraintMultiplier<double>>();
  //problem.addConstraint("Equality Constraint",econ,emul);
  
  //ROL::Ptr<ROL::Constraint<real_t>> ineq_constraint = ROL::makePtr<MassConstraint_TopOpt>(this, nodal_density_flag);
  //ROL::Ptr<ROL::Constraint<double>> lin_econ = ROL::makePtr<MyLinearEqualityConstraint<double>>();
  //ROL::Ptr<ROL::Vector<double>      lin_emul = ROL::makePtr<MyLinearEqualityConstraintMultiplier<double>>();
  //problem.addLinearConstraint("Linear Equality Constraint",lin_econ,lin_mul);

  //ROL::Ptr<ROL::Constraint<real_t>> ineq_constraint = ROL::makePtr<MassConstraint_TopOpt>(fea_elasticity, nodal_density_flag);
  //problem->addConstraint("Inequality Constraint",ineq_constraint,constraint_mul,constraint_bnd);
  //problem->addConstraint("equality Constraint 2",eq_constraint2,constraint_mul2);
  //problem->addConstraint("equality Constraint 3",eq_constraint3,constraint_mul3);
  //problem->addLinearConstraint("Equality Constraint",eq_constraint,constraint_mul);
  
  for(int imodule = 0; imodule < nTO_modules; imodule++){
    number_union.str(constraint_base);
    number_union << imodule + 1;
    constraint_name = number_union.str();
    ROL::Ptr<std::vector<real_t> > li_ptr = ROL::makePtr<std::vector<real_t>>(1,0.0);
    ROL::Ptr<ROL::Vector<real_t> > constraint_mul = ROL::makePtr<ROL::StdVector<real_t>>(li_ptr);
    if(TO_Function_Type[imodule] == Simulation_Parameters_Topology_Optimization::EQUALITY_CONSTRAINT){
      //pointers are reference counting
      ROL::Ptr<ROL::Constraint<real_t>> eq_constraint;
      if(TO_Module_List[imodule]=="Mass_Constraint"){
        
        *fos << " MASS CONSTRAINT EXPECTS FEA MODULE INDEX " <<TO_Module_My_FEA_Module[imodule] << std::endl;
        eq_constraint = ROL::makePtr<MassConstraint_TopOpt>(fea_modules[TO_Module_My_FEA_Module[imodule]], nodal_density_flag, Function_Arguments[imodule][0], false);
      }
      else if(TO_Module_List[imodule]=="Moment_of_Inertia_Constraint"){
        *fos << " MOMENT OF INERTIA CONSTRAINT EXPECTS FEA MODULE INDEX " <<TO_Module_My_FEA_Module[imodule] << std::endl;
        eq_constraint = ROL::makePtr<MomentOfInertiaConstraint_TopOpt>(fea_modules[TO_Module_My_FEA_Module[imodule]], nodal_density_flag, Function_Arguments[imodule][1], Function_Arguments[imodule][0], false);
      }
      else if(TO_Module_List[imodule]=="Strain_Energy_Constraint"){    
        *fos << " STRAIN ENERGY CONSTRAINT EXPECTS FEA MODULE INDEX " <<TO_Module_My_FEA_Module[imodule] << std::endl;
        eq_constraint = ROL::makePtr<StrainEnergyConstraint_TopOpt>(fea_modules[TO_Module_My_FEA_Module[imodule]], nodal_density_flag, Function_Arguments[imodule][0], false);
      }
      else if(TO_Module_List[imodule]=="Heat_Capacity_Potential_Constraint"){
        *fos << " HEAT CAPACITY POTENTIAL CONSTRAINT EXPECTS FEA MODULE INDEX " <<TO_Module_My_FEA_Module[imodule] << std::endl;
        eq_constraint = ROL::makePtr<HeatCapacityPotentialConstraint_TopOpt>(fea_modules[TO_Module_My_FEA_Module[imodule]], nodal_density_flag, Function_Arguments[imodule][0], false);
      }
      else{
        *fos << "PROGRAM IS ENDING DUE TO ERROR; UNDEFINED EQUALITY CONSTRAINT FUNCTION REQUESTED WITH NAME \"" <<TO_Module_List[imodule] <<"\"" << std::endl;
        exit_solver(0);
      }
      *fos << " ADDING CONSTRAINT " << constraint_name << std::endl;
      problem->addConstraint(constraint_name, eq_constraint, constraint_mul);
    }

    if(TO_Function_Type[imodule] == Simulation_Parameters_Topology_Optimization::INEQUALITY_CONSTRAINT){
      //pointers are reference counting
      ROL::Ptr<ROL::Constraint<real_t>> ineq_constraint;
      ROL::Ptr<std::vector<real_t> > ll_ptr = ROL::makePtr<std::vector<real_t>>(1,Function_Arguments[imodule][0]);
      ROL::Ptr<std::vector<real_t> > lu_ptr = ROL::makePtr<std::vector<real_t>>(1,Function_Arguments[imodule][1]);    
      ROL::Ptr<ROL::Vector<real_t> > ll = ROL::makePtr<ROL::StdVector<real_t>>(ll_ptr);
      ROL::Ptr<ROL::Vector<real_t> > lu = ROL::makePtr<ROL::StdVector<real_t>>(lu_ptr);
      ROL::Ptr<ROL::BoundConstraint<real_t>> constraint_bnd = ROL::makePtr<ROL::Bounds<real_t>>(ll,lu);
      if(TO_Module_List[imodule]=="Mass_Constraint"){
        *fos << " MASS CONSTRAINT EXPECTS FEA MODULE INDEX " <<TO_Module_My_FEA_Module[imodule] << std::endl;
        ineq_constraint = ROL::makePtr<MassConstraint_TopOpt>(fea_modules[TO_Module_My_FEA_Module[imodule]], nodal_density_flag);
      }
      else if(TO_Module_List[imodule]=="Moment_of_Inertia_Constraint"){
        *fos << " MOMENT OF INERTIA CONSTRAINT EXPECTS FEA MODULE INDEX " <<TO_Module_My_FEA_Module[imodule] << std::endl;
        ineq_constraint = ROL::makePtr<MassConstraint_TopOpt>(fea_modules[TO_Module_My_FEA_Module[imodule]], nodal_density_flag, Function_Arguments[imodule][2]);
      }
      else if(TO_Module_List[imodule]=="Strain_Energy_Constraint"){
        *fos << " STRAIN ENERGY CONSTRAINT EXPECTS FEA MODULE INDEX " <<TO_Module_My_FEA_Module[imodule] << std::endl;
        ineq_constraint = ROL::makePtr<StrainEnergyConstraint_TopOpt>(fea_modules[TO_Module_My_FEA_Module[imodule]], nodal_density_flag, Function_Arguments[imodule][2]);
      }
      else if(TO_Module_List[imodule]=="Heat_Capacity_Potential_Constraint"){
        *fos << " HEAT CAPACITY POTENTIAL CONSTRAINT EXPECTS FEA MODULE INDEX " <<TO_Module_My_FEA_Module[imodule] << std::endl;
        ineq_constraint = ROL::makePtr<HeatCapacityPotentialConstraint_TopOpt>(fea_modules[TO_Module_My_FEA_Module[imodule]], nodal_density_flag, Function_Arguments[imodule][2]);
      }
      else{
        *fos << "PROGRAM IS ENDING DUE TO ERROR; UNDEFINED INEQUALITY CONSTRAINT FUNCTION REQUESTED WITH NAME \"" <<TO_Module_List[imodule] <<"\"" << std::endl;
        exit_solver(0);
      }
      *fos << " ADDING CONSTRAINT " << constraint_name << std::endl;
      problem->addConstraint(constraint_name, ineq_constraint, constraint_mul, constraint_bnd);
    }
    number_union.clear();
  }

  //set bounds on design variables
  if(nodal_density_flag){
    dual_vec_array dual_node_densities_upper_bound = dual_vec_array("dual_node_densities_upper_bound", nlocal_nodes, 1);
    dual_vec_array dual_node_densities_lower_bound = dual_vec_array("dual_node_densities_lower_bound", nlocal_nodes, 1);
    host_vec_array node_densities_upper_bound = dual_node_densities_upper_bound.view_host();
    host_vec_array node_densities_lower_bound = dual_node_densities_lower_bound.view_host();
    //notify that the host view is going to be modified in the file readin
    dual_node_densities_upper_bound.modify_host();
    dual_node_densities_lower_bound.modify_host();

    //initialize densities to 1 for now; in the future there might be an option to read in an initial condition for each node
    for(int inode = 0; inode < nlocal_nodes; inode++){
      node_densities_upper_bound(inode,0) = 1;
      node_densities_lower_bound(inode,0) = DENSITY_EPSILON;
    }

    //set lower bounds for nodes on surfaces with boundary and loading conditions
    for(int imodule = 0; imodule < nfea_modules; imodule++){
      num_boundary_sets = fea_modules[imodule]->num_boundary_conditions;
      for(int iboundary = 0; iboundary < num_boundary_sets; iboundary++){

        num_bdy_patches_in_set = fea_modules[imodule]->NBoundary_Condition_Patches(iboundary);

        //loop over boundary patches for this boundary set
        for (int bdy_patch_gid = 0; bdy_patch_gid < num_bdy_patches_in_set; bdy_patch_gid++){
                
          // get the global id for this boundary patch
          patch_id = fea_modules[imodule]->Boundary_Condition_Patches(iboundary, bdy_patch_gid);
          Surface_Nodes = Boundary_Patches(patch_id).node_set;
          local_surface_id = Boundary_Patches(patch_id).local_patch_id;
          //debug print of local surface ids
          //std::cout << " LOCAL SURFACE IDS " << std::endl;
          //std::cout << local_surface_id << std::endl;
          //acquire set of nodes for this face
          for(int node_loop=0; node_loop < Surface_Nodes.size(); node_loop++){
            current_node_index = Surface_Nodes(node_loop);
            if(map->isNodeGlobalElement(current_node_index)){
              local_node_index = map->getLocalElement(current_node_index);
              node_densities_lower_bound(local_node_index,0) = 1;
            }
          }// node loop for
        }//boundary patch for
      }//boundary set for
    }//module for
  
    //sync device view
    dual_node_densities_upper_bound.sync_device();
    dual_node_densities_lower_bound.sync_device();
    
    //allocate global vector information
    upper_bound_node_densities_distributed = Teuchos::rcp(new MV(map, dual_node_densities_upper_bound));
    lower_bound_node_densities_distributed = Teuchos::rcp(new MV(map, dual_node_densities_lower_bound));
  
  }
  else{
    //initialize memory for volume storage
    vec_array Element_Densities_Upper_Bound("Element Densities_Upper_Bound", rnum_elem, 1);
    vec_array Element_Densities_Lower_Bound("Element Densities_Lower_Bound", rnum_elem, 1);
    for(int ielem = 0; ielem < rnum_elem; ielem++){
      Element_Densities_Upper_Bound(ielem,0) = 1;
      Element_Densities_Lower_Bound(ielem,0) = DENSITY_EPSILON;
    }

    //create global vector
    Global_Element_Densities_Upper_Bound = Teuchos::rcp(new MV(element_map, Element_Densities_Upper_Bound));
    Global_Element_Densities_Lower_Bound = Teuchos::rcp(new MV(element_map, Element_Densities_Lower_Bound));
  }
    
  // Bound constraint defining the possible range of design density variables
  ROL::Ptr<ROL::Vector<real_t> > lower_bounds;
  ROL::Ptr<ROL::Vector<real_t> > upper_bounds;
  if(nodal_density_flag){
    lower_bounds = ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO>>(lower_bound_node_densities_distributed);
    upper_bounds = ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO>>(upper_bound_node_densities_distributed);
  }
  else{
    lower_bounds = ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO>>(Global_Element_Densities_Lower_Bound);
    upper_bounds = ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO>>(Global_Element_Densities_Upper_Bound);
  }
  ROL::Ptr<ROL::BoundConstraint<real_t> > bnd = ROL::makePtr<ROL::Bounds<real_t>>(lower_bounds, upper_bounds);
  problem->addBoundConstraint(bnd);

  //compute initial constraint satisfaction
  //ROL::Ptr<ROL_MV> ROL_Element_Masses = ROL::makePtr<ROL_MV>(fea_elasticity->Global_Element_Masses);
  ROL::Elementwise::ReductionSum<real_t> sumreduc;
  if(nodal_density_flag)
    design_densities = design_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  else
    design_densities = Global_Element_Densities->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  //fea_elasticity->compute_element_masses(design_densities,true);
  
  //real_t initial_mass = ROL_Element_Masses->reduce(sumreduc);

  problem->setProjectionAlgorithm(*parlist);
  //finalize problem
  problem->finalize(false,true,*fos);
  //problem->check(true,std::cout);

  //debug checks
  ROL::Ptr<ROL::TpetraMultiVector<real_t,LO,GO,node_type>> rol_x =
   ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO,node_type>>(design_node_densities_distributed);
  //construct direction vector for check
  Teuchos::RCP<MV> directions_distributed = Teuchos::rcp(new MV(map, 1));
  //directions_distributed->putScalar(1);
  directions_distributed->randomize(-1,1);
  //real_t normd = directions_distributed->norm2();
  //directions_distributed->scale(normd);
  //set all but first component to 0 for debug
  host_vec_array directions = directions_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  //for(int init = 1; init < nlocal_nodes; init++)
  //directions(4,0) = -0.3;
  ROL::Ptr<ROL::TpetraMultiVector<real_t,LO,GO,node_type>> rol_d =
  ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO,node_type>>(directions_distributed);
  //obj->checkGradient(*rol_x, *rol_d);
  //obj->checkHessVec(*rol_x, *rol_d);
  //directions_distributed->putScalar(-0.000001);
  //obj->checkGradient(*rol_x, *rol_d);
  //directions_distributed->putScalar(-0.0000001);
  //obj->checkGradient(*rol_x, *rol_d);
  

  // Instantiate Solver.
  ROL::Solver<real_t> solver(problem,*parlist);
    
  // Solve optimization problem.
  //std::ostream outStream;
  solver.solve(*fos);

  //print final constraint satisfaction
  //fea_elasticity->compute_element_masses(design_densities,false);
  //real_t final_mass = ROL_Element_Masses->reduce(sumreduc);
  //if(myrank==0)
    //std::cout << "Final Mass Constraint is " << final_mass/initial_mass << std::endl;
}

/* ----------------------------------------------------------------------
   Find boundary surface segments that belong to this MPI rank
------------------------------------------------------------------------- */

void Implicit_Solver::Get_Boundary_Patches(){
  size_t npatches_repeat, npatches, element_npatches, num_nodes_in_patch, node_gid;
  int local_node_id;
  int num_dim = simparam->num_dim;
  CArrayKokkos<GO, array_layout, device_type, memory_traits> Surface_Nodes;
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  //Surface_Nodes = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(4, "Surface_Nodes");
  
  std::set<Node_Combination> my_patches;
  //inititializes type for the pair variable (finding the iterator type is annoying)
  std::pair<std::set<Node_Combination>::iterator, bool> current_combination;
  std::set<Node_Combination>::iterator it;
  
  //compute the number of patches in this MPI rank with repeats for adjacent cells
  npatches_repeat = 0;

  if(num_dim==2)
  for(int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
    element_npatches = elem2D->nsurfaces;
    npatches_repeat += element_npatches;
  }
  
  else if(num_dim==3)
  for(int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_3Delem_type(Element_Types(ielem), elem);
    element_npatches = elem->nsurfaces;
    npatches_repeat += element_npatches;
  }
  std::cout << "Starting boundary patch allocation of size " << npatches_repeat << std::endl <<std::flush;
  //information for all patches on this rank
  CArrayKokkos<Node_Combination,array_layout, device_type, memory_traits> Patch_Nodes(npatches_repeat, "Patch_Nodes");
  CArrayKokkos<size_t,array_layout, device_type, memory_traits> Patch_Boundary_Flags(npatches_repeat, "Patch_Boundary_Flags");
  if(myrank == 0)
    std::cout << "Done with boundary patch allocation" << std::endl <<std::flush;
  //initialize boundary patch flags
  for(int init = 0; init < npatches_repeat; init++)
    Patch_Boundary_Flags(init) = 1;
  
  if(myrank == 0)
    std::cout << "Done with boundary patch flags init" << std::endl <<std::flush;
  //use set of nodal combinations to find boundary set
  //boundary patches will not try to add nodal combinations twice
  //loop through elements in this rank to find boundary patches
  npatches_repeat = 0;
  if(num_dim==2)
  for(int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
    element_npatches = elem2D->nsurfaces;
    //loop through local surfaces
    for(int isurface = 0; isurface < element_npatches; isurface++){
      num_nodes_in_patch = elem2D->surface_to_dof_lid.stride(isurface);
      Surface_Nodes = CArrayKokkos<GO, array_layout, device_type, memory_traits>(num_nodes_in_patch, "Surface_Nodes");
      for(int inode = 0; inode < num_nodes_in_patch; inode++){
        local_node_id = elem2D->surface_to_dof_lid(isurface,inode);
        Surface_Nodes(inode) = nodes_in_elem(ielem, local_node_id);
      }
      Node_Combination temp(Surface_Nodes);
      //construct Node Combination object for this surface
      Patch_Nodes(npatches_repeat) = temp;
      Patch_Nodes(npatches_repeat).patch_id = npatches_repeat;
      Patch_Nodes(npatches_repeat).element_id = ielem;
      Patch_Nodes(npatches_repeat).local_patch_id = isurface;
      //test if this patch has already been added; if yes set boundary flags to 0
      current_combination = my_patches.insert(Patch_Nodes(npatches_repeat));
      //if the set determines this is a duplicate, access the original element's patch id and set flag to 0
      if(current_combination.second==false){
      //set original element flag to 0
      Patch_Boundary_Flags((*current_combination.first).patch_id) = 0;
      //set this current flag to 0 for the duplicate as well
      Patch_Boundary_Flags(npatches_repeat) = 0;
      npatches_repeat++;
      }

    }
  }

  if(num_dim==3)
  for(int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_3Delem_type(Element_Types(ielem), elem);
    element_npatches = elem->nsurfaces;
    //loop through local surfaces
    for(int isurface = 0; isurface < element_npatches; isurface++){
      num_nodes_in_patch = elem->surface_to_dof_lid.stride(isurface);
      //debug print
      //std::cout << "NUMBER OF PATCH NODES FOR ELEMENT " << ielem+1 << " ON LOCAL SURFACE " << isurface+1 << " IS " << elem->surface_to_dof_lid.start_index_[isurface+1] << std::endl;
      Surface_Nodes = CArrayKokkos<GO, array_layout, device_type, memory_traits>(num_nodes_in_patch, "Surface_Nodes");
      for(int inode = 0; inode < num_nodes_in_patch; inode++){
        local_node_id = elem->surface_to_dof_lid(isurface,inode);
        Surface_Nodes(inode) = nodes_in_elem(ielem, local_node_id);
      }
      Node_Combination temp(Surface_Nodes);
      //construct Node Combination object for this surface
      Patch_Nodes(npatches_repeat) = temp;
      Patch_Nodes(npatches_repeat).patch_id = npatches_repeat;
      Patch_Nodes(npatches_repeat).element_id = ielem;
      Patch_Nodes(npatches_repeat).local_patch_id = isurface;
      //test if this patch has already been added; if yes set boundary flags to 0
      current_combination = my_patches.insert(Patch_Nodes(npatches_repeat));
      //if the set determines this is a duplicate access the original element's patch id and set flag to 0
      if(current_combination.second==false){
        //set original element flag to 0
        Patch_Boundary_Flags((*current_combination.first).patch_id) = 0;
        //set this current flag to 0 for the duplicate as well
        Patch_Boundary_Flags(npatches_repeat) = 0;
      }
      npatches_repeat++;
    }
  }

  //debug print of all patches
  /*
  std::cout << " ALL PATCHES " << npatches_repeat <<std::endl;
  for(int iprint = 0; iprint < npatches_repeat; iprint++){
    std::cout << "Patch " << iprint + 1 << " ";
    for(int j = 0; j < Patch_Nodes(iprint).node_set.size(); j++)
      std::cout << Patch_Nodes(iprint).node_set(j) << " ";
    std::cout << std::endl;
  }
  */
  if(myrank == 0)
  std::cout << "Done with boundary patch loop" << std::endl <<std::flush;
  //loop through patch boundary flags to isolate boundary patches
  nboundary_patches = 0;
  for(int iflags = 0 ; iflags < npatches_repeat; iflags++){
    if(Patch_Boundary_Flags(iflags)) nboundary_patches++;
  }
  //upper bound that is not much larger
  Boundary_Patches = CArrayKokkos<Node_Combination, array_layout, device_type, memory_traits>(nboundary_patches, "Boundary_Patches");
  nboundary_patches = 0;
  bool my_rank_flag;
  size_t remote_count;
  for(int ipatch = 0 ; ipatch < npatches_repeat; ipatch++){
    if(Patch_Boundary_Flags(ipatch)){
      /*check if Nodes in the combination for this patch belong to this MPI rank.
        If all are local then this is a boundary patch belonging to this rank.
        If all nodes are remote then another rank must decide if that patch is a boundary.
        If only a subset of the nodes are local it must be a boundary patch; this
        case assigns the patch to the lowest mpi rank index the nodes in this patch belong to */
      num_nodes_in_patch = Patch_Nodes(ipatch).node_set.size();
      my_rank_flag = true;
      remote_count = 0;

      //assign as a local boundary patch if any of the nodes on the patch are local
      //only the local nodes on the patch will contribute to the equation assembly on this rank
      for(int inode = 0; inode < num_nodes_in_patch; inode++){
        node_gid = Patch_Nodes(ipatch).node_set(inode);
        if(!map->isNodeGlobalElement(node_gid)){
          //if(ghost_node_ranks(global2local_map.get(node_gid))<myrank)
          //my_rank_flag = false;
          remote_count++;
          //test
        } 
      }
      //all nodes were remote
      //if(remote_count == num_nodes_in_patch) my_rank_flag = false;

      //if all nodes were not local
      if(my_rank_flag){
        Boundary_Patches(nboundary_patches++) = Patch_Nodes(ipatch);
        boundary_patch_to_index[Patch_Nodes(ipatch)] = nboundary_patches-1;
      }
    }
  }

  //debug print of boundary patches
  /*std::cout << " BOUNDARY PATCHES ON TASK " << myrank << " = " << nboundary_patches <<std::endl;
  for(int iprint = 0; iprint < nboundary_patches; iprint++){
    std::cout << "Patch " << iprint + 1 << " ";
    for(int j = 0; j < Boundary_Patches(iprint).node_set.size(); j++)
      std::cout << Boundary_Patches(iprint).node_set(j) << " ";
    std::cout << std::endl;
  }
  std::fflush(stdout);
  */
}

/* ----------------------------------------------------------------------------
   Initialize sets of element boundary surfaces and arrays for input conditions
------------------------------------------------------------------------------- */

void Implicit_Solver::init_boundaries(){
  int num_boundary_sets = simparam->NB;
  int num_dim = simparam->num_dim;

  // build boundary mesh patches
  if(myrank == 0)
    std::cout << "Starting boundary patch setup" << std::endl <<std::flush;
  Get_Boundary_Patches();
  //std::cout << "Done with boundary patch setup" << std::endl <<std::flush;
  std::cout << "number of boundary patches on task " << myrank << " = " << nboundary_patches << std::endl;
  
  //disable for now
  if(0)
  init_topology_conditions(num_boundary_sets);
}

/* ----------------------------------------------------------------------
   Assign sets of element boundary surfaces corresponding to user BCs
------------------------------------------------------------------------- */

void Implicit_Solver::generate_tcs(){
  int num_dim = simparam->num_dim;
  int bdy_set_id;
  int tc_tag;
  int num_topology_conditions = 0;
  real_t value;
  real_t fix_limits[4];
  
  if(num_topology_conditions>0)
    init_topology_conditions(num_topology_conditions);

} // end generate_tcs

/* ----------------------------------------------------------------------
   initialize storage for DOF corresponding to user TCs
------------------------------------------------------------------------- */

void Implicit_Solver::init_topology_conditions (int num_sets){
 
  //surface conditions
  num_boundary_conditions = num_sets;
  if(num_sets == 0){
    std::cout << " Warning: number of boundary conditions = 0";
    return;
  }
  NTopology_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_sets, "NBoundary_Condition_Patches");
  Topology_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_sets, nboundary_patches, "Boundary_Condition_Patches");

  //initialize data
  for(int iset = 0; iset < num_sets; iset++) NTopology_Condition_Patches(iset) = 0;
}

/* ----------------------------------------------------------------------
   find which boundary patches correspond to the given BC.
   bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
   val = plane value, cylinder radius, shell radius
------------------------------------------------------------------------- */

void Implicit_Solver::tag_boundaries(int bc_tag, real_t val, int bdy_set, real_t *patch_limits){
  
  int num_boundary_sets = simparam->NB;
  int is_on_set;
  /*
  if (bdy_set == num_bdy_sets_){
    std::cout << " ERROR: number of boundary sets must be increased by "
      << bdy_set-num_bdy_sets_+1 << std::endl;
    exit(0);
  }
  */

  //test patch limits for feasibility
  if(patch_limits != NULL){
    //test for upper bounds being greater than lower bounds
    if(patch_limits[1] <= patch_limits[0]) std::cout << " Warning: patch limits for boundary condition are infeasible";
    if(patch_limits[2] <= patch_limits[3]) std::cout << " Warning: patch limits for boundary condition are infeasible";
  }
    
  // save the boundary vertices to this set that are on the plane
  int counter = 0;
  for (int iboundary_patch = 0; iboundary_patch < nboundary_patches; iboundary_patch++) {

    // check to see if this patch is on the specified plane
    is_on_set = check_boundary(Boundary_Patches(iboundary_patch), bc_tag, val, patch_limits); // no=0, yes=1
        
    if (is_on_set == 1){
      Topology_Condition_Patches(bdy_set,counter) = iboundary_patch;
      counter ++;
    }
  } // end for bdy_patch
    
  // save the number of bdy patches in the set
  NTopology_Condition_Patches(bdy_set) = counter;
    
  *fos << " tagged boundary patches " << std::endl;
}

/* ----------------------------------------------------------------------
   routine for checking to see if a patch is on a boundary set
   bc_tag = 0 xplane, 1 yplane, 3 zplane, 4 cylinder, 5 is shell
   val = plane value, radius, radius
------------------------------------------------------------------------- */


int Implicit_Solver::check_boundary(Node_Combination &Patch_Nodes, int bc_tag, real_t val, real_t *patch_limits){
  
  int is_on_set = 1;
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);

  //Nodes on the Patch
  auto node_list = Patch_Nodes.node_set;
  int num_dim = simparam->num_dim;
  size_t nnodes = node_list.size();
  size_t node_rid;
  real_t node_coord[num_dim];
  int dim_other1, dim_other2;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> node_on_flags(nnodes, "node_on_flags");

  //initialize
  for(int inode = 0; inode < nnodes; inode++) node_on_flags(inode) = 0;

  if(bc_tag==0){
    dim_other1 = 1;
    dim_other2 = 2;
  }
  else if(bc_tag==1){
    dim_other1 = 0;
    dim_other2 = 2;
  }
  else if(bc_tag==2){
    dim_other1 = 0;
    dim_other2 = 1;
  }
  
  
  //test for planes
  if(bc_tag < 3)
  for(int inode = 0; inode < nnodes; inode++){

    node_rid = all_node_map->getLocalElement(node_list(inode));
    for(int init=0; init < num_dim; init++){
      node_coord[init] = all_node_coords(node_rid,init);
    }
    if ( fabs(node_coord[bc_tag] - val) <= BC_EPSILON){ node_on_flags(inode) = 1;

      //test if within patch segment if user specified
      if(patch_limits!=NULL){
        if (node_coord[dim_other1] - patch_limits[0] <= -BC_EPSILON) node_on_flags(inode) = 0;
        if (node_coord[dim_other1] - patch_limits[1] >= BC_EPSILON) node_on_flags(inode) = 0;
        if (node_coord[dim_other2] - patch_limits[2] <= -BC_EPSILON) node_on_flags(inode) = 0;
        if (node_coord[dim_other2] - patch_limits[3] >= BC_EPSILON) node_on_flags(inode) = 0;
      }
    }
    //debug print of node id and node coord
    //std::cout << "node coords on task " << myrank << " for node " << node_rid << std::endl;
    //std::cout << "coord " <<node_coord << " flag " << node_on_flags(inode) << " bc_tag " << bc_tag << std::endl;
  }
    
    /*
    // cylinderical shell where radius = sqrt(x^2 + y^2)
    else if (this_bc_tag == 3){
        
        real_t R = sqrt(these_patch_coords[0]*these_patch_coords[0] +
                        these_patch_coords[1]*these_patch_coords[1]);
        
        if ( fabs(R - val) <= 1.0e-8 ) is_on_bdy = 1;

        
    }// end if on type
    
    // spherical shell where radius = sqrt(x^2 + y^2 + z^2)
    else if (this_bc_tag == 4){
        
        real_t R = sqrt(these_patch_coords[0]*these_patch_coords[0] +
                        these_patch_coords[1]*these_patch_coords[1] +
                        these_patch_coords[2]*these_patch_coords[2]);
        
        if ( fabs(R - val) <= 1.0e-8 ) is_on_bdy = 1;
        
    } // end if on type
    */
    //check if all nodes lie on the boundary set
  for(int inode = 0; inode < nnodes; inode++)
    if(!node_on_flags(inode)) is_on_set = 0;
  
  //debug print of return flag
  //std::cout << "patch flag on task " << myrank << " is " << is_on_set << std::endl;
  return is_on_set;
    
} // end method to check bdy

// Output writers

/* ----------------------------------------------------------------------
   Collect Nodal Information on Rank 0
------------------------------------------------------------------------- */

void Implicit_Solver::collect_information(){
  GO nreduce_nodes = 0;
  GO nreduce_elem = 0;
  int num_dim = simparam->num_dim;

  //collect nodal coordinate information
  if(myrank==0) nreduce_nodes = num_nodes;
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > global_reduce_map =
    Teuchos::rcp(new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),nreduce_nodes,0,comm));
  
  //importer from local node distribution to collected distribution
  Tpetra::Import<LO, GO> node_collection_importer(map, global_reduce_map);

  Teuchos::RCP<MV> collected_node_coords_distributed = Teuchos::rcp(new MV(global_reduce_map, num_dim));

  //comms to collect
  collected_node_coords_distributed->doImport(*node_coords_distributed, node_collection_importer, Tpetra::INSERT);

  //comms to collect FEA module related vector data
  for (int imodule = 0; imodule < nfea_modules; imodule++){
    fea_modules[imodule]->collect_output(global_reduce_map);
    //collected_node_displacements_distributed->doImport(*(fea_elasticity->node_displacements_distributed), dof_collection_importer, Tpetra::INSERT);
  }
  

  //collected nodal density information
  Teuchos::RCP<MV> collected_node_densities_distributed = Teuchos::rcp(new MV(global_reduce_map, 1));

  //comms to collect
  collected_node_densities_distributed->doImport(*design_node_densities_distributed, node_collection_importer, Tpetra::INSERT);

  //collect element connectivity data
  if(myrank==0) nreduce_elem = num_elem;
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > global_reduce_element_map =
    Teuchos::rcp(new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),nreduce_elem,0,comm));

  //importer from all element map to collected distribution
  Tpetra::Import<LO, GO> element_collection_importer(all_element_map, global_reduce_element_map);
  
  Teuchos::RCP<MCONN> collected_nodes_in_elem_distributed = Teuchos::rcp(new MCONN(global_reduce_element_map, max_nodes_per_element));

  //comms to collect
  collected_nodes_in_elem_distributed->doImport(*nodes_in_elem_distributed, element_collection_importer, Tpetra::INSERT);

  //collect element type data

  //set host views of the collected data to print out from
  if(myrank==0){
    collected_node_coords = collected_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    collected_node_densities = collected_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    collected_nodes_in_elem = collected_nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  }
}

/* ----------------------------------------------------------------------
   Output Model Information in tecplot format
------------------------------------------------------------------------- */

void Implicit_Solver::tecplot_writer(){
  
  int num_dim = simparam->num_dim;
	std::string current_file_name;
	std::string base_file_name= "TecplotTO";
  std::string base_file_name_undeformed= "TecplotTO_undeformed";
	std::stringstream ss;
	std::string file_extension= ".dat";
  std::string file_count;
	std::stringstream count_temp;
  int time_step = 0;
  int temp_convert;
  int noutput, nvector;
  bool displace_geometry;
  if(displacement_module!=-1)
    displace_geometry = fea_modules[displacement_module]->displaced_mesh_flag;
  const_host_vec_array current_collected_output;
  for (int imodule = 0; imodule < nfea_modules; imodule++){
    fea_modules[imodule]->compute_output();
  }
  
  collect_information();
  // Convert ijk index system to the finite element numbering convention
  // for vertices in cell
  CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> convert_ijk_to_ensight(max_nodes_per_element);
  CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> tmp_ijk_indx(max_nodes_per_element);
  convert_ijk_to_ensight(0) = 0;
  convert_ijk_to_ensight(1) = 1;
  convert_ijk_to_ensight(2) = 3;
  convert_ijk_to_ensight(3) = 2;
  convert_ijk_to_ensight(4) = 4;
  convert_ijk_to_ensight(5) = 5;
  convert_ijk_to_ensight(6) = 7;
  convert_ijk_to_ensight(7) = 6;

  //compared to primitive unit cell, assumes orthogonal primitive unit cell
    if(myrank==0){
      //initial undeformed geometry
      count_temp.str("");
      count_temp << file_index;
      //file_index++;
	    file_count = count_temp.str();
      if(displace_geometry&&displacement_module>=0)
        current_file_name = base_file_name_undeformed + file_count + file_extension;
      else
        current_file_name = base_file_name + file_count + file_extension;
      std::ofstream myfile (current_file_name.c_str()); //output filestream object for file output
	    //read in position data
	    myfile << std::fixed << std::setprecision(8);
		
		  //output header of the tecplot file

		  myfile << "TITLE=\"results for TO code\"" "\n";
      //myfile << "VARIABLES = \"x\", \"y\", \"z\", \"density\", \"sigmaxx\", \"sigmayy\", \"sigmazz\", \"sigmaxy\", \"sigmaxz\", \"sigmayz\"" "\n";
      //else
		  myfile << "VARIABLES = \"x\", \"y\", \"z\", \"density\"";
      for (int imodule = 0; imodule < nfea_modules; imodule++){
        for(int ioutput = 0; ioutput < fea_modules[imodule]->noutput; ioutput++){
          nvector = fea_modules[imodule]->output_vector_sizes[ioutput];
          for(int ivector = 0; ivector < nvector; ivector++){
            myfile << ", \"" << fea_modules[imodule]->output_dof_names[ioutput][ivector] << "\"";
          }
        }
      }
      myfile << "\n";

		  myfile << "ZONE T=\"load step " << time_step << "\", NODES= " << num_nodes
			  << ", ELEMENTS= " << num_elem << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" "\n";

		  for (int nodeline = 0; nodeline < num_nodes; nodeline++) {
			  myfile << std::setw(25) << collected_node_coords(nodeline,0) << " ";
			  myfile << std::setw(25) << collected_node_coords(nodeline,1) << " ";
        if(num_dim==3)
			  myfile << std::setw(25) << collected_node_coords(nodeline,2) << " ";
        myfile << std::setw(25) << collected_node_densities(nodeline,0) << " ";
        for (int imodule = 0; imodule < nfea_modules; imodule++){
          noutput = fea_modules[imodule]->noutput;
          for(int ioutput = 0; ioutput < noutput; ioutput++){
            current_collected_output = fea_modules[imodule]->collected_module_output[ioutput];
            if(fea_modules[imodule]->vector_style[ioutput] == FEA_Module::DOF){
              nvector = fea_modules[imodule]->output_vector_sizes[ioutput];
              for(int ivector = 0; ivector < nvector; ivector++){
                myfile << std::setw(25) << current_collected_output(nodeline*nvector + ivector,0) << " ";
              }
            }
            if(fea_modules[imodule]->vector_style[ioutput] == FEA_Module::NODAL){
              nvector = fea_modules[imodule]->output_vector_sizes[ioutput];
              for(int ivector = 0; ivector < nvector; ivector++){
                myfile << std::setw(25) << current_collected_output(nodeline,ivector) << " ";
              }
            }
          }
        }
        myfile << std::endl;
		  }
		  for (int elementline = 0; elementline < num_elem; elementline++) {
        //convert node ordering
			  for (int ii = 0; ii < max_nodes_per_element; ii++) {
          temp_convert = convert_ijk_to_ensight(ii);
				  myfile << std::setw(10) << collected_nodes_in_elem(elementline, temp_convert) + 1 << " ";
			  }
			  myfile << " \n";
		  }
      myfile.close();
    }
    if(myrank==0&&(displacement_module>=0&&displace_geometry)){
      //deformed geometry
      count_temp.str("");
      count_temp << file_index;
      file_index++;
	    file_count = count_temp.str();
    
      current_file_name = base_file_name + file_count + file_extension;
      std::ofstream myfile (current_file_name.c_str()); //output filestream object for file output
	    //read in position data
	    myfile << std::fixed << std::setprecision(8);
		
		  //output header of the tecplot file

		  myfile << "TITLE=\"results for TO code\" \n";
		  myfile << "VARIABLES = \"x\", \"y\", \"z\", \"density\"";
      for (int imodule = 0; imodule < nfea_modules; imodule++){
        for(int ioutput = 0; ioutput < fea_modules[imodule]->noutput; ioutput++){
          nvector = fea_modules[imodule]->output_vector_sizes[ioutput];
          for(int ivector = 0; ivector < nvector; ivector++){
            myfile << ", \"" << fea_modules[imodule]->output_dof_names[ioutput][ivector] << "\"";
          }
        }
      }
      myfile << "\n";

		  myfile << "ZONE T=\"load step " << time_step + 1 << "\", NODES= " << num_nodes
			<< ", ELEMENTS= " << num_elem << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" "\n";

		  for (int nodeline = 0; nodeline < num_nodes; nodeline++) {
			  myfile << std::setw(25) << collected_node_coords(nodeline,0) + fea_modules[displacement_module]->collected_displacement_output(nodeline*num_dim,0) << " ";
			  myfile << std::setw(25) << collected_node_coords(nodeline,1) + fea_modules[displacement_module]->collected_displacement_output(nodeline*num_dim + 1,0) << " ";
        if(num_dim==3)
			  myfile << std::setw(25) << collected_node_coords(nodeline,2) + fea_modules[displacement_module]->collected_displacement_output(nodeline*num_dim + 2,0) << " ";
        myfile << std::setw(25) << collected_node_densities(nodeline,0) << " ";
        for (int imodule = 0; imodule < nfea_modules; imodule++){
          noutput = fea_modules[imodule]->noutput;
          for(int ioutput = 0; ioutput < noutput; ioutput++){
            current_collected_output = fea_modules[imodule]->collected_module_output[ioutput];
            if(fea_modules[imodule]->vector_style[ioutput] == FEA_Module::DOF){
              nvector = fea_modules[imodule]->output_vector_sizes[ioutput];
              for(int ivector = 0; ivector < nvector; ivector++){
                myfile << std::setw(25) << current_collected_output(nodeline*nvector + ivector,0) << " ";
              }
            }
            if(fea_modules[imodule]->vector_style[ioutput] == FEA_Module::NODAL){
              nvector = fea_modules[imodule]->output_vector_sizes[ioutput];
              for(int ivector = 0; ivector < nvector; ivector++){
                myfile << std::setw(25) << current_collected_output(nodeline,ivector) << " ";
              }
            }
          }
        }
        myfile << std::endl;
		  }
		  for (int elementline = 0; elementline < num_elem; elementline++) {
        //convert node ordering
			  for (int ii = 0; ii < max_nodes_per_element; ii++) {
          temp_convert = convert_ijk_to_ensight(ii);
				  myfile << std::setw(10) << collected_nodes_in_elem(elementline, temp_convert) + 1 << " ";
			  }
			  myfile << " \n";
		  }
      myfile.close();
    }
}
/* ----------------------------------------------------------------------
   Output Model Information in vtk format
------------------------------------------------------------------------- */
/*
void Implicit_Solver::vtk_writer(){
    //local variable for host view in the dual view
    host_vec_array node_coords = dual_node_coords.view_host();
    const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    int graphics_id = simparam->graphics_id;
    int num_dim = simparam->num_dim;

    const int num_scalar_vars = 2;
    const int num_vec_vars = 1;
    const int num_cell_vars = 1;

    const char name[10] = {"Testing"};
    const char scalar_var_names[num_scalar_vars][15] = {
        "fake1", "elem"
    };
    const char vec_var_names[num_vec_vars][15] = {
        "velocity"
    };

    const char cell_var_names[num_vec_vars][15] = {
        "cell_rid"
    };
    //  ---------------------------------------------------------------------------
    //  Setup of file and directoring for exporting
    //  ---------------------------------------------------------------------------

    
    FILE *out[20];   // the output files that are written to
    char filename[128];
    
    struct stat st;
    
    if(stat("vtk",&st) != 0)
        system("mkdir vtk");
    
    if(stat("vtk/data",&st) != 0)
        system("mkdir vtk/data");


    //  ---------------------------------------------------------------------------
    //  Write the Geometry file
    //  ---------------------------------------------------------------------------
    

    
    sprintf(filename, "vtk/data/%s_%05d_%i.vtu", name, graphics_id, 0);
    // filename has the full string
    
    out[0] = fopen(filename, "w");
    
    int num_nodes = mesh->num_nodes();
    int num_cells = mesh->num_cells();


    fprintf(out[0],"<?xml version=\"1.0\"?> \n");
    fprintf(out[0],"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(out[0],"<UnstructuredGrid>\n");
    fprintf(out[0],"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", num_nodes, num_cells);

    

    //  ---------------------------------------------------------------------------
    //  Write point data
    //  ---------------------------------------------------------------------------


    fprintf(out[0],"<PointData> \n");


    fprintf(out[0],"</PointData> \n");
    fprintf(out[0],"\n");



    //  ---------------------------------------------------------------------------
    //  Write cell data
    //  ---------------------------------------------------------------------------

    fprintf(out[0],"<CellData> \n");

    for(int cell_var = 0; cell_var < num_cell_vars; cell_var++){
        
        fprintf(out[0],"<DataArray type=\"Float64\" Name=\"%s\" Format=\"ascii\">\n", cell_var_names[cell_var]);
        
        for (int cell_rid = 0; cell_rid < num_cells; cell_rid++){
            
            fprintf(out[0],"%f\n",(float) cell_rid);

        } // end for k over cells
    }
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"</CellData> \n");
    fprintf(out[0],"\n");


    //  ---------------------------------------------------------------------------
    //  Write node positions
    //  ---------------------------------------------------------------------------

    real_t min_coord = 0;
    real_t max_coord = 2.0;
    fprintf(out[0],"<Points> \n");

        
    fprintf(out[0],"<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"%i\" format=\"ascii\">\n", num_dim);
    
    for (int node_gid = 0; node_gid < num_nodes; node_gid++){
        
        fprintf(out[0],"%f   %f   %f   \n",node_coords(node_gid, 0),
                                           node_coords(node_gid, 1),
                                           node_coords(node_gid, 2));

    } 

    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"</Points> \n");
    fprintf(out[0],"\n");


    //  ---------------------------------------------------------------------------
    //  Write cell type definitions
    //  ---------------------------------------------------------------------------

    fprintf(out[0],"<Cells> \n");
    fprintf(out[0],"\n");


    // Cell connectivity

    fprintf(out[0],"<DataArray type=\"Int64\" Name=\"connectivity\">\n");

    // write nodes in a cell
    for (int cell_rid = 0; cell_rid < num_cells; cell_rid++){
        
        for(int node = 0; node < 8; node++){
            fprintf(out[0],"%i  ", nodes_in_elem(cell_rid, node));
        }

        fprintf(out[0],"\n");

    } // end for k over cells
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    fprintf(out[0],"<DataArray type=\"Int64\" Name=\"offsets\">\n");
    for (int cell_rid = 0; cell_rid < num_cells; cell_rid++){
        fprintf(out[0],"%i  \n", 8*(cell_rid+1));
    } // end for k over cells
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    fprintf(out[0],"<DataArray type=\"UInt64\" Name=\"types\">\n");
    for (int cell_rid = 0; cell_rid < num_cells; cell_rid++){
        fprintf(out[0],"%i  \n", 42);
    } // end for k over cells
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    fprintf(out[0],"<DataArray type=\"Int64\" Name=\"faces\">\n");
    
    for (int cell_rid = 0; cell_rid < num_cells; cell_rid++){
        fprintf(out[0],"%i  \n", 6);

        for(int patch_lid = 0; patch_lid < 6; patch_lid++){

            fprintf(out[0],"4  ");
            for(int node_lid = 0; node_lid < 4; node_lid++){
                fprintf(out[0],"%i  ", mesh->node_in_patch_in_cell(cell_rid, patch_lid, node_lid));
            }

            fprintf(out[0],"\n");
        }


    } // end for k over cells
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    int faceoffsets = 31;
    fprintf(out[0],"<DataArray type=\"Int64\" Name=\"faceoffsets\">\n");
    for (int cell_rid = 0; cell_rid < num_cells; cell_rid++){
        fprintf(out[0],"%i  \n", faceoffsets*(cell_rid+1));
    } // end for k over cells
    
    fprintf(out[0],"</DataArray> \n");
    fprintf(out[0],"\n");


    fprintf(out[0],"</Cells> \n");
    fprintf(out[0],"\n");


    fprintf(out[0],"\n");
    fprintf(out[0],"</Piece> \n");
    fprintf(out[0],"</UnstructuredGrid> \n");
    fprintf(out[0],"</VTKFile> \n");

    fclose(out[0]);
} // end vtk_writer
*/
/* ----------------------------------------------------------------------
   Output Model Information in ensight format
------------------------------------------------------------------------- */
/*
void Implicit_Solver::ensight_writer(){
    //local variable for host view in the dual view
    host_vec_array node_coords = dual_node_coords.view_host();
    const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    mat_pt_t *mat_pt = simparam->mat_pt;
    int &graphics_id = simparam->graphics_id;
    real_t *graphics_times = simparam->graphics_times;
    real_t &TIME = simparam->TIME;

    collect_information();

    auto convert_vert_list_ord_Ensight = CArray <int> (8);
    convert_vert_list_ord_Ensight(0) = 1;
    convert_vert_list_ord_Ensight(1) = 0;
    convert_vert_list_ord_Ensight(2) = 2;
    convert_vert_list_ord_Ensight(3) = 3;
    convert_vert_list_ord_Ensight(4) = 5;
    convert_vert_list_ord_Ensight(5) = 4;
    convert_vert_list_ord_Ensight(6) = 6;
    convert_vert_list_ord_Ensight(7) = 7;


    const int num_scalar_vars = 4;
    const int num_vec_vars = 1;

    const char name[10] = {"Testing"};
    const char scalar_var_names[num_scalar_vars][15] = {
       "cell_field1", "elem", "elem_id", "cell_field2"
    };
    const char vec_var_names[num_vec_vars][15] = {
        "position"
    };


    int num_nodes = mesh->num_nodes();
    int num_cells = mesh->num_cells();

    // save the cell state to an array for exporting to graphics files
    auto cell_fields = CArray <real_t> (num_cells, num_scalar_vars);

    int cell_cnt = 0;
    int c_in_e = mesh->num_cells_in_elem();
    int elem_val = 1;


    for (int cell_rid=0; cell_rid<num_cells; cell_rid++){
        cell_fields(cell_rid, 0) = mat_pt->field(cell_rid);    
    } // end for k over cells

    int num_elem = rnum_elem;

    int num_sub_1d;

    if(mesh->elem_order() == 0){
        num_sub_1d = 1;
    }
    
    else{
        num_sub_1d = mesh->elem_order()*2;
    }



    for (int elem_gid = 0; elem_gid < num_elem; elem_gid++){

        for(int k = 0; k < num_sub_1d; k++){
            for(int j = 0; j < num_sub_1d; j++){
                for(int i = 0; i < num_sub_1d; i++){

                    int cell_index = i + j*num_sub_1d + k*num_sub_1d*num_sub_1d;
                    int cell_mesh_index = cell_index + num_sub_1d*num_sub_1d*num_sub_1d*(elem_gid);
                    cell_fields(cell_mesh_index, 1) = elem_val;

                }

            }
        }
        elem_val *= -1;
    }

    for (int elem_gid = 0; elem_gid < num_elem; elem_gid++){

        for(int k = 0; k < num_sub_1d; k++){
            for(int j = 0; j < num_sub_1d; j++){
                for(int i = 0; i < num_sub_1d; i++){

                    int cell_index = i + j*num_sub_1d + k*num_sub_1d*num_sub_1d;
                    int cell_mesh_index = cell_index + num_sub_1d*num_sub_1d*num_sub_1d*(elem_gid);
                    cell_fields(cell_mesh_index, 2) = elem_gid;

                }

            }
        }
    }

    // Use average temp from each node as cell temp
    for (int cell_rid = 0; cell_rid < num_cells; cell_rid++){

        cell_fields(cell_rid, 3) = mat_pt->field(cell_rid);    

    } // end for k over cells


    // save the vertex vector fields to an array for exporting to graphics files
    auto vec_fields = CArray <real_t> (num_nodes, num_vec_vars, 3);

    for (int node_gid = 0; node_gid < num_nodes; node_gid++){
        
        vec_fields(node_gid, 0, 0) = node_coords(node_gid, 0); 
        vec_fields(node_gid, 0, 1) = node_coords(node_gid, 1);
        vec_fields(node_gid, 0, 2) = node_coords(node_gid, 2);

    } // end for loop over vertices

    //  ---------------------------------------------------------------------------
    //  Setup of file and directoring for exporting
    //  ---------------------------------------------------------------------------

    
    FILE *out[20];   // the output files that are written to
    char filename[128];
    
    struct stat st;
    
    if(stat("ensight",&st) != 0)
        system("mkdir ensight");
    
    if(stat("ensight/data",&st) != 0)
        system("mkdir ensight/data");


    //  ---------------------------------------------------------------------------
    //  Write the Geometry file
    //  ---------------------------------------------------------------------------
    
    
    sprintf(filename, "ensight/data/%s.%05d.geo", name, graphics_id);
    // filename has the full string
    
    out[0] = fopen(filename, "w");
    
    
    fprintf(out[0],"A graphics dump by Cercion \n");
    fprintf(out[0],"%s","EnSight Gold geometry\n");
    fprintf(out[0],"%s","node id assign\n");
    fprintf(out[0],"%s","element id assign\n");
    
    fprintf(out[0],"part\n");
    fprintf(out[0],"%10d\n",1);
    fprintf(out[0],"Mesh\n");
    
    
    // --- vertices ---
    fprintf(out[0],"coordinates\n");
    fprintf(out[0],"%10d\n",mesh->num_nodes());
    
    // write all components of the point coordinates
    for (int node_gid=0; node_gid<num_nodes; node_gid++){
        fprintf(out[0],"%12.5e\n",node_coords(node_gid, 0));
    }
    
    for (int node_gid=0; node_gid<num_nodes; node_gid++){
        fprintf(out[0],"%12.5e\n",node_coords(node_gid, 1));
    }
    
    for (int node_gid=0; node_gid<num_nodes; node_gid++){
        fprintf(out[0],"%12.5e\n",node_coords(node_gid, 2));
    }
    
    // convert_vert_list_ord_Ensight
    // --- cells ---
    fprintf(out[0],"hexa8\n");
    fprintf(out[0],"%10d\n",num_cells);
    int this_index;
    
    // write all global point numbers for this cell
    for (int cell_rid = 0; cell_rid<num_cells; cell_rid++) {

        for (int j=0; j<8; j++){
            this_index = convert_vert_list_ord_Ensight(j);
            fprintf(out[0],"%10d\t",nodes_in_elem(cell_rid, this_index)+1); // note node_gid starts at 1

        }
        fprintf(out[0],"\n");
    }
  
    fclose(out[0]);
    
 
    // ---------------------------------------------------------------------------
    // Write the Scalar variable files
    // ---------------------------------------------------------------------------
    // ensight_vars = (den, pres,...)
    for (int var = 0; var < num_scalar_vars; var++){
        
        // write a scalar value
        sprintf(filename,"ensight/data/%s.%05d.%s", name, graphics_id, scalar_var_names[var]);

        out[0]=fopen(filename,"w");
        
        fprintf(out[0],"Per_elem scalar values\n");
        fprintf(out[0],"part\n");
        fprintf(out[0],"%10d\n",1);
        
        fprintf(out[0],"hexa8\n");  // e.g., hexa8
        
        for (int cell_rid=0; cell_rid<num_cells; cell_rid++) {
            fprintf(out[0],"%12.5e\n", cell_fields(cell_rid, var));
        }
        
        fclose(out[0]);
        
    } // end for var

    //  ---------------------------------------------------------------------------
    //  Write the Vector variable files
    //  ---------------------------------------------------------------------------
    
    // ensight vector vars = (position, velocity, force)
    for (int var=0; var<num_vec_vars; var++){
        
        sprintf(filename,"ensight/data/%s.%05d.%s", name, graphics_id, vec_var_names[var]);
        
        out[0]=fopen(filename,"w");
        fprintf(out[0],"Per_node vector values\n");
        fprintf(out[0],"part\n");
        fprintf(out[0],"%10d\n",1);
        fprintf(out[0],"coordinates\n");
        
        for (int node_gid=0; node_gid<num_nodes; node_gid++){
            fprintf(out[0],"%12.5e\n",vec_fields(node_gid, var, 0));
        }
        
        for (int node_gid=0; node_gid<num_nodes; node_gid++){
            fprintf(out[0],"%12.5e\n",vec_fields(node_gid, var, 1));
        }
        
        for (int node_gid=0; node_gid<num_nodes; node_gid++){
            fprintf(out[0],"%12.5e\n",vec_fields(node_gid, var, 2));
        }
        
        fclose(out[0]);
    } // end for var

    // ---------------------------------------------------------------------------
    // Write the case file
    // ---------------------------------------------------------------------------
    
    sprintf(filename,"ensight/%s.case",name);
    out[0]=fopen(filename,"w");
    
    fprintf(out[0],"FORMAT\n");
    fprintf(out[0],"type: ensight gold\n");
    fprintf(out[0],"GEOMETRY\n");
    
    sprintf(filename,"model: data/%s.*****.geo\n",name);
    fprintf(out[0],"%s",filename);
    fprintf(out[0],"VARIABLE\n");
    
    for (int var=0; var<num_scalar_vars; var++){
        sprintf(filename,
                "scalar per element: %s data/%s.*****.%s\n",
                scalar_var_names[var], name, scalar_var_names[var]);
        fprintf(out[0],"%s",filename);
    }
    
    for (int var=0; var<num_vec_vars; var++){
        
        sprintf(filename,
                "vector per node: %s data/%s.*****.%s\n",
                vec_var_names[var], name, vec_var_names[var]);
        fprintf(out[0],"%s",filename);
    }
    
    fprintf(out[0],"TIME\n");
    fprintf(out[0],"time set: 1\n");
    fprintf(out[0],"number of steps: %4d\n",graphics_id+1);
    fprintf(out[0],"filename start number: 0\n");
    fprintf(out[0],"filename increment: 1\n");
    fprintf(out[0],"time values: \n");

    
    graphics_times[graphics_id]=TIME;
    
    for (int i=0;i<=graphics_id;i++) {
        fprintf(out[0],"%12.5e\n",graphics_times[i]);
    }
    fclose(out[0]);
    
    
    // ---------------------------------------------------------------------------
    // Done writing the graphics dump
    // ---------------------------------------------------------------------------

    // increment graphics id counter
    graphics_id++;
} // end ensight
*/

/* -----------------------------------------------------------------------------
   Initialize local views and global vectors needed to describe the design
-------------------------------------------------------------------------------- */

void Implicit_Solver::init_design(){
  int num_dim = simparam->num_dim;
  bool nodal_density_flag = simparam->nodal_density_flag;

  //set densities
  if(nodal_density_flag){
    if(!simparam->restart_file){
      dual_node_densities = dual_vec_array("dual_node_densities", nlocal_nodes, 1);
      host_vec_array node_densities = dual_node_densities.view_host();
      //notify that the host view is going to be modified in the file readin
      dual_node_densities.modify_host();
    
      //debug Tecplot file readin of initial densities (only works on runs with 1 MPI rank this way)
      std::string skip_line, read_line, substring;
      std::stringstream line_parse;
      real_t read_density;
      //in = new std::ifstream();
      //in->open("TecplotDensity.dat");
      //skip 3 lines
      //for (int j = 1; j <= 3; j++) {
      //getline(*in, skip_line);
      // std::cout << skip_line << std::endl;
      //}
  
    
      for(int inode = 0; inode < nlocal_nodes; inode++){
      //getline(*in,read_line);
      //line_parse.clear();
      //line_parse.str(read_line);
      //for(int iword = 0; iword < 10; iword++){
        //read portions of the line into the substring variable
        
       //if(iword==3){ line_parse >> read_density;}
        //else {line_parse >> substring;}
      //}
      //initialize densities to 1 for now; in the future there might be an option to read in an initial condition for each node
      //if(read_density < 0.3) read_density = 0.1;
      //node_densities(inode,0) = read_density;
      
      node_densities(inode,0) = 1;
    }

    //sync device view
    dual_node_densities.sync_device();
    }
    //allocate global vector information
    design_node_densities_distributed = Teuchos::rcp(new MV(map, dual_node_densities));
    all_node_densities_distributed = Teuchos::rcp(new MV(all_node_map, 1));

    //communicate ghost information to the all vector
    //create import object using local node indices map and all indices map
    Tpetra::Import<LO, GO> importer(map, all_node_map);

    //comms to get ghosts
    all_node_densities_distributed->doImport(*design_node_densities_distributed, importer, Tpetra::INSERT);

    //debug print
    //std::ostream &out = std::cout;
    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    //if(myrank==0)
    //*fos << "Node Densities with Ghosts :" << std::endl;
    //all_node_densities_distributed->describe(*fos,Teuchos::VERB_EXTREME);
    //*fos << std::endl;
    //std::fflush(stdout);
  
  }
  else{
    //initialize memory for volume storage
    vec_array Element_Densities("Element Densities", rnum_elem, 1);
    for(int ielem = 0; ielem < rnum_elem; ielem++)
      Element_Densities(ielem,0) = 1;

    //create global vector
    Global_Element_Densities = Teuchos::rcp(new MV(all_element_map, Element_Densities));

    //if(myrank==0)
    //*fos << "Global Element Densities:" << std::endl;
    //Global_Element_Densities->describe(*fos,Teuchos::VERB_EXTREME);
    //*fos << std::endl;
  }

}

/* ----------------------------------------------------------------------
   Return the CPU time for the current process in seconds very
   much in the same way as MPI_Wtime() returns the wall time.
------------------------------------------------------------------------- */

double Implicit_Solver::CPU_Time()
{
  double rv = 0.0;
/*
#ifdef _WIN32

  // from MSD docs.
  FILETIME ct,et,kt,ut;
  union { FILETIME ft; uint64_t ui; } cpu;
  if (GetProcessTimes(GetCurrentProcess(),&ct,&et,&kt,&ut)) {
    cpu.ft = ut;
    rv = cpu.ui * 0.0000001;
  }
*/

  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) == 0) {
    rv = (double) ru.ru_utime.tv_sec;
    rv += (double) ru.ru_utime.tv_usec * 0.000001;
  }


  return rv;
}

/* ----------------------------------------------------------------------
   Clock variable initialization
------------------------------------------------------------------------- */

void Implicit_Solver::init_clock(){
  double current_cpu = 0;
  initial_CPU_time = CPU_Time();
}
