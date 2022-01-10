/* 

Example code for smoothing some field on the mesh. 


A representative mesh is shown:

p
*---------*---------*
|         |         |
|         |         |
|    *z   |    *    |
|         |         |
|         |         |
*---------*---------*
|         |         |
|         |         |
|    *    |    *    |
|         |         |
|         |         |
*---------*---------*

The smoothing operation follows a two step process:

1. ) Loop over all the nodes (p) in a cell and 
average the field to the cell center materhial
point (z). 

2.) Loop over all of the cells (z) connected to a node (p)
and average values to the nodal field.


Each cell is within an element, and the number of cells is 
defined by the user using the p_order variable in the input

num_cells in element = (p_order*2)^3

*/

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
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

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
#include "MatrixMarket_Tpetra.hpp"
#include <set>

#include "elements.h"
#include "swage.h"
#include "matar.h"
#include "utilities.h"
#include "node_combination.h"
#include "Simulation_Parameters.h"
#include "Amesos2_Version.hpp"
#include "Amesos2.hpp"
#include "Parallel_Nonlinear_Solver.h"

//Optimization Package
#include "ROL_Algorithm.hpp"
#include "ROL_Solver.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Types.hpp"
#include "ROL_Stream.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_StdBoundConstraint.hpp"
#include "ROL_ParameterList.hpp"
#include <ROL_TpetraMultiVector.hpp>

//Objective Functions and Constraint Functions
#include "Mass_Objective.h"
#include "Mass_Constraint.h"
#include "Center_of_Mass_Constraint.h"
#include "Moment_of_Inertia_Constraint.h"
#include "Bounded_Strain_Constraint.h"
#include "Strain_Energy_Constraint.h"
#include "Strain_Energy_Minimize.h"
#include "Strain_Energy_Objective.h"

//Multigrid Solver
#include <Xpetra_Operator.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <MueLu.hpp>

#include <MueLu_BaseClass.hpp>
#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION
#include <MueLu_ExplicitInstantiation.hpp>
#endif
#include <MueLu_Level.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>
#include <DriverCore.hpp>

//debug and performance includes
#include <sys/time.h>
#include <sys/resource.h>
#include <ctime>

#define BUFFER_LINES 1000
#define MAX_WORD 30
#define MAX_ELEM_NODES 8
#define STRAIN_EPSILON 0.000000001
#define DENSITY_EPSILON 0.0001
#define BC_EPSILON 1.0e-8

using namespace utils;

/*

Swage is a reference to a swage block used in blacksmithing.  
Its a large metal block that has multiple shaps carved into 
each surface to use for hammering metal into to form it. 

*/

Parallel_Nonlinear_Solver::Parallel_Nonlinear_Solver() : Solver(){
  //create parameter object
  simparam = new Simulation_Parameters();
  // ---- Read input file, define state and boundary conditions ---- //
  simparam->input();
  //create ref element object
  ref_elem = new elements::ref_element();
  //create mesh objects
  //init_mesh = new swage::mesh_t(simparam);
  //mesh = new swage::mesh_t(simparam);

  element_select = new elements::element_selector();
  num_nodes = 0;
  update_count = 0;
  file_index = 0;

  Matrix_alloc=0;
  gradient_print_sync = 0;
  //RCP pointer to *this (Parallel Nonlinear Solver Object)
  //FEM_pass = Teuchos::rcp(this);

  //property initialization flags
  mass_init = false;
  com_init[0] = com_init[1] = com_init[2] = false;

  //property update counters
  mass_update = com_update[0] = com_update[1] = com_update[2] = -1;

  //RCP initialization
  mass_gradients_distributed = Teuchos::null;
  center_of_mass_gradients_distributed = Teuchos::null;
}

Parallel_Nonlinear_Solver::~Parallel_Nonlinear_Solver(){
   delete simparam;
   delete ref_elem;
   delete element_select;
   if(myrank==0)
   delete in;
}

//==============================================================================
//    Primary simulation runtime routine
//==============================================================================


void Parallel_Nonlinear_Solver::run(int argc, char *argv[]){
    
    std::cout << "Running TO Solver" << std::endl;
     // check to see of a mesh was supplied when running the code
    if (argc == 1) {
        std::cout << "\n\n**********************************\n\n";
        std::cout << " ERROR:\n";
        std::cout << " Please supply a mesh file as the second command line argument \n";
        std::cout << "**********************************\n\n" << std::endl;
        return;
    }
    //MPI info
    world = MPI_COMM_WORLD; //used for convenience to represent all the ranks in the job
    MPI_Comm_rank(world,&myrank);
    MPI_Comm_size(world,&nranks);

    //initialize Trilinos communicator class
    comm = Tpetra::getDefaultComm();

    //error handle for file input name
    //if(argc < 2)

    // ---- Read intial mesh, refine, and build connectivity ---- //
    if(simparam->tecplot_input)
      read_mesh_tecplot(argv[1]);
    else
      read_mesh_ensight(argv[1]);

    init_maps();
    
    std::cout << "Num elements = " << rnum_elem << std::endl;
    
    //initialize timing
    if(simparam->report_runtime_flag)
    init_clock();
    
    // ---- Find Boundaries on mesh ---- //
    generate_bcs();
    if(myrank == 0)
    std::cout << "Starting init assembly" << std::endl <<std::flush;
    //allocate and fill sparse structures needed for global solution
    init_assembly();

    //initialize and initialize TO design variable storage
    init_design();
    
    //assemble the global solution (stiffness matrix etc. and nodal forces)
    assemble_matrix();

    if(myrank == 0)
    std::cout << "Finished matrix assembly" << std::endl <<std::flush;
    
    assemble_vector();
    //return;
    //find each element's volume
    compute_element_volumes();

    linear_solver_parameters();
    
    if(myrank == 0)
    std::cout << "Starting First Solve" << std::endl <<std::flush;
    
    int solver_exit = solve();
    if(solver_exit == EXIT_SUCCESS){
      std::cout << "Linear Solver Error" << std::endl <<std::flush;
      return;
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
    std::cout << " RUNTIME OF CODE ON TASK " << myrank << " is "<< current_cpu-initial_CPU_time <<std::endl;
    //debug return to avoid printing further

    real_t dt = simparam->dt;
    int cycle_stop = simparam->cycle_stop;
    real_t &TIME = simparam->TIME;
    real_t TFINAL = simparam->TFINAL;
    int &cycle = simparam->cycle;
    int graphics_cyc_ival = simparam->graphics_cyc_ival;
    real_t graphics_dt_ival = simparam->graphics_dt_ival;
    real_t graphics_time = simparam->graphics_dt_ival;
    real_t &percent_comp = simparam->percent_comp;

    // Data writers
    tecplot_writer();
    // vtk_writer();
    if(myrank==0){
      std::cout << "Total number of solves and assembly " << update_count <<std::endl;
      std::cout << "End of Optimization" << std::endl;
    }
}

/* ----------------------------------------------------------------------
   Read Ensight format mesh file
------------------------------------------------------------------------- */
void Parallel_Nonlinear_Solver::read_mesh_ensight(char *MESH){

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
  dual_nodal_forces = dual_vec_array("dual_nodal_forces", nlocal_nodes*num_dim,1);

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
  
  std::cout << "RNUM ELEMENTS IS: " << rnum_elem << std::endl;
  //copy temporary element storage to multivector storage
  max_nodes_per_element = MAX_ELEM_NODES;

  dual_nodes_in_elem = dual_elem_conn_array("dual_nodes_in_elem", rnum_elem, max_nodes_per_element);
  nodes_in_elem = dual_nodes_in_elem.view_host();
  dual_nodes_in_elem.modify_host();

  Element_Types = CArrayKokkos<elements::elem_types::elem_type, array_layout, HostSpace, memory_traits>(rnum_elem);
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

  //simplified for now
  for(int ielem = 0; ielem < rnum_elem; ielem++)
    Element_Types(ielem) = elements::elem_types::Hex8;

  //element type selection (subject to change)
  // ---- Set Element Type ---- //
  // allocate element type memory
  //elements::elem_type_t* elem_choice;

  int NE = 1; // number of element types in problem
    
  //set base type pointer to one of the existing derived type object references
  if(simparam->num_dim==2)
  element_select->choose_2Delem_type(Element_Types(0), elem2D);
  else if(simparam->num_dim==3)
  element_select->choose_3Delem_type(Element_Types(0), elem);

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
void Parallel_Nonlinear_Solver::read_mesh_tecplot(char *MESH){

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
  dual_nodal_forces = dual_vec_array("dual_nodal_forces", nlocal_nodes*num_dim,1);

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
  std::cout<<"before initial mesh initialization"<<std::endl;
  
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
  std::cout << "BUFFER ITERATIONS IS: " << buffer_iterations << std::endl;
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
  max_nodes_per_element = MAX_ELEM_NODES;

  dual_nodes_in_elem = dual_elem_conn_array("dual_nodes_in_elem", rnum_elem, max_nodes_per_element);
  nodes_in_elem = dual_nodes_in_elem.view_host();
  dual_nodes_in_elem.modify_host();

  Element_Types = CArrayKokkos<elements::elem_types::elem_type, array_layout, HostSpace, memory_traits>(rnum_elem, array_layout, HostSpace, memory_traits);
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

  //simplified for now
  for(int ielem = 0; ielem < rnum_elem; ielem++)
    Element_Types(ielem) = elements::elem_types::Hex8;

  //element type selection (subject to change)
  // ---- Set Element Type ---- //
  // allocate element type memory
  //elements::elem_type_t* elem_choice;

  int NE = 1; // number of element types in problem
    
  //set base type pointer to one of the existing derived type object references
  if(simparam->num_dim==2)
  element_select->choose_2Delem_type(Element_Types(0), elem2D);
  else if(simparam->num_dim==3)
  element_select->choose_3Delem_type(Element_Types(0), elem);

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
   Initialize Ghost and Non-Overlapping Element Maps
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::init_maps(){
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
  node_displacements_distributed = Teuchos::rcp(new MV(local_dof_map, 1));
  all_node_displacements_distributed = Teuchos::rcp(new MV(all_dof_map, 1));
  //all_node_nconn_distributed = Teuchos::rcp(new MCONN(all_node_map, 1));
  if(num_dim==3) strain_count = 6;
  else strain_count = 3;
  node_strains_distributed = Teuchos::rcp(new MV(map, strain_count));
  all_node_strains_distributed = Teuchos::rcp(new MV(all_node_map, strain_count));
  Global_Nodal_Forces = Teuchos::rcp(new MV(local_dof_map, dual_nodal_forces));

  //initialize displacements to 0
  //local variable for host view in the dual view
  host_vec_array all_node_displacements = all_node_displacements_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  host_vec_array node_displacements = node_displacements_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  for(int init = 0; init < local_dof_map->getNodeNumElements(); init++)
    node_displacements(init,0) = 0;
  for(int init = 0; init < all_dof_map->getNodeNumElements(); init++)
    all_node_displacements(init,0) = 0;
  
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
  std::cout << "End of map setup " << std::endl;
}

/* ----------------------------------------------------------------------
   Setup Optimization Problem Object, Relevant Objective, and Constraints
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::setup_optimization_problem(){
  int num_dim = simparam->num_dim;
  bool nodal_density_flag = simparam->nodal_density_flag;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Surface_Nodes;
  GO current_node_index;
  LO local_node_index;
  int num_bdy_patches_in_set;
  size_t node_id, patch_id;
  int num_boundary_sets = num_boundary_conditions;
  int current_element_index, local_surface_id;
  const_host_vec_array design_densities;
  typedef ROL::TpetraMultiVector<real_t,LO,GO,node_type> ROL_MV;
  
  // fill parameter list with desired algorithmic options or leave as default
  // Read optimization input parameter list.
  std::string filename = "input_ex01.xml";
  auto parlist = ROL::getParametersFromXmlFile( filename );
  //ROL::ParameterList parlist;

  // Objective function
  ROL::Ptr<ROL::Objective<real_t>> obj = ROL::makePtr<StrainEnergyMinimize_TopOpt>(this, nodal_density_flag);
  //Design variables to optimize
  ROL::Ptr<ROL::Vector<real_t>> x;
  if(nodal_density_flag)
    x = ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO>>(design_node_densities_distributed);
  else
    x = ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO>>(Global_Element_Densities);
  //optimization problem interface that can have constraints added to it before passing to solver object
   ROL::Ptr<ROL::Problem<real_t>> problem = ROL::makePtr<ROL::Problem<real_t>>(obj,x);

  //ROL::Ptr<ROL::Constraint<double>> lin_econ = ROL::makePtr<MyLinearEqualityConstraint<double>>();
  //ROL::Ptr<ROL::Vector<double>      lin_emul = ROL::makePtr<MyLinearEqualityConstraintMultiplier<double>>();
  //problem.addLinearConstraint("Linear Equality Constraint",lin_econ,lin_mul);

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
    for(int iboundary = 0; iboundary < num_boundary_sets; iboundary++){

      num_bdy_patches_in_set = NBoundary_Condition_Patches(iboundary);

      //loop over boundary patches for this boundary set
      for (int bdy_patch_gid = 0; bdy_patch_gid < num_bdy_patches_in_set; bdy_patch_gid++){
                
      // get the global id for this boundary patch
      patch_id = Boundary_Condition_Patches(iboundary, bdy_patch_gid);
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
      }
      }
    }
  
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

  //ROL::Ptr<ROL::Constraint<double>>     lin_icon = ROL::makePtr<MyLinearInequalityConstraint<double>>();
  //ROL::Ptr<ROL::Vector<double>>         lin_imul = ROL::makePtr<MyLinearInequalityConstraintMultiplier<double>>();
  //ROL::Ptr<ROL:BoundConstraint<double>> lin_ibnd = ROL::makePtr<MyLinearInequalityConstraintBound<double>>();
  //problem.addLinearConstraint("Linear Inequality Constraint",lin_icon,lin_imul,lin_ibnd);
    
  // TypeG (generally constrained) specification
  //ROL::Ptr<ROL::Constraint<double>> econ = ROL::makePtr<MyEqualityConstraint<double>>();
  //ROL::Ptr<ROL::Vector<double>>     emul = ROL::makePtr<MyEqualityConstraintMultiplier<double>>();
  //problem.addConstraint("Equality Constraint",econ,emul);
  ROL::Ptr<std::vector<real_t> > li_ptr = ROL::makePtr<std::vector<real_t>>(1,0.0);
  ROL::Ptr<std::vector<real_t> > li_ptr2 = ROL::makePtr<std::vector<real_t>>(1,0.0);
  ROL::Ptr<std::vector<real_t> > li_ptr3 = ROL::makePtr<std::vector<real_t>>(1,0.0);
  ROL::Ptr<std::vector<real_t> > ll_ptr = ROL::makePtr<std::vector<real_t>>(1,0.0);
  ROL::Ptr<std::vector<real_t> > lu_ptr = ROL::makePtr<std::vector<real_t>>(1,0.15);

  ROL::Ptr<ROL::Vector<real_t> > constraint_mul = ROL::makePtr<ROL::StdVector<real_t>>(li_ptr);
  ROL::Ptr<ROL::Vector<real_t> > constraint_mul2 = ROL::makePtr<ROL::StdVector<real_t>>(li_ptr2);
  ROL::Ptr<ROL::Vector<real_t> > constraint_mul3 = ROL::makePtr<ROL::StdVector<real_t>>(li_ptr3);
  ROL::Ptr<ROL::Vector<real_t> > ll = ROL::makePtr<ROL::StdVector<real_t>>(ll_ptr);
  ROL::Ptr<ROL::Vector<real_t> > lu = ROL::makePtr<ROL::StdVector<real_t>>(lu_ptr);
  
  //ROL::Ptr<ROL::Constraint<real_t>> ineq_constraint = ROL::makePtr<MassConstraint_TopOpt>(this, nodal_density_flag);
  //compute initial mass
  ROL::Ptr<ROL_MV> ROL_Element_Masses = ROL::makePtr<ROL_MV>(Global_Element_Masses);
  ROL::Elementwise::ReductionSum<real_t> sumreduc;
  if(nodal_density_flag)
    design_densities = design_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  else
    design_densities = Global_Element_Densities->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  compute_element_masses(design_densities,true);
  
  real_t initial_mass = ROL_Element_Masses->reduce(sumreduc);

  //define constraint objects
  ROL::Ptr<ROL::Constraint<real_t>> eq_constraint = ROL::makePtr<MassConstraint_TopOpt>(this, nodal_density_flag, false, 0.2);
  ROL::Ptr<ROL::Constraint<real_t>> eq_constraint2 = ROL::makePtr<MomentOfInertiaConstraint_TopOpt>(this, nodal_density_flag, 2, false, 0.4);
  ROL::Ptr<ROL::Constraint<real_t>> eq_constraint3 = ROL::makePtr<MomentOfInertiaConstraint_TopOpt>(this, nodal_density_flag, 3, false);
  //ROL::Ptr<ROL::Constraint<real_t>> ineq_constraint = ROL::makePtr<MassConstraint_TopOpt>(this, nodal_density_flag);
  ROL::Ptr<ROL::BoundConstraint<real_t>> constraint_bnd = ROL::makePtr<ROL::Bounds<real_t>>(ll,lu);
  //problem->addConstraint("Inequality Constraint",ineq_constraint,constraint_mul,constraint_bnd);
  problem->addConstraint("equality Constraint 1",eq_constraint,constraint_mul);
  problem->addConstraint("equality Constraint 2",eq_constraint2,constraint_mul2);
  problem->addConstraint("equality Constraint 3",eq_constraint3,constraint_mul3);
  //problem->addLinearConstraint("Equality Constraint",eq_constraint,constraint_mul);
  problem->setProjectionAlgorithm(*parlist);
  //finalize problem
  problem->finalize(false,true,std::cout);
  //problem->check(true,std::cout);

  //debug checks
  ROL::Ptr<ROL::TpetraMultiVector<real_t,LO,GO,node_type>> rol_x =
   ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO,node_type>>(design_node_densities_distributed);
  //construct direction vector for check
  Teuchos::RCP<MV> directions_distributed = Teuchos::rcp(new MV(map, 1));
  directions_distributed->putScalar(0.1);
  ROL::Ptr<ROL::TpetraMultiVector<real_t,LO,GO,node_type>> rol_d =
  ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO,node_type>>(directions_distributed);
  //obj->checkGradient(*rol_x, *rol_d);
  //directions_distributed->putScalar(-0.000001);
  //obj->checkGradient(*rol_x, *rol_d);
  //directions_distributed->putScalar(-0.0000001);
  //obj->checkGradient(*rol_x, *rol_d);
  

  // Instantiate Solver.
  ROL::Solver<real_t> solver(problem,*parlist);
    
  // Solve optimization problem.
  //std::ostream outStream;
  solver.solve(std::cout);

  //print mass constraint for final design vector
  compute_element_masses(design_densities,false);
  real_t final_mass = ROL_Element_Masses->reduce(sumreduc);
  if(myrank==0)
    std::cout << "Final Mass Constraint is " << final_mass/initial_mass << std::endl;
}

/* ----------------------------------------------------------------------
   Find boundary surface segments that belong to this MPI rank
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::Get_Boundary_Patches(){
  size_t npatches_repeat, npatches, element_npatches, num_nodes_in_patch, node_gid;
  int local_node_id;
  int num_dim = simparam->num_dim;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Surface_Nodes;
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
  std::cout << "Done with boundary patch allocation" << std::endl <<std::flush;
  //initialize boundary patch flags
  for(int init = 0; init < npatches_repeat; init++)
    Patch_Boundary_Flags(init) = 1;

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
      Surface_Nodes = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_nodes_in_patch, "Surface_Nodes");
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
      //if the set determines this is a duplicate access the original element's patch id and set flag to 0
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
      Surface_Nodes = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_nodes_in_patch, "Surface_Nodes");
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
      if(my_rank_flag)
        Boundary_Patches(nboundary_patches++) = Patch_Nodes(ipatch);
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

/* ----------------------------------------------------------------------
   Assign sets of element boundary surfaces corresponding to user BCs
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::generate_bcs(){
    
  // build boundary mesh patches
  //mesh->build_bdy_patches();
  //std::cout << "Starting boundary patch setup" << std::endl <<std::flush;
  Get_Boundary_Patches();
  //std::cout << "Done with boundary patch setup" << std::endl <<std::flush;
  std::cout << "number of boundary patches on task " << myrank << " = " << nboundary_patches << std::endl;
  std::cout << "building boundary sets " << std::endl;
  // set the number of boundary sets
    
  int num_boundary_sets = simparam->NB;
  int num_surface_force_sets = simparam->NBSF;
  int num_surface_disp_sets = simparam->NBD;
  int num_dim = simparam->num_dim;
  int current_bdy_id = 0;
  int bdy_set_id;
  int surf_force_set_id = 0;
  int surf_disp_set_id = 0;

  init_boundary_sets(num_boundary_sets);
  Boundary_Condition_Type_List = CArrayKokkos<int, array_layout, HostSpace, memory_traits>(num_boundary_sets); 
  Boundary_Surface_Force_Densities = CArrayKokkos<real_t, array_layout, HostSpace, memory_traits>(num_surface_force_sets,3);
  Boundary_Surface_Displacements = CArrayKokkos<real_t, array_layout, HostSpace, memory_traits>(num_surface_disp_sets,3);
  //initialize
  for(int ibdy=0; ibdy < num_boundary_sets; ibdy++) Boundary_Condition_Type_List(ibdy) = NONE;
    
  // tag the z=0 plane,  (Direction, value, bdy_set)
  std::cout << "tagging z = 0 " << std::endl;
  int bc_tag = 2;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  real_t value = 0.0 * simparam->unit_scaling;
  real_t fix_limits[4];
  fix_limits[0] = fix_limits[2] = 4;
  fix_limits[1] = fix_limits[3] = 6;
  bdy_set_id = current_bdy_id++;
  tag_boundaries(bc_tag, value, bdy_set_id, fix_limits);
  Boundary_Condition_Type_List(bdy_set_id) = DISPLACEMENT_CONDITION;
  Boundary_Surface_Displacements(surf_disp_set_id,0) = 0;
  Boundary_Surface_Displacements(surf_disp_set_id,1) = 0;
  Boundary_Surface_Displacements(surf_disp_set_id,2) = 0;
  surf_disp_set_id++;
    
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
  std::cout << std::endl;
  /*
  //This part should be changed so it interfaces with simparam to handle multiple input cases
  // tag the y=0 plane,  (Direction, value, bdy_set)
  std::cout << "tagging y = 0 " << std::endl;
  bc_tag = 1;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 0.0;
  bdy_set_id = 1;
  mesh->tag_bdys(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = DISPLACEMENT_CONDITION;
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << mesh->num_bdy_patches_in_set(bdy_set_id) << std::endl;
  std::cout << std::endl;
    

  // tag the z=0 plane,  (Direction, value, bdy_set)
  std::cout << "tagging z = 0 " << std::endl;
  bc_tag = 2;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 0.0;
  bdy_set_id = 2;
  mesh->tag_bdys(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = DISPLACEMENT_CONDITION;
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << mesh->num_bdy_patches_in_set(bdy_set_id) << std::endl;
  std::cout << std::endl;
  
  */
  /*
  std::cout << "tagging z = 2 Force " << std::endl;
  bc_tag = 2;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 2 * simparam->unit_scaling;
  //value = 2;
  bdy_set_id = current_bdy_id++;
  //find boundary patches this BC corresponds to
  tag_boundaries(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = LOADING_CONDITION;
  Boundary_Surface_Force_Densities(surf_force_set_id,0) = 0;
  Boundary_Surface_Force_Densities(surf_force_set_id,1) = 0;
  Boundary_Surface_Force_Densities(surf_force_set_id,2) = 10/simparam->unit_scaling/simparam->unit_scaling;
  surf_force_set_id++;
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
  std::cout << std::endl;
  
  
  std::cout << "tagging z = 1 Force " << std::endl;
  bc_tag = 2;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 1 * simparam->unit_scaling;
  //value = 2;
  bdy_set_id = current_bdy_id++;
  //find boundary patches this BC corresponds to
  tag_boundaries(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = LOADING_CONDITION;
  Boundary_Surface_Force_Densities(surf_force_set_id,0) = 0;
  Boundary_Surface_Force_Densities(surf_force_set_id,1) = 0;
  Boundary_Surface_Force_Densities(surf_force_set_id,2) = 10/simparam->unit_scaling/simparam->unit_scaling;
  surf_force_set_id++;
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
  std::cout << std::endl;
  
  
  std::cout << "tagging beam x = 0 " << std::endl;
  bc_tag = 2;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 2 * simparam->unit_scaling;
  //value = 2;
  bdy_set_id = current_bdy_id++;
  //find boundary patches this BC corresponds to
  tag_boundaries(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = LOADING_CONDITION;
  Boundary_Surface_Force_Densities(surf_force_set_id,0) = 0;
  Boundary_Surface_Force_Densities(surf_force_set_id,1) = 0;
  Boundary_Surface_Force_Densities(surf_force_set_id,2) = 1/simparam->unit_scaling/simparam->unit_scaling;
  surf_force_set_id++;
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
  std::cout << std::endl;
  */
  /*
  std::cout << "tagging beam -x " << std::endl;
  bc_tag = 0;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 0 * simparam->unit_scaling;
  //value = 2;
  bdy_set_id = current_bdy_id++;
  //find boundary patches this BC corresponds to
  tag_boundaries(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = LOADING_CONDITION;
  Boundary_Surface_Force_Densities(surf_force_set_id,0) = 0;
  Boundary_Surface_Force_Densities(surf_force_set_id,1) = -1/simparam->unit_scaling/simparam->unit_scaling;
  Boundary_Surface_Force_Densities(surf_force_set_id,2) = 0;
  surf_force_set_id++;
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
  std::cout << std::endl;

  std::cout << "tagging beam +x " << std::endl;
  bc_tag = 0;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 10 * simparam->unit_scaling;
  //value = 2;
  bdy_set_id = current_bdy_id++;
  //find boundary patches this BC corresponds to
  tag_boundaries(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = LOADING_CONDITION;
  Boundary_Surface_Force_Densities(surf_force_set_id,0) = 0;
  Boundary_Surface_Force_Densities(surf_force_set_id,1) = 1/simparam->unit_scaling/simparam->unit_scaling;
  Boundary_Surface_Force_Densities(surf_force_set_id,2) = 0;
  surf_force_set_id++;
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
  std::cout << std::endl;
  */
  std::cout << "tagging beam -y " << std::endl;
  bc_tag = 1;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 0 * simparam->unit_scaling;
  real_t load_limits_left[4];
  load_limits_left[0] = load_limits_left[2] = 4;
  load_limits_left[1] = load_limits_left[3] = 6;
  //value = 2;
  bdy_set_id = current_bdy_id++;
  //find boundary patches this BC corresponds to
  tag_boundaries(bc_tag, value, bdy_set_id ,load_limits_left);
  Boundary_Condition_Type_List(bdy_set_id) = LOADING_CONDITION;
  Boundary_Surface_Force_Densities(surf_force_set_id,0) = 1/simparam->unit_scaling/simparam->unit_scaling;
  Boundary_Surface_Force_Densities(surf_force_set_id,1) = 0;
  Boundary_Surface_Force_Densities(surf_force_set_id,2) = 0;
  surf_force_set_id++;
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
  std::cout << std::endl;

  std::cout << "tagging beam +y " << std::endl;
  bc_tag = 1;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 10 * simparam->unit_scaling;
  real_t load_limits_right[4];
  load_limits_right[0] = load_limits_right[2] = 4;
  load_limits_right[1] = load_limits_right[3] = 6;
  //value = 2;
  bdy_set_id = current_bdy_id++;
  //find boundary patches this BC corresponds to
  tag_boundaries(bc_tag, value, bdy_set_id, load_limits_right);
  Boundary_Condition_Type_List(bdy_set_id) = LOADING_CONDITION;
  Boundary_Surface_Force_Densities(surf_force_set_id,0) = -1/simparam->unit_scaling/simparam->unit_scaling;
  Boundary_Surface_Force_Densities(surf_force_set_id,1) = 0;
  Boundary_Surface_Force_Densities(surf_force_set_id,2) = 0;
  surf_force_set_id++;
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
  std::cout << std::endl;
  
  
  
  /*
  std::cout << "tagging beam +z force " << std::endl;
  bc_tag = 2;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  //value = 0;
  value = 100;
  bdy_set_id = current_bdy_id++;
  //find boundary patches this BC corresponds to
  tag_boundaries(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = LOADING_CONDITION;
  Boundary_Surface_Force_Densities(surf_force_set_id,0) = 1/simparam->unit_scaling/simparam->unit_scaling;
  Boundary_Surface_Force_Densities(surf_force_set_id,1) = 0;
  Boundary_Surface_Force_Densities(surf_force_set_id,2) = 0;
  surf_force_set_id++;
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
  std::cout << std::endl;
  
  
  std::cout << "tagging y = 2 " << std::endl;
  bc_tag = 1;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 2.0;
  bdy_set_id = 4;
  mesh->tag_bdys(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = LOADING_CONDITION;
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << mesh->num_bdy_patches_in_set(bdy_set_id) << std::endl;
  std::cout << std::endl;

  std::cout << "tagging z = 2 " << std::endl;
  bc_tag = 2;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 2.0;
  bdy_set_id = 5;
  mesh->tag_bdys(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = LOADING_CONDITION;
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << mesh->num_bdy_patches_in_set(bdy_set_id) << std::endl;
  std::cout << std::endl;
  */

  //allocate nodal data
  Node_DOF_Boundary_Condition_Type = CArrayKokkos<int, array_layout, device_type, memory_traits>(nall_nodes*num_dim, "Node_DOF_Boundary_Condition_Type");
  Node_DOF_Displacement_Boundary_Conditions = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(nall_nodes*num_dim);
  Node_DOF_Force_Boundary_Conditions = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(nall_nodes*num_dim);

  //initialize
  for(int init=0; init < nall_nodes*num_dim; init++)
    Node_DOF_Boundary_Condition_Type(init) = NONE;

  //Tag nodes for Boundary conditions such as displacements
  Displacement_Boundary_Conditions();
} // end generate_bcs

/* ----------------------------------------------------------------------
   initialize storage for element boundary surfaces corresponding to user BCs
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::init_boundary_sets (int num_sets){
    
  num_boundary_conditions = num_sets;
  if(num_sets == 0){
    std::cout << " Warning: number of boundary conditions = 0";
    return;
  }
  Boundary_Condition_Patches_strides = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_sets, "Boundary_Condition_Patches_strides");
  NBoundary_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_sets, "NBoundary_Condition_Patches");
  Boundary_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_sets, nboundary_patches, "Boundary_Condition_Patches");

  //initialize data
  for(int iset = 0; iset < num_sets; iset++) NBoundary_Condition_Patches(iset) = 0;
}

/* ----------------------------------------------------------------------
   find which boundary patches correspond to the given BC.
   bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
   val = plane value, cylinder radius, shell radius
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::tag_boundaries(int bc_tag, real_t val, int bdy_set, real_t *patch_limits = NULL){
  
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
      Boundary_Condition_Patches(bdy_set,counter) = iboundary_patch;
      counter ++;
    }
  } // end for bdy_patch
    
  // save the number of bdy patches in the set
  NBoundary_Condition_Patches(bdy_set) = counter;
    
  std::cout << " tagged boundary patches " << std::endl;
}

/* ----------------------------------------------------------------------
   routine for checking to see if a patch is on a boundary set
   bc_tag = 0 xplane, 1 yplane, 3 zplane, 4 cylinder, 5 is shell
   val = plane value, radius, radius
------------------------------------------------------------------------- */


int Parallel_Nonlinear_Solver::check_boundary(Node_Combination &Patch_Nodes, int bc_tag, real_t val, real_t *patch_limits){
  
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

void Parallel_Nonlinear_Solver::collect_information(){
  size_t nreduce_dof = 0;
  size_t nreduce_nodes = 0;
  size_t nreduce_elem = 0;
  size_t num_dim = simparam->num_dim;
  int strain_count;
  int output_strain_flag = simparam->output_strain_flag;

  //collect nodal coordinate information
  if(myrank==0) nreduce_nodes = num_nodes;
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > global_reduce_map =
    Teuchos::rcp(new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),nreduce_nodes,0,comm));
  
  //importer from local node distribution to collected distribution
  Tpetra::Import<LO, GO> node_collection_importer(map, global_reduce_map);

  Teuchos::RCP<MV> collected_node_coords_distributed = Teuchos::rcp(new MV(global_reduce_map, num_dim));

  //comms to collect
  collected_node_coords_distributed->doImport(*node_coords_distributed, node_collection_importer, Tpetra::INSERT);
  
  //collect nodal displacement information
  if(myrank==0) nreduce_dof = num_nodes*num_dim;
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > global_reduce_dof_map =
    Teuchos::rcp(new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),nreduce_dof,0,comm));
  
  //importer from local node distribution to collected distribution
  Tpetra::Import<LO, GO> dof_collection_importer(local_dof_map, global_reduce_dof_map);

  Teuchos::RCP<MV> collected_node_displacements_distributed = Teuchos::rcp(new MV(global_reduce_dof_map, 1));

  //comms to collect
  collected_node_displacements_distributed->doImport(*node_displacements_distributed, dof_collection_importer, Tpetra::INSERT);

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
    collected_node_displacements = collected_node_displacements_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    collected_node_densities = collected_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    collected_nodes_in_elem = collected_nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  }

  //collect strain data
  if(output_strain_flag){
    if(num_dim==3) strain_count = 6;
    else strain_count = 3;

    //importer for strains, all nodes to global node set on rank 0
    //Tpetra::Import<LO, GO> strain_collection_importer(all_node_map, global_reduce_map);

    //collected nodal density information
    Teuchos::RCP<MV> collected_node_strains_distributed = Teuchos::rcp(new MV(global_reduce_map, strain_count));

    //comms to collect
    collected_node_strains_distributed->doImport(*node_strains_distributed, node_collection_importer, Tpetra::INSERT);

    //debug print
    //std::ostream &out = std::cout;
    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    //if(myrank==0)
    //*fos << "Collected nodal displacements :" << std::endl;
    //collected_node_strains_distributed->describe(*fos,Teuchos::VERB_EXTREME);
    //*fos << std::endl;
    //std::fflush(stdout);

    //host view to print from
    if(myrank==0)
    collected_node_strains = collected_node_strains_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  }
  
}

/* ----------------------------------------------------------------------
   Output Model Information in tecplot format
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::tecplot_writer(){
  
  size_t num_dim = simparam->num_dim;
	std::string current_file_name;
	std::string base_file_name= "TecplotTO";
  std::string base_file_name_undeformed= "TecplotTO_undeformed";
	std::stringstream ss;
	std::string file_extension= ".dat";
  std::string file_count;
	std::stringstream count_temp;
  int time_step = 0;
  int temp_convert;
  int displace_geometry = 1;
  int output_strain_flag = simparam->output_strain_flag;
  if(output_strain_flag) compute_nodal_strains();
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
  if(displace_geometry){
    if(myrank==0){
      //initial undeformed geometry
      count_temp.str("");
      count_temp << file_index;
      //file_index++;
	    file_count = count_temp.str();
    
      current_file_name = base_file_name_undeformed + file_count + file_extension;
      std::ofstream myfile (current_file_name.c_str()); //output filestream object for file output
	    //read in position data
	    myfile << std::fixed << std::setprecision(8);
		
		  //output header of the tecplot file

		  myfile << "TITLE=\"results for TO code\"" "\n";
      if(output_strain_flag)
      myfile << "VARIABLES = \"x\", \"y\", \"z\", \"density\", \"sigmaxx\", \"sigmayy\", \"sigmazz\", \"sigmaxy\", \"sigmaxz\", \"sigmayz\"" "\n";
      else
		  myfile << "VARIABLES = \"x\", \"y\", \"z\", \"density\"" "\n";

		  myfile << "ZONE T=\"load step " << time_step << "\", NODES= " << num_nodes
			  << ", ELEMENTS= " << num_elem << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" "\n";

		  for (int nodeline = 0; nodeline < num_nodes; nodeline++) {
			  myfile << std::setw(25) << collected_node_coords(nodeline,0) << " ";
			  myfile << std::setw(25) << collected_node_coords(nodeline,1) << " ";
        if(num_dim==3)
			  myfile << std::setw(25) << collected_node_coords(nodeline,2) << " ";
        if(output_strain_flag){
          myfile << std::setw(25) << collected_node_densities(nodeline,0) << " ";
          myfile << std::setw(25) << collected_node_strains(nodeline,0) << " ";
          myfile << std::setw(25) << collected_node_strains(nodeline,1) << " ";
          myfile << std::setw(25) << collected_node_strains(nodeline,2) << " ";
          if(num_dim==3){
            myfile << std::setw(25) << collected_node_strains(nodeline,3) << " ";
            myfile << std::setw(25) << collected_node_strains(nodeline,4) << " ";
            myfile << std::setw(25) << collected_node_strains(nodeline,5) << std::endl;
          }
        }
        else
          myfile << std::setw(25) << collected_node_densities(nodeline,0) << std::endl;
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
    if(myrank==0){
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

		  myfile << "TITLE=\"results for TO code\"" "\n";
		  if(output_strain_flag)
      myfile << "VARIABLES = \"x\", \"y\", \"z\", \"density\", \"sigmaxx\", \"sigmayy\", \"sigmazz\", \"sigmaxy\", \"sigmaxz\", \"sigmayz\"" "\n";
      else
		  myfile << "VARIABLES = \"x\", \"y\", \"z\", \"density\"" "\n";

		  myfile << "ZONE T=\"load step " << time_step + 1 << "\", NODES= " << num_nodes
			<< ", ELEMENTS= " << num_elem << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" "\n";

		  for (int nodeline = 0; nodeline < num_nodes; nodeline++) {
			  myfile << std::setw(25) << collected_node_coords(nodeline,0) + collected_node_displacements(nodeline*num_dim,0) << " ";
			  myfile << std::setw(25) << collected_node_coords(nodeline,1) + collected_node_displacements(nodeline*num_dim + 1,0) << " ";
        if(num_dim==3)
			  myfile << std::setw(25) << collected_node_coords(nodeline,2) + collected_node_displacements(nodeline*num_dim + 2,0) << " ";
			  if(output_strain_flag){
          myfile << std::setw(25) << collected_node_densities(nodeline,0) << " ";
          myfile << std::setw(25) << collected_node_strains(nodeline,0) << " ";
          myfile << std::setw(25) << collected_node_strains(nodeline,1) << " ";
          myfile << std::setw(25) << collected_node_strains(nodeline,2) << " ";
          if(num_dim==3){
            myfile << std::setw(25) << collected_node_strains(nodeline,3) << " ";
            myfile << std::setw(25) << collected_node_strains(nodeline,4) << " ";
            myfile << std::setw(25) << collected_node_strains(nodeline,5) << std::endl;
          }
        }
        else
          myfile << std::setw(25) << collected_node_densities(nodeline,0) << std::endl;
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
  else{
    if(myrank==0){
      count_temp << file_index;
      file_index++;
	    file_count = count_temp.str();
    
      current_file_name = base_file_name + file_count + file_extension;
      std::ofstream myfile (current_file_name.c_str()); //output filestream object for file output
	    //read in position data
	    myfile << std::fixed << std::setprecision(8);
		
		  //output header of the tecplot file

		  myfile << "TITLE=\"results for TO code\"" "\n";
		  myfile << "VARIABLES = \"x\", \"y\", \"z\", \"disp1\", \"disp2\", \"disp3\", \"density\"" "\n";

		  myfile << "ZONE T=\"load step " << time_step << "\", NODES= " << num_nodes
			  << ", ELEMENTS= " << num_elem << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" "\n";

		  for (int nodeline = 0; nodeline < num_nodes; nodeline++) {
			  myfile << std::setw(25) << collected_node_coords(nodeline,0) << " ";
			  myfile << std::setw(25) << collected_node_coords(nodeline,1) << " ";
        if(num_dim==3)
			  myfile << std::setw(25) << collected_node_coords(nodeline,2) << " ";
			  myfile << std::setw(25) << collected_node_displacements(nodeline*num_dim,0) << " ";
			  myfile << std::setw(25) << collected_node_displacements(nodeline*num_dim + 1,0) << " ";
        if(num_dim==3)
			  myfile << std::setw(25) << collected_node_displacements(nodeline*num_dim + 2,0) << " ";
			  myfile << std::setw(25) << collected_node_densities(nodeline,0) << std::endl;
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
  

}
/* ----------------------------------------------------------------------
   Output Model Information in vtk format
------------------------------------------------------------------------- */
/*
void Parallel_Nonlinear_Solver::vtk_writer(){
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
void Parallel_Nonlinear_Solver::ensight_writer(){
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

/* ----------------------------------------------------------------------
   Initialize global vectors and array maps needed for matrix assembly
------------------------------------------------------------------------- */
void Parallel_Nonlinear_Solver::init_assembly(){
  int num_dim = simparam->num_dim;
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  Stiffness_Matrix_Strides = CArrayKokkos<size_t, array_layout, device_type, memory_traits> (nlocal_nodes*num_dim, "Stiffness_Matrix_Strides");
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Graph_Fill(nall_nodes, "nall_nodes");
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> current_row_nodes_scanned;
  int current_row_n_nodes_scanned;
  int local_node_index, global_node_index, current_column_index;
  int max_stride = 0;
  size_t nodes_per_element;
  
  //allocate stride arrays
  CArrayKokkos <size_t, array_layout, device_type, memory_traits> Graph_Matrix_Strides_initial(nlocal_nodes, "Graph_Matrix_Strides_initial");
  Graph_Matrix_Strides = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(nlocal_nodes, "Graph_Matrix_Strides");

  //allocate storage for the sparse stiffness matrix map used in the assembly process
  Global_Stiffness_Matrix_Assembly_Map = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(rnum_elem,
                                         max_nodes_per_element,max_nodes_per_element, "Global_Stiffness_Matrix_Assembly_Map");

  //allocate array used to determine global node repeats in the sparse graph later
  CArrayKokkos <int, array_layout, device_type, memory_traits> node_indices_used(nall_nodes, "node_indices_used");

  /*allocate array that stores which column the node index occured on for the current row
    when removing repeats*/
  CArrayKokkos <size_t, array_layout, device_type, memory_traits> column_index(nall_nodes, "column_index");
  
  //initialize nlocal arrays
  for(int inode = 0; inode < nlocal_nodes; inode++){
    Graph_Matrix_Strides_initial(inode) = 0;
    Graph_Matrix_Strides(inode) = 0;
    Graph_Fill(inode) = 0;
  }

  //initialize nall arrays
  for(int inode = 0; inode < nall_nodes; inode++){
    node_indices_used(inode) = 0;
    column_index(inode) = 0;
  }
  
  //count upper bound of strides for Sparse Pattern Graph by allowing repeats due to connectivity
  if(num_dim == 2)
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
    nodes_per_element = elem2D->num_nodes();
    for (int lnode = 0; lnode < nodes_per_element; lnode++){
      global_node_index = nodes_in_elem(ielem, lnode);
      if(map->isNodeGlobalElement(global_node_index)){
        local_node_index = map->getLocalElement(global_node_index);
        Graph_Matrix_Strides_initial(local_node_index) += nodes_per_element;
      }
    }
  }

  if(num_dim == 3)
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_3Delem_type(Element_Types(ielem), elem);
    nodes_per_element = elem->num_nodes();
    for (int lnode = 0; lnode < nodes_per_element; lnode++){
      global_node_index = nodes_in_elem(ielem, lnode);
      if(map->isNodeGlobalElement(global_node_index)){
        local_node_index = map->getLocalElement(global_node_index);
        Graph_Matrix_Strides_initial(local_node_index) += nodes_per_element;
      }
    }
  }
  
  //equate strides for later
  for(int inode = 0; inode < nlocal_nodes; inode++)
    Graph_Matrix_Strides(inode) = Graph_Matrix_Strides_initial(inode);
  
  //compute maximum stride
  for(int inode = 0; inode < nlocal_nodes; inode++)
    if(Graph_Matrix_Strides_initial(inode) > max_stride) max_stride = Graph_Matrix_Strides_initial(inode);
  
  //allocate array used in the repeat removal process
  current_row_nodes_scanned = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(max_stride, "current_row_nodes_scanned");

  //allocate sparse graph with node repeats
  RaggedRightArrayKokkos<size_t, array_layout, device_type, memory_traits> Repeat_Graph_Matrix(Graph_Matrix_Strides_initial);
  RaggedRightArrayofVectorsKokkos<size_t, array_layout, device_type, memory_traits> Element_local_indices(Graph_Matrix_Strides_initial,num_dim);
  
  //Fill the initial Graph with repeats
  if(num_dim == 2)
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
    nodes_per_element = elem2D->num_nodes();
    for (int lnode = 0; lnode < nodes_per_element; lnode++){
      global_node_index = nodes_in_elem(ielem, lnode);
      if(map->isNodeGlobalElement(global_node_index)){
        local_node_index = map->getLocalElement(global_node_index);
        for (int jnode = 0; jnode < nodes_per_element; jnode++){
          current_column_index = Graph_Fill(local_node_index)+jnode;
          Repeat_Graph_Matrix(local_node_index, current_column_index) = nodes_in_elem(ielem,jnode);

          //fill inverse map
          Element_local_indices(local_node_index,current_column_index,0) = ielem;
          Element_local_indices(local_node_index,current_column_index,1) = lnode;
          Element_local_indices(local_node_index,current_column_index,2) = jnode;

          //fill forward map
          Global_Stiffness_Matrix_Assembly_Map(ielem,lnode,jnode) = current_column_index;
        }
        Graph_Fill(local_node_index) += nodes_per_element;
      }
    }
  }
  
  if(num_dim == 3)
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_3Delem_type(Element_Types(ielem), elem);
    nodes_per_element = elem->num_nodes();
    for (int lnode = 0; lnode < nodes_per_element; lnode++){
      global_node_index = nodes_in_elem(ielem, lnode);
      if(map->isNodeGlobalElement(global_node_index)){
        local_node_index = map->getLocalElement(global_node_index);
        for (int jnode = 0; jnode < nodes_per_element; jnode++){
          current_column_index = Graph_Fill(local_node_index)+jnode;
          Repeat_Graph_Matrix(local_node_index, current_column_index) = nodes_in_elem(ielem,jnode);

          //fill inverse map
          Element_local_indices(local_node_index,current_column_index,0) = ielem;
          Element_local_indices(local_node_index,current_column_index,1) = lnode;
          Element_local_indices(local_node_index,current_column_index,2) = jnode;

          //fill forward map
          Global_Stiffness_Matrix_Assembly_Map(ielem,lnode,jnode) = current_column_index;
        }
        Graph_Fill(local_node_index) += nodes_per_element;
      }
    }
  }
  
  //debug statement
  //std::cout << "started run" << std::endl;
  //std::cout << "Graph Matrix Strides Repeat on task " << myrank << std::endl;
  //for (int inode = 0; inode < nlocal_nodes; inode++)
    //std::cout << Graph_Matrix_Strides(inode) << std::endl;
  
  //remove repeats from the inital graph setup
  int current_node, current_element_index, element_row_index, element_column_index, current_stride;
  for (int inode = 0; inode < nlocal_nodes; inode++){
    current_row_n_nodes_scanned = 0;
    for (int istride = 0; istride < Graph_Matrix_Strides(inode); istride++){
      //convert global index in graph to its local index for the flagging array
      current_node = all_node_map->getLocalElement(Repeat_Graph_Matrix(inode,istride));
      //debug
      //if(current_node==-1)
      //std::cout << "Graph Matrix node access on task " << myrank << std::endl;
      //std::cout << Repeat_Graph_Matrix(inode,istride) << std::endl;
      if(node_indices_used(current_node)){
        //set global assembly map index to the location in the graph matrix where this global node was first found
        current_element_index = Element_local_indices(inode,istride,0);
        element_row_index = Element_local_indices(inode,istride,1);
        element_column_index = Element_local_indices(inode,istride,2);
        Global_Stiffness_Matrix_Assembly_Map(current_element_index,element_row_index, element_column_index) 
            = column_index(current_node);   

        
        //swap current node with the end of the current row and shorten the stride of the row
        //first swap information about the inverse and forward maps

        current_stride = Graph_Matrix_Strides(inode);
        if(istride!=current_stride-1){
        Element_local_indices(inode,istride,0) = Element_local_indices(inode,current_stride-1,0);
        Element_local_indices(inode,istride,1) = Element_local_indices(inode,current_stride-1,1);
        Element_local_indices(inode,istride,2) = Element_local_indices(inode,current_stride-1,2);
        current_element_index = Element_local_indices(inode,istride,0);
        element_row_index = Element_local_indices(inode,istride,1);
        element_column_index = Element_local_indices(inode,istride,2);

        Global_Stiffness_Matrix_Assembly_Map(current_element_index,element_row_index, element_column_index) 
            = istride;

        //now that the element map information has been copied, copy the global node index and delete the last index

        Repeat_Graph_Matrix(inode,istride) = Repeat_Graph_Matrix(inode,current_stride-1);
        }
        istride--;
        Graph_Matrix_Strides(inode)--;
      }
      else{
        /*this node hasn't shown up in the row before; add it to the list of nodes
          that have been scanned uniquely. Use this list to reset the flag array
          afterwards without having to loop over all the nodes in the system*/
        node_indices_used(current_node) = 1;
        column_index(current_node) = istride;
        current_row_nodes_scanned(current_row_n_nodes_scanned) = current_node;
        current_row_n_nodes_scanned++;
      }
    }
    //reset nodes used list for the next row of the sparse list
    for(int node_reset = 0; node_reset < current_row_n_nodes_scanned; node_reset++)
      node_indices_used(current_row_nodes_scanned(node_reset)) = 0;

  }

  //copy reduced content to non_repeat storage
  Graph_Matrix = RaggedRightArrayKokkos<size_t, array_layout, device_type, memory_traits>(Graph_Matrix_Strides);
  for(int inode = 0; inode < nlocal_nodes; inode++)
    for(int istride = 0; istride < Graph_Matrix_Strides(inode); istride++)
      Graph_Matrix(inode,istride) = Repeat_Graph_Matrix(inode,istride);

  //deallocate repeat matrix
  
  /*At this stage the sparse graph should have unique global indices on each row.
    The constructed Assembly map (to the global sparse matrix)
    is used to loop over each element's local stiffness matrix in the assembly process.*/
  
  //expand strides for stiffness matrix by multipling by dim
  for(int inode = 0; inode < nlocal_nodes; inode++){
    for (int idim = 0; idim < num_dim; idim++)
    Stiffness_Matrix_Strides(num_dim*inode + idim) = num_dim*Graph_Matrix_Strides(inode);
  }

  Stiffness_Matrix = RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout>(Stiffness_Matrix_Strides);
  DOF_Graph_Matrix = RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits> (Stiffness_Matrix_Strides);

  //set stiffness Matrix Graph
  //debug print
    //std::cout << "DOF GRAPH MATRIX ENTRIES ON TASK " << myrank << std::endl;
  for (int idof = 0; idof < num_dim*nlocal_nodes; idof++){
    for (int istride = 0; istride < Stiffness_Matrix_Strides(idof); istride++){
      DOF_Graph_Matrix(idof,istride) = Graph_Matrix(idof/num_dim,istride/num_dim)*num_dim + istride%num_dim;
      //debug print
      //std::cout << "{" <<istride + 1 << "," << DOF_Graph_Matrix(idof,istride) << "} ";
    }
    //debug print
    //std::cout << std::endl;
  }
  
  /*
  //debug print nodal positions and indices
  std::cout << " ------------NODAL POSITIONS--------------"<<std::endl;
  for (int inode = 0; inode < num_nodes; inode++){
      std::cout << "node: " << inode + 1 << " { ";
    for (int istride = 0; istride < num_dim; istride++){
        std::cout << node_coords(inode,istride) << " , ";
    }
    std::cout << " }"<< std::endl;
  }
  //debug print element edof
  
  std::cout << " ------------ELEMENT EDOF--------------"<<std::endl;

  for (int ielem = 0; ielem < rnum_elem; ielem++){
    std::cout << "elem:  " << ielem+1 << std::endl;
    for (int lnode = 0; lnode < nodes_per_elem; lnode++){
        std::cout << "{ ";
          std::cout << lnode+1 << " = " << nodes_in_elem(ielem,lnode) + 1 << " ";
        
        std::cout << " }"<< std::endl;
    }
    std::cout << std::endl;
  }
  

  //debug section; print stiffness matrix graph and per element map
  std::cout << " ------------SPARSE GRAPH MATRIX--------------"<<std::endl;
  for (int inode = 0; inode < num_nodes; inode++){
      std::cout << "row: " << inode + 1 << " { ";
    for (int istride = 0; istride < Graph_Matrix_Strides(inode); istride++){
        std::cout << istride + 1 << " = " << Repeat_Graph_Matrix(inode,istride) + 1 << " , " ;
    }
    std::cout << " }"<< std::endl;
  }

  std::cout << " ------------ELEMENT ASSEMBLY MAP--------------"<<std::endl;

  for (int ielem = 0; ielem < rnum_elem; ielem++){
    std::cout << "elem:  " << ielem+1 << std::endl;
    for (int lnode = 0; lnode < nodes_per_elem; lnode++){
        std::cout << "{ "<< std::endl;
        for (int jnode = 0; jnode < nodes_per_elem; jnode++){
          std::cout <<"(" << lnode+1 << "," << jnode+1 << "," << nodes_in_elem(ielem,lnode)+1 << ")"<< " = " << Global_Stiffness_Matrix_Assembly_Map(ielem,lnode, jnode) + 1 << " ";
        }
        std::cout << " }"<< std::endl;
    }
    std::cout << std::endl;
  }
  */
  
}

/* -----------------------------------------------------------------------------
   Initialize local views and global vectors needed to describe the design
-------------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::init_design(){
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

    std::ostream &out = std::cout;
    Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    //if(myrank==0)
    //*fos << "Global Element Densities:" << std::endl;
    //Global_Element_Densities->describe(*fos,Teuchos::VERB_EXTREME);
    //*fos << std::endl;
  }

  //create global vectors for mass and moment of inertia
  Global_Element_Masses = Teuchos::rcp(new MV(element_map, 1));
  Global_Element_Moments_x = Teuchos::rcp(new MV(element_map, 1));
  Global_Element_Moments_y = Teuchos::rcp(new MV(element_map, 1));
  Global_Element_Moments_z = Teuchos::rcp(new MV(element_map, 1));
  Global_Element_Moments_of_Inertia_xx = Teuchos::rcp(new MV(element_map, 1));
  Global_Element_Moments_of_Inertia_yy = Teuchos::rcp(new MV(element_map, 1));
  Global_Element_Moments_of_Inertia_zz = Teuchos::rcp(new MV(element_map, 1));
  Global_Element_Moments_of_Inertia_xy = Teuchos::rcp(new MV(element_map, 1));
  Global_Element_Moments_of_Inertia_xz = Teuchos::rcp(new MV(element_map, 1));
  Global_Element_Moments_of_Inertia_yz = Teuchos::rcp(new MV(element_map, 1));

}

/* ----------------------------------------------------------------------
   Assemble the Sparse Stiffness Matrix
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::assemble_matrix(){
  int num_dim = simparam->num_dim;
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  int nodes_per_element;
  int current_row_n_nodes_scanned;
  int local_dof_index, global_node_index, current_row, current_column;
  int max_stride = 0;
  
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> Local_Stiffness_Matrix(num_dim*max_nodes_per_element,num_dim*max_nodes_per_element);

  //initialize stiffness Matrix entries to 0
  //debug print
    //std::cout << "DOF GRAPH MATRIX ENTRIES ON TASK " << myrank << std::endl;
  for (int idof = 0; idof < num_dim*nlocal_nodes; idof++){
    for (int istride = 0; istride < Stiffness_Matrix_Strides(idof); istride++){
      Stiffness_Matrix(idof,istride) = 0;
      //debug print
      //std::cout << "{" <<istride + 1 << "," << DOF_Graph_Matrix(idof,istride) << "} ";
    }
    //debug print
    //std::cout << std::endl;
  }

  //assemble the global stiffness matrix
  if(num_dim==2)
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
    nodes_per_element = elem2D->num_nodes();
    //construct local stiffness matrix for this element
    local_matrix_multiply(ielem, Local_Stiffness_Matrix);
    //assign entries of this local matrix to the sparse global matrix storage;
    for (int inode = 0; inode < nodes_per_element; inode++){
      //see if this node is local
      global_node_index = nodes_in_elem(ielem,inode);
      if(!map->isNodeGlobalElement(global_node_index)) continue;
      //set dof row start index
      current_row = num_dim*map->getLocalElement(global_node_index);
      for(int jnode = 0; jnode < nodes_per_element; jnode++){
        
        current_column = num_dim*Global_Stiffness_Matrix_Assembly_Map(ielem,inode,jnode);
        for (int idim = 0; idim < num_dim; idim++){
          for (int jdim = 0; jdim < num_dim; jdim++){

            //debug print
            //if(current_row + idim==15&&current_column + jdim==4)
            //std::cout << " Local stiffness matrix contribution for row " << current_row + idim +1 << " and column " << current_column + jdim + 1 << " : " <<
            //Local_Stiffness_Matrix(num_dim*inode + idim,num_dim*jnode + jdim) << " from " << ielem +1 << " i: " << num_dim*inode+idim+1 << " j: " << num_dim*jnode + jdim +1 << std::endl << std::endl;
            //end debug

            Stiffness_Matrix(current_row + idim, current_column + jdim) += Local_Stiffness_Matrix(num_dim*inode + idim,num_dim*jnode + jdim);
          }
        }
      }
    }
  }

  if(num_dim==3)
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_3Delem_type(Element_Types(ielem), elem);
    nodes_per_element = elem->num_nodes();
    //construct local stiffness matrix for this element
    local_matrix_multiply(ielem, Local_Stiffness_Matrix);
    //assign entries of this local matrix to the sparse global matrix storage;
    for (int inode = 0; inode < nodes_per_element; inode++){
      //see if this node is local
      global_node_index = nodes_in_elem(ielem,inode);
      if(!map->isNodeGlobalElement(global_node_index)) continue;
      //set dof row start index
      current_row = num_dim*map->getLocalElement(global_node_index);
      for(int jnode = 0; jnode < nodes_per_element; jnode++){
        
        current_column = num_dim*Global_Stiffness_Matrix_Assembly_Map(ielem,inode,jnode);
        for (int idim = 0; idim < num_dim; idim++){
          for (int jdim = 0; jdim < num_dim; jdim++){

            //debug print
            //if(current_row + idim==15&&current_column + jdim==4)
            //std::cout << " Local stiffness matrix contribution for row " << current_row + idim +1 << " and column " << current_column + jdim + 1 << " : " <<
            //Local_Stiffness_Matrix(num_dim*inode + idim,num_dim*jnode + jdim) << " from " << ielem +1 << " i: " << num_dim*inode+idim+1 << " j: " << num_dim*jnode + jdim +1 << std::endl << std::endl;
            //end debug

            Stiffness_Matrix(current_row + idim, current_column + jdim) += Local_Stiffness_Matrix(num_dim*inode + idim,num_dim*jnode + jdim);
          }
        }
      }
    }
  }

  //construct distributed stiffness matrix and force vector from local kokkos data
  
  //build column map for the global stiffness matrix
  Teuchos::RCP<const Tpetra::Map<LO,GO,node_type> > colmap;
  const Teuchos::RCP<const Tpetra::Map<LO,GO,node_type> > dommap = local_dof_map;
  
  //debug print
  /*
    std::cout << "DOF GRAPH MATRIX ENTRIES ON TASK " << myrank << std::endl;
  for (int idof = 0; idof < num_dim*nlocal_nodes; idof++){
    for (int istride = 0; istride < Stiffness_Matrix_Strides(idof); istride++){
      //debug print
      std::cout << "{" <<istride + 1 << "," << DOF_Graph_Matrix(idof,istride) << "} ";
    }
    //debug print
    std::cout << std::endl;
  } */

  //debug print of stiffness matrix
  /*
  std::cout << " ------------SPARSE STIFFNESS MATRIX ON TASK"<< myrank << std::endl;
  for (int idof = 0; idof < num_dim*nlocal_nodes; idof++){
      std::cout << "row: " << idof + 1 << " { ";
    for (int istride = 0; istride < Stiffness_Matrix_Strides(idof); istride++){
        std::cout << istride + 1 << " = " << Stiffness_Matrix(idof,istride) << " , " ;
    }
    std::cout << " }"<< std::endl;
  }
  */
  Tpetra::Details::makeColMap<LO,GO,node_type>(colmap,dommap,DOF_Graph_Matrix.get_kokkos_view(), nullptr);

  size_t nnz = DOF_Graph_Matrix.size();

  //debug print
  //std::cout << "DOF GRAPH SIZE ON RANK " << myrank << " IS " << nnz << std::endl;
  
  //local indices in the graph using the constructed column map
  CArrayKokkos<LO, array_layout, device_type, memory_traits> stiffness_local_indices(nnz, "stiffness_local_indices");
  
  //row offsets with compatible template arguments
    row_pointers row_offsets = DOF_Graph_Matrix.start_index_;
    row_pointers row_offsets_pass("row_offsets", nlocal_nodes*num_dim+1);
    for(int ipass = 0; ipass < nlocal_nodes*num_dim + 1; ipass++){
      row_offsets_pass(ipass) = row_offsets(ipass);
    }

  size_t entrycount = 0;
  for(int irow = 0; irow < nlocal_nodes*num_dim; irow++){
    for(int istride = 0; istride < Stiffness_Matrix_Strides(irow); istride++){
      stiffness_local_indices(entrycount) = colmap->getLocalElement(DOF_Graph_Matrix(irow,istride));
      entrycount++;
    }
  }
  
  if(!Matrix_alloc){
  Global_Stiffness_Matrix = Teuchos::rcp(new MAT(local_dof_map, colmap, row_offsets_pass, stiffness_local_indices.get_kokkos_view(), Stiffness_Matrix.get_kokkos_view()));
  Global_Stiffness_Matrix->fillComplete();
  Matrix_alloc = 1;
  }

  //filter small negative numbers from floating point error
  for (int idof = 0; idof < num_dim*nlocal_nodes; idof++){
    for (int istride = 0; istride < Stiffness_Matrix_Strides(idof); istride++){
      if(Stiffness_Matrix(idof,istride)<0&&Stiffness_Matrix(idof,istride)>-0.00000000001)
      Stiffness_Matrix(idof,istride) = 0;
      //debug print
      //std::cout << "{" <<istride + 1 << "," << DOF_Graph_Matrix(idof,istride) << "} ";
    }
    //debug print
    //std::cout << std::endl;
  }
  //This completes the setup for A matrix of the linear system
  
  //file to debug print
  //std::ostream &out = std::cout;
  //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

  //debug print of A matrix
  //*fos << "Global Stiffness Matrix :" << std::endl;
  //Global_Stiffness_Matrix->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;

  //Print solution vector
  //*fos << "Global Nodal Forces :" << std::endl;
  //Global_Nodal_Forces->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
}

/* ----------------------------------------------------------------------
   Retrieve material properties associated with a finite element
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::Element_Material_Properties(size_t ielem, real_t &Element_Modulus, real_t &Poisson_Ratio, real_t density){
  real_t unit_scaling = simparam->unit_scaling;
  real_t penalty_product = 1;
  if(density < 0) density = 0;
  for(int i = 0; i < simparam->penalty_power; i++)
    penalty_product *= density;
  //relationship between density and stiffness
  Element_Modulus = (DENSITY_EPSILON + (1 - DENSITY_EPSILON)*penalty_product)*simparam->Elastic_Modulus/unit_scaling/unit_scaling;
  //Element_Modulus = density*simparam->Elastic_Modulus/unit_scaling/unit_scaling;
  Poisson_Ratio = simparam->Poisson_Ratio;

}

/* ----------------------------------------------------------------------
   Retrieve derivative of material properties with respect to local density
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::Gradient_Element_Material_Properties(size_t ielem, real_t &Element_Modulus_Derivative, real_t &Poisson_Ratio, real_t density){
  real_t unit_scaling = simparam->unit_scaling;
  real_t penalty_product = 1;
  if(density < 0) density = 0;
  for(int i = 0; i < simparam->penalty_power - 1; i++)
    penalty_product *= density;
  //relationship between density and stiffness
  Element_Modulus_Derivative = simparam->penalty_power*(1 - DENSITY_EPSILON)*penalty_product*simparam->Elastic_Modulus/unit_scaling/unit_scaling;
  //Element_Modulus_Derivative = simparam->Elastic_Modulus/unit_scaling/unit_scaling;
  Poisson_Ratio = simparam->Poisson_Ratio;

}

/* ----------------------------------------------------------------------
   Construct the local stiffness matrix
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::local_matrix(int ielem, CArrayKokkos<real_t, array_layout, device_type, memory_traits> &Local_Matrix){
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_vec_array Element_Densities;
  //local variable for host view of densities from the dual view
  bool nodal_density_flag = simparam->nodal_density_flag;
  const_host_vec_array all_node_densities;
  if(nodal_density_flag)
  all_node_densities = all_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  else
  Element_Densities = Global_Element_Densities->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
  int num_dim = simparam->num_dim;
  int nodes_per_elem = max_nodes_per_element;
  int num_gauss_points = simparam->num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  size_t local_node_id;

  direct_product_count = std::pow(num_gauss_points,num_dim);
  real_t Elastic_Constant, Shear_Term, Pressure_Term, matrix_term;
  real_t matrix_subterm1, matrix_subterm2, matrix_subterm3, Jacobian;
  real_t Element_Modulus, Poisson_Ratio;
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
  CArray<real_t> legendre_nodes_1D(num_gauss_points);
  CArray<real_t> legendre_weights_1D(num_gauss_points);
  real_t pointer_quad_coordinate[num_dim];
  real_t pointer_quad_coordinate_weight[num_dim];
  real_t pointer_interpolated_point[num_dim];
  real_t pointer_JT_row1[num_dim];
  real_t pointer_JT_row2[num_dim];
  real_t pointer_JT_row3[num_dim];
  ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate,num_dim);
  ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight,num_dim);
  ViewCArray<real_t> interpolated_point(pointer_interpolated_point,num_dim);
  ViewCArray<real_t> JT_row1(pointer_JT_row1,num_dim);
  ViewCArray<real_t> JT_row2(pointer_JT_row2,num_dim);
  ViewCArray<real_t> JT_row3(pointer_JT_row3,num_dim);

  real_t pointer_basis_values[elem->num_basis()];
  real_t pointer_basis_derivative_s1[elem->num_basis()];
  real_t pointer_basis_derivative_s2[elem->num_basis()];
  real_t pointer_basis_derivative_s3[elem->num_basis()];
  ViewCArray<real_t> basis_values(pointer_basis_values,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(),num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());

  //initialize weights
  elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
  elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);

  real_t current_density = 1;

  //acquire set of nodes for this local element
  for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
    local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, node_loop));
    nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
    nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
    nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
    if(nodal_density_flag) nodal_density(node_loop) = all_node_densities(local_node_id,0);
  }

  //initialize local stiffness matrix storage
  for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++)
      for(int jfill=0; jfill < num_dim*nodes_per_elem; jfill++)
      Local_Matrix(ifill,jfill) = 0;

  //loop over quadrature points
  for(int iquad=0; iquad < direct_product_count; iquad++){

    //set current quadrature point
    if(num_dim==3) z_quad = iquad/(num_gauss_points*num_gauss_points);
    y_quad = (iquad % (num_gauss_points*num_gauss_points))/num_gauss_points;
    x_quad = iquad % num_gauss_points;
    quad_coordinate(0) = legendre_nodes_1D(x_quad);
    quad_coordinate(1) = legendre_nodes_1D(y_quad);
    if(num_dim==3)
    quad_coordinate(2) = legendre_nodes_1D(z_quad);

    //set current quadrature weight
    quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
    quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
    if(num_dim==3)
    quad_coordinate_weight(2) = legendre_weights_1D(z_quad);

    //compute shape functions at this point for the element type
    elem->basis(basis_values,quad_coordinate);

    //compute density
    current_density = 0;
    if(nodal_density_flag)
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      current_density += nodal_density(node_loop)*basis_values(node_loop);
    }
    //default constant element density
    else{
      current_density = Element_Densities(ielem,0);
    }

    //look up element material properties as a function of density at the point
    Element_Material_Properties((size_t) ielem,Element_Modulus,Poisson_Ratio, current_density);
    Elastic_Constant = Element_Modulus/((1 + Poisson_Ratio)*(1 - 2*Poisson_Ratio));
    Shear_Term = 0.5-Poisson_Ratio;
    Pressure_Term = 1 - Poisson_Ratio;

    //compute all the necessary coordinates and derivatives at this point

    //compute shape function derivatives
    elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
    elem->partial_eta_basis(basis_derivative_s2,quad_coordinate);
    elem->partial_mu_basis(basis_derivative_s3,quad_coordinate);

    //compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
    //derivative of x,y,z w.r.t s
    JT_row1(0) = 0;
    JT_row1(1) = 0;
    JT_row1(2) = 0;
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      JT_row1(0) += nodal_positions(node_loop,0)*basis_derivative_s1(node_loop);
      JT_row1(1) += nodal_positions(node_loop,1)*basis_derivative_s1(node_loop);
      JT_row1(2) += nodal_positions(node_loop,2)*basis_derivative_s1(node_loop);
    }

    //derivative of x,y,z w.r.t t
    JT_row2(0) = 0;
    JT_row2(1) = 0;
    JT_row2(2) = 0;
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      JT_row2(0) += nodal_positions(node_loop,0)*basis_derivative_s2(node_loop);
      JT_row2(1) += nodal_positions(node_loop,1)*basis_derivative_s2(node_loop);
      JT_row2(2) += nodal_positions(node_loop,2)*basis_derivative_s2(node_loop);
    }

    //derivative of x,y,z w.r.t w
    JT_row3(0) = 0;
    JT_row3(1) = 0;
    JT_row3(2) = 0;
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      JT_row3(0) += nodal_positions(node_loop,0)*basis_derivative_s3(node_loop);
      JT_row3(1) += nodal_positions(node_loop,1)*basis_derivative_s3(node_loop);
      JT_row3(2) += nodal_positions(node_loop,2)*basis_derivative_s3(node_loop);
    }

    //compute the determinant of the Jacobian
    Jacobian = JT_row1(0)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
               JT_row1(1)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
               JT_row1(2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1));
    if(Jacobian<0) Jacobian = -Jacobian;

    //compute the contributions of this quadrature point to all the matrix elements
    int index_x,index_y,basis_index_x,basis_index_y,swap1,swap2;
    for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++)
      for(int jfill=0; jfill < num_dim*nodes_per_elem; jfill++){
        index_x = ifill%num_dim;
        index_y = jfill%num_dim;
        basis_index_x = ifill/num_dim;
        basis_index_y = jfill/num_dim;

        //compute stiffness matrix terms involving derivatives of the shape function and cofactor determinants from cramers rule
        if(index_x==0&&index_y==0){
          matrix_subterm1 = Pressure_Term*(basis_derivative_s1(basis_index_x)*
          (JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(basis_index_x)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(basis_index_x)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)))*
          (basis_derivative_s1(basis_index_y)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(basis_index_y)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(basis_index_y)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));

          matrix_subterm2 = Shear_Term*(-basis_derivative_s1(basis_index_x)*
          (JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(basis_index_x)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(basis_index_x)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)))*
          (-basis_derivative_s1(basis_index_y)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(basis_index_y)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(basis_index_y)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));

          matrix_subterm3 = Shear_Term*(basis_derivative_s1(basis_index_x)*
          (JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(basis_index_x)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(basis_index_x)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)))*
          (basis_derivative_s1(basis_index_y)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(basis_index_y)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(basis_index_y)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));

          matrix_term = matrix_subterm1 + matrix_subterm2 + matrix_subterm3;

           //debug print block
           /*
          if(ielem==0&&((jfill==3&&ifill==3))){
           std::cout << " ------------quadrature point "<< iquad + 1 <<"--------------"<<std::endl;
           std::cout << "row: " << ifill + 1 << " { ";
           std::cout << jfill + 1 << " = " << matrix_subterm1 << " , " << " bi:JT11 "<<JT_row1(0) << " bi:JT12 " <<  JT_row1(1) << " bi:JT13 " << JT_row1(2)
           <<  " bi:JT21 "<<JT_row2(0) << " bi:JT22 " <<  JT_row2(1) << " bi:JT23 " << JT_row2(2) <<  " bi:JT31 "<<JT_row3(0) << " bi:JT32 " <<  JT_row3(1) << " bi:JT33 " << JT_row3(2);
           std::cout << " }"<< std::endl;
          }
          
          if(ielem==0&&((jfill==3&&ifill==3))){
           std::cout << " ------------quadrature point "<< iquad + 1 <<"--------------"<<std::endl;
           std::cout << "row: " << ifill + 1 << " { ";
           std::cout << jfill + 1 << " = " << matrix_term*Elastic_Constant/Jacobian << " , " << " basis index x s1 "<< basis_derivative_s1(basis_index_x) << " quad x " <<  quad_coordinate(0)
           <<  " quad y "<< quad_coordinate(1) << " quad z " <<  quad_coordinate(2) << "Force Vector " << Local_Matrix(3,3);
           std::cout << " }"<< std::endl;
          }
          */
          
        }

        if(index_x==1&&index_y==1){
          
          matrix_subterm1 = Pressure_Term*(-basis_derivative_s1(basis_index_x)*
          (JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(basis_index_x)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(basis_index_x)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)))*
          (-basis_derivative_s1(basis_index_y)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(basis_index_y)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(basis_index_y)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));

          matrix_subterm2 = Shear_Term*(basis_derivative_s1(basis_index_x)*
          (JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(basis_index_x)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(basis_index_x)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)))*
          (basis_derivative_s1(basis_index_y)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(basis_index_y)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(basis_index_y)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));

          matrix_subterm3 = Shear_Term*(basis_derivative_s1(basis_index_x)*
          (JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(basis_index_x)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(basis_index_x)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)))*
          (basis_derivative_s1(basis_index_y)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(basis_index_y)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(basis_index_y)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));

          matrix_term = matrix_subterm1 + matrix_subterm2 + matrix_subterm3;
        }

        if(index_x==2&&index_y==2){
          
          matrix_subterm1 = Pressure_Term*(basis_derivative_s1(basis_index_x)*
          (JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(basis_index_x)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(basis_index_x)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)))*
          (basis_derivative_s1(basis_index_y)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(basis_index_y)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(basis_index_y)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));

          matrix_subterm2 = Shear_Term*(basis_derivative_s1(basis_index_x)*
          (JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(basis_index_x)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(basis_index_x)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)))*
          (basis_derivative_s1(basis_index_y)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(basis_index_y)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(basis_index_y)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));

          matrix_subterm3 = Shear_Term*(-basis_derivative_s1(basis_index_x)*
          (JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(basis_index_x)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(basis_index_x)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)))*
          (-basis_derivative_s1(basis_index_y)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(basis_index_y)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(basis_index_y)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));

          matrix_term = matrix_subterm1 + matrix_subterm2 + matrix_subterm3;
        }

        if((index_x==0&&index_y==1)||(index_x==1&&index_y==0)){
          if(index_x==1&&index_y==0){
            swap1 = basis_index_x;
            swap2 = basis_index_y;
          }
          else{
            swap1 = basis_index_y;
            swap2 = basis_index_x;
          }
          
          matrix_subterm1 = Poisson_Ratio*(-basis_derivative_s1(swap1)*
          (JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(swap1)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(swap1)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)))*
          (basis_derivative_s1(swap2)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(swap2)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(swap2)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
          
          
          matrix_subterm2 = Shear_Term*(basis_derivative_s1(swap1)*
          (JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(swap1)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(swap1)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)))*
          (-basis_derivative_s1(swap2)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(swap2)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(swap2)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));

          matrix_term = matrix_subterm1 + matrix_subterm2;

          /* debug print block
          if(iquad==0&&((jfill==4&&ifill==0)||(jfill==0&&ifill==4))){
           std::cout << " ------------LOCAL STIFFNESS MATRIX "<< ielem + 1 <<"--------------"<<std::endl;
           std::cout << "row: " << ifill + 1 << " { ";
           std::cout << jfill + 1 << " = " << matrix_subterm2 << " , " << " bi:JT11 "<<JT_row1(0) << " bi:JT12 " <<  JT_row1(1) << " bi:JT13 " << JT_row1(2)
           <<  " bi:JT21 "<<JT_row2(0) << " bi:JT22 " <<  JT_row2(1) << " bi:JT23 " << JT_row2(2) <<  " bi:JT31 "<<JT_row3(0) << " bi:JT32 " <<  JT_row3(1) << " bi:JT33 " << JT_row3(2);
           std::cout << " }"<< std::endl;
          }
          */
        }

        if((index_x==0&&index_y==2)||(index_x==2&&index_y==0)){
          if(index_x==2&index_y==0){
            swap1 = basis_index_x;
            swap2 = basis_index_y;
          }
          else{
            swap1 = basis_index_y;
            swap2 = basis_index_x;
          }
          matrix_subterm1 = Poisson_Ratio*(basis_derivative_s1(swap1)*
          (JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(swap1)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(swap1)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)))*
          (basis_derivative_s1(swap2)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(swap2)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(swap2)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));

          matrix_subterm2 = Shear_Term*(basis_derivative_s1(swap1)*
          (JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(swap1)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(swap1)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)))*
          (basis_derivative_s1(swap2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(swap2)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(swap2)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));

          matrix_term = matrix_subterm1 + matrix_subterm2;
        }

        if((index_x==1&&index_y==2)||(index_x==2&&index_y==1)){
          if(index_x==2&&index_y==1){
            swap1 = basis_index_x;
            swap2 = basis_index_y;
          }
          else{
            swap1 = basis_index_y;
            swap2 = basis_index_x;
          }
          matrix_subterm1 = Poisson_Ratio*(basis_derivative_s1(swap1)*
          (JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(swap1)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(swap1)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)))*
          (-basis_derivative_s1(swap2)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(swap2)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(swap2)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));

          matrix_subterm2 = Shear_Term*(-basis_derivative_s1(swap1)*
          (JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(swap1)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(swap1)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)))*
          (basis_derivative_s1(swap2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(swap2)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(swap2)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));

          matrix_term = matrix_subterm1 + matrix_subterm2;
        }
        
        Local_Matrix(ifill,jfill) += Elastic_Constant*quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*matrix_term/Jacobian;
      }
      
    }

    //debug print of local stiffness matrix
      /*
      std::cout << " ------------LOCAL STIFFNESS MATRIX "<< ielem + 1 <<"--------------"<<std::endl;
      for (int idof = 0; idof < num_dim*nodes_per_elem; idof++){
        std::cout << "row: " << idof + 1 << " { ";
        for (int istride = 0; istride < num_dim*nodes_per_elem; istride++){
          std::cout << istride + 1 << " = " << Local_Matrix(idof,istride) << " , " ;
        }
        std::cout << " }"<< std::endl;
        }
      */
}

/* ----------------------------------------------------------------------
   Construct the local stiffness matrix
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::local_matrix_multiply(int ielem, CArrayKokkos<real_t, array_layout, device_type, memory_traits> &Local_Matrix){
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_vec_array Element_Densities;
  //local variable for host view of densities from the dual view
  bool nodal_density_flag = simparam->nodal_density_flag;
  const_host_vec_array all_node_densities;
  if(nodal_density_flag)
  all_node_densities = all_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  else
  Element_Densities = Global_Element_Densities->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
  int num_dim = simparam->num_dim;
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam->num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  size_t local_node_id;

  direct_product_count = std::pow(num_gauss_points,num_dim);
  real_t Elastic_Constant, Shear_Term, Pressure_Term, matrix_term;
  real_t matrix_subterm1, matrix_subterm2, matrix_subterm3, Jacobian;
  real_t Element_Modulus, Poisson_Ratio;
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
  CArray<real_t> legendre_nodes_1D(num_gauss_points);
  CArray<real_t> legendre_weights_1D(num_gauss_points);
  real_t pointer_quad_coordinate[num_dim];
  real_t pointer_quad_coordinate_weight[num_dim];
  real_t pointer_interpolated_point[num_dim];
  real_t pointer_JT_row1[num_dim];
  real_t pointer_JT_row2[num_dim];
  real_t pointer_JT_row3[num_dim];
  ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate,num_dim);
  ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight,num_dim);
  ViewCArray<real_t> interpolated_point(pointer_interpolated_point,num_dim);
  ViewCArray<real_t> JT_row1(pointer_JT_row1,num_dim);
  ViewCArray<real_t> JT_row2(pointer_JT_row2,num_dim);
  ViewCArray<real_t> JT_row3(pointer_JT_row3,num_dim);

  real_t pointer_basis_values[elem->num_basis()];
  real_t pointer_basis_derivative_s1[elem->num_basis()];
  real_t pointer_basis_derivative_s2[elem->num_basis()];
  real_t pointer_basis_derivative_s3[elem->num_basis()];
  ViewCArray<real_t> basis_values(pointer_basis_values,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(),num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());

  size_t Brows;
  if(num_dim==2) Brows = 3;
  if(num_dim==3) Brows = 6;
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix_contribution(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix_contribution(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> C_matrix(Brows,Brows);

  //initialize weights
  elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
  elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);

  real_t current_density = 1;

  //acquire set of nodes for this local element
  for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
    local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, node_loop));
    nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
    nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
    nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
    if(nodal_density_flag) nodal_density(node_loop) = all_node_densities(local_node_id,0);
    /*
    if(myrank==1&&nodal_positions(node_loop,2)>10000000){
      std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
      std::cout << node_loop+1 <<" " << local_node_id <<" "<< nodes_in_elem(ielem, node_loop) << " "<< nodal_positions(node_loop,2) << std::endl;
      std::fflush(stdout);
    }
    */
    //std::cout << local_node_id << " " << nodes_in_elem(ielem, node_loop) << " " << nodal_positions(node_loop,0) << " " << nodal_positions(node_loop,1) << " "<< nodal_positions(node_loop,2) <<std::endl;
  }

  //initialize C matrix
  for(int irow = 0; irow < Brows; irow++)
    for(int icol = 0; icol < Brows; icol++)
      C_matrix(irow,icol) = 0;

  //initialize local stiffness matrix storage
  for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++)
      for(int jfill=0; jfill < num_dim*nodes_per_elem; jfill++)
      Local_Matrix(ifill,jfill) = 0;

  //B matrix initialization
  for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < num_dim*nodes_per_elem; icol++){
        CB_matrix(irow,icol) = B_matrix(irow,icol) = 0;
      }

  //loop over quadrature points
  for(int iquad=0; iquad < direct_product_count; iquad++){

    //set current quadrature point
    if(num_dim==3) z_quad = iquad/(num_gauss_points*num_gauss_points);
    y_quad = (iquad % (num_gauss_points*num_gauss_points))/num_gauss_points;
    x_quad = iquad % num_gauss_points;
    quad_coordinate(0) = legendre_nodes_1D(x_quad);
    quad_coordinate(1) = legendre_nodes_1D(y_quad);
    if(num_dim==3)
    quad_coordinate(2) = legendre_nodes_1D(z_quad);

    //set current quadrature weight
    quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
    quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
    if(num_dim==3)
    quad_coordinate_weight(2) = legendre_weights_1D(z_quad);

    //compute shape functions at this point for the element type
    elem->basis(basis_values,quad_coordinate);
    
    //compute density
    current_density = 0;
    if(nodal_density_flag)
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      current_density += nodal_density(node_loop)*basis_values(node_loop);
    }
    //default constant element density
    else{
      current_density = Element_Densities(ielem,0);
    }

    //debug print
    //std::cout << "Current Density " << current_density << std::endl;

    //look up element material properties at this point as a function of density
    Element_Material_Properties((size_t) ielem,Element_Modulus,Poisson_Ratio, current_density);
    Elastic_Constant = Element_Modulus/((1 + Poisson_Ratio)*(1 - 2*Poisson_Ratio));
    Shear_Term = 0.5-Poisson_Ratio;
    Pressure_Term = 1 - Poisson_Ratio;

    //debug print
    //std::cout << "Element Material Params " << Elastic_Constant << std::endl;

    //compute Elastic (C) matrix
    if(num_dim==2){
      C_matrix(0,0) = Pressure_Term;
      C_matrix(1,1) = Pressure_Term;
      C_matrix(0,1) = Poisson_Ratio;
      C_matrix(1,0) = Poisson_Ratio;
      C_matrix(2,2) = Shear_Term;
    }
    if(num_dim==3){
      C_matrix(0,0) = Pressure_Term;
      C_matrix(1,1) = Pressure_Term;
      C_matrix(2,2) = Pressure_Term;
      C_matrix(0,1) = Poisson_Ratio;
      C_matrix(0,2) = Poisson_Ratio;
      C_matrix(1,0) = Poisson_Ratio;
      C_matrix(1,2) = Poisson_Ratio;
      C_matrix(2,0) = Poisson_Ratio;
      C_matrix(2,1) = Poisson_Ratio;
      C_matrix(3,3) = Shear_Term;
      C_matrix(4,4) = Shear_Term;
      C_matrix(5,5) = Shear_Term;
    }
  
  /*
  //debug print of elasticity matrix
  std::cout << " ------------ELASTICITY MATRIX "<< ielem + 1 <<"--------------"<<std::endl;
  for (int idof = 0; idof < Brows; idof++){
    std::cout << "row: " << idof + 1 << " { ";
    for (int istride = 0; istride < Brows; istride++){
      std::cout << istride + 1 << " = " << C_matrix(idof,istride) << " , " ;
    }
    std::cout << " }"<< std::endl;
  }
  //end debug block
  */

    //compute all the necessary coordinates and derivatives at this point
    //compute shape function derivatives
    elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
    elem->partial_eta_basis(basis_derivative_s2,quad_coordinate);
    elem->partial_mu_basis(basis_derivative_s3,quad_coordinate);

    //compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
    //derivative of x,y,z w.r.t s
    JT_row1(0) = 0;
    JT_row1(1) = 0;
    JT_row1(2) = 0;
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      JT_row1(0) += nodal_positions(node_loop,0)*basis_derivative_s1(node_loop);
      JT_row1(1) += nodal_positions(node_loop,1)*basis_derivative_s1(node_loop);
      JT_row1(2) += nodal_positions(node_loop,2)*basis_derivative_s1(node_loop);
    }

    //derivative of x,y,z w.r.t t
    JT_row2(0) = 0;
    JT_row2(1) = 0;
    JT_row2(2) = 0;
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      JT_row2(0) += nodal_positions(node_loop,0)*basis_derivative_s2(node_loop);
      JT_row2(1) += nodal_positions(node_loop,1)*basis_derivative_s2(node_loop);
      JT_row2(2) += nodal_positions(node_loop,2)*basis_derivative_s2(node_loop);
    }

    //derivative of x,y,z w.r.t w
    JT_row3(0) = 0;
    JT_row3(1) = 0;
    JT_row3(2) = 0;
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      JT_row3(0) += nodal_positions(node_loop,0)*basis_derivative_s3(node_loop);
      JT_row3(1) += nodal_positions(node_loop,1)*basis_derivative_s3(node_loop);
      JT_row3(2) += nodal_positions(node_loop,2)*basis_derivative_s3(node_loop);
      //debug print
    /*if(myrank==1&&nodal_positions(node_loop,2)*basis_derivative_s3(node_loop)<-10000000){
      std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
      std::cout << node_loop+1 << " " << JT_row3(2) << " "<< nodal_positions(node_loop,2) <<" "<< basis_derivative_s3(node_loop) << std::endl;
      std::fflush(stdout);
    }*/
    }
    
    
    //compute the determinant of the Jacobian
    Jacobian = JT_row1(0)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
               JT_row1(1)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
               JT_row1(2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1));
    if(Jacobian<0) Jacobian = -Jacobian;

    //compute the contributions of this quadrature point to the B matrix
    if(num_dim==2)
    for(int ishape=0; ishape < nodes_per_elem; ishape++){
      B_matrix_contribution(0,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(1,ishape*num_dim) = 0;
      B_matrix_contribution(2,ishape*num_dim) = 0;
      B_matrix_contribution(3,ishape*num_dim) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(5,ishape*num_dim) = 0;
      B_matrix_contribution(0,ishape*num_dim+1) = 0;
      B_matrix_contribution(1,ishape*num_dim+1) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(2,ishape*num_dim+1) = 0;
      B_matrix_contribution(3,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim+1) = 0;
      B_matrix_contribution(5,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(0,ishape*num_dim+2) = 0;
      B_matrix_contribution(1,ishape*num_dim+2) = 0;
      B_matrix_contribution(2,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(3,ishape*num_dim+2) = 0;
      B_matrix_contribution(4,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(5,ishape*num_dim+2) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
    }
    if(num_dim==3)
    for(int ishape=0; ishape < nodes_per_elem; ishape++){
      B_matrix_contribution(0,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(1,ishape*num_dim) = 0;
      B_matrix_contribution(2,ishape*num_dim) = 0;
      B_matrix_contribution(3,ishape*num_dim) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(5,ishape*num_dim) = 0;
      B_matrix_contribution(0,ishape*num_dim+1) = 0;
      B_matrix_contribution(1,ishape*num_dim+1) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(2,ishape*num_dim+1) = 0;
      B_matrix_contribution(3,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim+1) = 0;
      B_matrix_contribution(5,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(0,ishape*num_dim+2) = 0;
      B_matrix_contribution(1,ishape*num_dim+2) = 0;
      B_matrix_contribution(2,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(3,ishape*num_dim+2) = 0;
      B_matrix_contribution(4,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(5,ishape*num_dim+2) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
    }
    /*
    //debug print of B matrix per quadrature point
    std::cout << " ------------B MATRIX QUADRATURE CONTRIBUTION"<< ielem + 1 <<"--------------"<<std::endl;
    for (int idof = 0; idof < Brows; idof++){
      std::cout << "row: " << idof + 1 << " { ";
      for (int istride = 0; istride < nodes_per_elem*num_dim; istride++){
        std::cout << istride + 1 << " = " << B_matrix_contribution(idof,istride) << " , " ;
      }
      std::cout << " }"<< std::endl;
    }
    //end debug block
    */
    //accumulate B matrix
    for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < num_dim*nodes_per_elem; icol++)
      B_matrix(irow,icol) += B_matrix_contribution(irow,icol);

    //compute the previous multiplied by the Elastic (C) Matrix
    for(int irow=0; irow < Brows; irow++){
      for(int icol=0; icol < num_dim*nodes_per_elem; icol++){
        CB_matrix_contribution(irow,icol) = 0;
        for(int span=0; span < Brows; span++){
          CB_matrix_contribution(irow,icol) += C_matrix(irow,span)*B_matrix_contribution(span,icol);
        }
      }
    }

    //accumulate CB matrix
    for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < num_dim*nodes_per_elem; icol++)
      CB_matrix(irow,icol) += CB_matrix_contribution(irow,icol);

    //compute the contributions of this quadrature point to all the local stiffness matrix elements
    for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++)
      for(int jfill=0; jfill < num_dim*nodes_per_elem; jfill++){
        matrix_term = 0;
        for(int span = 0; span < Brows; span++){
          matrix_term += B_matrix_contribution(span,ifill)*CB_matrix_contribution(span,jfill);
        }
        Local_Matrix(ifill,jfill) += Elastic_Constant*quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*matrix_term/Jacobian;
      }
    
    }
    
    /*
    //debug print of B matrix per quadrature point
    std::cout << " ------------B MATRIX"<< ielem + 1 <<"--------------"<<std::endl;
    for (int idof = 0; idof < Brows; idof++){
      std::cout << "row: " << idof + 1 << " { ";
      for (int istride = 0; istride < nodes_per_elem*num_dim; istride++){
        std::cout << istride + 1 << " = " << B_matrix(idof,istride) << " , " ;
      }
      std::cout << " }"<< std::endl;
    }
    //end debug block

    //debug print of B matrix per quadrature point
    std::cout << " ------------CB MATRIX"<< ielem + 1 <<"--------------"<<std::endl;
    for (int idof = 0; idof < Brows; idof++){
      std::cout << "row: " << idof + 1 << " { ";
      for (int istride = 0; istride < nodes_per_elem*num_dim; istride++){
        std::cout << istride + 1 << " = " << CB_matrix(idof,istride) << " , " ;
      }
      std::cout << " }"<< std::endl;
    }
    //end debug block
    */

    //debug print of local stiffness matrix
      /*
      std::cout << " ------------LOCAL STIFFNESS MATRIX "<< ielem + 1 <<"--------------"<<std::endl;
      for (int idof = 0; idof < num_dim*nodes_per_elem; idof++){
        std::cout << "row: " << idof + 1 << " { ";
        for (int istride = 0; istride < num_dim*nodes_per_elem; istride++){
          std::cout << istride + 1 << " = " << Local_Matrix(idof,istride) << " , " ;
        }
        std::cout << " }"<< std::endl;
        }
      */
}

/* ----------------------------------------------------------------------
   Loop through applied boundary conditions and tag node ids to redecule 
   necessary rows and columns from the assembled linear system
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::Displacement_Boundary_Conditions(){
  int num_bdy_patches_in_set, patch_id;
  int warning_flag = 0;
  int local_flag;
  int current_node_index, current_node_id;
  int num_boundary_sets = num_boundary_conditions;
  int surface_disp_set_id = 0;
  int num_dim = simparam->num_dim;
  int bc_option, bc_dim_set[3];
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> displacement(num_dim);
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Displacement_Conditions(num_dim);
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> first_condition_per_node(nall_nodes*num_dim);
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Surface_Nodes;
  Number_DOF_BCS = 0;
  Displacement_Conditions(0) = X_DISPLACEMENT_CONDITION;
  Displacement_Conditions(1) = Y_DISPLACEMENT_CONDITION;
  Displacement_Conditions(2) = Z_DISPLACEMENT_CONDITION;

  //host view of local nodal displacements
  host_vec_array node_displacements_host = node_displacements_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);

  //initialize to -1 (DO NOT make -1 an index for bdy sets)
  for(int inode = 0 ; inode < nlocal_nodes*num_dim; inode++)
    first_condition_per_node(inode) = -1;
  
  //scan for surface method of setting fixed nodal displacements
  for(int iboundary = 0; iboundary < num_boundary_sets; iboundary++){
    
    if(Boundary_Condition_Type_List(iboundary)==DISPLACEMENT_CONDITION){bc_option=3;}
    else if(Boundary_Condition_Type_List(iboundary)==X_DISPLACEMENT_CONDITION){bc_option=0;}
    else if(Boundary_Condition_Type_List(iboundary)==Y_DISPLACEMENT_CONDITION){bc_option=1;}
    else if(Boundary_Condition_Type_List(iboundary)==Z_DISPLACEMENT_CONDITION){bc_option=2;}
    else{
      continue;
    }
      
      //debug print of surface conditions
      //std::cout << " Surface BC types " << Boundary_Condition_Type_List(iboundary) <<std::endl;

      num_bdy_patches_in_set = NBoundary_Condition_Patches(iboundary);
      if(bc_option==0) {
        bc_dim_set[0]=1;
        displacement(0) = Boundary_Surface_Displacements(surface_disp_set_id,0);
      }
      else if(bc_option==1) {
        bc_dim_set[1]=1;
        displacement(1) = Boundary_Surface_Displacements(surface_disp_set_id,1);
      }
      else if(bc_option==2) {
        bc_dim_set[2]=1;
        displacement(2) = Boundary_Surface_Displacements(surface_disp_set_id,2);
      }
      else if(bc_option==3) {
        bc_dim_set[0]=1;
        bc_dim_set[1]=1;
        bc_dim_set[2]=1;
        displacement(0) = Boundary_Surface_Displacements(surface_disp_set_id,0);
        displacement(1) = Boundary_Surface_Displacements(surface_disp_set_id,1);
        displacement(2) = Boundary_Surface_Displacements(surface_disp_set_id,2);
      }
      surface_disp_set_id++;
      
      //loop over boundary set patches for their respective nodes
      for (int bdy_patch_gid = 0; bdy_patch_gid < num_bdy_patches_in_set; bdy_patch_gid++){
        //get number
        patch_id = Boundary_Condition_Patches(iboundary, bdy_patch_gid);
        Surface_Nodes = Boundary_Patches(patch_id).node_set;
        for(int inode = 0; inode < Surface_Nodes.size(); inode++){
          local_flag = 0;
          current_node_id = Surface_Nodes(inode);
          if(map->isNodeGlobalElement(current_node_id)) local_flag = 1;
          current_node_id = all_node_map->getLocalElement(current_node_id);
          
          /*
          //debug print of nodal bc settings
          std::cout << " Node BC types " << Node_DOF_Boundary_Condition_Type(current_node_id);
          std::cout << " node: " << inode << " { ";
          for (int istride = 0; istride < num_dim; istride++){
            std::cout << node_coords(current_node_id,istride) << " , ";
          }
          std::cout << " }"<< std::endl;
          */

          for(int idim=0; idim < num_dim; idim++){
          //warning for reapplied a displacement boundary condition (For now there is an output flag being set that triggers output later)
          if(Node_DOF_Boundary_Condition_Type(current_node_id*num_dim + idim)==DISPLACEMENT_CONDITION||
          Node_DOF_Boundary_Condition_Type(current_node_id*num_dim + idim)==Displacement_Conditions(idim)){
            //if overlap is just due to the loop over patches, a warning is not needed
            if(first_condition_per_node(current_node_id*num_dim + idim)!=iboundary) warning_flag = 1;
          }
          else{
            if(bc_dim_set[idim]){
              first_condition_per_node(current_node_id*num_dim + idim) = iboundary;
              Node_DOF_Boundary_Condition_Type(current_node_id*num_dim+idim) = Boundary_Condition_Type_List(iboundary);
              Node_DOF_Displacement_Boundary_Conditions(current_node_id*num_dim+idim) = displacement(idim);
              //counts local DOF being constrained
              if(local_flag){
              Number_DOF_BCS++;
              node_displacements_host(current_node_id*num_dim+idim,0) = displacement(idim);
              }
            }
          }
          }
        }
      }
  }

  //scan for direct setting of nodal displacements from input
  //indices for nodal BC settings referred to here start at num_boundary_sets

  //debug print of nodal bc settings
  /*
  std::cout << " ------------NODE BC SETTINGS--------------"<<std::endl;
  for(int inode=0; inode < num_nodes*num_dim; inode++)
  std::cout << " Node BC types " << Node_DOF_Boundary_Condition_Type(inode) <<std::endl;
  //end debug block
  */

  //print warning for overlapping boundary conditions
  if(warning_flag)
  std::cout << std::endl << "One or more displacement boundary conditions overlap on a subset of nodes; please revise input" << std::endl << std::endl;

}

/* ----------------------------------------------------------------------
   Construct the global applied force vector
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::assemble_vector(){
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  //local variable for host view in the dual view
  host_vec_array Nodal_Forces = Global_Nodal_Forces->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  int num_bdy_patches_in_set;
  size_t node_id, patch_id;
  GO current_node_index;
  LO local_node_index;
  int num_boundary_sets = num_boundary_conditions;
  int surface_force_set_id = 0;
  int num_dim = simparam->num_dim;
  int nodes_per_elem = max_nodes_per_element;
  int num_gauss_points = simparam->num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  int current_element_index, local_surface_id, surf_dim1, surf_dim2;
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
  CArray<real_t> legendre_nodes_1D(num_gauss_points);
  CArray<real_t> legendre_weights_1D(num_gauss_points);
  real_t pointer_quad_coordinate[num_dim];
  real_t pointer_quad_coordinate_weight[num_dim];
  real_t pointer_interpolated_point[num_dim];
  ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate,num_dim);
  ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight,num_dim);
  ViewCArray<real_t> interpolated_point(pointer_interpolated_point,num_dim);
  real_t force_density[3], wedge_product;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Surface_Nodes;
  
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> JT_row1(num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> JT_row2(num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> JT_row3(num_dim);

   //force vector initialization
  for(int i=0; i < num_dim*nlocal_nodes; i++)
    Nodal_Forces(i,0) = 0;

  /*Loop through boundary sets and check if they apply surface forces.
  These sets can have overlapping nodes since applied loading conditions
  are assumed to be additive*/
  for(int iboundary = 0; iboundary < num_boundary_sets; iboundary++){
    if(Boundary_Condition_Type_List(iboundary)!=LOADING_CONDITION) continue;
    //std::cout << "I REACHED THE LOADING BOUNDARY CONDITION" <<std::endl;
    num_bdy_patches_in_set = NBoundary_Condition_Patches(iboundary);
    
    force_density[0] = Boundary_Surface_Force_Densities(surface_force_set_id,0);
    force_density[1] = Boundary_Surface_Force_Densities(surface_force_set_id,1);
    force_density[2] = Boundary_Surface_Force_Densities(surface_force_set_id,2);
    surface_force_set_id++;

    real_t pointer_basis_values[elem->num_basis()];
    ViewCArray<real_t> basis_values(pointer_basis_values,elem->num_basis());

    //initialize weights
    elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
    elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);

    direct_product_count = std::pow(num_gauss_points,num_dim-1);
    //loop over boundary sets for their applied forces; use quadrature for distributed forces
    for (int bdy_patch_gid = 0; bdy_patch_gid < num_bdy_patches_in_set; bdy_patch_gid++){
                
    // get the global id for this boundary patch
    patch_id = Boundary_Condition_Patches(iboundary, bdy_patch_gid);
    Surface_Nodes = Boundary_Patches(patch_id).node_set;
    //find element index this boundary patch is on
    current_element_index = Boundary_Patches(patch_id).element_id;
    local_surface_id = Boundary_Patches(patch_id).local_patch_id;
    //debug print of local surface ids
    //std::cout << " LOCAL SURFACE IDS " << std::endl;
    //std::cout << local_surface_id << std::endl;

    //loop over quadrature points if this is a distributed force
    for(int iquad=0; iquad < direct_product_count; iquad++){
      
      if(Element_Types(current_element_index)==elements::elem_types::Hex8){

      real_t pointer_basis_values[8];
      real_t pointer_basis_derivative_s1[8];
      real_t pointer_basis_derivative_s2[8];
      ViewCArray<real_t> basis_values(pointer_basis_values,8);
      ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1,8);
      ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2,8);
      CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(4,num_dim);
      CArrayKokkos<real_t, array_layout, device_type, memory_traits> surf_basis_derivative_s1(4,num_dim);
      CArrayKokkos<real_t, array_layout, device_type, memory_traits> surf_basis_derivative_s2(4,num_dim);
      CArrayKokkos<real_t, array_layout, device_type, memory_traits> surf_basis_values(4,num_dim);
      int local_nodes[4];
      //set current quadrature point
      y_quad = iquad / num_gauss_points;
      x_quad = iquad % num_gauss_points;
      
      if(local_surface_id<2){
      surf_dim1 = 0;
      surf_dim2 = 1;
      quad_coordinate(0) = legendre_nodes_1D(x_quad);
      quad_coordinate(1) = legendre_nodes_1D(y_quad);
      //set to -1 or 1 for an isoparametric space
        if(local_surface_id%2==0)
        quad_coordinate(2) = -1;
        else
        quad_coordinate(2) = 1;
      }
      else if(local_surface_id<4){
      surf_dim1 = 0;
      surf_dim2 = 2;
      quad_coordinate(0) = legendre_nodes_1D(x_quad);
      quad_coordinate(2) = legendre_nodes_1D(y_quad);
      //set to -1 or 1 for an isoparametric space
        if(local_surface_id%2==0)
        quad_coordinate(1) = -1;
        else
        quad_coordinate(1) = 1;
      }
      else if(local_surface_id<6){
      surf_dim1 = 1;
      surf_dim2 = 2;
      quad_coordinate(1) = legendre_nodes_1D(x_quad);
      quad_coordinate(2) = legendre_nodes_1D(y_quad);
      //set to -1 or 1 for an isoparametric space
        if(local_surface_id%2==0)
        quad_coordinate(0) = -1;
        else
        quad_coordinate(0) = 1;
      }
      
      //set current quadrature weight
      quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
      quad_coordinate_weight(1) = legendre_weights_1D(y_quad);

      //find local dof set for this surface
      local_nodes[0] = elem->surface_to_dof_lid(local_surface_id,0);
      local_nodes[1] = elem->surface_to_dof_lid(local_surface_id,1);
      local_nodes[2] = elem->surface_to_dof_lid(local_surface_id,2);
      local_nodes[3] = elem->surface_to_dof_lid(local_surface_id,3);

      //acquire set of nodes for this face
      for(int node_loop=0; node_loop < 4; node_loop++){
        current_node_index = Surface_Nodes(node_loop);
        local_node_index = all_node_map->getLocalElement(current_node_index);
        nodal_positions(node_loop,0) = all_node_coords(local_node_index,0);
        nodal_positions(node_loop,1) = all_node_coords(local_node_index,1);
        nodal_positions(node_loop,2) = all_node_coords(local_node_index,2);
      }

      if(local_surface_id<2){
        //compute shape function derivatives
        elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
        elem->partial_eta_basis(basis_derivative_s2,quad_coordinate);
      }
      else if(local_surface_id<4){
        //compute shape function derivatives
        elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
        elem->partial_mu_basis(basis_derivative_s2,quad_coordinate);
      }
      else if(local_surface_id<6){
        //compute shape function derivatives
        elem->partial_eta_basis(basis_derivative_s1,quad_coordinate);
        elem->partial_mu_basis(basis_derivative_s2,quad_coordinate);
      }

      //set values relevant to this surface
      for(int node_loop=0; node_loop < 4; node_loop++){
        surf_basis_derivative_s1(node_loop) = basis_derivative_s1(local_nodes[node_loop]);
        surf_basis_derivative_s2(node_loop) = basis_derivative_s2(local_nodes[node_loop]);
      }

      //compute derivatives of x,y,z w.r.t the s,t coordinates of this surface; needed to compute dA in surface integral
      //derivative of x,y,z w.r.t s
      JT_row1(0) = 0;
      JT_row1(1) = 0;
      JT_row1(2) = 0;
      for(int node_loop=0; node_loop < 4; node_loop++){
        JT_row1(0) += nodal_positions(node_loop,0)*surf_basis_derivative_s1(node_loop);
        JT_row1(1) += nodal_positions(node_loop,1)*surf_basis_derivative_s1(node_loop);
        JT_row1(2) += nodal_positions(node_loop,2)*surf_basis_derivative_s1(node_loop);
      }

      //derivative of x,y,z w.r.t t
      JT_row2(0) = 0;
      JT_row2(1) = 0;
      JT_row2(2) = 0;
      for(int node_loop=0; node_loop < 4; node_loop++){
        JT_row2(0) += nodal_positions(node_loop,0)*surf_basis_derivative_s2(node_loop);
        JT_row2(1) += nodal_positions(node_loop,1)*surf_basis_derivative_s2(node_loop);
        JT_row2(2) += nodal_positions(node_loop,2)*surf_basis_derivative_s2(node_loop);
      }
      

      //compute jacobian for this surface
      //compute the determinant of the Jacobian
      wedge_product = sqrt(pow(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2),2)+
               pow(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2),2)+
               pow(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1),2));

      //compute shape functions at this point for the element type
      elem->basis(basis_values,quad_coordinate);

      // loop over nodes of this face and 
      for(int node_count = 0; node_count < 4; node_count++){
            
        node_id = nodes_in_elem(current_element_index, local_nodes[node_count]);
        //check if node is local to alter Nodal Forces vector
        if(!map->isNodeGlobalElement(node_id)) continue;
        node_id = map->getLocalElement(node_id);
        
        /*
        //debug print block
        std::cout << " ------------Element "<< current_element_index + 1 <<"--------------"<<std::endl;
        std::cout <<  " = , " << " Wedge Product: " << wedge_product << " local node " << local_nodes[node_count] << " node " << node_gid + 1<< " : s " 
        << quad_coordinate(0) << " t " << quad_coordinate(1) << " w " << quad_coordinate(2) << " basis value  "<< basis_values(local_nodes[node_count])
        << " Nodal Force value"<< Nodal_Forces(num_dim*node_gid)+wedge_product*quad_coordinate_weight(0)*quad_coordinate_weight(1)*force_density[0]*basis_values(local_nodes[node_count]);
        
        std::cout << " }"<< std::endl;
        //end debug print block
        */

        // Accumulate force vector contribution from this quadrature point
        for(int idim = 0; idim < num_dim; idim++){
          if(force_density[idim]!=0)
          //Nodal_Forces(num_dim*node_gid + idim) += wedge_product*quad_coordinate_weight(0)*quad_coordinate_weight(1)*force_density[idim]*basis_values(local_nodes[node_count]);
          Nodal_Forces(num_dim*node_id + idim,0) += wedge_product*quad_coordinate_weight(0)*quad_coordinate_weight(1)*force_density[idim]*basis_values(local_nodes[node_count]);
        }
      }
      }
    }
  }
  }

    //apply line distribution of forces

    //apply point forces

    //apply contribution from non-zero displacement boundary conditions

    //debug print of force vector
    /*
    std::cout << "---------FORCE VECTOR-------------" << std::endl;
    for(int iforce=0; iforce < num_nodes*num_dim; iforce++)
      std::cout << " DOF: "<< iforce+1 << ", "<< Nodal_Forces(iforce) << std::endl;
    */

}

/* ----------------------------------------------------------------------
   Compute the mass of each element; estimated with quadrature
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::compute_element_masses(const_host_vec_array design_densities, bool max_flag){
  //local number of uniquely assigned elements
  size_t nonoverlap_nelements = element_map->getNodeNumElements();
  //initialize memory for volume storage
  host_vec_array Element_Masses = Global_Element_Masses->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
  const_host_vec_array Element_Volumes = Global_Element_Volumes->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_vec_array all_design_densities;
  bool nodal_density_flag = simparam->nodal_density_flag;
  if(nodal_density_flag)
  all_design_densities = all_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  int num_dim = simparam->num_dim;
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam->num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  size_t local_node_id;
  LO ielem;
  GO global_element_index;

  real_t Jacobian, current_density;
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
  CArray<real_t> legendre_nodes_1D(num_gauss_points);
  CArray<real_t> legendre_weights_1D(num_gauss_points);
  real_t pointer_quad_coordinate[num_dim];
  real_t pointer_quad_coordinate_weight[num_dim];
  real_t pointer_interpolated_point[num_dim];
  real_t pointer_JT_row1[num_dim];
  real_t pointer_JT_row2[num_dim];
  real_t pointer_JT_row3[num_dim];
  ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate,num_dim);
  ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight,num_dim);
  ViewCArray<real_t> interpolated_point(pointer_interpolated_point,num_dim);
  ViewCArray<real_t> JT_row1(pointer_JT_row1,num_dim);
  ViewCArray<real_t> JT_row2(pointer_JT_row2,num_dim);
  ViewCArray<real_t> JT_row3(pointer_JT_row3,num_dim);

  real_t pointer_basis_values[elem->num_basis()];
  real_t pointer_basis_derivative_s1[elem->num_basis()];
  real_t pointer_basis_derivative_s2[elem->num_basis()];
  real_t pointer_basis_derivative_s3[elem->num_basis()];
  ViewCArray<real_t> basis_values(pointer_basis_values,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(),num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());

  //initialize weights
  elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
  elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);

  //loop over elements and use quadrature rule to compute volume from Jacobian determinant
  for(int nonoverlapping_ielem = 0; nonoverlapping_ielem < nonoverlap_nelements; nonoverlapping_ielem++){
    global_element_index = element_map->getGlobalElement(nonoverlapping_ielem);
    ielem = all_element_map->getLocalElement(global_element_index);
    if(nodal_density_flag){
    //acquire set of nodes for this local element
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, node_loop));
      nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
      nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
      nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
      if(nodal_density_flag) nodal_density(node_loop) = all_design_densities(local_node_id,0);
      /*
      if(myrank==1&&nodal_positions(node_loop,2)>10000000){
        std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
        std::cout << node_loop+1 <<" " << local_node_id <<" "<< nodes_in_elem(ielem, node_loop) << " "<< nodal_positions(node_loop,2) << std::endl;
        std::fflush(stdout);
      }
      */
      //std::cout << local_node_id << " " << nodes_in_elem(ielem, node_loop) << " "
       //<< nodal_positions(node_loop,0) << " " << nodal_positions(node_loop,1) << " "<< nodal_positions(node_loop,2) << " " << nodal_density(node_loop) <<std::endl;
    }
    
    //debug print of index
    //std::cout << "nonoverlap element id on TASK " << myrank << " is " << nonoverlapping_ielem << std::endl;
    //std::fflush(stdout);

    //initialize element mass
    Element_Masses(nonoverlapping_ielem,0) = 0;
    
    if(Element_Types(ielem)==elements::elem_types::Hex8){
      direct_product_count = std::pow(num_gauss_points,num_dim);
    }
  
    //loop over quadrature points
    for(int iquad=0; iquad < direct_product_count; iquad++){

      //set current quadrature point
      if(num_dim==3) z_quad = iquad/(num_gauss_points*num_gauss_points);
      y_quad = (iquad % (num_gauss_points*num_gauss_points))/num_gauss_points;
      x_quad = iquad % num_gauss_points;
      quad_coordinate(0) = legendre_nodes_1D(x_quad);
      quad_coordinate(1) = legendre_nodes_1D(y_quad);
      if(num_dim==3)
      quad_coordinate(2) = legendre_nodes_1D(z_quad);

      //set current quadrature weight
      quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
      quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
      if(num_dim==3)
      quad_coordinate_weight(2) = legendre_weights_1D(z_quad);

      //compute shape functions at this point for the element type
      elem->basis(basis_values,quad_coordinate);

      //compute all the necessary coordinates and derivatives at this point

      //compute shape function derivatives
      elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
      elem->partial_eta_basis(basis_derivative_s2,quad_coordinate);
      elem->partial_mu_basis(basis_derivative_s3,quad_coordinate);

      //compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
      //derivative of x,y,z w.r.t s
      JT_row1(0) = 0;
      JT_row1(1) = 0;
      JT_row1(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row1(0) += nodal_positions(node_loop,0)*basis_derivative_s1(node_loop);
        JT_row1(1) += nodal_positions(node_loop,1)*basis_derivative_s1(node_loop);
        JT_row1(2) += nodal_positions(node_loop,2)*basis_derivative_s1(node_loop);
      }

      //derivative of x,y,z w.r.t t
      JT_row2(0) = 0;
      JT_row2(1) = 0;
      JT_row2(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row2(0) += nodal_positions(node_loop,0)*basis_derivative_s2(node_loop);
        JT_row2(1) += nodal_positions(node_loop,1)*basis_derivative_s2(node_loop);
        JT_row2(2) += nodal_positions(node_loop,2)*basis_derivative_s2(node_loop);
      }

      //derivative of x,y,z w.r.t w
      JT_row3(0) = 0;
      JT_row3(1) = 0;
      JT_row3(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row3(0) += nodal_positions(node_loop,0)*basis_derivative_s3(node_loop);
        JT_row3(1) += nodal_positions(node_loop,1)*basis_derivative_s3(node_loop);
        JT_row3(2) += nodal_positions(node_loop,2)*basis_derivative_s3(node_loop);
        //debug print
        /*if(myrank==1&&nodal_positions(node_loop,2)*basis_derivative_s3(node_loop)<-10000000){
        std::cout << " ELEMENT VOLUME JACOBIAN DEBUG ON TASK " << myrank << std::endl;
        std::cout << node_loop+1 << " " << JT_row3(2) << " "<< nodal_positions(node_loop,2) <<" "<< basis_derivative_s3(node_loop) << std::endl;
        std::fflush(stdout);
        }*/
      }
    
    
      //compute the determinant of the Jacobian
      Jacobian = JT_row1(0)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                 JT_row1(1)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                 JT_row1(2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1));
      if(Jacobian<0) Jacobian = -Jacobian;

      //compute density
      current_density = 0;
      if(max_flag){
        current_density = 1;
      }
      else{
        for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
          current_density += nodal_density(node_loop)*basis_values(node_loop);
        }
      }
    
      Element_Masses(nonoverlapping_ielem,0) += current_density*quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*Jacobian;
    }
    }
    else{
      Element_Masses(nonoverlapping_ielem,0) = Element_Volumes(nonoverlapping_ielem,0)*design_densities(nonoverlapping_ielem,0);
    }
  }

  //std::ostream &out = std::cout;
  //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  //if(myrank==0)
  //*fos << "Global Element Masses:" << std::endl;
  //Global_Element_Masses->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
}

/* ----------------------------------------------------------------------
   Compute the gradients of mass function with respect to nodal densities
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::compute_nodal_gradients(const_host_vec_array design_variables, host_vec_array design_gradients){
  //local number of uniquely assigned elements
  size_t nonoverlap_nelements = element_map->getNodeNumElements();
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  int num_dim = simparam->num_dim;
  const_host_vec_array all_node_densities;
  bool nodal_density_flag = simparam->nodal_density_flag;
  if(nodal_density_flag)
  all_node_densities = all_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam->num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  size_t local_node_id;
  LO ielem;
  GO global_element_index;
  
  real_t Jacobian;
  #CArrayKokkos<real_t> legendre_nodes_1D(num_gauss_points);
  #CArrayKokkos<real_t> legendre_weights_1D(num_gauss_points);
  CArray<real_t> legendre_nodes_1D(num_gauss_points);
  CArray<real_t> legendre_weights_1D(num_gauss_points);
  real_t pointer_quad_coordinate[num_dim];
  real_t pointer_quad_coordinate_weight[num_dim];
  real_t pointer_interpolated_point[num_dim];
  real_t pointer_JT_row1[num_dim];
  real_t pointer_JT_row2[num_dim];
  real_t pointer_JT_row3[num_dim];
  ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate,num_dim);
  ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight,num_dim);
  ViewCArray<real_t> interpolated_point(pointer_interpolated_point,num_dim);
  ViewCArray<real_t> JT_row1(pointer_JT_row1,num_dim);
  ViewCArray<real_t> JT_row2(pointer_JT_row2,num_dim);
  ViewCArray<real_t> JT_row3(pointer_JT_row3,num_dim);

  real_t pointer_basis_values[elem->num_basis()];
  real_t pointer_basis_derivative_s1[elem->num_basis()];
  real_t pointer_basis_derivative_s2[elem->num_basis()];
  real_t pointer_basis_derivative_s3[elem->num_basis()];
  ViewCArray<real_t> basis_values(pointer_basis_values,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(),num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());

  //initialize weights
  elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
  elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);

  //initialize design gradients to 0
  for(int init = 0; init < nlocal_nodes; init++)
    design_gradients(init,0) = 0;

  //loop over elements and use quadrature rule to compute volume from Jacobian determinant
  for(int ielem = 0; ielem < rnum_elem; ielem++){
    //acquire set of nodes for this local element
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, node_loop));
      nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
      nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
      nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
      if(nodal_density_flag) nodal_density(node_loop) = all_node_densities(local_node_id,0);
      /*
      if(myrank==1&&nodal_positions(node_loop,2)>10000000){
        std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
        std::cout << node_loop+1 <<" " << local_node_id <<" "<< nodes_in_elem(ielem, node_loop) << " "<< nodal_positions(node_loop,2) << std::endl;
        std::fflush(stdout);
      }
      */
      //std::cout << local_node_id << " " << nodes_in_elem(ielem, node_loop) << " " << nodal_positions(node_loop,0) << " " << nodal_positions(node_loop,1) << " "<< nodal_positions(node_loop,2) <<std::endl;
    }
    
    if(Element_Types(ielem)==elements::elem_types::Hex8){
      direct_product_count = std::pow(num_gauss_points,num_dim);
    }
  
    //loop over quadrature points
    for(int iquad=0; iquad < direct_product_count; iquad++){

      //set current quadrature point
      if(num_dim==3) z_quad = iquad/(num_gauss_points*num_gauss_points);
      y_quad = (iquad % (num_gauss_points*num_gauss_points))/num_gauss_points;
      x_quad = iquad % num_gauss_points;
      quad_coordinate(0) = legendre_nodes_1D(x_quad);
      quad_coordinate(1) = legendre_nodes_1D(y_quad);
      if(num_dim==3)
      quad_coordinate(2) = legendre_nodes_1D(z_quad);

      //set current quadrature weight
      quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
      quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
      if(num_dim==3)
      quad_coordinate_weight(2) = legendre_weights_1D(z_quad);

      //compute shape functions at this point for the element type
      elem->basis(basis_values,quad_coordinate);

      //compute all the necessary coordinates and derivatives at this point

      //compute shape function derivatives
      elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
      elem->partial_eta_basis(basis_derivative_s2,quad_coordinate);
      elem->partial_mu_basis(basis_derivative_s3,quad_coordinate);

      //compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
      //derivative of x,y,z w.r.t s
      JT_row1(0) = 0;
      JT_row1(1) = 0;
      JT_row1(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row1(0) += nodal_positions(node_loop,0)*basis_derivative_s1(node_loop);
        JT_row1(1) += nodal_positions(node_loop,1)*basis_derivative_s1(node_loop);
        JT_row1(2) += nodal_positions(node_loop,2)*basis_derivative_s1(node_loop);
      }

      //derivative of x,y,z w.r.t t
      JT_row2(0) = 0;
      JT_row2(1) = 0;
      JT_row2(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row2(0) += nodal_positions(node_loop,0)*basis_derivative_s2(node_loop);
        JT_row2(1) += nodal_positions(node_loop,1)*basis_derivative_s2(node_loop);
        JT_row2(2) += nodal_positions(node_loop,2)*basis_derivative_s2(node_loop);
      }

      //derivative of x,y,z w.r.t w
      JT_row3(0) = 0;
      JT_row3(1) = 0;
      JT_row3(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row3(0) += nodal_positions(node_loop,0)*basis_derivative_s3(node_loop);
        JT_row3(1) += nodal_positions(node_loop,1)*basis_derivative_s3(node_loop);
        JT_row3(2) += nodal_positions(node_loop,2)*basis_derivative_s3(node_loop);
        //debug print
        /*if(myrank==1&&nodal_positions(node_loop,2)*basis_derivative_s3(node_loop)<-10000000){
        std::cout << " ELEMENT VOLUME JACOBIAN DEBUG ON TASK " << myrank << std::endl;
        std::cout << node_loop+1 << " " << JT_row3(2) << " "<< nodal_positions(node_loop,2) <<" "<< basis_derivative_s3(node_loop) << std::endl;
        std::fflush(stdout);
        }*/
      }
    
    
      //compute the determinant of the Jacobian
      Jacobian = JT_row1(0)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                 JT_row1(1)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                 JT_row1(2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1));
      if(Jacobian<0) Jacobian = -Jacobian;
      
      //assign contribution to every local node this element has
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        if(map->isNodeGlobalElement(nodes_in_elem(ielem, node_loop))){
          local_node_id = map->getLocalElement(nodes_in_elem(ielem, node_loop));
          design_gradients(local_node_id,0)+=quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*basis_values(node_loop)*Jacobian;
        }
      }
    }
    
  }

}

/* ------------------------------------------------------------------------------------------------------------------------
   Compute the moment of inertia of each element for a specified component of the inertia tensor; estimated with quadrature
--------------------------------------------------------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::compute_element_moments(const_host_vec_array design_densities, bool max_flag, int moment_component){
  //local number of uniquely assigned elements
  size_t nonoverlap_nelements = element_map->getNodeNumElements();
  //initialize memory for volume storage
  host_vec_array Element_Masses = Global_Element_Masses->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
  host_vec_array Element_Moments;

  if(moment_component==0) Element_Moments = Global_Element_Moments_x->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
  if(moment_component==1) Element_Moments = Global_Element_Moments_y->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
  if(moment_component==2) Element_Moments = Global_Element_Moments_z->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);

  const_host_vec_array Element_Volumes = Global_Element_Volumes->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_vec_array all_design_densities;
  bool nodal_density_flag = simparam->nodal_density_flag;
  if(nodal_density_flag)
  all_design_densities = all_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  int num_dim = simparam->num_dim;
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam->num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  size_t local_node_id;
  LO ielem;
  GO global_element_index;

  real_t Jacobian, current_density;
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
  CArray<real_t> legendre_nodes_1D(num_gauss_points);
  CArray<real_t> legendre_weights_1D(num_gauss_points);
  real_t pointer_quad_coordinate[num_dim];
  real_t pointer_quad_coordinate_weight[num_dim];
  real_t pointer_interpolated_point[num_dim];
  real_t pointer_JT_row1[num_dim];
  real_t pointer_JT_row2[num_dim];
  real_t pointer_JT_row3[num_dim];
  ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate,num_dim);
  ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight,num_dim);
  ViewCArray<real_t> interpolated_point(pointer_interpolated_point,num_dim);
  ViewCArray<real_t> JT_row1(pointer_JT_row1,num_dim);
  ViewCArray<real_t> JT_row2(pointer_JT_row2,num_dim);
  ViewCArray<real_t> JT_row3(pointer_JT_row3,num_dim);

  real_t pointer_basis_values[elem->num_basis()];
  real_t pointer_basis_derivative_s1[elem->num_basis()];
  real_t pointer_basis_derivative_s2[elem->num_basis()];
  real_t pointer_basis_derivative_s3[elem->num_basis()];
  ViewCArray<real_t> basis_values(pointer_basis_values,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(),num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_position(num_dim);

  //initialize weights
  elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
  elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);

  //loop over elements and use quadrature rule to compute volume from Jacobian determinant
  for(int nonoverlapping_ielem = 0; nonoverlapping_ielem < nonoverlap_nelements; nonoverlapping_ielem++){
    global_element_index = element_map->getGlobalElement(nonoverlapping_ielem);
    ielem = all_element_map->getLocalElement(global_element_index);
    //acquire set of nodes for this local element
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, node_loop));
      nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
      nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
      nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
      if(nodal_density_flag) nodal_density(node_loop) = all_design_densities(local_node_id,0);
      /*
      if(myrank==1&&nodal_positions(node_loop,2)>10000000){
        std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
        std::cout << node_loop+1 <<" " << local_node_id <<" "<< nodes_in_elem(ielem, node_loop) << " "<< nodal_positions(node_loop,2) << std::endl;
        std::fflush(stdout);
      }
      */
      //std::cout << local_node_id << " " << nodes_in_elem(ielem, node_loop) << " "
       //<< nodal_positions(node_loop,0) << " " << nodal_positions(node_loop,1) << " "<< nodal_positions(node_loop,2) << " " << nodal_density(node_loop) <<std::endl;
    }
    
    //debug print of index
    //std::cout << "nonoverlap element id on TASK " << myrank << " is " << nonoverlapping_ielem << std::endl;
    //std::fflush(stdout);

    //initialize element mass
    Element_Moments(nonoverlapping_ielem,0) = 0;
    
    if(Element_Types(ielem)==elements::elem_types::Hex8){
      direct_product_count = std::pow(num_gauss_points,num_dim);
    }
  
    //loop over quadrature points
    for(int iquad=0; iquad < direct_product_count; iquad++){

      //set current quadrature point
      if(num_dim==3) z_quad = iquad/(num_gauss_points*num_gauss_points);
      y_quad = (iquad % (num_gauss_points*num_gauss_points))/num_gauss_points;
      x_quad = iquad % num_gauss_points;
      quad_coordinate(0) = legendre_nodes_1D(x_quad);
      quad_coordinate(1) = legendre_nodes_1D(y_quad);
      if(num_dim==3)
      quad_coordinate(2) = legendre_nodes_1D(z_quad);

      //set current quadrature weight
      quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
      quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
      if(num_dim==3)
      quad_coordinate_weight(2) = legendre_weights_1D(z_quad);

      //compute shape functions at this point for the element type
      elem->basis(basis_values,quad_coordinate);

      //compute all the necessary coordinates and derivatives at this point

      //compute shape function derivatives
      elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
      elem->partial_eta_basis(basis_derivative_s2,quad_coordinate);
      elem->partial_mu_basis(basis_derivative_s3,quad_coordinate);

      //compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
      //derivative of x,y,z w.r.t s
      JT_row1(0) = 0;
      JT_row1(1) = 0;
      JT_row1(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row1(0) += nodal_positions(node_loop,0)*basis_derivative_s1(node_loop);
        JT_row1(1) += nodal_positions(node_loop,1)*basis_derivative_s1(node_loop);
        JT_row1(2) += nodal_positions(node_loop,2)*basis_derivative_s1(node_loop);
      }

      //derivative of x,y,z w.r.t t
      JT_row2(0) = 0;
      JT_row2(1) = 0;
      JT_row2(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row2(0) += nodal_positions(node_loop,0)*basis_derivative_s2(node_loop);
        JT_row2(1) += nodal_positions(node_loop,1)*basis_derivative_s2(node_loop);
        JT_row2(2) += nodal_positions(node_loop,2)*basis_derivative_s2(node_loop);
      }

      //derivative of x,y,z w.r.t w
      JT_row3(0) = 0;
      JT_row3(1) = 0;
      JT_row3(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row3(0) += nodal_positions(node_loop,0)*basis_derivative_s3(node_loop);
        JT_row3(1) += nodal_positions(node_loop,1)*basis_derivative_s3(node_loop);
        JT_row3(2) += nodal_positions(node_loop,2)*basis_derivative_s3(node_loop);
        //debug print
        /*if(myrank==1&&nodal_positions(node_loop,2)*basis_derivative_s3(node_loop)<-10000000){
        std::cout << " ELEMENT VOLUME JACOBIAN DEBUG ON TASK " << myrank << std::endl;
        std::cout << node_loop+1 << " " << JT_row3(2) << " "<< nodal_positions(node_loop,2) <<" "<< basis_derivative_s3(node_loop) << std::endl;
        std::fflush(stdout);
        }*/
      }
    
    
      //compute the determinant of the Jacobian
      Jacobian = JT_row1(0)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                 JT_row1(1)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                 JT_row1(2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1));
      if(Jacobian<0) Jacobian = -Jacobian;

      //compute density
      if(max_flag){
        current_density = 1;
      }
      else{
        if(nodal_density_flag){
          current_density = 0;
          for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
            current_density += nodal_density(node_loop)*basis_values(node_loop);
          }
        }// if
        else{
          current_density = design_densities(nonoverlapping_ielem,0);
        }
      }

      //compute current position
      current_position(0) = current_position(1) = current_position(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        current_position(moment_component) += nodal_positions(node_loop,moment_component)*basis_values(node_loop);
      }

      Element_Moments(nonoverlapping_ielem,0) += current_density*current_position(moment_component)*quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*Jacobian;
    }
  }

  //std::ostream &out = std::cout;
  //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  //if(myrank==0)
  //*fos << "Global Element Masses:" << std::endl;
  //Global_Element_Masses->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
}

/* ---------------------------------------------------------------------------------------------------
   Compute the gradients of the specified moment of inertia component with respect to design densities
------------------------------------------------------------------------------------------------------ */

void Parallel_Nonlinear_Solver::compute_moment_gradients(const_host_vec_array design_variables, host_vec_array design_gradients, int moment_component){
  //local number of uniquely assigned elements
  size_t nonoverlap_nelements = element_map->getNodeNumElements();
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  int num_dim = simparam->num_dim;
  const_host_vec_array all_node_densities;
  bool nodal_density_flag = simparam->nodal_density_flag;
  if(nodal_density_flag)
  all_node_densities = all_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam->num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  size_t local_node_id;
  LO ielem;
  GO global_element_index;
  
  real_t Jacobian;
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
  CArray<real_t> legendre_nodes_1D(num_gauss_points);
  CArray<real_t> legendre_weights_1D(num_gauss_points);
  real_t pointer_quad_coordinate[num_dim];
  real_t pointer_quad_coordinate_weight[num_dim];
  real_t pointer_interpolated_point[num_dim];
  real_t pointer_JT_row1[num_dim];
  real_t pointer_JT_row2[num_dim];
  real_t pointer_JT_row3[num_dim];
  ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate,num_dim);
  ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight,num_dim);
  ViewCArray<real_t> interpolated_point(pointer_interpolated_point,num_dim);
  ViewCArray<real_t> JT_row1(pointer_JT_row1,num_dim);
  ViewCArray<real_t> JT_row2(pointer_JT_row2,num_dim);
  ViewCArray<real_t> JT_row3(pointer_JT_row3,num_dim);

  real_t pointer_basis_values[elem->num_basis()];
  real_t pointer_basis_derivative_s1[elem->num_basis()];
  real_t pointer_basis_derivative_s2[elem->num_basis()];
  real_t pointer_basis_derivative_s3[elem->num_basis()];
  ViewCArray<real_t> basis_values(pointer_basis_values,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(),num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_position(num_dim);

  //initialize weights
  elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
  elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);

  //initialize design gradients to 0
  for(int init = 0; init < nlocal_nodes; init++)
    design_gradients(init,moment_component) = 0;

  //loop over elements and use quadrature rule to compute volume from Jacobian determinant
  for(int ielem = 0; ielem < rnum_elem; ielem++){
    //acquire set of nodes for this local element
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, node_loop));
      nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
      nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
      nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
      if(nodal_density_flag) nodal_density(node_loop) = all_node_densities(local_node_id,0);
      /*
      if(myrank==1&&nodal_positions(node_loop,2)>10000000){
        std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
        std::cout << node_loop+1 <<" " << local_node_id <<" "<< nodes_in_elem(ielem, node_loop) << " "<< nodal_positions(node_loop,2) << std::endl;
        std::fflush(stdout);
      }
      */
      //std::cout << local_node_id << " " << nodes_in_elem(ielem, node_loop) << " " << nodal_positions(node_loop,0) << " " << nodal_positions(node_loop,1) << " "<< nodal_positions(node_loop,2) <<std::endl;
    }
    
    if(Element_Types(ielem)==elements::elem_types::Hex8){
      direct_product_count = std::pow(num_gauss_points,num_dim);
    }
  
    //loop over quadrature points
    for(int iquad=0; iquad < direct_product_count; iquad++){

      //set current quadrature point
      if(num_dim==3) z_quad = iquad/(num_gauss_points*num_gauss_points);
      y_quad = (iquad % (num_gauss_points*num_gauss_points))/num_gauss_points;
      x_quad = iquad % num_gauss_points;
      quad_coordinate(0) = legendre_nodes_1D(x_quad);
      quad_coordinate(1) = legendre_nodes_1D(y_quad);
      if(num_dim==3)
      quad_coordinate(2) = legendre_nodes_1D(z_quad);

      //set current quadrature weight
      quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
      quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
      if(num_dim==3)
      quad_coordinate_weight(2) = legendre_weights_1D(z_quad);

      //compute shape functions at this point for the element type
      elem->basis(basis_values,quad_coordinate);

      //compute all the necessary coordinates and derivatives at this point

      //compute shape function derivatives
      elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
      elem->partial_eta_basis(basis_derivative_s2,quad_coordinate);
      elem->partial_mu_basis(basis_derivative_s3,quad_coordinate);

      //compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
      //derivative of x,y,z w.r.t s
      JT_row1(0) = 0;
      JT_row1(1) = 0;
      JT_row1(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row1(0) += nodal_positions(node_loop,0)*basis_derivative_s1(node_loop);
        JT_row1(1) += nodal_positions(node_loop,1)*basis_derivative_s1(node_loop);
        JT_row1(2) += nodal_positions(node_loop,2)*basis_derivative_s1(node_loop);
      }

      //derivative of x,y,z w.r.t t
      JT_row2(0) = 0;
      JT_row2(1) = 0;
      JT_row2(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row2(0) += nodal_positions(node_loop,0)*basis_derivative_s2(node_loop);
        JT_row2(1) += nodal_positions(node_loop,1)*basis_derivative_s2(node_loop);
        JT_row2(2) += nodal_positions(node_loop,2)*basis_derivative_s2(node_loop);
      }

      //derivative of x,y,z w.r.t w
      JT_row3(0) = 0;
      JT_row3(1) = 0;
      JT_row3(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row3(0) += nodal_positions(node_loop,0)*basis_derivative_s3(node_loop);
        JT_row3(1) += nodal_positions(node_loop,1)*basis_derivative_s3(node_loop);
        JT_row3(2) += nodal_positions(node_loop,2)*basis_derivative_s3(node_loop);
        //debug print
        /*if(myrank==1&&nodal_positions(node_loop,2)*basis_derivative_s3(node_loop)<-10000000){
        std::cout << " ELEMENT VOLUME JACOBIAN DEBUG ON TASK " << myrank << std::endl;
        std::cout << node_loop+1 << " " << JT_row3(2) << " "<< nodal_positions(node_loop,2) <<" "<< basis_derivative_s3(node_loop) << std::endl;
        std::fflush(stdout);
        }*/
      }
    
    
      //compute the determinant of the Jacobian
      Jacobian = JT_row1(0)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                 JT_row1(1)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                 JT_row1(2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1));
      if(Jacobian<0) Jacobian = -Jacobian;

      //compute current position
      current_position(0) = current_position(1) = current_position(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        current_position(moment_component) += nodal_positions(node_loop,moment_component)*basis_values(node_loop);
      }
      
      //assign contribution to every local node this element has
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        if(map->isNodeGlobalElement(nodes_in_elem(ielem, node_loop))){
          local_node_id = map->getLocalElement(nodes_in_elem(ielem, node_loop));
            design_gradients(local_node_id,moment_component)+=quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*basis_values(node_loop)*current_position(moment_component)*Jacobian;
        }
      }
    }
    
  }

}


/* ------------------------------------------------------------------------------------------------------------------------
   Compute the moment of inertia of each element for a specified component of the inertia tensor; estimated with quadrature
--------------------------------------------------------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::compute_element_moments_of_inertia(const_host_vec_array design_densities, bool max_flag, int inertia_component){
  //local number of uniquely assigned elements
  size_t nonoverlap_nelements = element_map->getNodeNumElements();
  //initialize memory for volume storage
  host_vec_array Element_Masses = Global_Element_Masses->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
  host_vec_array Element_Moments_of_Inertia;

  if(inertia_component==0) Element_Moments_of_Inertia = Global_Element_Moments_of_Inertia_xx->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
  if(inertia_component==1) Element_Moments_of_Inertia = Global_Element_Moments_of_Inertia_yy->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
  if(inertia_component==2) Element_Moments_of_Inertia = Global_Element_Moments_of_Inertia_zz->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
  if(inertia_component==3) Element_Moments_of_Inertia = Global_Element_Moments_of_Inertia_xy->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
  if(inertia_component==4) Element_Moments_of_Inertia = Global_Element_Moments_of_Inertia_xz->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
  if(inertia_component==5) Element_Moments_of_Inertia = Global_Element_Moments_of_Inertia_yz->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);

  const_host_vec_array Element_Volumes = Global_Element_Volumes->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_vec_array all_design_densities;
  bool nodal_density_flag = simparam->nodal_density_flag;
  if(nodal_density_flag)
  all_design_densities = all_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  int num_dim = simparam->num_dim;
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam->num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  size_t local_node_id;
  LO ielem;
  GO global_element_index;
  real_t delx1, delx2;

  real_t Jacobian, current_density;
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
  CArray<real_t> legendre_nodes_1D(num_gauss_points);
  CArray<real_t> legendre_weights_1D(num_gauss_points);
  real_t pointer_quad_coordinate[num_dim];
  real_t pointer_quad_coordinate_weight[num_dim];
  real_t pointer_interpolated_point[num_dim];
  real_t pointer_JT_row1[num_dim];
  real_t pointer_JT_row2[num_dim];
  real_t pointer_JT_row3[num_dim];
  ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate,num_dim);
  ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight,num_dim);
  ViewCArray<real_t> interpolated_point(pointer_interpolated_point,num_dim);
  ViewCArray<real_t> JT_row1(pointer_JT_row1,num_dim);
  ViewCArray<real_t> JT_row2(pointer_JT_row2,num_dim);
  ViewCArray<real_t> JT_row3(pointer_JT_row3,num_dim);

  real_t pointer_basis_values[elem->num_basis()];
  real_t pointer_basis_derivative_s1[elem->num_basis()];
  real_t pointer_basis_derivative_s2[elem->num_basis()];
  real_t pointer_basis_derivative_s3[elem->num_basis()];
  ViewCArray<real_t> basis_values(pointer_basis_values,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(),num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_position(num_dim);

  //initialize weights
  elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
  elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);

  //loop over elements and use quadrature rule to compute volume from Jacobian determinant
  for(int nonoverlapping_ielem = 0; nonoverlapping_ielem < nonoverlap_nelements; nonoverlapping_ielem++){
    global_element_index = element_map->getGlobalElement(nonoverlapping_ielem);
    ielem = all_element_map->getLocalElement(global_element_index);
    //acquire set of nodes for this local element
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, node_loop));
      nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
      nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
      nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
      if(nodal_density_flag) nodal_density(node_loop) = all_design_densities(local_node_id,0);
      /*
      if(myrank==1&&nodal_positions(node_loop,2)>10000000){
        std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
        std::cout << node_loop+1 <<" " << local_node_id <<" "<< nodes_in_elem(ielem, node_loop) << " "<< nodal_positions(node_loop,2) << std::endl;
        std::fflush(stdout);
      }
      */
      //std::cout << local_node_id << " " << nodes_in_elem(ielem, node_loop) << " "
       //<< nodal_positions(node_loop,0) << " " << nodal_positions(node_loop,1) << " "<< nodal_positions(node_loop,2) << " " << nodal_density(node_loop) <<std::endl;
    }
    
    //debug print of index
    //std::cout << "nonoverlap element id on TASK " << myrank << " is " << nonoverlapping_ielem << std::endl;
    //std::fflush(stdout);

    //initialize element mass
    Element_Moments_of_Inertia(nonoverlapping_ielem,0) = 0;
    
    if(Element_Types(ielem)==elements::elem_types::Hex8){
      direct_product_count = std::pow(num_gauss_points,num_dim);
    }
  
    //loop over quadrature points
    for(int iquad=0; iquad < direct_product_count; iquad++){

      //set current quadrature point
      if(num_dim==3) z_quad = iquad/(num_gauss_points*num_gauss_points);
      y_quad = (iquad % (num_gauss_points*num_gauss_points))/num_gauss_points;
      x_quad = iquad % num_gauss_points;
      quad_coordinate(0) = legendre_nodes_1D(x_quad);
      quad_coordinate(1) = legendre_nodes_1D(y_quad);
      if(num_dim==3)
      quad_coordinate(2) = legendre_nodes_1D(z_quad);

      //set current quadrature weight
      quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
      quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
      if(num_dim==3)
      quad_coordinate_weight(2) = legendre_weights_1D(z_quad);

      //compute shape functions at this point for the element type
      elem->basis(basis_values,quad_coordinate);

      //compute all the necessary coordinates and derivatives at this point

      //compute shape function derivatives
      elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
      elem->partial_eta_basis(basis_derivative_s2,quad_coordinate);
      elem->partial_mu_basis(basis_derivative_s3,quad_coordinate);

      //compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
      //derivative of x,y,z w.r.t s
      JT_row1(0) = 0;
      JT_row1(1) = 0;
      JT_row1(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row1(0) += nodal_positions(node_loop,0)*basis_derivative_s1(node_loop);
        JT_row1(1) += nodal_positions(node_loop,1)*basis_derivative_s1(node_loop);
        JT_row1(2) += nodal_positions(node_loop,2)*basis_derivative_s1(node_loop);
      }

      //derivative of x,y,z w.r.t t
      JT_row2(0) = 0;
      JT_row2(1) = 0;
      JT_row2(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row2(0) += nodal_positions(node_loop,0)*basis_derivative_s2(node_loop);
        JT_row2(1) += nodal_positions(node_loop,1)*basis_derivative_s2(node_loop);
        JT_row2(2) += nodal_positions(node_loop,2)*basis_derivative_s2(node_loop);
      }

      //derivative of x,y,z w.r.t w
      JT_row3(0) = 0;
      JT_row3(1) = 0;
      JT_row3(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row3(0) += nodal_positions(node_loop,0)*basis_derivative_s3(node_loop);
        JT_row3(1) += nodal_positions(node_loop,1)*basis_derivative_s3(node_loop);
        JT_row3(2) += nodal_positions(node_loop,2)*basis_derivative_s3(node_loop);
        //debug print
        /*if(myrank==1&&nodal_positions(node_loop,2)*basis_derivative_s3(node_loop)<-10000000){
        std::cout << " ELEMENT VOLUME JACOBIAN DEBUG ON TASK " << myrank << std::endl;
        std::cout << node_loop+1 << " " << JT_row3(2) << " "<< nodal_positions(node_loop,2) <<" "<< basis_derivative_s3(node_loop) << std::endl;
        std::fflush(stdout);
        }*/
      }
    
    
      //compute the determinant of the Jacobian
      Jacobian = JT_row1(0)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                 JT_row1(1)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                 JT_row1(2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1));
      if(Jacobian<0) Jacobian = -Jacobian;

      //compute density
      if(max_flag){
        current_density = 1;
      }
      else{
        if(nodal_density_flag){
          current_density = 0;
          for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
            current_density += nodal_density(node_loop)*basis_values(node_loop);
          }
        }// if
        else{
          current_density = design_densities(nonoverlapping_ielem,0);
        }
      }

      //compute current position
      current_position(0) = current_position(1) = current_position(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        current_position(0) += nodal_positions(node_loop,0)*basis_values(node_loop);
        current_position(1) += nodal_positions(node_loop,1)*basis_values(node_loop);
        current_position(2) += nodal_positions(node_loop,2)*basis_values(node_loop);
      }

      if(inertia_component==0){
        delx1 = current_position(1) - center_of_mass[1];
        delx2 = current_position(2) - center_of_mass[2];
        Element_Moments_of_Inertia(nonoverlapping_ielem,0) += current_density*(delx1*delx1 + delx2*delx2)*quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*Jacobian;
      }
      if(inertia_component==1){
        delx1 = current_position(0) - center_of_mass[0];
        delx2 = current_position(2) - center_of_mass[2];
        Element_Moments_of_Inertia(nonoverlapping_ielem,0) += current_density*(delx1*delx1 + delx2*delx2)*quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*Jacobian;
      }
      if(inertia_component==2){
        delx1 = current_position(0) - center_of_mass[0];
        delx2 = current_position(1) - center_of_mass[1];
        Element_Moments_of_Inertia(nonoverlapping_ielem,0) += current_density*(delx1*delx1 + delx2*delx2)*quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*Jacobian;
      }
      if(inertia_component==3){
        delx1 = current_position(0) - center_of_mass[0];
        delx2 = current_position(1) - center_of_mass[1];
        Element_Moments_of_Inertia(nonoverlapping_ielem,0) -= current_density*(delx1*delx2)*quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*Jacobian;
      }
      if(inertia_component==4){
        delx1 = current_position(0) - center_of_mass[0];
        delx2 = current_position(2) - center_of_mass[2];
        Element_Moments_of_Inertia(nonoverlapping_ielem,0) -= current_density*(delx1*delx2)*quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*Jacobian;
      }
      if(inertia_component==5){
        delx1 = current_position(1) - center_of_mass[1];
        delx2 = current_position(2) - center_of_mass[2];
        Element_Moments_of_Inertia(nonoverlapping_ielem,0) -= current_density*(delx1*delx2)*quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*Jacobian;
      }
    }
  }

  //std::ostream &out = std::cout;
  //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  //if(myrank==0)
  //*fos << "Global Element Masses:" << std::endl;
  //Global_Element_Masses->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
}

/* ---------------------------------------------------------------------------------------------------
   Compute the gradients of the specified moment of inertia component with respect to design densities
------------------------------------------------------------------------------------------------------ */

void Parallel_Nonlinear_Solver::compute_moment_of_inertia_gradients(const_host_vec_array design_variables, host_vec_array design_gradients, int inertia_component){
  //local number of uniquely assigned elements
  size_t nonoverlap_nelements = element_map->getNodeNumElements();
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  int num_dim = simparam->num_dim;
  const_host_vec_array all_node_densities;
  bool nodal_density_flag = simparam->nodal_density_flag;
  if(nodal_density_flag)
  all_node_densities = all_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam->num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  size_t local_node_id;
  LO ielem;
  GO global_element_index;
  real_t delx1, delx2;
  
  real_t Jacobian;
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
  CArray<real_t> legendre_nodes_1D(num_gauss_points);
  CArray<real_t> legendre_weights_1D(num_gauss_points);
  real_t pointer_quad_coordinate[num_dim];
  real_t pointer_quad_coordinate_weight[num_dim];
  real_t pointer_interpolated_point[num_dim];
  real_t pointer_JT_row1[num_dim];
  real_t pointer_JT_row2[num_dim];
  real_t pointer_JT_row3[num_dim];
  ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate,num_dim);
  ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight,num_dim);
  ViewCArray<real_t> interpolated_point(pointer_interpolated_point,num_dim);
  ViewCArray<real_t> JT_row1(pointer_JT_row1,num_dim);
  ViewCArray<real_t> JT_row2(pointer_JT_row2,num_dim);
  ViewCArray<real_t> JT_row3(pointer_JT_row3,num_dim);

  real_t pointer_basis_values[elem->num_basis()];
  real_t pointer_basis_derivative_s1[elem->num_basis()];
  real_t pointer_basis_derivative_s2[elem->num_basis()];
  real_t pointer_basis_derivative_s3[elem->num_basis()];
  ViewCArray<real_t> basis_values(pointer_basis_values,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(),num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_position(num_dim);

  //initialize weights
  elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
  elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);

  //initialize design gradients to 0
  for(int init = 0; init < nlocal_nodes; init++)
    design_gradients(init,0) = 0;

  //loop over elements and use quadrature rule to compute volume from Jacobian determinant
  for(int ielem = 0; ielem < rnum_elem; ielem++){
    //acquire set of nodes for this local element
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, node_loop));
      nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
      nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
      nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
      if(nodal_density_flag) nodal_density(node_loop) = all_node_densities(local_node_id,0);
      /*
      if(myrank==1&&nodal_positions(node_loop,2)>10000000){
        std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
        std::cout << node_loop+1 <<" " << local_node_id <<" "<< nodes_in_elem(ielem, node_loop) << " "<< nodal_positions(node_loop,2) << std::endl;
        std::fflush(stdout);
      }
      */
      //std::cout << local_node_id << " " << nodes_in_elem(ielem, node_loop) << " " << nodal_positions(node_loop,0) << " " << nodal_positions(node_loop,1) << " "<< nodal_positions(node_loop,2) <<std::endl;
    }
    
    if(Element_Types(ielem)==elements::elem_types::Hex8){
      direct_product_count = std::pow(num_gauss_points,num_dim);
    }
  
    //loop over quadrature points
    for(int iquad=0; iquad < direct_product_count; iquad++){

      //set current quadrature point
      if(num_dim==3) z_quad = iquad/(num_gauss_points*num_gauss_points);
      y_quad = (iquad % (num_gauss_points*num_gauss_points))/num_gauss_points;
      x_quad = iquad % num_gauss_points;
      quad_coordinate(0) = legendre_nodes_1D(x_quad);
      quad_coordinate(1) = legendre_nodes_1D(y_quad);
      if(num_dim==3)
      quad_coordinate(2) = legendre_nodes_1D(z_quad);

      //set current quadrature weight
      quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
      quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
      if(num_dim==3)
      quad_coordinate_weight(2) = legendre_weights_1D(z_quad);

      //compute shape functions at this point for the element type
      elem->basis(basis_values,quad_coordinate);

      //compute all the necessary coordinates and derivatives at this point

      //compute shape function derivatives
      elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
      elem->partial_eta_basis(basis_derivative_s2,quad_coordinate);
      elem->partial_mu_basis(basis_derivative_s3,quad_coordinate);

      //compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
      //derivative of x,y,z w.r.t s
      JT_row1(0) = 0;
      JT_row1(1) = 0;
      JT_row1(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row1(0) += nodal_positions(node_loop,0)*basis_derivative_s1(node_loop);
        JT_row1(1) += nodal_positions(node_loop,1)*basis_derivative_s1(node_loop);
        JT_row1(2) += nodal_positions(node_loop,2)*basis_derivative_s1(node_loop);
      }

      //derivative of x,y,z w.r.t t
      JT_row2(0) = 0;
      JT_row2(1) = 0;
      JT_row2(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row2(0) += nodal_positions(node_loop,0)*basis_derivative_s2(node_loop);
        JT_row2(1) += nodal_positions(node_loop,1)*basis_derivative_s2(node_loop);
        JT_row2(2) += nodal_positions(node_loop,2)*basis_derivative_s2(node_loop);
      }

      //derivative of x,y,z w.r.t w
      JT_row3(0) = 0;
      JT_row3(1) = 0;
      JT_row3(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row3(0) += nodal_positions(node_loop,0)*basis_derivative_s3(node_loop);
        JT_row3(1) += nodal_positions(node_loop,1)*basis_derivative_s3(node_loop);
        JT_row3(2) += nodal_positions(node_loop,2)*basis_derivative_s3(node_loop);
        //debug print
        /*if(myrank==1&&nodal_positions(node_loop,2)*basis_derivative_s3(node_loop)<-10000000){
        std::cout << " ELEMENT VOLUME JACOBIAN DEBUG ON TASK " << myrank << std::endl;
        std::cout << node_loop+1 << " " << JT_row3(2) << " "<< nodal_positions(node_loop,2) <<" "<< basis_derivative_s3(node_loop) << std::endl;
        std::fflush(stdout);
        }*/
      }
    
    
      //compute the determinant of the Jacobian
      Jacobian = JT_row1(0)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                 JT_row1(1)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                 JT_row1(2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1));
      if(Jacobian<0) Jacobian = -Jacobian;

      //compute current position
      current_position(0) = current_position(1) = current_position(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        current_position(0) += nodal_positions(node_loop,0)*basis_values(node_loop);
        current_position(1) += nodal_positions(node_loop,1)*basis_values(node_loop);
        current_position(2) += nodal_positions(node_loop,2)*basis_values(node_loop);
      }
      
      //assign contribution to every local node this element has
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        if(map->isNodeGlobalElement(nodes_in_elem(ielem, node_loop))){
          local_node_id = map->getLocalElement(nodes_in_elem(ielem, node_loop));
            if(inertia_component==0){
              delx1 = current_position(1) - center_of_mass[1];
              delx2 = current_position(2) - center_of_mass[2];
              design_gradients(local_node_id,0)+=quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*basis_values(node_loop)*(delx1*delx1 + delx2*delx2)*Jacobian;
            }
            if(inertia_component==1){
              delx1 = current_position(0) - center_of_mass[0];
              delx2 = current_position(2) - center_of_mass[2];
              design_gradients(local_node_id,0)+=quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*basis_values(node_loop)*(delx1*delx1 + delx2*delx2)*Jacobian;
            }
            if(inertia_component==2){
              delx1 = current_position(0) - center_of_mass[0];
              delx2 = current_position(1) - center_of_mass[1];
              design_gradients(local_node_id,0)+=quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*basis_values(node_loop)*(delx1*delx1 + delx2*delx2)*Jacobian;
            }
            if(inertia_component==3){
              delx1 = current_position(0) - center_of_mass[0];
              delx2 = current_position(1) - center_of_mass[1];
              design_gradients(local_node_id,0)-=quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*basis_values(node_loop)*(delx1*delx2)*Jacobian;
            }
            if(inertia_component==4){
              delx1 = current_position(0) - center_of_mass[0];
              delx2 = current_position(2) - center_of_mass[2];
              design_gradients(local_node_id,0)-=quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*basis_values(node_loop)*(delx1*delx2)*Jacobian;
            }
            if(inertia_component==5){
              delx1 = current_position(1) - center_of_mass[1];
              delx2 = current_position(2) - center_of_mass[2];
              design_gradients(local_node_id,0)-=quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*basis_values(node_loop)*(delx1*delx2)*Jacobian;
            }
        }
      }
    }
    
  }

}

/* ----------------------------------------------------------------------
   Compute the gradient of strain energy with respect to nodal densities
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::compute_adjoint_gradients(const_host_vec_array design_variables, host_vec_array design_gradients){
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_vec_array all_node_displacements = all_node_displacements_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  
  const_host_vec_array Element_Densities;
  //local variable for host view of densities from the dual view
  bool nodal_density_flag = simparam->nodal_density_flag;
  const_host_vec_array all_node_densities;
  if(nodal_density_flag)
  all_node_densities = all_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  else
  Element_Densities = Global_Element_Densities->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
  int num_dim = simparam->num_dim;
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam->num_gauss_points;
  int strain_max_flag = simparam->strain_max_flag;
  int z_quad,y_quad,x_quad, direct_product_count;
  int solve_flag, zero_strain_flag;
  size_t local_node_id, local_dof_idx, local_dof_idy, local_dof_idz;
  GO current_global_index;

  direct_product_count = std::pow(num_gauss_points,num_dim);
  real_t Element_Modulus_Gradient, Poisson_Ratio;
  real_t Elastic_Constant, Shear_Term, Pressure_Term;
  real_t inner_product, matrix_term, Jacobian;
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
  CArray<real_t> legendre_nodes_1D(num_gauss_points);
  CArray<real_t> legendre_weights_1D(num_gauss_points);
  real_t pointer_quad_coordinate[num_dim];
  real_t pointer_quad_coordinate_weight[num_dim];
  real_t pointer_interpolated_point[num_dim];
  real_t pointer_JT_row1[num_dim];
  real_t pointer_JT_row2[num_dim];
  real_t pointer_JT_row3[num_dim];
  ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate,num_dim);
  ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight,num_dim);
  ViewCArray<real_t> interpolated_point(pointer_interpolated_point,num_dim);
  ViewCArray<real_t> JT_row1(pointer_JT_row1,num_dim);
  ViewCArray<real_t> JT_row2(pointer_JT_row2,num_dim);
  ViewCArray<real_t> JT_row3(pointer_JT_row3,num_dim);

  real_t pointer_basis_values[elem->num_basis()];
  real_t pointer_basis_derivative_s1[elem->num_basis()];
  real_t pointer_basis_derivative_s2[elem->num_basis()];
  real_t pointer_basis_derivative_s3[elem->num_basis()];
  ViewCArray<real_t> basis_values(pointer_basis_values,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(),num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_nodal_displacements(elem->num_basis()*num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());

  size_t Brows;
  if(num_dim==2) Brows = 3;
  if(num_dim==3) Brows = 6;
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix_contribution(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix_contribution(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> C_matrix(Brows,Brows);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> Local_Matrix_Contribution(num_dim*nodes_per_elem,num_dim*nodes_per_elem);

  //initialize weights
  elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
  elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);
  
  real_t current_density = 1;

  //initialize gradient value to zero
  for(size_t inode = 0; inode < nlocal_nodes; inode++)
    design_gradients(inode,0) = 0;

  //loop through each element and assign the contribution to compliance gradient for each of its local nodes
  for(size_t ielem = 0; ielem < rnum_elem; ielem++){
    nodes_per_elem = elem->num_basis();

    //initialize C matrix
    for(int irow = 0; irow < Brows; irow++)
      for(int icol = 0; icol < Brows; icol++)
        C_matrix(irow,icol) = 0;

    //B matrix initialization
    for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < num_dim*nodes_per_elem; icol++){
        CB_matrix(irow,icol) = 0;
      }

    //acquire set of nodes and nodal displacements for this local element
    for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
      local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, node_loop));
      local_dof_idx = all_dof_map->getLocalElement(nodes_in_elem(ielem, node_loop)*num_dim);
      local_dof_idy = local_dof_idx + 1;
      local_dof_idz = local_dof_idx + 2;
      nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
      nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
      nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
      current_nodal_displacements(node_loop*num_dim) = all_node_displacements(local_dof_idx,0);
      current_nodal_displacements(node_loop*num_dim+1) = all_node_displacements(local_dof_idy,0);
      current_nodal_displacements(node_loop*num_dim+2) = all_node_displacements(local_dof_idz,0);
      
      if(nodal_density_flag) nodal_density(node_loop) = all_node_densities(local_node_id,0);
      //debug print
      /*
      std::cout << "node index access x "<< local_node_id << std::endl;
      std::cout << "local index access x "<< local_dof_idx << " displacement x " << current_nodal_displacements(node_loop*num_dim) <<std::endl;
      std::cout << "local index access y "<< local_dof_idy << " displacement y " << current_nodal_displacements(node_loop*num_dim + 1) << std::endl;
      std::cout << "local index access z "<< local_dof_idz << " displacement z " << current_nodal_displacements(node_loop*num_dim + 2) << std::endl; 
      */
    }

    //debug print of current_nodal_displacements
    /*
    std::cout << " ------------nodal displacements for Element "<< ielem + 1 <<"--------------"<<std::endl;
    std::cout << " { ";
    for (int idof = 0; idof < num_dim*nodes_per_elem; idof++){
      std::cout << idof + 1 << " = " << current_nodal_displacements(idof) << " , " ;
    }
    std::cout << " }"<< std::endl;
    */

    //loop over quadrature points
    for(int iquad=0; iquad < direct_product_count; iquad++){

    //set current quadrature point
    if(num_dim==3) z_quad = iquad/(num_gauss_points*num_gauss_points);
    y_quad = (iquad % (num_gauss_points*num_gauss_points))/num_gauss_points;
    x_quad = iquad % num_gauss_points;
    quad_coordinate(0) = legendre_nodes_1D(x_quad);
    quad_coordinate(1) = legendre_nodes_1D(y_quad);
    if(num_dim==3)
    quad_coordinate(2) = legendre_nodes_1D(z_quad);

    //set current quadrature weight
    quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
    quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
    if(num_dim==3)
    quad_coordinate_weight(2) = legendre_weights_1D(z_quad);

    //compute shape functions at this point for the element type
    elem->basis(basis_values,quad_coordinate);

    //compute all the necessary coordinates and derivatives at this point
    //compute shape function derivatives
    elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
    elem->partial_eta_basis(basis_derivative_s2,quad_coordinate);
    elem->partial_mu_basis(basis_derivative_s3,quad_coordinate);

    //compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
    //derivative of x,y,z w.r.t s
    JT_row1(0) = 0;
    JT_row1(1) = 0;
    JT_row1(2) = 0;
    for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
      JT_row1(0) += nodal_positions(node_loop,0)*basis_derivative_s1(node_loop);
      JT_row1(1) += nodal_positions(node_loop,1)*basis_derivative_s1(node_loop);
      JT_row1(2) += nodal_positions(node_loop,2)*basis_derivative_s1(node_loop);
    }

    //derivative of x,y,z w.r.t t
    JT_row2(0) = 0;
    JT_row2(1) = 0;
    JT_row2(2) = 0;
    for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
      JT_row2(0) += nodal_positions(node_loop,0)*basis_derivative_s2(node_loop);
      JT_row2(1) += nodal_positions(node_loop,1)*basis_derivative_s2(node_loop);
      JT_row2(2) += nodal_positions(node_loop,2)*basis_derivative_s2(node_loop);
    }

    //derivative of x,y,z w.r.t w
    JT_row3(0) = 0;
    JT_row3(1) = 0;
    JT_row3(2) = 0;
    for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
      JT_row3(0) += nodal_positions(node_loop,0)*basis_derivative_s3(node_loop);
      JT_row3(1) += nodal_positions(node_loop,1)*basis_derivative_s3(node_loop);
      JT_row3(2) += nodal_positions(node_loop,2)*basis_derivative_s3(node_loop);
    }
    
    //compute the determinant of the Jacobian
    Jacobian = JT_row1(0)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
               JT_row1(1)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
               JT_row1(2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1));
    if(Jacobian<0) Jacobian = -Jacobian;

    //compute density
    current_density = 0;
    if(nodal_density_flag)
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      current_density += nodal_density(node_loop)*basis_values(node_loop);
    }
    //default constant element density
    else{
      current_density = Element_Densities(ielem,0);
    }

    //debug print
    //std::cout << "Current Density " << current_density << std::endl;

    //compute the contributions of this quadrature point to the B matrix
    if(num_dim==2)
    for(int ishape=0; ishape < nodes_per_elem; ishape++){
      B_matrix_contribution(0,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(1,ishape*num_dim) = 0;
      B_matrix_contribution(2,ishape*num_dim) = 0;
      B_matrix_contribution(3,ishape*num_dim) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(5,ishape*num_dim) = 0;
      B_matrix_contribution(0,ishape*num_dim+1) = 0;
      B_matrix_contribution(1,ishape*num_dim+1) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(2,ishape*num_dim+1) = 0;
      B_matrix_contribution(3,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim+1) = 0;
      B_matrix_contribution(5,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(0,ishape*num_dim+2) = 0;
      B_matrix_contribution(1,ishape*num_dim+2) = 0;
      B_matrix_contribution(2,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(3,ishape*num_dim+2) = 0;
      B_matrix_contribution(4,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(5,ishape*num_dim+2) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
    }
    if(num_dim==3)
    for(int ishape=0; ishape < nodes_per_elem; ishape++){
      B_matrix_contribution(0,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(1,ishape*num_dim) = 0;
      B_matrix_contribution(2,ishape*num_dim) = 0;
      B_matrix_contribution(3,ishape*num_dim) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(5,ishape*num_dim) = 0;
      B_matrix_contribution(0,ishape*num_dim+1) = 0;
      B_matrix_contribution(1,ishape*num_dim+1) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(2,ishape*num_dim+1) = 0;
      B_matrix_contribution(3,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim+1) = 0;
      B_matrix_contribution(5,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(0,ishape*num_dim+2) = 0;
      B_matrix_contribution(1,ishape*num_dim+2) = 0;
      B_matrix_contribution(2,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(3,ishape*num_dim+2) = 0;
      B_matrix_contribution(4,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(5,ishape*num_dim+2) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
    }
    
    
    //evaluate local stiffness matrix gradient with respect to igradient
    for(int igradient=0; igradient < nodes_per_elem; igradient++){
      if(!map->isNodeGlobalElement(nodes_in_elem(ielem, igradient))) continue;
      local_node_id = map->getLocalElement(nodes_in_elem(ielem, igradient));
      //look up element material properties at this point as a function of density
      Gradient_Element_Material_Properties(ielem, Element_Modulus_Gradient, Poisson_Ratio, current_density);
      Elastic_Constant = Element_Modulus_Gradient/((1 + Poisson_Ratio)*(1 - 2*Poisson_Ratio));
      Shear_Term = 0.5 - Poisson_Ratio;
      Pressure_Term = 1 - Poisson_Ratio;

      //debug print
      //std::cout << "Element Material Params " << Elastic_Constant << std::endl;

      //compute Elastic (C) matrix
      if(num_dim==2){
        C_matrix(0,0) = Pressure_Term;
        C_matrix(1,1) = Pressure_Term;
        C_matrix(0,1) = Poisson_Ratio;
        C_matrix(1,0) = Poisson_Ratio;
        C_matrix(2,2) = Shear_Term;
      }
      if(num_dim==3){
        C_matrix(0,0) = Pressure_Term;
        C_matrix(1,1) = Pressure_Term;
        C_matrix(2,2) = Pressure_Term;
        C_matrix(0,1) = Poisson_Ratio;
        C_matrix(0,2) = Poisson_Ratio;
        C_matrix(1,0) = Poisson_Ratio;
        C_matrix(1,2) = Poisson_Ratio;
        C_matrix(2,0) = Poisson_Ratio;
        C_matrix(2,1) = Poisson_Ratio;
        C_matrix(3,3) = Shear_Term;
        C_matrix(4,4) = Shear_Term;
        C_matrix(5,5) = Shear_Term;
      }

      //compute the previous multiplied by the Elastic (C) Matrix
      for(int irow=0; irow < Brows; irow++){
        for(int icol=0; icol < num_dim*nodes_per_elem; icol++){
          CB_matrix_contribution(irow,icol) = 0;
          for(int span=0; span < Brows; span++){
            CB_matrix_contribution(irow,icol) += C_matrix(irow,span)*B_matrix_contribution(span,icol);
          }
        }
      }

      //compute the contributions of this quadrature point to all the local stiffness matrix elements
      for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++){
        for(int jfill=0; jfill < num_dim*nodes_per_elem; jfill++){
          matrix_term = 0;
          for(int span = 0; span < Brows; span++){
            matrix_term += B_matrix_contribution(span,ifill)*CB_matrix_contribution(span,jfill);
          }
          Local_Matrix_Contribution(ifill,jfill) = Elastic_Constant*basis_values(igradient)*quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*matrix_term/Jacobian;
        }
      }
      
      //compute inner product for this quadrature point contribution
      inner_product = 0;
      for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++){
        for(int jfill=0; jfill < num_dim*nodes_per_elem; jfill++){
          inner_product += Local_Matrix_Contribution(ifill, jfill)*current_nodal_displacements(ifill)*current_nodal_displacements(jfill);
          //debug
          //if(Local_Matrix_Contribution(ifill, jfill)<0) Local_Matrix_Contribution(ifill, jfill) = - Local_Matrix_Contribution(ifill, jfill);
          //inner_product += Local_Matrix_Contribution(ifill, jfill);
        }
      }
      
      //debug print
      //std::cout << "contribution for " << igradient + 1 << " is " << inner_product << std::endl;
      design_gradients(local_node_id,0) -= inner_product;
      }
    }
  }
  
    
  //debug print

}

/* -------------------------------------------------------------------------------------------
   Communicate ghosts using the current optimization design data
---------------------------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::comm_variables(Teuchos::RCP<const MV> zp){
  
  //set density vector to the current value chosen by the optimizer
  test_node_densities_distributed = zp;
  
  //debug print of design vector
      //std::ostream &out = std::cout;
      //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
      //if(myrank==0)
      //*fos << "Density data :" << std::endl;
      //node_densities_distributed->describe(*fos,Teuchos::VERB_EXTREME);
      //*fos << std::endl;
      //std::fflush(stdout);

  //communicate design densities
  //create import object using local node indices map and all indices map
  Tpetra::Import<LO, GO> importer(map, all_node_map);

  //comms to get ghosts
  all_node_densities_distributed->doImport(*test_node_densities_distributed, importer, Tpetra::INSERT);

  //update_count++;
  //if(update_count==1){
      //MPI_Barrier(world);
      //MPI_Abort(world,4);
  //}
}

/* -------------------------------------------------------------------------------------------
   update nodal displacement information in accordance with current optimization vector
---------------------------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::update_linear_solve(Teuchos::RCP<const MV> zp){
  
  //set density vector to the current value chosen by the optimizer
  test_node_densities_distributed = zp;

  assemble_matrix();
  
  //solve for new nodal displacements
  int solver_exit = solve();
  if(solver_exit == EXIT_SUCCESS){
    std::cout << "Linear Solver Error" << std::endl <<std::flush;
    return;
  }
  
  update_count++;
  //if(update_count==1){
      //MPI_Barrier(world);
      //MPI_Abort(world,4);
  //}
}

/* -------------------------------------------------------------------------------------------
   Compute the maximum nodal strains resulting from minimizing the L2 error
   between strain (subspace solution) and a nodal interpolation (nodal strains defined at each node)
   for each element. Mainly used for output and is approximate.
---------------------------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::compute_nodal_strains(){
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_vec_array all_node_displacements = all_node_displacements_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  host_vec_array all_node_strains = all_node_strains_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  host_vec_array node_strains = node_strains_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  const_host_elem_conn_array node_nconn = node_nconn_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  int num_dim = simparam->num_dim;
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam->num_gauss_points;
  int strain_max_flag = simparam->strain_max_flag;
  int z_quad,y_quad,x_quad, direct_product_count;
  int solve_flag, zero_strain_flag;
  size_t local_node_id, local_dof_idx, local_dof_idy, local_dof_idz;
  //real_t J_min = std::numeric_limits<real_t>::max();
  GO current_global_index;

  direct_product_count = std::pow(num_gauss_points,num_dim);
  real_t matrix_term, current_strain;
  real_t matrix_subterm1, matrix_subterm2, matrix_subterm3, Jacobian;
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
  CArray<real_t> legendre_nodes_1D(num_gauss_points);
  CArray<real_t> legendre_weights_1D(num_gauss_points);
  real_t pointer_quad_coordinate[num_dim];
  real_t pointer_quad_coordinate_weight[num_dim];
  real_t pointer_interpolated_point[num_dim];
  real_t pointer_JT_row1[num_dim];
  real_t pointer_JT_row2[num_dim];
  real_t pointer_JT_row3[num_dim];
  ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate,num_dim);
  ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight,num_dim);
  ViewCArray<real_t> interpolated_point(pointer_interpolated_point,num_dim);
  ViewCArray<real_t> JT_row1(pointer_JT_row1,num_dim);
  ViewCArray<real_t> JT_row2(pointer_JT_row2,num_dim);
  ViewCArray<real_t> JT_row3(pointer_JT_row3,num_dim);

  real_t pointer_basis_values[elem->num_basis()];
  real_t pointer_basis_derivative_s1[elem->num_basis()];
  real_t pointer_basis_derivative_s2[elem->num_basis()];
  real_t pointer_basis_derivative_s3[elem->num_basis()];
  ViewCArray<real_t> basis_values(pointer_basis_values,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(),num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_nodal_displacements(elem->num_basis()*num_dim);

  size_t Brows;
  if(num_dim==2) Brows = 3;
  if(num_dim==3) Brows = 6;
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix_contribution(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> quad_strain(Brows);
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> projection_matrix(max_nodes_per_element,max_nodes_per_element);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> projection_vector(Brows,max_nodes_per_element);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> strain_vector(max_nodes_per_element);
  //Teuchos::SerialSymDenseMatrix<LO,real_t> projection_matrix_pass;
  Teuchos::RCP<Teuchos::SerialDenseMatrix<LO,real_t>> projection_matrix_pass;
  //Teuchos::SerialDenseVector<LO,real_t> projection_vector_pass;
  Teuchos::RCP<Teuchos::SerialDenseVector<LO,real_t>> projection_vector_pass;
  //Teuchos::SerialDenseVector<LO,real_t> strain_vector_pass;
  Teuchos::RCP<Teuchos::SerialDenseVector<LO,real_t>> strain_vector_pass;
  //Teuchos::SerialSpdDenseSolver<LO,real_t> projection_solver;
  Teuchos::SerialDenseSolver<LO,real_t> projection_solver;

  //initialize weights
  elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
  elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);

  //initialize strains to 0
  //local variable for host view in the dual view
  for(int init = 0; init < map->getNodeNumElements(); init++)
    for(int istrain = 0; istrain < Brows; istrain++)
      node_strains(init,istrain) = 0;

  for(int init = 0; init < all_node_map->getNodeNumElements(); init++)
    for(int istrain = 0; istrain < Brows; istrain++)
      all_node_strains(init,istrain) = 0;
  

  real_t current_density = 1;

  //compute nodal strains as the result of minimizing the L2 error between the strain field and the nodal interpolant over each element
  //assign the maximum nodal strain computed from connected elements to nodes sharing elements (safer to know the largest positive or negative values)
  
  for(int ielem = 0; ielem < rnum_elem; ielem++){
    nodes_per_elem = elem->num_basis();
    //B matrix initialization
    for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < num_dim*nodes_per_elem; icol++){
        B_matrix(irow,icol) = 0;
      }

    //initialize projection matrix
    for(int irow=0; irow < max_nodes_per_element; irow++)
      for(int icol=0; icol < max_nodes_per_element; icol++){
        projection_matrix(irow,icol) = 0;
      }

    //initialize projection vector
    for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < max_nodes_per_element; icol++){
        projection_vector(irow,icol) = 0;
      }

    //acquire set of nodes and nodal displacements for this local element
    for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
      local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, node_loop));
      local_dof_idx = all_dof_map->getLocalElement(nodes_in_elem(ielem, node_loop)*num_dim);
      local_dof_idy = local_dof_idx + 1;
      local_dof_idz = local_dof_idx + 2;
      nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
      nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
      nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
      current_nodal_displacements(node_loop*num_dim) = all_node_displacements(local_dof_idx,0);
      current_nodal_displacements(node_loop*num_dim+1) = all_node_displacements(local_dof_idy,0);
      current_nodal_displacements(node_loop*num_dim+2) = all_node_displacements(local_dof_idz,0);

      //debug print
      /*
      std::cout << "node index access x "<< local_node_id << std::endl;
      std::cout << "local index access x "<< local_dof_idx << " displacement x " << current_nodal_displacements(node_loop*num_dim) <<std::endl;
      std::cout << "local index access y "<< local_dof_idy << " displacement y " << current_nodal_displacements(node_loop*num_dim + 1) << std::endl;
      std::cout << "local index access z "<< local_dof_idz << " displacement z " << current_nodal_displacements(node_loop*num_dim + 2) << std::endl; 
      */
    }

    //debug print of current_nodal_displacements
    /*
    std::cout << " ------------nodal displacements for Element "<< ielem + 1 <<"--------------"<<std::endl;
    std::cout << " { ";
    for (int idof = 0; idof < num_dim*nodes_per_elem; idof++){
      std::cout << idof + 1 << " = " << current_nodal_displacements(idof) << " , " ;
    }
    std::cout << " }"<< std::endl;
    */
    //loop over quadrature points
    for(int iquad=0; iquad < direct_product_count; iquad++){

    //set current quadrature point
    if(num_dim==3) z_quad = iquad/(num_gauss_points*num_gauss_points);
    y_quad = (iquad % (num_gauss_points*num_gauss_points))/num_gauss_points;
    x_quad = iquad % num_gauss_points;
    quad_coordinate(0) = legendre_nodes_1D(x_quad);
    quad_coordinate(1) = legendre_nodes_1D(y_quad);
    if(num_dim==3)
    quad_coordinate(2) = legendre_nodes_1D(z_quad);

    //set current quadrature weight
    quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
    quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
    if(num_dim==3)
    quad_coordinate_weight(2) = legendre_weights_1D(z_quad);

    //compute shape functions at this point for the element type
    elem->basis(basis_values,quad_coordinate);

    //compute all the necessary coordinates and derivatives at this point
    //compute shape function derivatives
    elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
    elem->partial_eta_basis(basis_derivative_s2,quad_coordinate);
    elem->partial_mu_basis(basis_derivative_s3,quad_coordinate);

    //compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
    //derivative of x,y,z w.r.t s
    JT_row1(0) = 0;
    JT_row1(1) = 0;
    JT_row1(2) = 0;
    for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
      JT_row1(0) += nodal_positions(node_loop,0)*basis_derivative_s1(node_loop);
      JT_row1(1) += nodal_positions(node_loop,1)*basis_derivative_s1(node_loop);
      JT_row1(2) += nodal_positions(node_loop,2)*basis_derivative_s1(node_loop);
    }

    //derivative of x,y,z w.r.t t
    JT_row2(0) = 0;
    JT_row2(1) = 0;
    JT_row2(2) = 0;
    for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
      JT_row2(0) += nodal_positions(node_loop,0)*basis_derivative_s2(node_loop);
      JT_row2(1) += nodal_positions(node_loop,1)*basis_derivative_s2(node_loop);
      JT_row2(2) += nodal_positions(node_loop,2)*basis_derivative_s2(node_loop);
    }

    //derivative of x,y,z w.r.t w
    JT_row3(0) = 0;
    JT_row3(1) = 0;
    JT_row3(2) = 0;
    for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
      JT_row3(0) += nodal_positions(node_loop,0)*basis_derivative_s3(node_loop);
      JT_row3(1) += nodal_positions(node_loop,1)*basis_derivative_s3(node_loop);
      JT_row3(2) += nodal_positions(node_loop,2)*basis_derivative_s3(node_loop);
    }
    
    
    //compute the determinant of the Jacobian
    Jacobian = JT_row1(0)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
               JT_row1(1)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
               JT_row1(2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1));
    if(Jacobian<0) Jacobian = -Jacobian;

    //compute the contributions of this quadrature point to the B matrix
    if(num_dim==2)
    for(int ishape=0; ishape < nodes_per_elem; ishape++){
      B_matrix_contribution(0,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(1,ishape*num_dim) = 0;
      B_matrix_contribution(2,ishape*num_dim) = 0;
      B_matrix_contribution(3,ishape*num_dim) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(5,ishape*num_dim) = 0;
      B_matrix_contribution(0,ishape*num_dim+1) = 0;
      B_matrix_contribution(1,ishape*num_dim+1) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(2,ishape*num_dim+1) = 0;
      B_matrix_contribution(3,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim+1) = 0;
      B_matrix_contribution(5,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(0,ishape*num_dim+2) = 0;
      B_matrix_contribution(1,ishape*num_dim+2) = 0;
      B_matrix_contribution(2,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(3,ishape*num_dim+2) = 0;
      B_matrix_contribution(4,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(5,ishape*num_dim+2) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
    }
    if(num_dim==3)
    for(int ishape=0; ishape < nodes_per_elem; ishape++){
      B_matrix_contribution(0,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(1,ishape*num_dim) = 0;
      B_matrix_contribution(2,ishape*num_dim) = 0;
      B_matrix_contribution(3,ishape*num_dim) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(5,ishape*num_dim) = 0;
      B_matrix_contribution(0,ishape*num_dim+1) = 0;
      B_matrix_contribution(1,ishape*num_dim+1) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      B_matrix_contribution(2,ishape*num_dim+1) = 0;
      B_matrix_contribution(3,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(4,ishape*num_dim+1) = 0;
      B_matrix_contribution(5,ishape*num_dim+1) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(0,ishape*num_dim+2) = 0;
      B_matrix_contribution(1,ishape*num_dim+2) = 0;
      B_matrix_contribution(2,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      B_matrix_contribution(3,ishape*num_dim+2) = 0;
      B_matrix_contribution(4,ishape*num_dim+2) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
          basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
          basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(5,ishape*num_dim+2) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
          basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
          basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
    }
    
    //multiply by displacement vector to get strain at this quadrature point
    //division by J ommited since it is multiplied out later
    for(int irow=0; irow < Brows; irow++){
      quad_strain(irow) = 0;
      for(int icol=0; icol < num_dim*nodes_per_elem; icol++){
        quad_strain(irow) += B_matrix_contribution(irow,icol)*current_nodal_displacements(icol);
      }
    }

    //debug print of quad strain
    /*
    std::cout << " ------------Strain at quadrature point for Element "<< ielem + 1 <<"--------------"<<std::endl;
    std::cout << " { ";
    for (int idof = 0; idof < Brows; idof++){
      std::cout << idof + 1 << " = " << quad_strain(idof)/Jacobian*100 << " , " ;
    }
    std::cout << " }"<< std::endl;
    std::fflush(stdout);
    */

    //compute contribution to RHS projection vector
    for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < nodes_per_elem; icol++){
        //debug print
        /*
        if(irow==2){
        std::cout << " ------------Strain Projection Vector Contribution"<< ielem + 1 <<"--------------"<<std::endl;
        std::cout << icol + 1 << " = " << projection_vector(irow,icol) << " " << quad_strain(irow) << " " << basis_values(icol) << " " << quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2) << std::endl;
        std::fflush(stdout);
        }
        */
        projection_vector(irow,icol) += quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*quad_strain(irow)*basis_values(icol);
        
      }

    //compute contribution to projection matrix (only upper part is set)
    for(int irow=0; irow < nodes_per_elem; irow++)
      for(int icol=0; icol < nodes_per_elem; icol++){
        //if(irow<=icol)
        projection_matrix(irow,icol) += quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*basis_values(irow)*basis_values(icol)*Jacobian;
      }

    //accumulate B matrix
    //for(int irow=0; irow < Brows; irow++)
      //for(int icol=0; icol < num_dim*nodes_per_elem; icol++)
        //B_matrix(irow,icol) += quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*B_matrix_contribution(irow,icol)*Jacobian;

    
    //debug print of B matrix per quadrature point
    /*
    std::cout << " ------------B MATRIX QUADRATURE CONTRIBUTION"<< ielem + 1 <<"--------------"<<std::endl;
    for (int idof = 0; idof < Brows; idof++){
      std::cout << "row: " << idof + 1 << " { ";
      for (int istride = 0; istride < nodes_per_elem*num_dim; istride++){
        std::cout << istride + 1 << " = " << B_matrix_contribution(idof,istride) << " , " ;
      }
      std::cout << " }"<< std::endl;
    }
    */
    //end debug block
    
    
    }
    
    //debug print of projection vector across all components
    /*
    std::cout << " ------------Strain Projection Vector"<< ielem + 1 <<"--------------"<<std::endl;
    for (int idof = 0; idof < Brows; idof++){
      std::cout << "row: " << idof + 1 << " { ";
      for (int istride = 0; istride < nodes_per_elem; istride++){
        std::cout << istride + 1 << " = " << projection_vector(idof,istride) << " , " ;
      }
      std::cout << " }"<< std::endl;
    }
    std::fflush(stdout);
    */
    //end debug block
    
    //solve small linear system for nodal strain values
    //scale for conditioning
    /*
    for(int irow=0; irow < nodes_per_elem; irow++)
      for(int icol=0; icol < nodes_per_elem; icol++){
        //if(irow<=icol)
        projection_matrix(irow,icol) /= J_min;
      }
      */
    //use percentages
    for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < nodes_per_elem; icol++){
        projection_vector(irow,icol) *= 100;
      }
      
    //construct matrix and vector wrappers for dense solver
    //projection_matrix_pass = Teuchos::rcp( new Teuchos::SerialSymDenseMatrix<LO,real_t>(Teuchos::View, true, projection_matrix.pointer(), nodes_per_elem, nodes_per_elem));

    //debug print of matrix
    //projection_matrix_pass->print(std::cout);

    strain_vector_pass = Teuchos::rcp( new Teuchos::SerialDenseVector<LO,real_t>(Teuchos::View, strain_vector.pointer(), nodes_per_elem));
    //loop through strain components and solve for nodal values of that component
    for(int istrain = 0; istrain < Brows; istrain++){
      //check if projection vector is zero due to zero strains
      zero_strain_flag = 1;
      for(int icol=0; icol < nodes_per_elem; icol++){
        if(fabs(projection_vector(istrain,icol))>STRAIN_EPSILON) zero_strain_flag = 0;
      }
      if(!zero_strain_flag){
        projection_vector_pass = Teuchos::rcp( new Teuchos::SerialDenseVector<LO,real_t>(Teuchos::View, &projection_vector(istrain,0), nodes_per_elem));
        projection_matrix_pass = Teuchos::rcp( new Teuchos::SerialDenseMatrix<LO,real_t>(Teuchos::Copy, projection_matrix.pointer(), nodes_per_elem, nodes_per_elem, nodes_per_elem));
        //debug print of vectors
        //projection_vector_pass->print(std::cout);
        //projection_matrix_pass->print(std::cout);
        projection_solver.setMatrix(projection_matrix_pass);
        projection_solver.setVectors(strain_vector_pass, projection_vector_pass);
        projection_solver.factorWithEquilibration(true);
        solve_flag = projection_solver.solve();
        
        //debug print
        //std::cout << "STRAIN COMPONENT ON NODES " << istrain + 1 << std::endl;
        //strain_vector_pass->print(std::cout);
        if(solve_flag) std::cout << "Projection Solve Failed With: " << solve_flag << std::endl;

        //contribute equivalent nodal strain for this element to corresponding global nodes
        //replace if abs greater than abs of currently assigned strain; accumulate average if flagged for node average
        for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
          current_global_index = nodes_in_elem(ielem, node_loop);
          local_node_id = all_node_map->getLocalElement(current_global_index);
          //debug print
          //std::cout << "STRAIN CHECK " << (*strain_vector_pass)(node_loop) << " " << all_node_strains(local_node_id, istrain) << std::endl;
          current_strain = (*strain_vector_pass)(node_loop);
          if(strain_max_flag){
            if(fabs(current_strain) > all_node_strains(local_node_id, istrain)){
              all_node_strains(local_node_id, istrain) = current_strain;

              if(map->isNodeGlobalElement(current_global_index)){
                local_node_id = map->getLocalElement(current_global_index);
                node_strains(local_node_id, istrain) = current_strain;
              }
            }
          }
          else{
            //debug print
            //std::cout << current_strain/(double)all_node_nconn(local_node_id,0) << all_node_nconn(local_node_id,0) << " ";
            if(map->isNodeGlobalElement(current_global_index)){
              local_node_id = map->getLocalElement(current_global_index);
              node_strains(local_node_id, istrain) += current_strain/(double)node_nconn(local_node_id,0);
            }
          }
        }
      }
    }
    
  }
    
  //debug print
  
  //std::ostream &out = std::cout;
  //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  //if(myrank==0)
  //*fos << "Local Node Strains :" << std::endl;
  //all_node_strains_distributed->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  //std::fflush(stdout);
  

}

/* ----------------------------------------------------------------------
   Compute the volume of each element; estimated with quadrature
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::compute_element_volumes(){
  //local number of uniquely assigned elements
  size_t nonoverlap_nelements = element_map->getNodeNumElements();
  //initialize memory for volume storage
  vec_array Element_Volumes("Element Volumes", nonoverlap_nelements, 1);
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  int num_dim = simparam->num_dim;
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam->num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  size_t local_node_id;
  LO ielem;
  GO global_element_index;

  real_t Jacobian;
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
  #CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
  CArray<real_t> legendre_nodes_1D(num_gauss_points);
  CArray<real_t> legendre_weights_1D(num_gauss_points);
  real_t pointer_quad_coordinate[num_dim];
  real_t pointer_quad_coordinate_weight[num_dim];
  real_t pointer_interpolated_point[num_dim];
  real_t pointer_JT_row1[num_dim];
  real_t pointer_JT_row2[num_dim];
  real_t pointer_JT_row3[num_dim];
  ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate,num_dim);
  ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight,num_dim);
  ViewCArray<real_t> interpolated_point(pointer_interpolated_point,num_dim);
  ViewCArray<real_t> JT_row1(pointer_JT_row1,num_dim);
  ViewCArray<real_t> JT_row2(pointer_JT_row2,num_dim);
  ViewCArray<real_t> JT_row3(pointer_JT_row3,num_dim);

  real_t pointer_basis_values[elem->num_basis()];
  real_t pointer_basis_derivative_s1[elem->num_basis()];
  real_t pointer_basis_derivative_s2[elem->num_basis()];
  real_t pointer_basis_derivative_s3[elem->num_basis()];
  ViewCArray<real_t> basis_values(pointer_basis_values,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2,elem->num_basis());
  ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(elem->num_basis(),num_dim);

  //initialize weights
  elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
  elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);

  //loop over elements and use quadrature rule to compute volume from Jacobian determinant
  for(int nonoverlapping_ielem = 0; nonoverlapping_ielem < nonoverlap_nelements; nonoverlapping_ielem++){
    global_element_index = element_map->getGlobalElement(nonoverlapping_ielem);
    ielem = all_element_map->getLocalElement(global_element_index);
    //debug print
    //std::cout << "ELEMENT INDEX IS: " << ielem << " " <<global_element_index << std::endl;

    //acquire set of nodes for this local element
    for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
      local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, node_loop));
      nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
      nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
      nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
      /*
      if(myrank==1&&nodal_positions(node_loop,2)>10000000){
        std::cout << " LOCAL MATRIX DEBUG ON TASK " << myrank << std::endl;
        std::cout << node_loop+1 <<" " << local_node_id <<" "<< host_elem_conn_(ielem, node_loop) << " "<< nodal_positions(node_loop,2) << std::endl;
        std::fflush(stdout);
      }
      */
      //std::cout << local_node_id << " " << nodes_in_elem(ielem, node_loop) << " " << nodal_positions(node_loop,0) << " " << nodal_positions(node_loop,1) << " "<< nodal_positions(node_loop,2) <<std::endl;
    }
    
    //initialize element volume
    Element_Volumes(nonoverlapping_ielem,0) = 0;
    
    if(Element_Types(ielem)==elements::elem_types::Hex8){
      direct_product_count = std::pow(num_gauss_points,num_dim);
    }
  
    //loop over quadrature points
    for(int iquad=0; iquad < direct_product_count; iquad++){

      //set current quadrature point
      if(num_dim==3) z_quad = iquad/(num_gauss_points*num_gauss_points);
      y_quad = (iquad % (num_gauss_points*num_gauss_points))/num_gauss_points;
      x_quad = iquad % num_gauss_points;
      quad_coordinate(0) = legendre_nodes_1D(x_quad);
      quad_coordinate(1) = legendre_nodes_1D(y_quad);
      if(num_dim==3)
      quad_coordinate(2) = legendre_nodes_1D(z_quad);

      //set current quadrature weight
      quad_coordinate_weight(0) = legendre_weights_1D(x_quad);
      quad_coordinate_weight(1) = legendre_weights_1D(y_quad);
      if(num_dim==3)
      quad_coordinate_weight(2) = legendre_weights_1D(z_quad);

      //compute shape functions at this point for the element type
      elem->basis(basis_values,quad_coordinate);

      //compute all the necessary coordinates and derivatives at this point

      //compute shape function derivatives
      elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
      elem->partial_eta_basis(basis_derivative_s2,quad_coordinate);
      elem->partial_mu_basis(basis_derivative_s3,quad_coordinate);

      //compute derivatives of x,y,z w.r.t the s,t,w isoparametric space needed by JT (Transpose of the Jacobian)
      //derivative of x,y,z w.r.t s
      JT_row1(0) = 0;
      JT_row1(1) = 0;
      JT_row1(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row1(0) += nodal_positions(node_loop,0)*basis_derivative_s1(node_loop);
        JT_row1(1) += nodal_positions(node_loop,1)*basis_derivative_s1(node_loop);
        JT_row1(2) += nodal_positions(node_loop,2)*basis_derivative_s1(node_loop);
      }

      //derivative of x,y,z w.r.t t
      JT_row2(0) = 0;
      JT_row2(1) = 0;
      JT_row2(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row2(0) += nodal_positions(node_loop,0)*basis_derivative_s2(node_loop);
        JT_row2(1) += nodal_positions(node_loop,1)*basis_derivative_s2(node_loop);
        JT_row2(2) += nodal_positions(node_loop,2)*basis_derivative_s2(node_loop);
      }

      //derivative of x,y,z w.r.t w
      JT_row3(0) = 0;
      JT_row3(1) = 0;
      JT_row3(2) = 0;
      for(int node_loop=0; node_loop < elem->num_basis(); node_loop++){
        JT_row3(0) += nodal_positions(node_loop,0)*basis_derivative_s3(node_loop);
        JT_row3(1) += nodal_positions(node_loop,1)*basis_derivative_s3(node_loop);
        JT_row3(2) += nodal_positions(node_loop,2)*basis_derivative_s3(node_loop);
        //debug print
        /*if(myrank==1&&nodal_positions(node_loop,2)*basis_derivative_s3(node_loop)<-10000000){
        std::cout << " ELEMENT VOLUME JACOBIAN DEBUG ON TASK " << myrank << std::endl;
        std::cout << node_loop+1 << " " << JT_row3(2) << " "<< nodal_positions(node_loop,2) <<" "<< basis_derivative_s3(node_loop) << std::endl;
        std::fflush(stdout);
        }*/
      }
    
    
    //compute the determinant of the Jacobian
    Jacobian = JT_row1(0)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
               JT_row1(1)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
               JT_row1(2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1));
    if(Jacobian<0) Jacobian = -Jacobian;
    
    Element_Volumes(nonoverlapping_ielem,0) += quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2)*Jacobian;
    }
  }

  //create global vector
  Global_Element_Volumes = Teuchos::rcp(new MV(element_map, Element_Volumes));

  std::ostream &out = std::cout;
  Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  //if(myrank==0)
  //*fos << "Global Element Volumes:" << std::endl;
  //Global_Element_Volumes->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
}

/* ----------------------------------------------------------------------
   Solve the FEA linear system
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::linear_solver_parameters(){
  if(simparam->direct_solver_flag){
    Linear_Solve_Params = Teuchos::rcp(new Teuchos::ParameterList("Amesos2"));
    auto superlu_params = Teuchos::sublist(Teuchos::rcpFromRef(*Linear_Solve_Params), "SuperLU_DIST");
    superlu_params->set("Equil", true);
    //superlu_params.set("Trans","TRANS","Whether to solve with A^T");
    //superlu_params.set("ColPerm","NATURAL","Use 'natural' ordering of columns");
  
  }
  else{
    Linear_Solve_Params = Teuchos::rcp(new Teuchos::ParameterList("MueLu"));
    std::string xmlFileName = "elasticity3D.xml";
    //std::string xmlFileName = "simple_test.xml";
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&(*Linear_Solve_Params)), *comm);
  }
}

/* ----------------------------------------------------------------------
   Solve the FEA linear system
------------------------------------------------------------------------- */

int Parallel_Nonlinear_Solver::solve(){
  //local variable for host view in the dual view
  const_host_vec_array node_coords = node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  //local variable for host view in the dual view
  const_host_vec_array Nodal_Forces = Global_Nodal_Forces->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  int num_dim = simparam->num_dim;
  int nodes_per_elem = max_nodes_per_element;
  int local_node_index, current_row, current_column;
  size_t reduced_index;
  int max_stride = 0;
  size_t access_index, row_access_index, row_counter;
  global_size_t reduced_row_count;
  GO global_index, global_dof_index;
  LO reduced_local_dof_index;

  //ROL::Ptr<MV> > x_ptr  = ROL::makePtr<MV>(2);
  //assumes ROL::Ptr<T> was compiled as Teuchos::RCP<T>
  ROL::Ptr<ROL::Vector<real_t> > x = ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO>>(node_coords_distributed);
  
  std::ostream &out = std::cout;
  Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

  //*fos << Amesos2::version() << std::endl << std::endl;

  bool printTiming   = true;
  bool verbose       = false;
  std::string filename("arc130.mtx");
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("print-timing","no-print-timing",&printTiming,"Print solver timing statistics");

  const size_t numVectors = 1;
  
  //construct global data for example; normally this would be from file input and distributed
  //according to the row map at that point
  global_size_t nrows = num_nodes*num_dim;
  
 
  //number of boundary conditions on this mpi rank
  global_size_t local_nboundaries = Number_DOF_BCS;
  size_t local_nrows_reduced = nlocal_nodes*num_dim - local_nboundaries;

  //obtain total number of boundary conditions on all ranks
  global_size_t global_nboundaries = 0;
  MPI_Allreduce(&local_nboundaries,&global_nboundaries,1,MPI_INT,MPI_SUM,world);
  global_size_t nrows_reduced = nrows - global_nboundaries;

  //Rebalance distribution of the global stiffness matrix rows here later since
  //rows and columns are being removed.

  //global_size_t *entries_per_row = new global_size_t[local_nrows];
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Reduced_Stiffness_Matrix_Strides(local_nrows_reduced,"Reduced_Stiffness_Matrix_Strides");
  row_pointers reduced_row_offsets_pass = row_pointers("reduced_row_offsets_pass", local_nrows_reduced+1);

  //init row_offsets
  for(int i=0; i < local_nrows_reduced+1; i++){
    reduced_row_offsets_pass(i) = 0;
  }
  
  //stores global indices belonging to this MPI rank from the non-reduced map corresponding to the reduced system
  CArrayKokkos<GO, array_layout, device_type, memory_traits> Free_Indices(local_nrows_reduced,"Free_Indices");
  reduced_index = 0;
  for(LO i=0; i < nlocal_nodes*num_dim; i++)
    if((Node_DOF_Boundary_Condition_Type(i)!=DISPLACEMENT_CONDITION)&&
       (Node_DOF_Boundary_Condition_Type(i)!=X_DISPLACEMENT_CONDITION)&&
       (Node_DOF_Boundary_Condition_Type(i)!=Y_DISPLACEMENT_CONDITION)&&
       (Node_DOF_Boundary_Condition_Type(i)!=Z_DISPLACEMENT_CONDITION)){
        Free_Indices(reduced_index) = local_dof_map->getGlobalElement(i);
        reduced_index++;
      }
    
  
  //compute reduced matrix strides
  for(LO i=0; i < local_nrows_reduced; i++){
    access_index = local_dof_map->getLocalElement(Free_Indices(i));
    reduced_row_count = 0;
    for(LO j=0; j < Stiffness_Matrix_Strides(access_index); j++){
      global_dof_index = DOF_Graph_Matrix(access_index,j);
      row_access_index = all_dof_map->getLocalElement(global_dof_index);
      if((Node_DOF_Boundary_Condition_Type(row_access_index)!=DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(row_access_index)!=X_DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(row_access_index)!=Y_DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(row_access_index)!=Z_DISPLACEMENT_CONDITION)){
        reduced_row_count++;
      }
    }
    Reduced_Stiffness_Matrix_Strides(i) = reduced_row_count;
  }
  
  RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits> Reduced_DOF_Graph_Matrix(Reduced_Stiffness_Matrix_Strides); //stores global dof indices
  RaggedRightArrayKokkos<LO, array_layout, device_type, memory_traits> Reduced_Local_DOF_Graph_Matrix(Reduced_Stiffness_Matrix_Strides); //stores local dof indices
  RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout> Reduced_Stiffness_Matrix(Reduced_Stiffness_Matrix_Strides);

  //template compatible row offsets (may change to shallow copy later if it works on device types etc.)
  row_pointers reduced_row_offsets = Reduced_DOF_Graph_Matrix.start_index_;
    for(int ipass = 0; ipass < local_nrows_reduced+1; ipass++){
      reduced_row_offsets_pass(ipass) = reduced_row_offsets(ipass);
    }
  
  //construct maps to define set of global indices for the reduced node set
  //stores indices that aren't contiguous due to removal of BCS
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > local_reduced_dof_original_map =
    Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(nrows_reduced,Free_Indices.get_kokkos_view(),0,comm) );

  //stores contiguous indices with an unbalanced local distribution
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > local_reduced_dof_map = 
    Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(nrows_reduced,local_nrows_reduced,0,comm));
  
  //dual view of the local global index to reduced global index map
  dual_vec_array dual_reduced_dof_original("dual_reduced_dof_original",local_nrows_reduced,1);

  //local variable for host view in the dual view
  host_vec_array reduced_dof_original = dual_reduced_dof_original.view_host();
  //notify that the host view is going to be modified
  dual_reduced_dof_original.modify_host();

  //set contents
  for(LO i=0; i < local_nrows_reduced; i++){
    reduced_dof_original(i,0) = local_reduced_dof_map->getGlobalElement(i);
  }
  
  //create a multivector where each local index entry stores the new reduced global index associated with each old global index
  Teuchos::RCP<MV> local_reduced_global_indices = Teuchos::rcp(new MV(local_reduced_dof_original_map, dual_reduced_dof_original));
  
  //construct map of all indices including ghosts for the reduced system
  //stores global indices belonging to this MPI rank and ghosts from the non-reduced map corresponding to the reduced system
  size_t all_nrows_reduced = local_nrows_reduced + nghost_nodes*num_dim;
  for(LO i=nlocal_nodes*num_dim; i < nall_nodes*num_dim; i++){
      if((Node_DOF_Boundary_Condition_Type(i)==DISPLACEMENT_CONDITION)||
         (Node_DOF_Boundary_Condition_Type(i)==X_DISPLACEMENT_CONDITION)||
         (Node_DOF_Boundary_Condition_Type(i)==Y_DISPLACEMENT_CONDITION)||
         (Node_DOF_Boundary_Condition_Type(i)==Z_DISPLACEMENT_CONDITION))
      all_nrows_reduced--;
  }
  
  CArrayKokkos<GO, array_layout, device_type, memory_traits> All_Free_Indices(all_nrows_reduced,"All_Free_Indices");

  //debug print
  /*
  if(myrank==0||myrank==4){
  std::cout << "DOF flags global :" << std::endl;
  std::cout << "Reduced DOF Graph Matrix on Rank " << myrank << std::endl;
  for(LO i=0; i < nall_nodes*num_dim; i++){
    std::cout << all_dof_map->getGlobalElement(i) << " " << Node_DOF_Boundary_Condition_Type(i) <<" ";   
    std::cout << std::endl;
  }
  std::fflush(stdout);
  }
  */
  
  reduced_index = 0;
  for(LO i=0; i < nall_nodes*num_dim; i++)
    if((Node_DOF_Boundary_Condition_Type(i)!=DISPLACEMENT_CONDITION)&&
       (Node_DOF_Boundary_Condition_Type(i)!=X_DISPLACEMENT_CONDITION)&&
       (Node_DOF_Boundary_Condition_Type(i)!=Y_DISPLACEMENT_CONDITION)&&
       (Node_DOF_Boundary_Condition_Type(i)!=Z_DISPLACEMENT_CONDITION)){
        All_Free_Indices(reduced_index) = all_dof_map->getGlobalElement(i);
        reduced_index++;
      }
  
  //construct map to define set of global indices for the reduced node set including ghosts
  //passing invalid forces the map to count the global elements
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > all_reduced_dof_original_map =
    Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),All_Free_Indices.get_kokkos_view(),0,comm) );

  //debug print
  /*
  if(myrank==0)
  *fos << "All reduced dof original indices :" << std::endl;
  local_reduced_dof_original_map->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;
  std::fflush(stdout);

  //debug print
  if(myrank==0)
  *fos << "All reduced dof original indices :" << std::endl;
  all_reduced_dof_original_map->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;
  std::fflush(stdout);
  */

  //communicate the new reduced global indices for ghost dof indices using the local information on other ranks through import
  //create import object using local node indices map and all indices map
  Tpetra::Import<LO, GO> importer(local_reduced_dof_original_map, all_reduced_dof_original_map);

  Teuchos::RCP<MV> all_reduced_global_indices = Teuchos::rcp(new MV(all_reduced_dof_original_map, 1));

  //comms to get reduced global indices for ghosts that are still free of BCs
  all_reduced_global_indices->doImport(*local_reduced_global_indices, importer, Tpetra::INSERT);

  const_host_vec_array all_reduced_global_indices_host = all_reduced_global_indices->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);

  //debug print
  /*
  if(myrank==0)
  *fos << "All reduced global indices :" << std::endl;
  all_reduced_global_indices->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;
  std::fflush(stdout);
  */
  
  //store the new global indices for the reduced matrix graph
  for(LO i=0; i < local_nrows_reduced; i++){
    access_index = local_dof_map->getLocalElement(Free_Indices(i));
    row_counter = 0;
    for(LO j=0; j < Stiffness_Matrix_Strides(access_index); j++){
      global_dof_index = DOF_Graph_Matrix(access_index,j);
      row_access_index = all_dof_map->getLocalElement(global_dof_index);
      if((Node_DOF_Boundary_Condition_Type(row_access_index)!=DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(row_access_index)!=X_DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(row_access_index)!=Y_DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(row_access_index)!=Z_DISPLACEMENT_CONDITION)){
        reduced_local_dof_index = all_reduced_dof_original_map->getLocalElement(global_dof_index);
        //std::cout << "REDUCED LOCAL INDEX ON TASK " << myrank << " is " << Reduced_Stiffness_Matrix_Strides(i) << Reduced_DOF_Graph_Matrix(i,row_counter++) << std::endl;
        Reduced_DOF_Graph_Matrix(i,row_counter++) = all_reduced_global_indices_host(reduced_local_dof_index,0);
      }
    }
  }
  
  //store reduced stiffness matrix values
  for(LO i=0; i < local_nrows_reduced; i++){
    access_index = local_dof_map->getLocalElement(Free_Indices(i));
    row_counter = 0;
    for(LO j=0; j < Stiffness_Matrix_Strides(access_index); j++){
      global_dof_index = DOF_Graph_Matrix(access_index,j);
      row_access_index = all_dof_map->getLocalElement(global_dof_index);
      if((Node_DOF_Boundary_Condition_Type(row_access_index)!=DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(row_access_index)!=X_DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(row_access_index)!=Y_DISPLACEMENT_CONDITION)&&
         (Node_DOF_Boundary_Condition_Type(row_access_index)!=Z_DISPLACEMENT_CONDITION)){
        Reduced_Stiffness_Matrix(i,row_counter++) = Stiffness_Matrix(access_index,j);
      }
    }
  }
  
  // create a Map for the reduced global stiffness matrix that is evenly distributed amongst mpi ranks
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > local_balanced_reduced_dof_map = 
    Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(nrows_reduced,0,comm));

  //build column map
  Teuchos::RCP<const Tpetra::Map<LO,GO,node_type> > colmap;
  const Teuchos::RCP<const Tpetra::Map<LO,GO,node_type> > dommap = local_reduced_dof_map;

  //debug print
  /*
  if(myrank==4){
  std::cout << "Reduced DOF Graph Matrix on Rank " << myrank << std::endl;
  for(LO i=0; i < local_nrows_reduced; i++){
    for(LO j=0; j < Reduced_Stiffness_Matrix_Strides(i); j++){
      std::cout << Reduced_DOF_Graph_Matrix(i,j) <<" ";
    }
    std::cout << std::endl;
  }
  }
  */

  Tpetra::Details::makeColMap<LO,GO,node_type>(colmap,dommap,Reduced_DOF_Graph_Matrix.get_kokkos_view(), nullptr);

  /*//debug print of reduced row offsets
  std::cout << " DEBUG PRINT FOR ROW OFFSETS" << std::endl;
  for(int debug = 0; debug < local_nrows_reduced+1; debug++)
  std::cout <<  reduced_row_offsets_pass(debug) << " ";
  std::cout << std::endl;
  //end debug print
  */

  //convert global indices to local indices using column map
  for(LO i=0; i < local_nrows_reduced; i++)
    for(LO j=0; j < Reduced_Stiffness_Matrix_Strides(i); j++)
      Reduced_Local_DOF_Graph_Matrix(i,j) = colmap->getLocalElement(Reduced_DOF_Graph_Matrix(i,j));
  
  Teuchos::RCP<MAT> unbalanced_A = Teuchos::rcp(new MAT(local_reduced_dof_map, colmap, reduced_row_offsets_pass,
                   Reduced_Local_DOF_Graph_Matrix.get_kokkos_view(), Reduced_Stiffness_Matrix.get_kokkos_view()));
  unbalanced_A->fillComplete();
  Teuchos::RCP<const_MAT> const_unbalanced_A(new const_MAT(*unbalanced_A));
  
  //This completes the setup for A matrix of the linear system

  //debug print of A matrix before balancing
  //*fos << "Reduced Stiffness Matrix :" << std::endl;
  //const_unbalanced_A->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;

  //communicate reduced stiffness matrix entries for better load balancing
  //create import object using the unbalanced map and the balanced map
  Tpetra::Import<LO, GO> matrix_importer(local_reduced_dof_map, local_balanced_reduced_dof_map);
  Teuchos::RCP<MAT> balanced_A = Tpetra::importAndFillCompleteCrsMatrix(const_unbalanced_A, matrix_importer, local_balanced_reduced_dof_map, local_balanced_reduced_dof_map);

  //debug print of map
  //if(myrank==0)
  //*fos << "Reduced DOF Map :" << std::endl;
  //local_balanced_reduced_dof_map->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  //std::fflush(stdout);

  //debug print of A matrix after balancing
  //if(myrank==0)
  //*fos << "Reduced Stiffness Matrix :" << std::endl;
  //balanced_A->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  //std::fflush(stdout);
  
  // Create random X vector
  size_t balanced_local_nrows = local_balanced_reduced_dof_map->getNodeNumElements();
  //vec_array Xview_pass = vec_array("Xview_pass", balanced_local_nrows, 1);
  //Xview_pass.assign_data(Xview.pointer());
  Teuchos::RCP<MV> X = Teuchos::rcp(new MV(local_balanced_reduced_dof_map, 1));
  //return !EXIT_SUCCESS;
  //X->randomize();

  // Create Kokkos view of RHS B vector (Force Vector)  
  vec_array Bview_pass = vec_array("Bview_pass", local_nrows_reduced,1);

  //set bview to force vector data
  for(LO i=0; i < local_nrows_reduced; i++){
    access_index = local_dof_map->getLocalElement(Free_Indices(i));
    Bview_pass(i,0) = Nodal_Forces(access_index,0);
  }
  
  Teuchos::RCP<MV> unbalanced_B = Teuchos::rcp(new MV(local_reduced_dof_map, Bview_pass));
  
  //import object to rebalance force vector
  Tpetra::Import<LO, GO> Bvec_importer(local_reduced_dof_map, local_balanced_reduced_dof_map);

  Teuchos::RCP<MV> balanced_B = Teuchos::rcp(new MV(local_balanced_reduced_dof_map, 1));
  
  //comms to rebalance force vector
  balanced_B->doImport(*unbalanced_B, Bvec_importer, Tpetra::INSERT);
  
  //if(myrank==0)
  //*fos << "RHS :" << std::endl;
  //balanced_B->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;

  //debug print
  //if(update_count==42){
    //Tpetra::MatrixMarket::Writer<MAT> market_writer();
    //Tpetra::MatrixMarket::Writer<MAT>::writeSparseFile("A_matrix.txt", *balanced_A, "A_matrix", "Stores stiffness matrix values");
  //}
  //return !EXIT_SUCCESS;
  // Create solver interface to KLU2 with Amesos2 factory method
  //std::cout << "Creating solver" << std::endl <<std::flush;
  if(simparam->direct_solver_flag){
    // Before we do anything, check that KLU2 is enabled
    if( !Amesos2::query("SuperLUDist") ){
      std::cerr << "SuperLUDist not enabled in this run.  Exiting..." << std::endl;
      return EXIT_SUCCESS;        // Otherwise CTest will pick it up as
                                  // failure, which it isn't really
    }
    Teuchos::RCP<Amesos2::Solver<MAT,MV>> solver = Amesos2::create<MAT,MV>("SuperLUDist", balanced_A, X, balanced_B);
    //Teuchos::RCP<Amesos2::Solver<MAT,MV>> solver = Amesos2::create<MAT,MV>("SuperLUDist", balanced_A, X, balanced_B);
    //Teuchos::RCP<Amesos2::Solver<MAT,MV>> solver = Amesos2::create<MAT,MV>("KLU2", balanced_A, X, balanced_B);
  
    solver->setParameters( Teuchos::rcpFromRef(*Linear_Solve_Params) );
    //declare non-contiguous map
    //Create a Teuchos::ParameterList to hold solver parameters
    //Teuchos::ParameterList amesos2_params("Amesos2");
    //amesos2_params.sublist("KLU2").set("IsContiguous", false, "Are GIDs Contiguous");
    //solver->setParameters( Teuchos::rcpFromRef(amesos2_params) );
  
    //Solve the system
    //std::cout << "BEFORE LINEAR SOLVE" << std::endl << std::flush;
    solver->symbolicFactorization().numericFactorization().solve();
    //debug print of displacements
    //std::ostream &out = std::cout;
    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    //if(myrank==0)
    //*fos << "balanced_X: " << update_count << std::endl;
    //X->describe(*fos,Teuchos::VERB_EXTREME);
    //*fos << std::endl;
    //std::fflush(stdout);
    //std::cout << "AFTER LINEAR SOLVE" << std::endl<< std::flush;
  }
  else{
    //dimension of the nullspace for linear elasticity
    int nulldim = 6;
    if(num_dim == 2) nulldim = 3;
    using impl_scalar_type =
      typename Kokkos::Details::ArithTraits<real_t>::val_type;
    using mag_type = typename Kokkos::ArithTraits<impl_scalar_type>::mag_type;
    // Instead of checking each time for rank, create a rank 0 stream
    Teuchos::RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream& out = *fancy;

    Teuchos::RCP<Xpetra::MultiVector<real_t,LO,GO,node_type>> xbalanced_B = Teuchos::rcp(new Xpetra::TpetraMultiVector<real_t,LO,GO,node_type>(balanced_B));
    Teuchos::RCP<Xpetra::MultiVector<real_t,LO,GO,node_type>> xX = Teuchos::rcp(new Xpetra::TpetraMultiVector<real_t,LO,GO,node_type>(X));
    //Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > reduced_node_map = 
    //Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(nrows_reduced/num_dim,0,comm));

    //set coordinates vector
    Teuchos::RCP<MV> unbalanced_coordinates_distributed = Teuchos::rcp(new MV(local_reduced_dof_map, num_dim));
    //loop through dofs and set coordinates, duplicated for each dim to imitate MueLu example for now (no idea why this was done that way)

    host_vec_array unbalanced_coordinates_view = unbalanced_coordinates_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    int dim_index;
    real_t node_x, node_y, node_z;
    //set coordinates components
      for(LO i=0; i < local_nrows_reduced; i++){
        global_dof_index = Free_Indices(i);
        global_index = Free_Indices(i)/num_dim;
        dim_index = global_dof_index % num_dim;
        access_index = map->getLocalElement(global_index);
        node_x = node_coords(access_index, 0);
        node_y = node_coords(access_index, 1);
        node_z = node_coords(access_index, 2);

        unbalanced_coordinates_view(i,0) = node_x;
        unbalanced_coordinates_view(i,1) = node_y;
        unbalanced_coordinates_view(i,2) = node_z;
      }// for
    
    //balanced coordinates vector
    Teuchos::RCP<MV> tcoordinates = Teuchos::rcp(new MV(local_balanced_reduced_dof_map, num_dim));
    //rebalance coordinates vector
    tcoordinates->doImport(*unbalanced_coordinates_distributed, Bvec_importer, Tpetra::INSERT);
    Teuchos::RCP<Xpetra::MultiVector<real_t,LO,GO,node_type>> coordinates = Teuchos::rcp(new Xpetra::TpetraMultiVector<real_t,LO,GO,node_type>(tcoordinates));
    
    //nullspace vector
    Teuchos::RCP<MV> unbalanced_nullspace_distributed = Teuchos::rcp(new MV(local_reduced_dof_map, nulldim));
    //set nullspace components
    //init
    unbalanced_nullspace_distributed->putScalar(0);
    //loop through dofs and compute nullspace components for each
    host_vec_array unbalanced_nullspace_view = unbalanced_nullspace_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);

    //compute center
    // Calculate center
	  real_t cx = tcoordinates->getVector(0)->meanValue();
	  real_t cy = tcoordinates->getVector(1)->meanValue();
    real_t cz;
    if(num_dim==3)
	    cz = tcoordinates->getVector(2)->meanValue();

    if(num_dim==3){
      for(LO i=0; i < local_nrows_reduced; i++){
        global_dof_index = Free_Indices(i);
        global_index = Free_Indices(i)/num_dim;
        dim_index = global_dof_index % num_dim;
        access_index = map->getLocalElement(global_index);
        node_x = node_coords(access_index, 0);
        node_y = node_coords(access_index, 1);
        node_z = node_coords(access_index, 2);
        //set translational component
        unbalanced_nullspace_view(i,dim_index) = 1;
        //set rotational components
        if(dim_index==0){
          unbalanced_nullspace_view(i,3) = -node_y + cy;
          unbalanced_nullspace_view(i,5) = node_z - cz;
        }
        if(dim_index==1){
          unbalanced_nullspace_view(i,3) = node_x - cx;
          unbalanced_nullspace_view(i,4) = -node_z + cz;
        }
        if(dim_index==2){
          unbalanced_nullspace_view(i,4) = node_y - cy;
          unbalanced_nullspace_view(i,5) = -node_x + cx;
        }
      }// for
    }
    else{
      for(LO i=0; i < local_nrows_reduced; i++){
        global_dof_index = Free_Indices(i);
        global_index = Free_Indices(i)/num_dim;
        dim_index = global_dof_index % num_dim;
        access_index = map->getLocalElement(global_index);
        node_x = node_coords(access_index, 0);
        node_y = node_coords(access_index, 1);
        //set translational component
        unbalanced_nullspace_view(i,dim_index) = 1;
        //set rotational components
        if(dim_index==0){
          unbalanced_nullspace_view(i,3) = -node_y + cy;
        }
        if(dim_index==1){
          unbalanced_nullspace_view(i,3) = node_x - cx;
        }
      }// for
    }
    
    //balanced nullspace vector
    Teuchos::RCP<MV> tnullspace = Teuchos::rcp(new MV(local_balanced_reduced_dof_map, nulldim));
    //rebalance nullspace vector
    tnullspace->doImport(*unbalanced_nullspace_distributed, Bvec_importer, Tpetra::INSERT);
    Teuchos::RCP<Xpetra::MultiVector<real_t,LO,GO,node_type>> nullspace = Teuchos::rcp(new Xpetra::TpetraMultiVector<real_t,LO,GO,node_type>(tnullspace));

    //normalize components
    Kokkos::View<mag_type*, Kokkos::HostSpace> norms2("norms2", nulldim);
    tnullspace->norm2(norms2);
    Kokkos::View<impl_scalar_type*, device_type> scaling_values("scaling_values", nulldim);
    for (int i = 0; i < nulldim; i++)
        scaling_values(i) = norms2(0) / norms2(i);
    tnullspace->scale(scaling_values);

    Teuchos::RCP<Xpetra::MultiVector<real_t,LO,GO,node_type>> material = Teuchos::null;
    Teuchos::RCP<Xpetra::CrsMatrix<real_t,LO,GO,node_type>> xbalanced_A = Teuchos::rcp(new Xpetra::TpetraCrsMatrix<real_t,LO,GO,node_type>(balanced_A));
    Teuchos::RCP<Xpetra::Matrix<real_t,LO,GO,node_type>> xwrap_balanced_A = Teuchos::rcp(new Xpetra::CrsMatrixWrap<real_t,LO,GO,node_type>(xbalanced_A));
    //xwrap_balanced_A->SetFixedBlockSize(1);

    //randomize initial vector
    xX->setSeed(100);
    xX->randomize();
    
    
    //debug print
    //if(myrank==0)
    //*fos << "Xpetra A matrix :" << std::endl;
    //xX->describe(*fos,Teuchos::VERB_EXTREME);
    //*fos << std::endl;
    //std::fflush(stdout);
    
    int num_iter = 2000;
    double solve_tol = 1e-12;
    int cacheSize = 1000;
    std::string solveType         = "belos";
    std::string belosType         = "cg";
    // =========================================================================
    // Preconditioner construction
    // =========================================================================
    //bool useML   = Linear_Solve_Params->isParameter("use external multigrid package") && (Linear_Solve_Params->get<std::string>("use external multigrid package") == "ml");
    //out<<"*********** MueLu ParameterList ***********"<<std::endl;
    //out<<*Linear_Solve_Params;
    //out<<"*******************************************"<<std::endl;
    
    Teuchos::RCP<MueLu::Hierarchy<real_t,LO,GO,node_type>> H;
    Teuchos::RCP<Xpetra::Operator<real_t,LO,GO,node_type>> Prec;
    {
      comm->barrier();
      //PreconditionerSetup(A,coordinates,nullspace,material,paramList,false,false,useML,0,H,Prec);
      PreconditionerSetup(xwrap_balanced_A,coordinates,nullspace,material,*Linear_Solve_Params,false,false,false,0,H,Prec);
      comm->barrier();
      //H->Write(-1, -1);
      //H->describe(*fos,Teuchos::VERB_EXTREME);
    }

    // =========================================================================
    // System solution (Ax = b)
    // =========================================================================
    {
      comm->barrier();
      SystemSolve(xwrap_balanced_A,xX,xbalanced_B,H,Prec,out,solveType,belosType,false,false,false,cacheSize,0,true,true,num_iter,solve_tol);
      comm->barrier();
    }
  }
  //return !EXIT_SUCCESS;
  //timing statistics for LU solver
  //solver->printTiming(*fos);
  
  //Print solution vector
  //print allocation of the solution vector to check distribution
  
  //if(myrank==0)
  //*fos << "Solution:" << std::endl;
  //X->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  

  //communicate solution on reduced map to the all node map vector for post processing of strain etc.
  //intermediate storage on the unbalanced reduced system
  Teuchos::RCP<MV> reduced_node_displacements_distributed = Teuchos::rcp(new MV(local_reduced_dof_map, 1));
  //create import object using local node indices map and all indices map
  Tpetra::Import<LO, GO> reduced_displacement_importer(local_balanced_reduced_dof_map, local_reduced_dof_map);

  //comms to get displacements on reduced unbalanced displacement vector
  reduced_node_displacements_distributed->doImport(*X, reduced_displacement_importer, Tpetra::INSERT);

  //populate node displacement multivector on the local dof map
  const_host_vec_array reduced_node_displacements_host = reduced_node_displacements_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  host_vec_array node_displacements_host = node_displacements_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);

  for(int init = 0; init < local_dof_map->getNodeNumElements(); init++)
    node_displacements_host(init,0) = 0;

  for(LO i=0; i < local_nrows_reduced; i++){
    access_index = local_dof_map->getLocalElement(Free_Indices(i));
    node_displacements_host(access_index,0) = reduced_node_displacements_host(i,0);
  }
  
  //import for displacement of ghosts
  Tpetra::Import<LO, GO> ghost_displacement_importer(local_dof_map, all_dof_map);

  //comms to get displacements on all node map
  all_node_displacements_distributed->doImport(*node_displacements_distributed, ghost_displacement_importer, Tpetra::INSERT);

  //if(myrank==0)
  //*fos << "All displacements :" << std::endl;
  //all_node_displacements_distributed->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  
  return !EXIT_SUCCESS;
}

/* ----------------------------------------------------------------------
   Return the CPU time for the current process in seconds very
   much in the same way as MPI_Wtime() returns the wall time.
------------------------------------------------------------------------- */

double Parallel_Nonlinear_Solver::CPU_Time()
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

void Parallel_Nonlinear_Solver::init_clock(){
  double current_cpu = 0;
  initial_CPU_time = CPU_Time();
}
