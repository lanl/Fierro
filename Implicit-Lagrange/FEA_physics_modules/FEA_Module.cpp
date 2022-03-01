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
#include "ROL_Elementwise_Reduce.hpp"
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
  hessvec_count = update_count = 0;
  file_index = 0;
  linear_solve_time = hessvec_time = hessvec_linear_time = 0;

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

  //boundary condition data
  current_bdy_id = 0;

  //boundary condition flags
  body_force_flag = gravity_flag = thermal_flag = electric_flag = false;

  //preconditioner construction
  Hierarchy_Constructed = false;

  //Trilinos output stream
  std::ostream &out = std::cout;
  fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  (*fos).setOutputToRootOnly(0);
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
    else
      read_mesh_ensight(argv[1]);

    init_maps();
    
    std::cout << "Num elements on process " << myrank << " = " << rnum_elem << std::endl;
    
    //initialize timing
    if(simparam->report_runtime_flag)
    init_clock();
    
    // ---- Find Boundaries on mesh ---- //
    init_boundaries();

    //set boundary conditions
    generate_bcs();

    //set applied loading conditions
    generate_applied_loads();

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
    std::cout << " RUNTIME OF CODE ON TASK " << myrank << " is "<< current_cpu-initial_CPU_time << " update solve time " << linear_solve_time << " hess solve time " << hessvec_linear_time <<std::endl;
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
      std::cout << "Total number of hessvec counts " << hessvec_count <<std::endl;
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
  if(myrank == 0)
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
  //problem->addConstraint("equality Constraint 2",eq_constraint2,constraint_mul2);
  //problem->addConstraint("equality Constraint 3",eq_constraint3,constraint_mul3);
  //problem->addLinearConstraint("Equality Constraint",eq_constraint,constraint_mul);
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

/* ----------------------------------------------------------------------------
   Initialize sets of element boundary surfaces and arrays for input conditions
------------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::init_boundaries(){
  int num_boundary_sets = simparam->NB;
  int num_surface_force_sets = simparam->NBSF;
  int num_surface_disp_sets = simparam->NBD;
  int num_dim = simparam->num_dim;

  // build boundary mesh patches
  if(myrank == 0)
    std::cout << "Starting boundary patch setup" << std::endl <<std::flush;
  Get_Boundary_Patches();
  //std::cout << "Done with boundary patch setup" << std::endl <<std::flush;
  std::cout << "number of boundary patches on task " << myrank << " = " << nboundary_patches << std::endl;
  
  // set the number of boundary sets
  if(myrank == 0)
    std::cout << "building boundary sets " << std::endl;
  
  init_boundary_sets(num_boundary_sets);
  Boundary_Condition_Type_List = CArrayKokkos<int, array_layout, HostSpace, memory_traits>(num_boundary_sets); 
  Boundary_Surface_Force_Densities = CArrayKokkos<real_t, array_layout, HostSpace, memory_traits>(num_surface_force_sets,3);
  Boundary_Surface_Displacements = CArrayKokkos<real_t, array_layout, HostSpace, memory_traits>(num_surface_disp_sets,3);

  //initialize
  for(int ibdy=0; ibdy < num_boundary_sets; ibdy++) Boundary_Condition_Type_List(ibdy) = NONE;

  //allocate nodal data
  Node_DOF_Boundary_Condition_Type = CArrayKokkos<int, array_layout, device_type, memory_traits>(nall_nodes*num_dim, "Node_DOF_Boundary_Condition_Type");
  Node_DOF_Displacement_Boundary_Conditions = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(nall_nodes*num_dim);
  Node_DOF_Force_Boundary_Conditions = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(nall_nodes*num_dim);

  //initialize
  for(int init=0; init < nall_nodes*num_dim; init++)
    Node_DOF_Boundary_Condition_Type(init) = NONE;
}

/* ----------------------------------------------------------------------
   Assign sets of element boundary surfaces corresponding to user BCs
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::generate_bcs(){
  int num_dim = simparam->num_dim;
  int bdy_set_id;
  int surf_disp_set_id = 0;
  int bc_tag;
  real_t value;
  real_t fix_limits[4];

  // tag the z=0 plane,  (Direction, value, bdy_set)
  *fos << "tagging z = 0 " << std::endl;
  bc_tag = 2;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 0.0 * simparam->unit_scaling;
  fix_limits[0] = fix_limits[2] = 4;
  fix_limits[1] = fix_limits[3] = 6;
  bdy_set_id = current_bdy_id++;
  //tag_boundaries(bc_tag, value, bdy_set_id, fix_limits);
  tag_boundaries(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = DISPLACEMENT_CONDITION;
  Boundary_Surface_Displacements(surf_disp_set_id,0) = 0;
  Boundary_Surface_Displacements(surf_disp_set_id,1) = 0;
  Boundary_Surface_Displacements(surf_disp_set_id,2) = 0;
  surf_disp_set_id++;
    
  *fos << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
  *fos << std::endl;
 /*
  // tag the y=10 plane,  (Direction, value, bdy_set)
  std::cout << "tagging y = 10 " << std::endl;
  bc_tag = 1;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 10.0 * simparam->unit_scaling;
  fix_limits[0] = fix_limits[2] = 4;
  fix_limits[1] = fix_limits[3] = 6;
  bdy_set_id = current_bdy_id++;
  //tag_boundaries(bc_tag, value, bdy_set_id, fix_limits);
  tag_boundaries(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = DISPLACEMENT_CONDITION;
  Boundary_Surface_Displacements(surf_disp_set_id,0) = 0;
  Boundary_Surface_Displacements(surf_disp_set_id,1) = 0;
  Boundary_Surface_Displacements(surf_disp_set_id,2) = 0;
  surf_disp_set_id++;
    
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
  std::cout << std::endl;

  // tag the x=10 plane,  (Direction, value, bdy_set)
  std::cout << "tagging y = 10 " << std::endl;
  bc_tag = 0;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 10.0 * simparam->unit_scaling;
  fix_limits[0] = fix_limits[2] = 4;
  fix_limits[1] = fix_limits[3] = 6;
  bdy_set_id = current_bdy_id++;
  //tag_boundaries(bc_tag, value, bdy_set_id, fix_limits);
  tag_boundaries(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = DISPLACEMENT_CONDITION;
  Boundary_Surface_Displacements(surf_disp_set_id,0) = 0;
  Boundary_Surface_Displacements(surf_disp_set_id,1) = 0;
  Boundary_Surface_Displacements(surf_disp_set_id,2) = 0;
  surf_disp_set_id++;
    
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
  std::cout << std::endl;
 
  // tag the +z beam plane,  (Direction, value, bdy_set)
  std::cout << "tagging z = 100 " << std::endl;
  bc_tag = 2;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 100.0 * simparam->unit_scaling;
  //real_t fix_limits[4];
  fix_limits[0] = fix_limits[2] = 4;
  fix_limits[1] = fix_limits[3] = 6;
  bdy_set_id = current_bdy_id++;
  //tag_boundaries(bc_tag, value, bdy_set_id, fix_limits);
  tag_boundaries(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = DISPLACEMENT_CONDITION;
  Boundary_Surface_Displacements(surf_disp_set_id,0) = 0;
  Boundary_Surface_Displacements(surf_disp_set_id,1) = 0;
  Boundary_Surface_Displacements(surf_disp_set_id,2) = 0;
  surf_disp_set_id++;
    
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
  std::cout << std::endl;
  
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

  //Tag nodes for Boundary conditions such as displacements
  Displacement_Boundary_Conditions();
} // end generate_bcs

/* ----------------------------------------------------------------------
   Assign sets of element boundary surfaces corresponding to user BCs
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::generate_applied_loads(){
  int num_dim = simparam->num_dim;
  int bdy_set_id;
  int surf_force_set_id = 0;
  int bc_tag;
  real_t value;
  
  //Surface Forces Section

  /*
  std::cout << "tagging z = 2 Force " << std::endl;
  bc_tag = 2;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 2 * simparam->unit_scaling;
  //value = 2;
  bdy_set_id = current_bdy_id++;
  //find boundary patches this BC corresponds to
  tag_boundaries(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = SURFACE_LOADING_CONDITION;
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
  Boundary_Condition_Type_List(bdy_set_id) = SURFACE_LOADING_CONDITION;
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
  Boundary_Condition_Type_List(bdy_set_id) = SURFACE_LOADING_CONDITION;
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
  Boundary_Condition_Type_List(bdy_set_id) = SURFACE_LOADING_CONDITION;
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
  Boundary_Condition_Type_List(bdy_set_id) = SURFACE_LOADING_CONDITION;
  Boundary_Surface_Force_Densities(surf_force_set_id,0) = 0;
  Boundary_Surface_Force_Densities(surf_force_set_id,1) = 1/simparam->unit_scaling/simparam->unit_scaling;
  Boundary_Surface_Force_Densities(surf_force_set_id,2) = 0;
  surf_force_set_id++;
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
  std::cout << std::endl;
  */

  /*
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
  Boundary_Condition_Type_List(bdy_set_id) = SURFACE_LOADING_CONDITION;
  Boundary_Surface_Force_Densities(surf_force_set_id,0) = 10/simparam->unit_scaling/simparam->unit_scaling;
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
  Boundary_Condition_Type_List(bdy_set_id) = SURFACE_LOADING_CONDITION;
  Boundary_Surface_Force_Densities(surf_force_set_id,0) = -10/simparam->unit_scaling/simparam->unit_scaling;
  Boundary_Surface_Force_Densities(surf_force_set_id,1) = 0;
  Boundary_Surface_Force_Densities(surf_force_set_id,2) = 0;
  surf_force_set_id++;
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
  std::cout << std::endl;
  */
  
  *fos << "tagging beam +z force " << std::endl;
  bc_tag = 2;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  //value = 0;
  value = 100;
  bdy_set_id = current_bdy_id++;
  //find boundary patches this BC corresponds to
  tag_boundaries(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = SURFACE_LOADING_CONDITION;
  Boundary_Surface_Force_Densities(surf_force_set_id,0) = 0.5/simparam->unit_scaling/simparam->unit_scaling;
  Boundary_Surface_Force_Densities(surf_force_set_id,1) = 0;
  Boundary_Surface_Force_Densities(surf_force_set_id,2) = 0;
  surf_force_set_id++;
  *fos << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
  *fos << std::endl;
  
  /*
  std::cout << "tagging y = 2 " << std::endl;
  bc_tag = 1;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 2.0;
  bdy_set_id = 4;
  mesh->tag_bdys(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = SURFACE_LOADING_CONDITION;
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << mesh->num_bdy_patches_in_set(bdy_set_id) << std::endl;
  std::cout << std::endl;

  std::cout << "tagging z = 2 " << std::endl;
  bc_tag = 2;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 2.0;
  bdy_set_id = 5;
  mesh->tag_bdys(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = SURFACE_LOADING_CONDITION;
  std::cout << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << mesh->num_bdy_patches_in_set(bdy_set_id) << std::endl;
  std::cout << std::endl;
  */
  
  //Body Forces Section

  //apply gravity
  gravity_flag = simparam->gravity_flag;
  gravity_vector = simparam->gravity_vector;

  if(electric_flag||gravity_flag||thermal_flag) body_force_flag = true;

}

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

void Parallel_Nonlinear_Solver::tag_boundaries(int bc_tag, real_t val, int bdy_set, real_t *patch_limits){
  
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
    
  *fos << " tagged boundary patches " << std::endl;
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

  //filter small negative numbers (that should have been 0 from cancellation) from floating point error
  /*
  for (int idof = 0; idof < num_dim*nlocal_nodes; idof++){
    for (int istride = 0; istride < Stiffness_Matrix_Strides(idof); istride++){
      if(Stiffness_Matrix(idof,istride)<0.000000001*simparam->Elastic_Modulus*DENSITY_EPSILON||Stiffness_Matrix(idof,istride)>-0.000000001*simparam->Elastic_Modulus*DENSITY_EPSILON)
      Stiffness_Matrix(idof,istride) = 0;
      //debug print
      //std::cout << "{" <<istride + 1 << "," << DOF_Graph_Matrix(idof,istride) << "} ";
    }
    //debug print
    //std::cout << std::endl;
  }
  */
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
   Construct the global applied force vector
------------------------------------------------------------------------- */

void Parallel_Nonlinear_Solver::assemble_vector(){
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  //local variable for host view in the dual view
  host_vec_array Nodal_Forces = Global_Nodal_Forces->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  const_host_vec_array Element_Densities;
  //local variable for host view of densities from the dual view
  bool nodal_density_flag = simparam->nodal_density_flag;
  const_host_vec_array all_node_densities;
  if(nodal_density_flag)
  all_node_densities = all_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  else
  Element_Densities = Global_Element_Densities->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
  int num_bdy_patches_in_set;
  size_t patch_id;
  GO current_node_index;
  LO local_node_id;
  LO node_id;
  int num_boundary_sets = num_boundary_conditions;
  int surface_force_set_id = 0;
  int num_dim = simparam->num_dim;
  int nodes_per_elem = max_nodes_per_element;
  int num_gauss_points = simparam->num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  int current_element_index, local_surface_id, surf_dim1, surf_dim2;
  //CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
  //CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
  CArray<real_t> legendre_nodes_1D(num_gauss_points);
  CArray<real_t> legendre_weights_1D(num_gauss_points);
  real_t pointer_quad_coordinate[num_dim];
  real_t pointer_quad_coordinate_weight[num_dim];
  real_t pointer_interpolated_point[num_dim];
  ViewCArray<real_t> quad_coordinate(pointer_quad_coordinate,num_dim);
  ViewCArray<real_t> quad_coordinate_weight(pointer_quad_coordinate_weight,num_dim);
  ViewCArray<real_t> interpolated_point(pointer_interpolated_point,num_dim);
  real_t force_density[3], wedge_product, Jacobian, current_density, weight_multiply;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Surface_Nodes;
  
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> JT_row1(num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> JT_row2(num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> JT_row3(num_dim);

  real_t pointer_basis_values[nodes_per_elem];
  real_t pointer_basis_derivative_s1[nodes_per_elem];
  real_t pointer_basis_derivative_s2[nodes_per_elem];
  real_t pointer_basis_derivative_s3[nodes_per_elem];
  ViewCArray<real_t> basis_values(pointer_basis_values,nodes_per_elem);
  ViewCArray<real_t> basis_derivative_s1(pointer_basis_derivative_s1,nodes_per_elem);
  ViewCArray<real_t> basis_derivative_s2(pointer_basis_derivative_s2,nodes_per_elem);
  ViewCArray<real_t> basis_derivative_s3(pointer_basis_derivative_s3,nodes_per_elem);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_positions(nodes_per_elem,num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(nodes_per_elem);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> surf_basis_derivative_s1(nodes_per_elem,num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> surf_basis_derivative_s2(nodes_per_elem,num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> surf_basis_values(nodes_per_elem,num_dim);

   //force vector initialization
  for(int i=0; i < num_dim*nlocal_nodes; i++)
    Nodal_Forces(i,0) = 0;

  /*Loop through boundary sets and check if they apply surface forces.
  These sets can have overlapping nodes since applied loading conditions
  are assumed to be additive*/
  for(int iboundary = 0; iboundary < num_boundary_sets; iboundary++){
    if(Boundary_Condition_Type_List(iboundary)!=SURFACE_LOADING_CONDITION) continue;
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
        local_node_id = all_node_map->getLocalElement(current_node_index);
        nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
        nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
        nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
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

    //apply body forces
    if(body_force_flag){
      //initialize quadrature data
      elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
      elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);
      direct_product_count = std::pow(num_gauss_points,num_dim);

      for(size_t ielem = 0; ielem < rnum_elem; ielem++){

      //acquire set of nodes and nodal displacements for this local element
      for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
        local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, node_loop));
        nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
        nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
        nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
        if(nodal_density_flag) nodal_density(node_loop) = all_node_densities(local_node_id,0);
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
      else
      quad_coordinate_weight(2) = 1;
      weight_multiply = quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2);

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
      //look up element material properties at this point as a function of density
      Body_Force(ielem, current_density, force_density);
    
      //evaluate contribution to force vector component
      for(int ibasis=0; ibasis < nodes_per_elem; ibasis++){
        if(!map->isNodeGlobalElement(nodes_in_elem(ielem, ibasis))) continue;
        local_node_id = map->getLocalElement(nodes_in_elem(ielem, ibasis));

        for(int idim = 0; idim < num_dim; idim++){
            if(force_density[idim]!=0)
            Nodal_Forces(num_dim*local_node_id + idim,0) += Jacobian*weight_multiply*force_density[idim]*basis_values(ibasis);
          }
        }
      }
    }
  }
    //debug print of force vector
    /*
    std::cout << "---------FORCE VECTOR-------------" << std::endl;
    for(int iforce=0; iforce < num_nodes*num_dim; iforce++)
      std::cout << " DOF: "<< iforce+1 << ", "<< Nodal_Forces(iforce) << std::endl;
    */

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

  if(body_force_flag)
    assemble_vector();
  
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
