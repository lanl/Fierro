/**********************************************************************************************
 © 2020. Triad National Security, LLC. All rights reserved.
 This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
 National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
 Department of Energy/National Nuclear Security Administration. All rights in the program are
 reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
 Security Administration. The Government is granted for itself and others acting on its behalf a
 nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
 derivative works, distribute copies to the public, perform publicly and display publicly, and
 to permit others to do so.
 This program is open source under the BSD-3 License.
 Redistribution and use in source and binary forms, with or without modification, are permitted
 provided that the following conditions are met:
 
 1.  Redistributions of source code must retain the above copyright notice, this list of
 conditions and the following disclaimer.
 
 2.  Redistributions in binary form must reproduce the above copyright notice, this list of
 conditions and the following disclaimer in the documentation and/or other materials
 provided with the distribution.
 
 3.  Neither the name of the copyright holder nor the names of its contributors may be used
 to endorse or promote products derived from this software without specific prior
 written permission.
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **********************************************************************************************/

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

#include "Tpetra_Import.hpp"

#include "elements.h"
#include "swage.h"
#include "matar.h"
#include "utilities.h"
#include "node_combination.h"
#include "Simulation_Parameters/Simulation_Parameters_Implicit.h"
#include "Simulation_Parameters/FEA_Module/FEA_Module_Headers.h"
#include "Implicit_Solver.h"
#include "FEA_Modules_Headers.h"

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
#include "Implicit_Optimization_Function_Headers.h"


#define BUFFER_LINES 20000
#define MAX_WORD 30
#define MAX_ELEM_NODES 32
#define STRAIN_EPSILON 0.000000001
#define BC_EPSILON 1.0e-8

using namespace utils;

/*

Swage is a reference to a swage block used in blacksmithing.  
Its a large metal block that has multiple shaps carved into 
each surface to use for hammering metal into to form it. 

*/

Implicit_Solver::Implicit_Solver(Simulation_Parameters_Implicit params) : Solver(params) {
  simparam = params;
  //create ref element object
  ref_elem = std::make_shared<elements::ref_element>(elements::ref_element());
  //create mesh objects
  //init_mesh = new swage::mesh_t(simparam);
  //mesh = new swage::mesh_t(simparam);

  element_select = std::make_shared<elements::element_selector>(elements::element_selector());
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

  //file readin parameter
  active_node_ordering_convention = IJK;
}

Implicit_Solver::~Implicit_Solver(){
  if(myrank==0 && in != NULL)
    delete in;
}

//==============================================================================
//    Primary simulation runtime routine
//==============================================================================


void Implicit_Solver::run(){
    
    //MPI info
    world = MPI_COMM_WORLD; //used for convenience to represent all the ranks in the job
    MPI_Comm_rank(world,&myrank);
    MPI_Comm_size(world,&nranks);
    
    if(myrank == 0){
      std::cout << "Running Implicit Solver" << std::endl;
    }
    //initialize Trilinos communicator class
    comm = Tpetra::getDefaultComm();
    if (simparam.input_options.has_value()) {
      const Input_Options& input_options = simparam.input_options.value();
      const char* mesh_file_name = input_options.mesh_file_name.c_str();
      switch (input_options.mesh_file_format) {
        case MESH_FORMAT::tecplot:
          read_mesh_tecplot(mesh_file_name);
          break;
        case MESH_FORMAT::vtk:
          read_mesh_vtk(mesh_file_name);
          break;
        case MESH_FORMAT::ansys_dat:
          read_mesh_ansys_dat(mesh_file_name);
          break;
        case MESH_FORMAT::ensight:
          read_mesh_ensight(mesh_file_name);
          break;
        case MESH_FORMAT::abaqus_inp:
          read_mesh_abaqus_inp(mesh_file_name);
          break;
        default:
          *fos << "ERROR: MESH FILE FORMAT NOT SUPPORTED BY IMPLICIT SOLVER" << std::endl;
          exit_solver(0);
      }
    } else {
      generate_mesh(simparam.mesh_generation_options.value());
    }

    //debug
    //return;
    init_maps();

    //equate pointers for this solver
    initial_node_coords_distributed = node_coords_distributed;
    all_initial_node_coords_distributed = all_node_coords_distributed;

    //print element imbalance stats
    long long int rnum_global_sum = 0;
    long long int temp_rnum_elem = rnum_elem;
    MPI_Allreduce(&temp_rnum_elem, &rnum_global_sum, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
    double local_imbalance, max_imbalance, avg_elem;
    max_imbalance = 0;
    avg_elem = rnum_global_sum/((double) nranks);
    local_imbalance = rnum_elem/avg_elem;
    MPI_Allreduce(&local_imbalance, &max_imbalance, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    *fos << "Element imbalance: " << max_imbalance << std::endl;
    
    //std::cout << "Num elements on process " << myrank << " = " << rnum_elem << std::endl;
    
    init_clock();

    //initialize runtime counters and timers
    int hessvec_count = 0;
    int update_count = 0;
    file_index = 0;
    real_t linear_solve_time = 0;
    real_t hessvec_time = 0;
    real_t hessvec_linear_time = 0;

    if(simparam.output_options.convert_to_tecplot){
      output_design(0,simparam.output_options.convert_to_tecplot);
    }
    
    // ---- Find Boundaries on mesh ---- //
    init_boundaries();

    //set boundary conditions
    generate_tcs();

    //initialize optmization design variable storage
    init_design();

    //process process list of requested FEA modules to construct list of objects
    FEA_module_setup();

    //Have modules read in boundary/loading conditions if file format provides it
    for(int imodule = 0; imodule < nfea_modules; imodule++){
      // TODO: This is almost certainly wrong given the changes to simparam

      //update set options for next FEA module in loop with synchronized set options in solver's simparam class
      FEA_MODULE_TYPE m_type = simparam.fea_module_parameters[imodule]->type;
      bool yaml_supplied_bcs = (bool)(simparam.fea_module_parameters[imodule]->boundary_conditions.size());
      bool yaml_supplied_lcs = (bool)(simparam.fea_module_parameters[imodule]->loading_conditions.size());
      bool yaml_supplied_no_conditions = !yaml_supplied_bcs && !yaml_supplied_lcs;
      bool replace_import_bcs = simparam.fea_module_parameters[imodule]->replace_import_bcs;
      bool imported_conditions = false;
      if(fea_module_must_read.find(m_type) != fea_module_must_read.end() && !replace_import_bcs){
        if (simparam.input_options.has_value()) {
          const Input_Options& input_options = simparam.input_options.value();
          switch (input_options.mesh_file_format) {
          case MESH_FORMAT::ansys_dat:
            fea_modules[imodule]->read_conditions_ansys_dat(in, before_condition_header);
            break;
          case MESH_FORMAT::abaqus_inp:
            fea_modules[imodule]->read_conditions_abaqus_inp(in, before_condition_header);
            break;
          default:
            *fos << "ERROR: MESH FILE FORMAT DOESN'T SUPPORT CONDITION READ" << std::endl;
            exit_solver(0);
          }
          imported_conditions = true;
          fea_modules[imodule]->bcs_initialized = true;
        }
      }
      if((yaml_supplied_bcs||yaml_supplied_lcs)&&!fea_modules[imodule]->bcs_initialized){
        fea_modules[imodule]->init_boundaries();
        fea_modules[imodule]->bcs_initialized = true;
      }
      if(yaml_supplied_bcs){
        //set boundary conditions for FEA modules
        fea_modules[imodule]->generate_bcs();
      }
      if(yaml_supplied_lcs){
        //set applied loading conditions for FEA modules
        fea_modules[imodule]->generate_applied_loads();
      }

      if(simparam.fea_module_parameters[imodule]->requires_conditions&&yaml_supplied_no_conditions&&!imported_conditions){
        *fos << "ERROR: FEA MODULE " << imodule << " WAS NOT ASSIGNED ANY CONDTIONS" << std::endl;
            exit_solver(0);
      }
    }

    //check for errors in the yaml input and exit if any found with an error message
    // TODO: Should figure out a way to add this functionality back
    // int map_size = simparam.unapplied_settings();
    // if(map_size) {
    //   *fos << "YAML input has encountered an error; please correct options that were not applied, or remove unnecessary options." << std::endl;
    //   exit_solver(0);
    // }

    for(int imodule = 0; imodule < nfea_modules; imodule++){
      if(myrank == 0)
        std::cout << "Starting init assembly for module " << imodule <<std::endl <<std::flush;
      //allocate and fill sparse structures needed for global solution in each FEA module
      fea_modules[imodule]->init_assembly();
    
      //assemble the global solution (stiffness matrix etc. and nodal forces)
      fea_modules[imodule]->assemble_matrix();
      std::cout << "Finished matrix assembly for module " << imodule << std::endl <<std::flush;
    
      fea_modules[imodule]->assemble_vector();

      fea_modules[imodule]->linear_solver_parameters();
    
      if(myrank == 0)
        std::cout << "Starting First Solve for module " << imodule << std::endl <<std::flush;

      int solver_exit = fea_modules[imodule]->solve();
      if(solver_exit != EXIT_SUCCESS){
        std::cout << "Linear Solver Error for module " << imodule << std::endl <<std::flush;
        return;
      }

      if(fea_modules_modal_analysis[imodule])
      {
        int eigensolver_exit = fea_modules[imodule]->eigensolve();
        if(solver_exit != 0){
          std::cout << "Eigen-Solver Error for module " << imodule << std::endl <<std::flush;
          return;
        }
      }
    }

    //std::cout << "FEA MODULES " << nfea_modules << " " << simparam.nfea_modules << std::endl;
    //call boundary routines on fea modules
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
    if(simparam.topology_optimization_on||simparam.shape_optimization_on)
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
    output_design(last_print_step+1);
    // vtk_writer();
    if(myrank==0){
      std::cout << "Total number of solves and assembly " << fea_modules[0]->update_count <<std::endl;
      std::cout << "Total number of hessvec counts " << fea_modules[0]->hessvec_count <<std::endl;
      std::cout << "End of Optimization" << std::endl;
    }
}

/* ----------------------------------------------------------------------
   Read ANSYS dat format mesh file
------------------------------------------------------------------------- */
void Implicit_Solver::read_mesh_ansys_dat(const char *MESH){
  Input_Options input_options = simparam.input_options.value();
  char ch;
  int num_dim = simparam.num_dims;
  int p_order = input_options.p_order;
  real_t unit_scaling = input_options.unit_scaling;
  bool restart_file = simparam.restart_file;
  int local_node_index, current_column_index;
  size_t strain_count;
  std::string skip_line, read_line, substring, token;
  std::stringstream line_parse;
  CArrayKokkos<char, array_layout, HostSpace, memory_traits> read_buffer;
  int buffer_loop, buffer_iteration, buffer_iterations, dof_limit, scan_loop, nodes_per_element;
  size_t read_index_start, node_rid, elem_gid, zone_num_elem;
  GO node_gid;
  real_t dof_value;
  host_vec_array node_densities;
  bool zero_index_base = input_options.zero_index_base;
  int negative_index_found = 0;
  int global_negative_index_found = 0;
  num_elem = 0;
  rnum_elem = 0;
  int etype_index = 0;
  bool searching_for_elements = true;
  bool no_elements_found = true;
  elements::elem_types::elem_type mesh_element_type;
  std::vector<size_t> element_temp;
  std::vector<size_t> global_indices_temp;

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
  global_size_t local_nrows = map->getLocalNumElements();
  nlocal_nodes = local_nrows;
  //populate local row offset data from global data
  global_size_t min_gid = map->getMinGlobalIndex();
  global_size_t max_gid = map->getMaxGlobalIndex();
  global_size_t index_base = map->getIndexBase();
  //debug print
  //std::cout << "local node count on task: " << " " << nlocal_nodes << std::endl;

  //allocate node storage with dual view
  //dual_node_coords = dual_vec_array("dual_node_coords", nlocal_nodes,num_dim);
  //if(restart_file)
    //dual_node_densities = dual_vec_array("dual_node_densities", nlocal_nodes,1);

  //local variable for host view in the dual view
  node_coords_distributed = Teuchos::rcp(new MV(map, num_dim));
  //view scope
  {
  host_vec_array node_coords = node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  //host_vec_array node_coords = dual_node_coords.view_host();
  if(restart_file){
    design_node_densities_distributed = Teuchos::rcp(new MV(map, 1));
    node_densities = design_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  }
  //notify that the host view is going to be modified in the file readin
  //dual_node_coords.modify_host();
  //if(restart_file)
    //dual_node_densities.modify_host();

  //old swage method
  //mesh->init_nodes(local_nrows); // add 1 for index starting at 1
    
  //std::cout << "Num nodes assigned to task " << myrank << " = " << nlocal_nodes << std::endl;

  // read the initial mesh coordinates
  // x-coords
  /*only task 0 reads in nodes and elements from the input file
  stores node data in a buffer and communicates once the buffer cap is reached
  or the data ends*/

  words_per_line = input_options.words_per_line;
  //if(restart_file) words_per_line++;
  elem_words_per_line = input_options.elem_words_per_line;

  //allocate read buffer
  read_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(BUFFER_LINES,words_per_line,MAX_WORD);

  dof_limit = num_nodes;
  buffer_iterations = dof_limit/BUFFER_LINES;
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
  } //end view scope
  //repartition node distribution
  repartition_nodes();

  //synchronize device data
  //dual_node_coords.sync_device();
  //dual_node_coords.modify_device();
  //if(restart_file){
    //dual_node_densities.sync_device();
    //dual_node_densities.modify_device();
  //}

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
  //seek element connectivity zone; ansys dat can store element data in multiple file zones for different declared bodies
  size_t num_elem_zone_readins = 0;
  while(searching_for_elements){
    num_elem_zone_readins++;
    if(myrank==0){
      //search for element zone header
      bool searching_for_element_zone = true;
      std::string substring2, substring3;
      while (in->good()&&searching_for_element_zone) {
        getline(*in, skip_line);
        //std::cout << skip_line << std::endl;
        // line_parse.clear();
        // line_parse.str(skip_line);
        // //stop when the NODES= string is reached
        // while (!line_parse.eof()){
        //   line_parse >> substring;
        //   //std::cout << substring << std::endl;
        //   if(!substring.compare("Elements")){
        //     no_elements_found = false;
        //     searching_for_element_zone = false;
        //     break;
        //   }
        // } //while
        if(skip_line.find("Elements for Body")!=std::string::npos){
            no_elements_found = false;
            searching_for_element_zone = false;
            break;
        }
        
      }

      //check if no more element zones could be found
      if(!in->good()){
        searching_for_elements = false;
      }

      if(no_elements_found){
        std::cout << "FILE FORMAT ERROR; NO ELEMENTS FOUND" << std::endl;
      }
    }

    //broadcast element zone find flag
    MPI_Bcast(&searching_for_elements,1,MPI_CXX_BOOL,0,world);
    if(!searching_for_elements){
      break;
    }
    //std::cout << "NO MORE ELEMENT READINS AFTER " << num_elem_zone_readins << std::endl;

    if(myrank==0){
      //read in element type from following line
      getline(*in, read_line);
      std::cout << read_line << std::endl;
      // line_parse.clear();
      // line_parse.str(read_line);
      // line_parse >> substring;
      //std::cout << substring << std::endl;
      if(read_line.find("185")!=std::string::npos){
        //Hex8 type
        etype_index = 1;

      }
      else if(read_line.find("186")!=std::string::npos){
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
      zone_num_elem = std::stoi(token);

      //skip line
      for (int j = 0; j < 1; j++) {
        getline(*in, skip_line);
        std::cout << skip_line << std::endl;
      }
    }

    //broadcast element type
    MPI_Bcast(&etype_index,1,MPI_INT,0,world);

    int elem_words_per_line_no_nodes = 11;
    if(etype_index==1){
      mesh_element_type = elements::elem_types::Hex8;
      nodes_per_element = 8;
      elem_words_per_line = elem_words_per_line_no_nodes + 8;
      max_nodes_per_patch = 4;
    }
    else if(etype_index==2){
      mesh_element_type = elements::elem_types::Hex20;
      nodes_per_element = 20;
      elem_words_per_line = elem_words_per_line_no_nodes + 20;
      max_nodes_per_patch = 8;
    }
    else if(etype_index==3){
      mesh_element_type = elements::elem_types::Hex32;
      nodes_per_element = 32;
      elem_words_per_line = elem_words_per_line_no_nodes + 32;
      max_nodes_per_patch = 12;
    }
    else{
      *fos << "ERROR: ANSYS ELEMENT TYPE NOT FOUND OR RECOGNIZED" << std::endl;
      exit_solver(0);
    }
    
    //broadcast number of elements
    MPI_Bcast(&zone_num_elem,1,MPI_LONG_LONG_INT,0,world);
    
    *fos << "declared element count in this zone: " << zone_num_elem << std::endl;
    num_elem += zone_num_elem;
    *fos << "total declared elements is now: " << num_elem << std::endl;

    //std::cout<<"before initial mesh initialization"<<std::endl;
    
    //read in element connectivity
    //we're gonna reallocate for the words per line expected for the element connectivity
    read_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(BUFFER_LINES,elem_words_per_line,MAX_WORD); 
    CArrayKokkos<int, array_layout, HostSpace, memory_traits> node_store(nodes_per_element);

    //calculate buffer iterations to read number of lines
    buffer_iterations = zone_num_elem/BUFFER_LINES;
    int assign_flag;

    //dynamic buffer used to store elements before we know how many this rank needs
    if(num_elem_zone_readins == 1){
      element_temp = std::vector<size_t>(BUFFER_LINES*elem_words_per_line);
      global_indices_temp = std::vector<size_t>(BUFFER_LINES);
    }
    size_t buffer_max = BUFFER_LINES*elem_words_per_line;
    size_t indices_buffer_max = BUFFER_LINES;

    if(zone_num_elem%BUFFER_LINES!=0) buffer_iterations++;
    read_index_start = 0;
    //std::cout << "ELEMENT BUFFER ITERATIONS: " << buffer_iterations << std::endl;
    //std::cout << "BUFFER ITERATIONS IS: " << buffer_iterations << std::endl;
    for(buffer_iteration = 0; buffer_iteration < buffer_iterations; buffer_iteration++){
      //pack buffer on rank 0
      if(myrank==0&&buffer_iteration<buffer_iterations-1){
        for (buffer_loop = 0; buffer_loop < BUFFER_LINES; buffer_loop++) {
          getline(*in,read_line);
          //std::cout << read_line << std::endl;
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
        while(buffer_iteration*BUFFER_LINES+buffer_loop < zone_num_elem) {
          getline(*in,read_line);
          //std::cout << read_line << std::endl;
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
        elem_gid = read_index_start + scan_loop + num_elem - zone_num_elem;
        //std::cout << elem_gid << std::endl;
        //add this element to the local list if any of its nodes belong to this rank according to the map
        //get list of nodes for each element line and check if they belong to the map
        assign_flag = 0;
        for(int inode = elem_words_per_line_no_nodes; inode < elem_words_per_line; inode++){
          //as we loop through the nodes belonging to this element we store them
          //if any of these nodes belongs to this rank this list is used to store the element locally
          node_gid = atoi(&read_buffer(scan_loop,inode,0));
          //std::cout << node_gid << " ";
          if(zero_index_base)
            node_store(inode-elem_words_per_line_no_nodes) = node_gid; //subtract 1 since file index start is 1 but code expects 0
          else
            node_store(inode-elem_words_per_line_no_nodes) = node_gid - 1; //subtract 1 since file index start is 1 but code expects 0
          if(node_store(inode-elem_words_per_line_no_nodes) < 0){
            negative_index_found = 1;
          }
          //first we add the elements to a dynamically allocated list
          if(zero_index_base){
            if(map->isNodeGlobalElement(node_gid)&&!assign_flag){
              assign_flag = 1;
              rnum_elem++;
            }
          }
          else{
            if(map->isNodeGlobalElement(node_gid-1)&&!assign_flag){
              assign_flag = 1;
              rnum_elem++;
            }
          }
        }
        //std::cout << std::endl;

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

    //if the next attempt to find an element zone failed; this will allow moving the stream position so we can find other data such as boundary conditions
    if(myrank==0){
      if(in->good())
        prev_read_elem_zone_end = in->tellg();
    }
  } //element read in while
  
  //move streamposition back to the end of the last element connectivity data zone
  if(myrank==0){
    in->clear();
    in->seekg(prev_read_elem_zone_end);
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
        if(!substring.compare("Supports")||!substring.compare("Pressure")||
           !substring.compare("Displacements")||!substring.compare("Displacement\"")){
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
    // check that the input file has configured some kind of acceptable module
    simparam.validate_module_is_specified(FEA_MODULE_TYPE::Elasticity);
    simparam.fea_module_must_read.insert(FEA_MODULE_TYPE::Elasticity);
  }

  // Close mesh input file if no further readin is done by FEA modules for conditions
  if(myrank==0&&No_Conditions){
    in->close();
  }

  std::cout << "RNUM ELEMENTS IS: " << rnum_elem << std::endl;
  //copy temporary element storage to multivector storage
  Element_Types = CArrayKokkos<elements::elem_types::elem_type, array_layout, HostSpace, memory_traits>(rnum_elem);

  //set element object pointer
  if(simparam.num_dims==2){
    element_select->choose_2Delem_type(mesh_element_type, elem2D);
     max_nodes_per_element = elem2D->num_nodes();
  }
  else if(simparam.num_dims==3){
    element_select->choose_3Delem_type(mesh_element_type, elem);
     max_nodes_per_element = elem->num_nodes();
  }

  //1 type per mesh for now
  for(int ielem = 0; ielem < rnum_elem; ielem++)
    Element_Types(ielem) = mesh_element_type;

  dual_nodes_in_elem = dual_elem_conn_array("dual_nodes_in_elem", rnum_elem, max_nodes_per_element);
  host_elem_conn_array nodes_in_elem = dual_nodes_in_elem.view_host();
  dual_nodes_in_elem.modify_host();

  for(int ielem = 0; ielem < rnum_elem; ielem++)
    for(int inode = 0; inode < nodes_per_element; inode++){
      nodes_in_elem(ielem, inode) = element_temp[ielem*nodes_per_element + inode];
    }

  //view storage for all local elements connected to local nodes on this rank
  Kokkos::DualView <GO*, array_layout, device_type, memory_traits> All_Element_Global_Indices("All_Element_Global_Indices",rnum_elem);

  //copy temporary global indices storage to view storage
  for(int ielem = 0; ielem < rnum_elem; ielem++){
    All_Element_Global_Indices.h_view(ielem) = global_indices_temp[ielem];
    if(global_indices_temp[ielem]<0){
      negative_index_found = 1;
    }
  }
  
  MPI_Allreduce(&negative_index_found,&global_negative_index_found,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  if(global_negative_index_found){
    if(myrank==0){
    std::cout << "Node index less than or equal to zero detected; set \"zero_index_base: true\" under \"input_options\" in your yaml file if indices start at 0" << std::endl;
    }
    exit_solver(0);
  }
  
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

  All_Element_Global_Indices.modify_host();
  All_Element_Global_Indices.sync_device();
  
  //construct overlapping element map (since different ranks can own the same elements due to the local node map)
  all_element_map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),All_Element_Global_Indices.d_view,0,comm));


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
} // end read_mesh

/* ----------------------------------------------------------------------
   Construct list of objects for FEA modules 
------------------------------------------------------------------------- */

void Implicit_Solver::FEA_module_setup(){
  nfea_modules = simparam.fea_module_parameters.size();
  fea_module_must_read = simparam.fea_module_must_read;
  //allocate lists to size
  fea_module_types = std::vector<FEA_MODULE_TYPE>();
  fea_modules = std::vector<FEA_Module*>();
  fea_modules_modal_analysis = std::vector<bool>();
  bool module_found = false;
  displacement_module = -1;
  for (auto& param : simparam.fea_module_parameters) {
    fea_modules_modal_analysis.push_back(false);
    fea_module_types.push_back(param->type);
    param->apply(
      [&](Elasticity_Parameters& param) {
        displacement_module = fea_modules.size();
        fea_modules.push_back(new FEA_Module_Elasticity(param, this, fea_modules.size()));
      },
      [&](Inertial_Parameters& param) {
        fea_modules.push_back(new FEA_Module_Inertial(param, this, fea_modules.size()));
      },
      [&](Heat_Conduction_Parameters& param) {
        fea_modules.push_back(new FEA_Module_Heat_Conduction(param, this, fea_modules.size()));
      },
      [&](Thermo_Elasticity_Parameters& param) {
        //ensure another momentum conservation module was not allocated
        if (displacement_module >= 0){
          *fos << "PROGRAM IS ENDING DUE TO ERROR; MORE THAN ONE ELASTIC MODULE ALLOCATED \"" 
                << to_string(param.type) << "\"" << std::endl;
          exit_solver(0);
        }
        displacement_module = fea_modules.size();
        fea_modules.push_back(new FEA_Module_Thermo_Elasticity(param, this, fea_modules.size()));
      },
      [&](const FEA_Module_Parameters& param) {
        *fos << "PROGRAM IS ENDING DUE TO ERROR; UNDEFINED FEA MODULE REQUESTED WITH NAME \"" << to_string(param.type) <<"\"" << std::endl;
        exit_solver(0);
      }
    );
    
    *fos << " " << fea_module_types.back() << " MODULE ALLOCATED AS " << fea_module_types.size() - 1 << std::endl;
  }
}

/* ----------------------------------------------------------------------
   Setup Optimization Problem Object, Relevant Objective, and Constraints
------------------------------------------------------------------------- */

void Implicit_Solver::setup_optimization_problem(){
  int num_dim = simparam.num_dims;
  bool nodal_density_flag = simparam.nodal_density_flag;
  int nOpt_modules = simparam.Optimization_Module_List.size();
  int nmulti_objective_modules = simparam.Multi_Objective_Modules.size();
  std::vector<OPTIMIZATION_MODULE_TYPE> Optimization_Module_List = simparam.Optimization_Module_List;
  std::vector<int> Optimization_Module_My_FEA_Module = simparam.Optimization_Module_My_FEA_Module;
  std::vector<int> Multi_Objective_Modules = simparam.Multi_Objective_Modules;
  std::vector<real_t> Multi_Objective_Weights = simparam.Multi_Objective_Weights;
  std::vector<std::vector<real_t>> Function_Arguments = simparam.Function_Arguments;
  std::vector<FUNCTION_TYPE> Optimization_Function_Type = simparam.Optimization_Function_Type;
  std::vector<ROL::Ptr<ROL::Objective<real_t>>> Multi_Objective_Terms;

  std::string constraint_base, constraint_name;
  std::stringstream number_union;
  CArray<GO> Surface_Nodes;
  GO current_node_index, current_element_index;
  LO local_node_index, local_element_id;
  int num_bdy_patches_in_set;
  size_t node_id, patch_id, module_id;
  int num_boundary_sets;
  int local_surface_id;
  const_host_vec_array design_densities;
  typedef ROL::TpetraMultiVector<real_t,LO,GO,node_type> ROL_MV;
  const_host_elem_conn_array nodes_in_elem = global_nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  
  // fill parameter list with desired algorithmic options or leave as default
  // Read optimization input parameter list.
  Teuchos::RCP<Teuchos::ParameterList> parlist;
  if(simparam.optimization_options.optimization_parameters_xml_file){
    std::string xmlFileName = simparam.optimization_options.xml_parameters_file_name;
    
    //check if parameter file exists
    std::ifstream param_file_check(xmlFileName);
    if (!param_file_check.is_open()) {
        *fos << "Unable to find xml parameter file required for optimization with the ROL library"  << std::endl;
        exit_solver(0);
    }
    param_file_check.close();

    parlist = ROL::getParametersFromXmlFile(xmlFileName);
  }
  else{
    parlist = Teuchos::rcp(new Teuchos::ParameterList("Inputs"));
    set_rol_params(parlist);
  }

  //ROL::ParameterList parlist;

  //Design variables to optimize
  ROL::Ptr<ROL::Vector<real_t>> x;
  if(nodal_density_flag)
    x = ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO>>(design_node_densities_distributed);
  else
    x = ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO>>(Global_Element_Densities);

  //set bounds on design variables
  ROL::Ptr<ROL::BoundConstraint<real_t> > bnd, mma_bnd;
  // Bound constraint defining the possible range of design density variables
  ROL::Ptr<ROL::Vector<real_t> > lower_bounds, mma_lower_bounds;
  ROL::Ptr<ROL::Vector<real_t> > upper_bounds, mma_upper_bounds;
  if(nodal_density_flag){
    //allocate global vector information
    upper_bound_node_densities_distributed = Teuchos::rcp(new MV(map, 1));
    all_lower_bound_node_densities_distributed = Teuchos::rcp(new MV(all_node_map, 1));
    lower_bound_node_densities_distributed = Teuchos::rcp(new MV(*all_lower_bound_node_densities_distributed, map));
    host_vec_array node_densities_upper_bound = upper_bound_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    host_vec_array node_densities_lower_bound = all_lower_bound_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);


    //initialize densities to 1 for now; in the future there might be an option to read in an initial condition for each node
    for(int inode = 0; inode < nlocal_nodes; inode++){
      if(simparam.optimization_options.maximum_density<1){
        node_densities_upper_bound(inode,0) = simparam.optimization_options.maximum_density;
      }
      else{
        node_densities_upper_bound(inode,0) = 1;
      }
    }
    for(int inode = 0; inode < nall_nodes; inode++){
      if(simparam.optimization_options.minimum_density>simparam.optimization_options.density_epsilon){
        node_densities_lower_bound(inode,0) = simparam.optimization_options.minimum_density;
      }
      else{
        node_densities_lower_bound(inode,0) = simparam.optimization_options.density_epsilon;
      }
    }

    if(simparam.optimization_options.method_of_moving_asymptotes){
      Teuchos::RCP<MV> mma_upper_bound_node_densities_distributed = Teuchos::rcp(new MV(map, 1));
      Teuchos::RCP<MV> mma_lower_bound_node_densities_distributed = Teuchos::rcp(new MV(map, 1));
      //mma_lower_bound_node_densities_distributed->assign(*lower_bound_node_densities_distributed);
      //mma_upper_bound_node_densities_distributed->assign(*upper_bound_node_densities_distributed);
      
      mma_lower_bound_node_densities_distributed->putScalar(-0.1);
      mma_upper_bound_node_densities_distributed->putScalar(0.1);
      mma_lower_bounds = ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO>>(mma_lower_bound_node_densities_distributed);
      mma_upper_bounds = ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO>>(mma_upper_bound_node_densities_distributed);
      
      mma_bnd = ROL::makePtr<ROL::Bounds<real_t>>(mma_lower_bounds, mma_upper_bounds);
    }

    //set lower bounds for nodes on surfaces with boundary and loading conditions
    if(!simparam.optimization_options.variable_outer_shell){
      for(int imodule = 0; imodule < nfea_modules; imodule++){
        num_boundary_sets = fea_modules[imodule]->num_boundary_conditions;
        for(int iboundary = 0; iboundary < num_boundary_sets; iboundary++){

          num_bdy_patches_in_set = fea_modules[imodule]->NBoundary_Condition_Patches(iboundary);

          //loop over boundary patches for this boundary set
          for (int bdy_patch_gid = 0; bdy_patch_gid < num_bdy_patches_in_set; bdy_patch_gid++){
                  
            // get the global id for this boundary patch
            patch_id = fea_modules[imodule]->Boundary_Condition_Patches(iboundary, bdy_patch_gid);
            if(simparam.optimization_options.thick_condition_boundary){
              Surface_Nodes = Boundary_Patches(patch_id).node_set;
              current_element_index = Boundary_Patches(patch_id).element_id;
              //debug print of local surface ids
              //std::cout << " LOCAL SURFACE IDS " << std::endl;
              //std::cout << local_surface_id << std::endl;
              //acquire set of nodes for this face
              for(int node_loop=0; node_loop < max_nodes_per_element; node_loop++){
                current_node_index = nodes_in_elem(current_element_index,node_loop);
                local_node_index = all_node_map->getLocalElement(current_node_index);
                node_densities_lower_bound(local_node_index,0) = simparam.optimization_options.shell_density;
              }// node loop for
            }//if
            else{
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
                  node_densities_lower_bound(local_node_index,0) = simparam.optimization_options.shell_density;
                }
              }// node loop for
            }//if
          }//boundary patch for
        }//boundary set for

        //set node conditions due to point BCS that might not show up in boundary sets
        //possible to have overlap in which nodes are set with the previous loop
        if(fea_modules[imodule]->node_specified_bcs)
          fea_modules[imodule]->node_density_constraints(node_densities_lower_bound);
      }//module for

      //enforce rho=1 on all boundary patches if user setting enabled
      if(simparam.optimization_options.retain_outer_shell){
        //loop over boundary patches for this boundary set
          for (int bdy_patch_gid = 0; bdy_patch_gid < nboundary_patches; bdy_patch_gid++){
            patch_id = bdy_patch_gid;
            if(simparam.optimization_options.thick_condition_boundary){
              Surface_Nodes = Boundary_Patches(patch_id).node_set;
              current_element_index = Boundary_Patches(patch_id).element_id;
              //debug print of local surface ids
              //std::cout << " LOCAL SURFACE IDS " << std::endl;
              //std::cout << local_surface_id << std::endl;
              //acquire set of nodes for this face
              for(int node_loop=0; node_loop < max_nodes_per_element; node_loop++){
                current_node_index = nodes_in_elem(current_element_index,node_loop);
                local_node_index = all_node_map->getLocalElement(current_node_index);
                node_densities_lower_bound(local_node_index,0) = simparam.optimization_options.shell_density;
              }// node loop for
            }//if
            else{
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
                  node_densities_lower_bound(local_node_index,0) = simparam.optimization_options.shell_density;
                }
              }// node loop for
            }//if
          }//boundary patch for
      }
    }//if to check if shell should have any constraints

    //constraints due to specified user regions
    const size_t num_fills = simparam.optimization_options.volume_bound_constraints.size();
    const_host_vec_array node_coords_view = node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    const DCArrayKokkos <Optimization_Bound_Constraint_Region> mat_fill = simparam.optimization_options.optimization_bound_constraint_volumes;
    
    for(int ifill = 0; ifill < num_fills; ifill++){
      for(int inode = 0; inode < nlocal_nodes; inode++){
        real_t node_coords[3];
        node_coords[0] = node_coords_view(inode,0);
        node_coords[1] = node_coords_view(inode,1);
        node_coords[2] = node_coords_view(inode,2);
        bool fill_this = mat_fill(ifill).volume.contains(node_coords);
        if(fill_this){
          node_densities_lower_bound(inode,0) = mat_fill(ifill).set_lower_density_bound;
          //make sure 0 setting is increased to the epsilon value setting
          if(node_densities_lower_bound(inode,0) < simparam.optimization_options.density_epsilon){
            node_densities_lower_bound(inode,0) = simparam.optimization_options.density_epsilon;
          }
          node_densities_upper_bound(inode,0) = mat_fill(ifill).set_upper_density_bound;
        }
      }//node for
    }//fill region for

    //communicate constraints vector
    lower_bound_node_densities_distributed->doExport(*all_lower_bound_node_densities_distributed, *exporter, Tpetra::ABSMAX, true);
  }
  else{
    //initialize memory for volume storage
    vec_array Element_Densities_Upper_Bound("Element Densities_Upper_Bound", rnum_elem, 1);
    vec_array Element_Densities_Lower_Bound("Element Densities_Lower_Bound", rnum_elem, 1);
    for(int ielem = 0; ielem < rnum_elem; ielem++){
      Element_Densities_Upper_Bound(ielem,0) = 1;
      Element_Densities_Lower_Bound(ielem,0) = simparam.optimization_options.density_epsilon;
    }

    //create global vector
    Global_Element_Densities_Upper_Bound = Teuchos::rcp(new MV(element_map, Element_Densities_Upper_Bound));
    Global_Element_Densities_Lower_Bound = Teuchos::rcp(new MV(element_map, Element_Densities_Lower_Bound));
  }
    
  if(nodal_density_flag){
    lower_bounds = ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO>>(lower_bound_node_densities_distributed);
    upper_bounds = ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO>>(upper_bound_node_densities_distributed);
  }
  else{
    lower_bounds = ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO>>(Global_Element_Densities_Lower_Bound);
    upper_bounds = ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO>>(Global_Element_Densities_Upper_Bound);
  }
  bnd = ROL::makePtr<ROL::Bounds<real_t>>(lower_bounds, upper_bounds);
  
  //Instantiate (the one) objective function for the problem
  ROL::Ptr<ROL::Objective<real_t>> obj, sub_obj;
  bool objective_declared = false;
  for(int imodule = 0; imodule < nOpt_modules; imodule++){
    if(Optimization_Function_Type[imodule] == FUNCTION_TYPE::OBJECTIVE){
      //check if previous module already defined an objective, there must be one objective module
      if(objective_declared){
        *fos << "PROGRAM IS ENDING DUE TO ERROR; ANOTHER OBJECTIVE FUNCTION WITH NAME \"" 
              << to_string(Optimization_Module_List[imodule]) <<"\" ATTEMPTED TO REPLACE A PREVIOUS OBJECTIVE; THERE MUST BE ONE OBJECTIVE." << std::endl;
          exit_solver(0);
      }
      if(Optimization_Module_List[imodule] == OPTIMIZATION_MODULE_TYPE::Strain_Energy_Minimize_TopOpt){
        //debug print
        *fos << " STRAIN ENERGY OBJECTIVE EXPECTS FEA MODULE INDEX " <<Optimization_Module_My_FEA_Module[imodule] << std::endl;
        if(simparam.optimization_options.method_of_moving_asymptotes){
          sub_obj = ROL::makePtr<StrainEnergyMinimize_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[imodule]], nodal_density_flag);
          obj = ROL::makePtr<ObjectiveMMA>(sub_obj, mma_bnd, x);
        }
        else{
          obj = ROL::makePtr<StrainEnergyMinimize_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[imodule]], nodal_density_flag);
        }
      }
      else if(Optimization_Module_List[imodule] == OPTIMIZATION_MODULE_TYPE::Heat_Capacity_Potential_Minimize_TopOpt){
        //debug print
        *fos << " HEAT CAPACITY POTENTIAL OBJECTIVE EXPECTS FEA MODULE INDEX " <<Optimization_Module_My_FEA_Module[imodule] << std::endl;
        if(simparam.optimization_options.method_of_moving_asymptotes){
          sub_obj = ROL::makePtr<HeatCapacityPotentialMinimize_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[imodule]], nodal_density_flag);
          obj = ROL::makePtr<ObjectiveMMA>(sub_obj, mma_bnd, x);
        }
        else{
          obj = ROL::makePtr<HeatCapacityPotentialMinimize_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[imodule]], nodal_density_flag);
        }
      }
      //Multi-Objective case
      else if(Optimization_Module_List[imodule] == OPTIMIZATION_MODULE_TYPE::Multi_Objective){
        //allocate vector of Objective Functions to pass
        Multi_Objective_Terms = std::vector<ROL::Ptr<ROL::Objective<real_t>>>(nmulti_objective_modules);
        for(int imulti = 0; imulti < nmulti_objective_modules; imulti++){
          //get module index for objective term
          module_id = Multi_Objective_Modules[imulti];
          if(Optimization_Module_List[module_id] == OPTIMIZATION_MODULE_TYPE::Strain_Energy_Minimize_TopOpt){
            //debug print
            *fos << " STRAIN ENERGY OBJECTIVE EXPECTS FEA MODULE INDEX " <<Optimization_Module_My_FEA_Module[module_id] << std::endl;
            if(simparam.optimization_options.method_of_moving_asymptotes){
              sub_obj = ROL::makePtr<StrainEnergyMinimize_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[module_id]], nodal_density_flag);
              Multi_Objective_Terms[imulti] = ROL::makePtr<ObjectiveMMA>(sub_obj, mma_bnd, x);
            }
            else{
              Multi_Objective_Terms[imulti] = ROL::makePtr<StrainEnergyMinimize_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[module_id]], nodal_density_flag);
            }
          }
          else if(Optimization_Module_List[module_id] == OPTIMIZATION_MODULE_TYPE::Heat_Capacity_Potential_Minimize_TopOpt){
            //debug print
            *fos << " HEAT CAPACITY POTENTIAL OBJECTIVE EXPECTS FEA MODULE INDEX " <<Optimization_Module_My_FEA_Module[module_id] << std::endl;
            if(simparam.optimization_options.method_of_moving_asymptotes){
              sub_obj = ROL::makePtr<HeatCapacityPotentialMinimize_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[module_id]], nodal_density_flag);
              Multi_Objective_Terms[imulti] = ROL::makePtr<ObjectiveMMA>(sub_obj, mma_bnd, x);}
            else{
              Multi_Objective_Terms[imulti] = ROL::makePtr<HeatCapacityPotentialMinimize_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[module_id]], nodal_density_flag);
            }
          }
          else if(Optimization_Module_List[module_id] == OPTIMIZATION_MODULE_TYPE::Thermo_Elastic_Strain_Energy_Minimize_TopOpt){
            //debug print
            *fos << " HEAT CAPACITY POTENTIAL OBJECTIVE EXPECTS FEA MODULE INDEX " <<Optimization_Module_My_FEA_Module[module_id] << std::endl;
            if(simparam.optimization_options.method_of_moving_asymptotes){
              sub_obj = ROL::makePtr<HeatCapacityPotentialMinimize_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[module_id]], nodal_density_flag);
              Multi_Objective_Terms[imulti] = ROL::makePtr<ObjectiveMMA>(sub_obj, mma_bnd, x);}
            else{
              Multi_Objective_Terms[imulti] = ROL::makePtr<HeatCapacityPotentialMinimize_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[module_id]], nodal_density_flag);
            }
          }
        }
        //allocate multi objective function
        obj = ROL::makePtr<MultiObjective_TopOpt>(Multi_Objective_Terms, Multi_Objective_Weights);
      }
      else{
        *fos << "PROGRAM IS ENDING DUE TO ERROR; UNDEFINED OBJECTIVE FUNCTION REQUESTED WITH NAME \"" 
              << Optimization_Module_List[imodule] <<"\"" << std::endl;
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
  
  for(int imodule = 0; imodule < nOpt_modules; imodule++){
    number_union.str(constraint_base);
    number_union << imodule + 1;
    constraint_name = number_union.str();
    ROL::Ptr<std::vector<real_t> > li_ptr = ROL::makePtr<std::vector<real_t>>(1,0.0);
    ROL::Ptr<ROL::Vector<real_t> > constraint_mul = ROL::makePtr<ROL::StdVector<real_t>>(li_ptr);
    if(Optimization_Function_Type[imodule] == FUNCTION_TYPE::EQUALITY_CONSTRAINT){
      //pointers are reference counting
      ROL::Ptr<ROL::Constraint<real_t>> eq_constraint;
      if(Optimization_Module_List[imodule] == OPTIMIZATION_MODULE_TYPE::Mass_Constraint_TopOpt){
        
        *fos << " MASS CONSTRAINT EXPECTS FEA MODULE INDEX " <<Optimization_Module_My_FEA_Module[imodule] << std::endl;
        eq_constraint = ROL::makePtr<MassConstraint_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[imodule]], nodal_density_flag, Function_Arguments[imodule][0], false);
      }
      else if(Optimization_Module_List[imodule] == OPTIMIZATION_MODULE_TYPE::Displacement_Constraint_TopOpt){
        //simple test of assignment to the vector for constraint dofs
        Teuchos::RCP<MV> target_displacements = Teuchos::rcp(new MV(local_dof_map, 1));
        Teuchos::RCP<Tpetra::MultiVector<int,LO,GO>> active_dofs = Teuchos::rcp(new Tpetra::MultiVector<int,LO,GO>(local_dof_map, 1));
        active_dofs->putScalar(0);
        host_vec_array target_displacements_view = target_displacements->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
        host_ivec_array active_dofs_view = active_dofs->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
        std::string constraint_filename = simparam.optimization_options.constraints[imodule-nmulti_objective_modules-1].argument_file_name;
        //std::cout << "DISPLACEMENT CONSTRAINT SETTING FILE NAME " << constraint_filename << std::endl;
        std::ifstream *disp_in = NULL;
        disp_in = new std::ifstream();
        disp_in->open(constraint_filename);
        if (!(*disp_in)) throw std::runtime_error(std::string("Can't open ") + constraint_filename);
        if (disp_in->is_open()) {
        std::string line;
          while (std::getline(*disp_in, line)) { // Read each line until the end of the file
            std::istringstream ss(line); // Create a string stream from the line
            GO tag_id;
            real_t target_displacement_value;
              // Process each number (replace this with your own logic)
              ss >> tag_id;
              ss >> target_displacement_value;
              //std::cout << "Read number: " << tag_id << " " << target_displacement_value << std::endl;
              if(map->isNodeGlobalElement(tag_id)){
                LO local_node_id = map->getLocalElement(tag_id);
                active_dofs_view(local_node_id*num_dim,0) = 1;
                target_displacements_view(local_node_id*num_dim,0) = target_displacement_value;
              }

          }
          disp_in->close(); // Close the file
        }
        *fos << " DISPLACEMENT CONSTRAINT EXPECTS FEA MODULE INDEX " <<Optimization_Module_My_FEA_Module[imodule] << std::endl;
        eq_constraint = ROL::makePtr<DisplacementConstraint_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[imodule]], nodal_density_flag, target_displacements, active_dofs, 0, false);

        //ROL::Ptr<ROL::TpetraMultiVector<real_t,LO,GO,node_type>> rol_x =
        //ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO,node_type>>(design_node_densities_distributed);
        //construct direction vector for check
        //Teuchos::RCP<MV> directions_distributed = Teuchos::rcp(new MV(map, 1));
        //directions_distributed->putScalar(1);
        //directions_distributed->randomize(-0.8,1);
        // ROL::Ptr<ROL::TpetraMultiVector<real_t,LO,GO,node_type>> rol_d =
        // ROL::makePtr<ROL::TpetraMultiVector<real_t,LO,GO,node_type>>(directions_distributed);
        
        // ROL::Ptr<std::vector<real_t> > c_ptr = ROL::makePtr<std::vector<real_t>>(1,1.0);
        // ROL::Ptr<ROL::Vector<real_t> > constraint_buf = ROL::makePtr<ROL::StdVector<real_t>>(c_ptr);
        //std::cout << " VALUE TEST " << obj_value << std::endl;
        //eq_constraint->checkApplyJacobian(*rol_x, *rol_d, *constraint_buf);
        //eq_constraint->checkApplyAdjointHessian(*rol_x, *constraint_buf, *rol_d, *rol_x);
      }
      else if(Optimization_Module_List[imodule] == OPTIMIZATION_MODULE_TYPE::Moment_of_Inertia_Constraint_TopOpt){
        *fos << " MOMENT OF INERTIA CONSTRAINT EXPECTS FEA MODULE INDEX " <<Optimization_Module_My_FEA_Module[imodule] << std::endl;
        eq_constraint = ROL::makePtr<MomentOfInertiaConstraint_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[imodule]], nodal_density_flag, Function_Arguments[imodule][1], Function_Arguments[imodule][0], false);
      }
      else if(Optimization_Module_List[imodule] == OPTIMIZATION_MODULE_TYPE::Center_of_Mass_Constraint_TopOpt){
        *fos << " CENTER OF MASS CONSTRAINT EXPECTS FEA MODULE INDEX " <<Optimization_Module_My_FEA_Module[imodule] << std::endl;
        eq_constraint = ROL::makePtr<CenterOfMassConstraint_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[imodule]], nodal_density_flag, Function_Arguments[imodule][1], Function_Arguments[imodule][0], false);
      }
      else if(Optimization_Module_List[imodule] == OPTIMIZATION_MODULE_TYPE::Strain_Energy_Constraint_TopOpt){    
        *fos << " STRAIN ENERGY CONSTRAINT EXPECTS FEA MODULE INDEX " <<Optimization_Module_My_FEA_Module[imodule] << std::endl;
        eq_constraint = ROL::makePtr<StrainEnergyConstraint_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[imodule]], nodal_density_flag, Function_Arguments[imodule][0], false);
      }
      else if(Optimization_Module_List[imodule] == OPTIMIZATION_MODULE_TYPE::Heat_Capacity_Potential_Constraint_TopOpt){
        *fos << " HEAT CAPACITY POTENTIAL CONSTRAINT EXPECTS FEA MODULE INDEX " <<Optimization_Module_My_FEA_Module[imodule] << std::endl;
        eq_constraint = ROL::makePtr<HeatCapacityPotentialConstraint_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[imodule]], nodal_density_flag, Function_Arguments[imodule][0], false);
      }
      else{
        *fos << "PROGRAM IS ENDING DUE TO ERROR; UNDEFINED EQUALITY CONSTRAINT FUNCTION REQUESTED WITH NAME \"" <<Optimization_Module_List[imodule] <<"\"" << std::endl;
        exit_solver(0);
      }
      *fos << " ADDING CONSTRAINT " << constraint_name << std::endl;
      problem->addConstraint(constraint_name, eq_constraint, constraint_mul);
    }

    if(Optimization_Function_Type[imodule] == FUNCTION_TYPE::INEQUALITY_CONSTRAINT){
      //pointers are reference counting
      ROL::Ptr<ROL::Constraint<real_t>> ineq_constraint;
      ROL::Ptr<std::vector<real_t> > ll_ptr = ROL::makePtr<std::vector<real_t>>(1,Function_Arguments[imodule][0]);
      ROL::Ptr<std::vector<real_t> > lu_ptr = ROL::makePtr<std::vector<real_t>>(1,Function_Arguments[imodule][1]);    
      ROL::Ptr<ROL::Vector<real_t> > ll = ROL::makePtr<ROL::StdVector<real_t>>(ll_ptr);
      ROL::Ptr<ROL::Vector<real_t> > lu = ROL::makePtr<ROL::StdVector<real_t>>(lu_ptr);
      ROL::Ptr<ROL::BoundConstraint<real_t>> constraint_bnd = ROL::makePtr<ROL::Bounds<real_t>>(ll,lu);
      if(Optimization_Module_List[imodule] == OPTIMIZATION_MODULE_TYPE::Mass_Constraint_TopOpt){
        *fos << " MASS CONSTRAINT EXPECTS FEA MODULE INDEX " <<Optimization_Module_My_FEA_Module[imodule] << std::endl;
        ineq_constraint = ROL::makePtr<MassConstraint_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[imodule]], nodal_density_flag);
      }
      else if(Optimization_Module_List[imodule] == OPTIMIZATION_MODULE_TYPE::Displacement_Constraint_TopOpt){
        //simple test of assignment to the vector for constraint dofs
        Teuchos::RCP<MV> target_displacements = Teuchos::rcp(new MV(local_dof_map, 1));
        Teuchos::RCP<Tpetra::MultiVector<int,LO,GO>> active_dofs = Teuchos::rcp(new Tpetra::MultiVector<int,LO,GO>(local_dof_map, 1));
        active_dofs->putScalar(0);
        host_vec_array target_displacements_view = target_displacements->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
        host_ivec_array active_dofs_view = active_dofs->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
        *fos << " DISPLACEMENT CONSTRAINT EXPECTS FEA MODULE INDEX " <<Optimization_Module_My_FEA_Module[imodule] << std::endl;
        ineq_constraint = ROL::makePtr<DisplacementConstraint_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[imodule]], nodal_density_flag, target_displacements, active_dofs);
      }
      else if(Optimization_Module_List[imodule] == OPTIMIZATION_MODULE_TYPE::Moment_of_Inertia_Constraint_TopOpt){
        *fos << " MOMENT OF INERTIA CONSTRAINT EXPECTS FEA MODULE INDEX " <<Optimization_Module_My_FEA_Module[imodule] << std::endl;
        ineq_constraint = ROL::makePtr<MomentOfInertiaConstraint_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[imodule]], nodal_density_flag, Function_Arguments[imodule][1], Function_Arguments[imodule][0]);
      }
      else if(Optimization_Module_List[imodule] == OPTIMIZATION_MODULE_TYPE::Center_of_Mass_Constraint_TopOpt){
        *fos << " CENTER OF MASS CONSTRAINT EXPECTS FEA MODULE INDEX " <<Optimization_Module_My_FEA_Module[imodule] << std::endl;
        ineq_constraint = ROL::makePtr<CenterOfMassConstraint_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[imodule]], nodal_density_flag, Function_Arguments[imodule][1], Function_Arguments[imodule][0]);
      }
      else if(Optimization_Module_List[imodule] == OPTIMIZATION_MODULE_TYPE::Strain_Energy_Constraint_TopOpt){
        *fos << " STRAIN ENERGY CONSTRAINT EXPECTS FEA MODULE INDEX " <<Optimization_Module_My_FEA_Module[imodule] << std::endl;
        ineq_constraint = ROL::makePtr<StrainEnergyConstraint_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[imodule]], nodal_density_flag, Function_Arguments[imodule][2]);
      }
      else if(Optimization_Module_List[imodule] == OPTIMIZATION_MODULE_TYPE::Heat_Capacity_Potential_Constraint_TopOpt){
        *fos << " HEAT CAPACITY POTENTIAL CONSTRAINT EXPECTS FEA MODULE INDEX " <<Optimization_Module_My_FEA_Module[imodule] << std::endl;
        ineq_constraint = ROL::makePtr<HeatCapacityPotentialConstraint_TopOpt>(fea_modules[Optimization_Module_My_FEA_Module[imodule]], nodal_density_flag, Function_Arguments[imodule][2]);
      }
      else{
        *fos << "PROGRAM IS ENDING DUE TO ERROR; UNDEFINED INEQUALITY CONSTRAINT FUNCTION REQUESTED WITH NAME \"" 
              << to_string(Optimization_Module_List[imodule]) <<"\"" << std::endl;
        exit_solver(0);
      }
      *fos << " ADDING CONSTRAINT " << constraint_name << std::endl;
      problem->addConstraint(constraint_name, ineq_constraint, constraint_mul, constraint_bnd);
    }
    number_union.clear();
  }

  
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
  //real_t tol = 0.0001;
  //obj->update(*rol_x,ROL::UpdateType::Initial);
  //real_t obj_value = obj->value(*rol_x,tol);
  //std::cout << " VALUE TEST " << obj_value << std::endl;
  if(simparam.optimization_options.check_objective_gradient){
    obj->checkGradient(*rol_x, *rol_d);
  }
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

/* ----------------------------------------------------------------------------
   Initialize sets of element boundary surfaces and arrays for input conditions
------------------------------------------------------------------------------- */

void Implicit_Solver::init_boundaries(){
  int num_boundary_sets = 0;
  int num_dim = simparam.num_dims;

  // build boundary mesh patches
  if(myrank == 0)
    std::cout << "Starting boundary patch setup" << std::endl <<std::flush;
  Get_Boundary_Patches();
  //std::cout << "Done with boundary patch setup" << std::endl <<std::flush;
  //std::cout << "number of boundary patches on task " << myrank << " = " << nboundary_patches << std::endl;
  
  //disable for now
  if(0) {
    throw std::runtime_error("simparam.NB is not set. Figure out the number of BCs some other way.");
    init_topology_conditions(num_boundary_sets);
  }
}

/* ----------------------------------------------------------------------
   Assign sets of element boundary surfaces corresponding to user BCs
------------------------------------------------------------------------- */

void Implicit_Solver::generate_tcs(){
  int num_dim = simparam.num_dims;
  int bdy_set_id;
  int tc_tag;
  int num_topology_conditions = 0;
  real_t value;
  real_t fix_limits[4];
  
  if(num_topology_conditions>0) {
    throw std::runtime_error("simparam.NB is not set. Figure out the number of BCs some other way.");
    init_topology_conditions(num_topology_conditions);
  }

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
  int num_dim = simparam.num_dims;
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
  int num_dim = simparam.num_dims;

  //collect nodal coordinate information
  if(myrank==0) nreduce_nodes = num_nodes;
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > global_reduce_map =
    Teuchos::rcp(new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),nreduce_nodes,0,comm));
  
  //importer from local node distribution to collected distribution
  Tpetra::Import<LO, GO> node_collection_importer(map, global_reduce_map);

  collected_node_coords_distributed = Teuchos::rcp(new MV(global_reduce_map, num_dim));

  //comms to collect
  collected_node_coords_distributed->doImport(*node_coords_distributed, node_collection_importer, Tpetra::INSERT);

  //comms to collect FEA module related vector data
  for (int imodule = 0; imodule < nfea_modules; imodule++){
    fea_modules[imodule]->collect_output(global_reduce_map);
    //collected_node_displacements_distributed->doImport(*(fea_elasticity->node_displacements_distributed), dof_collection_importer, Tpetra::INSERT);
  }
  

  //collected nodal density information
  collected_node_densities_distributed = Teuchos::rcp(new MV(global_reduce_map, 1));

  //comms to collect
  collected_node_densities_distributed->doImport(*design_node_densities_distributed, node_collection_importer, Tpetra::INSERT);

  //collect element connectivity data
  if(myrank==0) nreduce_elem = num_elem;
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > global_reduce_element_map =
    Teuchos::rcp(new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),nreduce_elem,0,comm));

  //importer from all element map to collected distribution
  Tpetra::Import<LO, GO> element_collection_importer(all_element_map, global_reduce_element_map);
  
  collected_nodes_in_elem_distributed = Teuchos::rcp(new MCONN(global_reduce_element_map, max_nodes_per_element));

  //comms to collect
  collected_nodes_in_elem_distributed->doImport(*global_nodes_in_elem_distributed, element_collection_importer, Tpetra::INSERT);

  //collect element type data

}

/* ----------------------------------------------------------------------
   Sort Information for parallel file output
------------------------------------------------------------------------- */

void Implicit_Solver::sort_information(bool mesh_conversion_flag){
  GO nreduce_nodes = 0;
  GO nreduce_elem = 0;
  int num_dim = simparam.num_dims;

  //map with sorted continuous global indices
  sorted_map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(num_nodes,0,comm));
  
  //importer from local node distribution to sorted distribution
  Tpetra::Import<LO, GO> node_sorting_importer(map, sorted_map);

  sorted_node_coords_distributed = Teuchos::rcp(new MV(sorted_map, num_dim));

  //comms to collect
  sorted_node_coords_distributed->doImport(*node_coords_distributed, node_sorting_importer, Tpetra::INSERT);

  if(!mesh_conversion_flag){
    //comms to collect FEA module related vector data
    for (int imodule = 0; imodule < nfea_modules; imodule++){
      fea_modules[imodule]->sort_output(sorted_map);
      //collected_node_displacements_distributed->doImport(*(fea_elasticity->node_displacements_distributed), dof_collection_importer, Tpetra::INSERT);
    }
    

    //collected nodal density information
    sorted_node_densities_distributed = Teuchos::rcp(new MV(sorted_map, 1));

    //comms to collect
    sorted_node_densities_distributed->doImport(*design_node_densities_distributed, node_sorting_importer, Tpetra::INSERT);
  }
  //sort element connectivity
  sorted_element_map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(num_elem,0,comm));
  Tpetra::Import<LO, GO> element_sorting_importer(all_element_map, sorted_element_map);
  
  sorted_nodes_in_elem_distributed = Teuchos::rcp(new MCONN(sorted_element_map, max_nodes_per_element));

  //comms
  sorted_nodes_in_elem_distributed->doImport(*global_nodes_in_elem_distributed, element_sorting_importer, Tpetra::INSERT);

}

/* ----------------------------------------------------------------------
   Output Model Information in tecplot format
------------------------------------------------------------------------- */

void Implicit_Solver::output_design(int current_step, bool mesh_conversion_flag){
  if(current_step!=last_print_step){
    last_print_step = current_step;
    parallel_tecplot_writer(mesh_conversion_flag);
  }

}

/* ----------------------------------------------------------------------
   Output Model Information in tecplot format
------------------------------------------------------------------------- */

void Implicit_Solver::parallel_tecplot_writer(bool mesh_conversion_flag){
  
  int num_dim = simparam.num_dims;
	std::string current_file_name;
	std::string base_file_name= "TecplotTO";
  std::string base_file_name_undeformed= "TecplotTO_undeformed";
  std::stringstream current_line_stream;
  std::string current_line;
	std::string file_extension= ".dat";
  std::string file_count;
	std::stringstream count_temp;
  int temp_convert;
  int noutput, nvector;
  bool displace_geometry = false;
  int displacement_index;
  if(displacement_module!=-1){
    displace_geometry = simparam.output(FIELD::displaced_mesh);
    displacement_index = fea_modules[displacement_module]->displacement_index;
  }
  const_host_vec_array current_sorted_output;
  for (int imodule = 0; imodule < nfea_modules; imodule++){
    fea_modules[imodule]->compute_output();
  }
  
  sort_information(mesh_conversion_flag);
  const_host_vec_array sorted_node_coords = sorted_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_vec_array sorted_node_densities;
  if(!mesh_conversion_flag)
    sorted_node_densities = sorted_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array sorted_nodes_in_elem = sorted_nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
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
  
  MPI_File myfile_parallel, myfile_parallel_deformed;
  MPI_Offset header_stream_offset = 0;
  //initial undeformed geometry
  count_temp.str("");
  count_temp << file_index;
  if(!displace_geometry) {
    file_index++;
  }
	file_count = count_temp.str();
  if(displace_geometry&&displacement_module>=0)
    current_file_name = base_file_name_undeformed + file_count + file_extension;
  else
    current_file_name = base_file_name + file_count + file_extension;
  
  if(mesh_conversion_flag) current_file_name = "Converted_Tecplot_Mesh.dat";

  MPI_File_open(MPI_COMM_WORLD, current_file_name.c_str(), 
                MPI_MODE_CREATE|MPI_MODE_WRONLY, 
                MPI_INFO_NULL, &myfile_parallel);
  
  int err = MPI_File_open(MPI_COMM_WORLD, current_file_name.c_str(), MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &myfile_parallel);
  //allows overwriting the file if it already existed in the directory
  if (err != MPI_SUCCESS)  {
    if (myrank == 0){
      MPI_File_delete(current_file_name.c_str(),MPI_INFO_NULL);
    }
    MPI_File_open(MPI_COMM_WORLD, current_file_name.c_str(), MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &myfile_parallel);
  }
		
	//output header of the tecplot file

  //std::cout << current_file_name << std::endl;
	current_line_stream << "TITLE=\"results for TO simulation\"" "\n";
  current_line = current_line_stream.str();
  if(myrank == 0)
    MPI_File_write(myfile_parallel,current_line.c_str(),current_line.length(), MPI_CHAR, MPI_STATUS_IGNORE);
  header_stream_offset += current_line.length();
  //myfile << "VARIABLES = \"x\", \"y\", \"z\", \"density\", \"sigmaxx\", \"sigmayy\", \"sigmazz\", \"sigmaxy\", \"sigmaxz\", \"sigmayz\"" "\n";
  //else
  current_line_stream.str("");
  if(mesh_conversion_flag){
    if(num_dim == 2)
      current_line_stream << "VARIABLES = \"x\", \"y\"";
    else if(num_dim == 3)
      current_line_stream << "VARIABLES = \"x\", \"y\", \"z\"";
  }
  else{
    if(num_dim == 2)
      current_line_stream << "VARIABLES = \"x\", \"y\", \"density\"";
    else if(num_dim == 3)
      current_line_stream << "VARIABLES = \"x\", \"y\", \"z\", \"density\"";

    for (int imodule = 0; imodule < nfea_modules; imodule++){
      for(int ioutput = 0; ioutput < fea_modules[imodule]->noutput; ioutput++){
        nvector = fea_modules[imodule]->output_vector_sizes[ioutput];
        for(int ivector = 0; ivector < nvector; ivector++){
          current_line_stream << ", \"" << fea_modules[imodule]->output_dof_names[ioutput][ivector] << "\"";
        }
      }
    }
  }
  
  current_line = current_line_stream.str();
  if(myrank == 0)
    MPI_File_write(myfile_parallel,current_line.c_str(),current_line.length(), MPI_CHAR, MPI_STATUS_IGNORE);
  header_stream_offset += current_line.length();
    /*
    for (int imodule = 0; imodule < nfea_modules; imodule++){
      for(int ioutput = 0; ioutput < fea_modules[imodule]->noutput; ioutput++){
        nvector = fea_modules[imodule]->output_vector_sizes[ioutput];
        for(int ivector = 0; ivector < nvector; ivector++){
          myfile << ", \"" << fea_modules[imodule]->output_dof_names[ioutput][ivector] << "\"";
        }
      }
    }
    */
  current_line_stream.str("");
	current_line_stream << "\n";
  current_line = current_line_stream.str();
  if(myrank == 0)
    MPI_File_write(myfile_parallel,current_line.c_str(),current_line.length(), MPI_CHAR, MPI_STATUS_IGNORE);
  header_stream_offset += current_line.length();

	current_line_stream.str("");
	if(num_dim==2){
	  current_line_stream << "ZONE T=\"load step " << file_index << "\", NODES= " << num_nodes
		  << ", ELEMENTS= " << num_elem << ", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" "\n";
  }
  else if(num_dim==3){
   	current_line_stream << "ZONE T=\"load step " << file_index << "\", NODES= " << num_nodes
		<< ", ELEMENTS= " << num_elem << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" "\n"; 
  }
  current_line = current_line_stream.str();
  if(myrank == 0)
    MPI_File_write(myfile_parallel,current_line.c_str(),current_line.length(), MPI_CHAR, MPI_STATUS_IGNORE);
  header_stream_offset += current_line.length();

  //output nodal data
  //compute buffer output size and file stream offset for this MPI rank
  int default_dof_count = num_dim;
  if(!mesh_conversion_flag) default_dof_count++;
  int buffer_size_per_node_line = 26*default_dof_count + 1; //25 width + 1 space per number plus line terminator
  if(!mesh_conversion_flag){
    for (int imodule = 0; imodule < nfea_modules; imodule++){
      noutput = fea_modules[imodule]->noutput;
      for(int ioutput = 0; ioutput < noutput; ioutput++){
        nvector = fea_modules[imodule]->output_vector_sizes[ioutput];  
        buffer_size_per_node_line += 26*nvector;
      }
    }
  }
  int nlocal_sorted_nodes = sorted_map->getLocalNumElements();
  GO first_node_global_id = sorted_map->getGlobalElement(0);
  CArrayKokkos<char, array_layout, HostSpace, memory_traits> print_buffer(buffer_size_per_node_line*nlocal_sorted_nodes);
  MPI_Offset file_stream_offset = buffer_size_per_node_line*first_node_global_id;

  //populate buffer
  long long current_buffer_position = 0;
  current_line_stream << std::fixed << std::setprecision(8);
  for (int nodeline = 0; nodeline < nlocal_sorted_nodes; nodeline++) {
    current_line_stream.str("");
		current_line_stream << std::setw(25) << sorted_node_coords(nodeline,0) << " ";
		current_line_stream << std::setw(25) << sorted_node_coords(nodeline,1) << " ";
    if(num_dim==3)
		current_line_stream << std::setw(25) << sorted_node_coords(nodeline,2) << " ";

    //velocity print
    if(!mesh_conversion_flag){
      current_line_stream << std::setw(25) << sorted_node_densities(nodeline,0) << " ";

      for (int imodule = 0; imodule < nfea_modules; imodule++){
        noutput = fea_modules[imodule]->noutput;
        for(int ioutput = 0; ioutput < noutput; ioutput++){
          current_sorted_output = fea_modules[imodule]->module_outputs[ioutput];
          nvector = fea_modules[imodule]->output_vector_sizes[ioutput];
          if(fea_modules[imodule]->vector_style[ioutput] == FEA_Module::DOF){
            for(int ivector = 0; ivector < nvector; ivector++){
            current_line_stream << std::setw(25) << current_sorted_output(nodeline*nvector + ivector,0) << " ";
            }
          }
          if(fea_modules[imodule]->vector_style[ioutput] == FEA_Module::NODAL){
            for(int ivector = 0; ivector < nvector; ivector++){
              current_line_stream << std::setw(25) << current_sorted_output(nodeline,ivector) << " ";
            }
          }
        }
      }
    }
    current_line_stream << std::endl;
    
    current_line = current_line_stream.str();

    //copy current line over to C style string buffer (wrapped by matar)
    strcpy(&print_buffer(current_buffer_position),current_line.c_str());

    current_buffer_position += current_line.length();
	}

	//print buffers at offsets with collective MPI write
  //MPI_Offset current_stream_position = MPI_File_get_position(myfile_parallel,0);
  MPI_File_write_at_all(myfile_parallel, file_stream_offset + header_stream_offset, print_buffer.get_kokkos_view().data(), buffer_size_per_node_line*nlocal_sorted_nodes, MPI_CHAR, MPI_STATUS_IGNORE);
  //MPI_File_close(&myfile_parallel);
  
  MPI_Offset current_stream_position;
  MPI_Barrier(world);
  MPI_File_sync(myfile_parallel);
  MPI_File_seek_shared(myfile_parallel, 0, MPI_SEEK_END);
  MPI_File_sync(myfile_parallel);
  MPI_File_get_position_shared(myfile_parallel, &current_stream_position);
  
  //debug check 
  //std::cout << "offset on rank " << myrank << " is " << file_stream_offset + header_stream_offset + current_buffer_position << std::endl;
  //std::cout << "get position on rank " << myrank << " is " << current_stream_position << std::endl;
  
  //expand print buffer if needed
  int buffer_size_per_element_line = 11*max_nodes_per_element + 1; //25 width per number plus 6 spaces plus line terminator
  int nlocal_elements = sorted_element_map->getLocalNumElements();
  GO first_element_global_id = sorted_element_map->getGlobalElement(0);
  if(buffer_size_per_element_line*nlocal_elements > print_buffer.size())
    print_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(buffer_size_per_element_line*nlocal_elements);
  file_stream_offset = buffer_size_per_element_line*first_element_global_id + current_stream_position;
  
  current_buffer_position = 0;
  for (int elementline = 0; elementline < nlocal_elements; elementline++) {
    current_line_stream.str("");
    //convert node ordering
		for (int ii = 0; ii < max_nodes_per_element; ii++) {
      if(active_node_ordering_convention == IJK)
        temp_convert = convert_ijk_to_ensight(ii);
      else
        temp_convert = ii;
				current_line_stream << std::setw(10) << sorted_nodes_in_elem(elementline, temp_convert) + 1 << " ";
		}
		current_line_stream << std::endl;
    current_line = current_line_stream.str();

    //copy current line over to C style string buffer (wrapped by matar)
    strcpy(&print_buffer(current_buffer_position),current_line.c_str());

    current_buffer_position += current_line.length();
	}

  MPI_Barrier(world);
  MPI_File_write_at_all(myfile_parallel, file_stream_offset, print_buffer.get_kokkos_view().data(), buffer_size_per_element_line*nlocal_elements, MPI_CHAR, MPI_STATUS_IGNORE);

  if(simparam.output_options.optimization_restart_file){
      // Write commented restart data to be used by Fierro
      MPI_Offset current_stream_position;
      MPI_Barrier(world);
      MPI_File_sync(myfile_parallel);
      MPI_File_seek_shared(myfile_parallel, 0, MPI_SEEK_END);
      MPI_File_sync(myfile_parallel);
      MPI_File_get_position_shared(myfile_parallel, &current_stream_position);
      current_line_stream.str("");
      current_line_stream << std::endl << "#RESTART DATA: Objective_Normalization_Constant " <<
                  simparam.optimization_options.objective_normalization_constant << std::endl;
      if (myrank == 0)
      {
          MPI_File_write_at(myfile_parallel, current_stream_position, current_line_stream.str().c_str(),
                      current_line_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
      }
  }
  
  MPI_File_close(&myfile_parallel);
  MPI_Barrier(world);

  //Displaced Geometry File option
  if(displacement_module>=0&&(displace_geometry&&!mesh_conversion_flag)){
    current_line_stream.str("");
    header_stream_offset = 0;
    //deformed geometry
    count_temp.str("");
    count_temp << file_index;
    file_index++;
	  file_count = count_temp.str();
    
    current_file_name = base_file_name + file_count + file_extension;
    MPI_File_open(MPI_COMM_WORLD, current_file_name.c_str(), 
              MPI_MODE_CREATE|MPI_MODE_WRONLY, 
              MPI_INFO_NULL, &myfile_parallel_deformed);
  
    int err = MPI_File_open(MPI_COMM_WORLD, current_file_name.c_str(), MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &myfile_parallel_deformed);
    //allows overwriting the file if it already existed in the directory
    if (err != MPI_SUCCESS)  {
      if (myrank == 0){
        MPI_File_delete(current_file_name.c_str(),MPI_INFO_NULL);
      }
      MPI_File_open(MPI_COMM_WORLD, current_file_name.c_str(), MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &myfile_parallel_deformed);
    }
		
	  //output header of the tecplot file
	  current_line_stream << "TITLE=\"results for TO simulation\"" "\n";
    current_line = current_line_stream.str();
    if(myrank == 0)
      MPI_File_write(myfile_parallel_deformed,current_line.c_str(),current_line.length(), MPI_CHAR, MPI_STATUS_IGNORE);
    header_stream_offset += current_line.length();
    //myfile << "VARIABLES = \"x\", \"y\", \"z\", \"density\", \"sigmaxx\", \"sigmayy\", \"sigmazz\", \"sigmaxy\", \"sigmaxz\", \"sigmayz\"" "\n";
    //else
    current_line_stream.str("");
	  current_line_stream << "VARIABLES = \"x\", \"y\", \"z\", \"density\"";
    for (int imodule = 0; imodule < nfea_modules; imodule++){
      for(int ioutput = 0; ioutput < fea_modules[imodule]->noutput; ioutput++){
        nvector = fea_modules[imodule]->output_vector_sizes[ioutput];
        for(int ivector = 0; ivector < nvector; ivector++){
          current_line_stream << ", \"" << fea_modules[imodule]->output_dof_names[ioutput][ivector] << "\"";
        }
      }
    }
    current_line = current_line_stream.str();
    if(myrank == 0)
      MPI_File_write(myfile_parallel_deformed,current_line.c_str(),current_line.length(), MPI_CHAR, MPI_STATUS_IGNORE);
    header_stream_offset += current_line.length();
    /*
    for (int imodule = 0; imodule < nfea_modules; imodule++){
      for(int ioutput = 0; ioutput < fea_modules[imodule]->noutput; ioutput++){
        nvector = fea_modules[imodule]->output_vector_sizes[ioutput];
        for(int ivector = 0; ivector < nvector; ivector++){
          myfile << ", \"" << fea_modules[imodule]->output_dof_names[ioutput][ivector] << "\"";
        }
      }
    }
    */
    current_line_stream.str("");
	  current_line_stream << "\n";
    current_line = current_line_stream.str();
    if(myrank == 0)
      MPI_File_write(myfile_parallel_deformed,current_line.c_str(),current_line.length(), MPI_CHAR, MPI_STATUS_IGNORE);
    header_stream_offset += current_line.length();

	  current_line_stream.str("");
	  if(num_dim==2){
	  current_line_stream << "ZONE T=\"load step " << file_index << "\", NODES= " << num_nodes
		  << ", ELEMENTS= " << num_elem << ", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" "\n";
  }
  else if(num_dim==3){
   	current_line_stream << "ZONE T=\"load step " << file_index << "\", NODES= " << num_nodes
		<< ", ELEMENTS= " << num_elem << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" "\n"; 
  }
    current_line = current_line_stream.str();
    if(myrank == 0)
      MPI_File_write(myfile_parallel_deformed,current_line.c_str(),current_line.length(), MPI_CHAR, MPI_STATUS_IGNORE);
    header_stream_offset += current_line.length();

		//output nodal data
    //compute buffer output size and file stream offset for this MPI rank
    int default_dof_count = num_dim;
    default_dof_count++;
    int buffer_size_per_node_line = 26*default_dof_count + 1; //25 width + 1 space per number plus line terminator
    for (int imodule = 0; imodule < nfea_modules; imodule++){
      noutput = fea_modules[imodule]->noutput;
      for(int ioutput = 0; ioutput < noutput; ioutput++){
        nvector = fea_modules[imodule]->output_vector_sizes[ioutput];  
        buffer_size_per_node_line += 26*nvector;
      }
    }
    int nlocal_sorted_nodes = sorted_map->getLocalNumElements();
    GO first_node_global_id = sorted_map->getGlobalElement(0);
    CArrayKokkos<char, array_layout, HostSpace, memory_traits> print_buffer(buffer_size_per_node_line*nlocal_sorted_nodes);
    MPI_Offset file_stream_offset = buffer_size_per_node_line*first_node_global_id;
    //populate buffer
    long long current_buffer_position = 0;
    current_line_stream << std::fixed << std::setprecision(8);
    for (int nodeline = 0; nodeline < nlocal_sorted_nodes; nodeline++) {
      current_sorted_output = fea_modules[displacement_module]->module_outputs[displacement_index];
      current_line_stream.str("");
		  current_line_stream << std::setw(25) << sorted_node_coords(nodeline,0) + current_sorted_output(nodeline*num_dim,0) << " ";
		  current_line_stream << std::setw(25) << sorted_node_coords(nodeline,1) + current_sorted_output(nodeline*num_dim + 1,0) << " ";
      if(num_dim==3)
		  current_line_stream << std::setw(25) << sorted_node_coords(nodeline,2) + current_sorted_output(nodeline*num_dim + 2,0) << " ";

      //velocity print
      current_line_stream << std::setw(25) << sorted_node_densities(nodeline,0) << " ";
    
      for (int imodule = 0; imodule < nfea_modules; imodule++){
        noutput = fea_modules[imodule]->noutput;
        for(int ioutput = 0; ioutput < noutput; ioutput++){
          current_sorted_output = fea_modules[imodule]->module_outputs[ioutput];
          nvector = fea_modules[imodule]->output_vector_sizes[ioutput];
          if(fea_modules[imodule]->vector_style[ioutput] == FEA_Module::DOF){
            for(int ivector = 0; ivector < nvector; ivector++){
             current_line_stream << std::setw(25) << current_sorted_output(nodeline*nvector + ivector,0) << " ";
            }
          }
          if(fea_modules[imodule]->vector_style[ioutput] == FEA_Module::NODAL){
            for(int ivector = 0; ivector < nvector; ivector++){
              current_line_stream << std::setw(25) << current_sorted_output(nodeline,ivector) << " ";
            }
          }
        }
      }
      current_line_stream << std::endl;
    
      current_line = current_line_stream.str();

      //copy current line over to C style string buffer (wrapped by matar)
      strcpy(&print_buffer(current_buffer_position),current_line.c_str());

      current_buffer_position += current_line.length();
	  }
		//print buffers at offsets with collective MPI write
    //MPI_Offset current_stream_position = MPI_File_get_position(myfile_parallel,0);
    MPI_Barrier(world);
    MPI_File_write_at_all(myfile_parallel_deformed, file_stream_offset + header_stream_offset, print_buffer.get_kokkos_view().data(), buffer_size_per_node_line*nlocal_sorted_nodes, MPI_CHAR, MPI_STATUS_IGNORE);
    //MPI_File_close(&myfile_parallel);
  
    MPI_Offset current_stream_position;
    MPI_Barrier(world);
    MPI_File_sync(myfile_parallel_deformed);
    MPI_File_seek_shared(myfile_parallel_deformed, 0, MPI_SEEK_END);
    MPI_File_sync(myfile_parallel_deformed);
    MPI_File_get_position_shared(myfile_parallel_deformed, &current_stream_position);
  
    //debug check 
    //std::cout << "offset on rank " << myrank << " is " << file_stream_offset + header_stream_offset + current_buffer_position << std::endl;
    //std::cout << "get position on rank " << myrank << " is " << current_stream_position << std::endl;
  
    //expand print buffer if needed
    int buffer_size_per_element_line = 11*max_nodes_per_element + 1; //25 width per number plus 6 spaces plus line terminator
    int nlocal_elements = sorted_element_map->getLocalNumElements();
    GO first_element_global_id = sorted_element_map->getGlobalElement(0);
    if(buffer_size_per_element_line*nlocal_elements > print_buffer.size())
      print_buffer = CArrayKokkos<char, array_layout, HostSpace, memory_traits>(buffer_size_per_element_line*nlocal_elements);
    file_stream_offset = buffer_size_per_element_line*first_element_global_id + current_stream_position;
  
    current_buffer_position = 0;
    for (int elementline = 0; elementline < nlocal_elements; elementline++) {
      current_line_stream.str("");
      //convert node ordering
		  for (int ii = 0; ii < max_nodes_per_element; ii++) {
        if(active_node_ordering_convention == IJK)
          temp_convert = convert_ijk_to_ensight(ii);
        else
          temp_convert = ii;
				  current_line_stream << std::setw(10) << sorted_nodes_in_elem(elementline, temp_convert) + 1 << " ";
		  }
		  current_line_stream << std::endl;
      current_line = current_line_stream.str();

      //copy current line over to C style string buffer (wrapped by matar)
      strcpy(&print_buffer(current_buffer_position),current_line.c_str());

      current_buffer_position += current_line.length();
	  }
    
    MPI_Barrier(world);
    MPI_File_write_at_all(myfile_parallel_deformed, file_stream_offset, print_buffer.get_kokkos_view().data(), buffer_size_per_element_line*nlocal_elements, MPI_CHAR, MPI_STATUS_IGNORE);
    
    MPI_Barrier(world);
    MPI_File_sync(myfile_parallel);
    MPI_File_close(&myfile_parallel_deformed);
  }
}

/* ----------------------------------------------------------------------
   Output Model Information in tecplot format
------------------------------------------------------------------------- */

void Implicit_Solver::tecplot_writer(){
  
  int num_dim = simparam.num_dims;
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
  bool displace_geometry = false;
  int displacement_index;
  if(displacement_module!=-1){
    displace_geometry = simparam.output(FIELD::displaced_mesh);
    displacement_index = fea_modules[displacement_module]->displacement_index;
  }
  const_host_vec_array current_collected_output;
  for (int imodule = 0; imodule < nfea_modules; imodule++){
    fea_modules[imodule]->compute_output();
  }
  
  collect_information();
  const_host_vec_array collected_node_coords = collected_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_vec_array collected_node_densities = collected_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array collected_nodes_in_elem = collected_nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  
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
      if(!displace_geometry) {
        file_index++;
      }
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

		  if(num_dim==2){
		    myfile << "ZONE T=\"load step " << file_index << "\", NODES= " << num_nodes
			    << ", ELEMENTS= " << num_elem << ", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" "\n";
      }
      else if(num_dim==3){
		    myfile << "ZONE T=\"load step " << file_index << "\", NODES= " << num_nodes
			    << ", ELEMENTS= " << num_elem << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" "\n";
      }

		  for (int nodeline = 0; nodeline < num_nodes; nodeline++) {
			  myfile << std::setw(25) << collected_node_coords(nodeline,0) << " ";
			  myfile << std::setw(25) << collected_node_coords(nodeline,1) << " ";
        if(num_dim==3)
			  myfile << std::setw(25) << collected_node_coords(nodeline,2) << " ";
        myfile << std::setw(25) << collected_node_densities(nodeline,0) << " ";
        for (int imodule = 0; imodule < nfea_modules; imodule++){
          noutput = fea_modules[imodule]->noutput;
          for(int ioutput = 0; ioutput < noutput; ioutput++){
            current_collected_output = fea_modules[imodule]->module_outputs[ioutput];
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

		  if(num_dim==2){
		    myfile << "ZONE T=\"load step " << file_index << "\", NODES= " << num_nodes
			    << ", ELEMENTS= " << num_elem << ", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" "\n";
      }
      else if(num_dim==3){
		    myfile << "ZONE T=\"load step " << file_index << "\", NODES= " << num_nodes
			    << ", ELEMENTS= " << num_elem << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" "\n";
      }

		  for (int nodeline = 0; nodeline < num_nodes; nodeline++) {
        current_collected_output = fea_modules[displacement_module]->module_outputs[displacement_index];
			  myfile << std::setw(25) << collected_node_coords(nodeline,0) + current_collected_output(nodeline*num_dim,0) << " ";
			  myfile << std::setw(25) << collected_node_coords(nodeline,1) + current_collected_output(nodeline*num_dim + 1,0) << " ";
        if(num_dim==3)
			  myfile << std::setw(25) << collected_node_coords(nodeline,2) + current_collected_output(nodeline*num_dim + 2,0) << " ";
        myfile << std::setw(25) << collected_node_densities(nodeline,0) << " ";
        for (int imodule = 0; imodule < nfea_modules; imodule++){
          noutput = fea_modules[imodule]->noutput;
          for(int ioutput = 0; ioutput < noutput; ioutput++){
            current_collected_output = fea_modules[imodule]->module_outputs[ioutput];
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
    const_host_elem_conn_array nodes_in_elem = global_nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    int graphics_id = simparam.graphics_id;
    int num_dim = simparam.num_dims;

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
    const_host_elem_conn_array nodes_in_elem = global_nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    mat_pt_t *mat_pt = simparam.mat_pt;
    int &graphics_id = simparam.graphics_id;
    real_t *graphics_times = simparam.graphics_times;
    real_t &TIME = simparam.TIME;

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
  int num_dim = simparam.num_dims;
  bool nodal_density_flag = simparam.nodal_density_flag;

  //set densities
  if(nodal_density_flag){
    if(!simparam.restart_file){
      design_node_densities_distributed = Teuchos::rcp(new MV(map, 1));
      if(simparam.optimization_options.density_filter == DENSITY_FILTER::helmholtz_filter)
        filtered_node_densities_distributed = Teuchos::rcp(new MV(map, 1));
      host_vec_array node_densities = design_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
      //notify that the host view is going to be modified in the file readin
      //dual_node_densities.modify_host();
    
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
    //dual_node_densities.sync_device();
    }
    //allocate global vector information
    all_node_densities_distributed = Teuchos::rcp(new MV(all_node_map, 1));

    //communicate ghost information to the all vector
    //create import object using local node indices map and all indices map
    //Tpetra::Import<LO, GO> importer(map, all_node_map);

    //comms to get ghosts
    all_node_densities_distributed->doImport(*design_node_densities_distributed, *importer, Tpetra::INSERT);

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
