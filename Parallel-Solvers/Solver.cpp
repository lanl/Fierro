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
#include <set>
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

#include <Xpetra_MultiVector.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_Parallel_Reduce.hpp>
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"
#include "Tpetra_Details_FixedHashTable.hpp"
#include "Tpetra_Import.hpp"

#include "elements.h"
#include "matar.h"
#include "utilities.h"
#include "node_combination.h"
#include "Simulation_Parameters.h"
#include "Solver.h"

//Repartition Package
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_InputTraits.hpp>

#define BUFFER_LINES 20000
#define MAX_WORD 30
#define MAX_ELEM_NODES 32
#define STRAIN_EPSILON 0.000000001
#define DENSITY_EPSILON 0.0001
#define BC_EPSILON 1.0e-8

Solver::Solver(){
  //default flags assume optional routines are off
  setup_flag = finalize_flag = 0;
}

void Solver::exit_solver(int status){
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  exit(status);
}

Solver::~Solver(){}

/* ----------------------------------------------------------------------
   Read Ensight format mesh file
------------------------------------------------------------------------- */

void Solver::read_mesh_ensight(char *MESH){

  char ch;
  int num_dim = simparam->num_dim;
  std::cout << " NUM DIM IS " << num_dim << std::endl;
  int p_order = simparam->p_order;
  real_t unit_scaling = simparam->unit_scaling;
  int local_node_index, current_column_index;
  size_t strain_count;
  std::string skip_line, read_line, substring;
  std::stringstream line_parse;
  CArrayKokkos<char, array_layout, HostSpace, memory_traits> read_buffer;
  int buffer_loop, buffer_iteration, buffer_iterations, dof_limit, scan_loop;
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
    //skip 8 lines
    for (int j = 1; j <= 8; j++) {
      getline(*in, skip_line);
      std::cout << skip_line << std::endl;
    } //for
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
  nlocal_nodes = map->getLocalNumElements();
  //populate local row offset data from global data
  global_size_t min_gid = map->getMinGlobalIndex();
  global_size_t max_gid = map->getMaxGlobalIndex();
  global_size_t index_base = map->getIndexBase();
  //debug print
  //std::cout << "local node count on task: " << " " << nlocal_nodes << std::endl;

  //allocate node storage with dual view
  //dual_node_coords = dual_vec_array("dual_node_coords", nlocal_nodes,num_dim);

  //local variable for host view in the dual view
  
  node_coords_distributed = Teuchos::rcp(new MV(map, num_dim));

  //scope ensures view is destroyed for now to avoid calling a device view with an active host view later
  {
  host_vec_array node_coords = node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  //host_vec_array node_coords = dual_node_coords.view_host();
  //notify that the host view is going to be modified in the file readin
  //dual_node_coords.modify_host();

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

  dof_limit = num_nodes;
  buffer_iterations = dof_limit/BUFFER_LINES;
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

    //loop through read buffer and store coords in node coords view
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
        if(num_dim==3)
          node_coords(node_rid, 2) = dof_value * unit_scaling;
      }
    }
    read_index_start+=BUFFER_LINES;
  }
  } //end active view scope
  //repartition node distribution
  repartition_nodes();

  //synchronize device data
  //dual_node_coords.sync_device();
  //dual_node_coords.modify_device();

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

  if(simparam->num_dim==2){
    if(simparam->element_type == "Quad4"){
      mesh_element_type = elements::elem_types::Quad4;
    }
    else if(simparam->element_type == "Quad8"){
      mesh_element_type = elements::elem_types::Quad8;
    }
    else if(simparam->element_type == "Quad12"){
      mesh_element_type = elements::elem_types::Quad12;
    }
    element_select->choose_2Delem_type(mesh_element_type, elem2D);
    max_nodes_per_element = elem2D->num_nodes();
  }

  if(simparam->num_dim==3){
    if(simparam->element_type == "Hex8"){
      mesh_element_type = elements::elem_types::Hex8;
    }
    else if(simparam->element_type == "Hex20"){
      mesh_element_type = elements::elem_types::Hex20;
    }
    else if(simparam->element_type == "Hex32"){
      mesh_element_type = elements::elem_types::Hex32;
    }
    element_select->choose_3Delem_type(mesh_element_type, elem);
    max_nodes_per_element = elem->num_nodes();
  }

  //1 type per mesh for now
  for(int ielem = 0; ielem < rnum_elem; ielem++)
    Element_Types(ielem) = mesh_element_type;
  
  //copy temporary element storage to multivector storage
  dual_nodes_in_elem = dual_elem_conn_array("dual_nodes_in_elem", rnum_elem, max_nodes_per_element);
  host_elem_conn_array nodes_in_elem = dual_nodes_in_elem.view_host();
  dual_nodes_in_elem.modify_host();

  for(int ielem = 0; ielem < rnum_elem; ielem++)
    for(int inode = 0; inode < elem_words_per_line; inode++)
      nodes_in_elem(ielem, inode) = element_temp[ielem*elem_words_per_line + inode];

  //view storage for all local elements connected to local nodes on this rank
  //DCArrayKokkos<GO, array_layout, device_type, memory_traits> All_Element_Global_Indices(rnum_elem);
  Kokkos::DualView <GO*, array_layout, device_type, memory_traits> All_Element_Global_Indices("All_Element_Global_Indices",rnum_elem);
  //copy temporary global indices storage to view storage
  for(int ielem = 0; ielem < rnum_elem; ielem++)
    All_Element_Global_Indices.h_view(ielem) = global_indices_temp[ielem];

  //delete temporary element connectivity and index storage
  std::vector<size_t>().swap(element_temp);
  std::vector<size_t>().swap(global_indices_temp);
  
  //construct overlapping element map (since different ranks can own the same elements due to the local node map)
  All_Element_Global_Indices.modify_host();
  All_Element_Global_Indices.sync_device();

  //debug print
  /*
  Kokkos::View <GO*, array_layout, device_type, memory_traits> All_Element_Global_Indices_pass("All_Element_Global_Indices_pass",rnum_elem);
  deep_copy(All_Element_Global_Indices_pass, All_Element_Global_Indices.h_view);
  std::cout << " ------------ELEMENT GLOBAL INDICES ON TASK " << myrank << " --------------"<<std::endl;
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    std::cout << "elem: " << All_Element_Global_Indices_pass(ielem) + 1;
    std::cout << std::endl;
  }
  */

  all_element_map = Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),All_Element_Global_Indices.d_view,0,comm));

  //element type selection (subject to change)
  // ---- Set Element Type ---- //
  // allocate element type memory
  //elements::elem_type_t* elem_choice;

  int NE = 1; // number of element types in problem
    

  // Convert ensight index system to the ijk finite element numbering convention
  // for vertices in cell
  if(active_node_ordering_convention == IJK){
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
  }
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
   Rebalance the initial node decomposition with Zoltan2
------------------------------------------------------------------------- */

void Solver::repartition_nodes(){
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
  
  //construct input adapted needed by Zoltan2 problem
  typedef Xpetra::MultiVector<real_t,LO,GO,node_type> xvector_t;
  typedef Zoltan2::XpetraMultiVectorAdapter<xvector_t> inputAdapter_t;
  typedef Zoltan2::EvaluatePartition<inputAdapter_t> quality_t;
  
  Teuchos::RCP<xvector_t> xpetra_node_coords = Teuchos::rcp(new Xpetra::TpetraMultiVector<real_t,LO,GO,node_type>(node_coords_distributed));
  Teuchos::RCP<inputAdapter_t> problem_adapter =  Teuchos::rcp(new inputAdapter_t(xpetra_node_coords));

  // Create parameters for an RCB problem

  double tolerance = 1.05;

  Teuchos::ParameterList params("Node Partition Params");
  params.set("debug_level", "basic_status");
  params.set("debug_procs", "0");
  params.set("error_check_level", "debug_mode_assertions");

  //params.set("algorithm", "rcb");
  params.set("algorithm", "multijagged");
  params.set("imbalance_tolerance", tolerance );
  params.set("num_global_parts", nranks);
  params.set("partitioning_objective", "minimize_cut_edge_count");
  
  Teuchos::RCP<Zoltan2::PartitioningProblem<inputAdapter_t> > problem =
           Teuchos::rcp(new Zoltan2::PartitioningProblem<inputAdapter_t>(&(*problem_adapter), &params));
   
  // Solve the problem

  problem->solve();

  // create metric object where communicator is Teuchos default

  quality_t *metricObject1 = new quality_t(&(*problem_adapter), &params, //problem1->getComm(),
					   &problem->getSolution());
  // Check the solution.

  if (myrank == 0) {
    metricObject1->printMetrics(std::cout);
  }

  if (myrank == 0){
    real_t imb = metricObject1->getObjectCountImbalance();
    if (imb <= tolerance)
      std::cout << "pass: " << imb << std::endl;
    else
      std::cout << "fail: " << imb << std::endl;
    std::cout << std::endl;
  }
  delete metricObject1;

  //migrate rows of the vector so they correspond to the partition recommended by Zoltan2
  Teuchos::RCP<MV> partitioned_node_coords_distributed = Teuchos::rcp(new MV(map,num_dim));
  Teuchos::RCP<xvector_t> xpartitioned_node_coords_distributed =
                          Teuchos::rcp(new Xpetra::TpetraMultiVector<real_t,LO,GO,node_type>(partitioned_node_coords_distributed));

  problem_adapter->applyPartitioningSolution(*xpetra_node_coords, xpartitioned_node_coords_distributed, problem->getSolution());
  *partitioned_node_coords_distributed = Xpetra::toTpetra<real_t,LO,GO,node_type>(*xpartitioned_node_coords_distributed);
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > partitioned_map = Teuchos::rcp(new Tpetra::Map<LO,GO,node_type>(*(partitioned_node_coords_distributed->getMap())));
  Teuchos::RCP<const Tpetra::Map<LO,GO,node_type> > partitioned_map_one_to_one;
  partitioned_map_one_to_one = Tpetra::createOneToOne<LO,GO,node_type>(partitioned_map);
  Teuchos::RCP<MV> partitioned_node_coords_one_to_one_distributed = Teuchos::rcp(new MV(partitioned_map_one_to_one,num_dim));

  Tpetra::Import<LO, GO> importer_one_to_one(partitioned_map, partitioned_map_one_to_one);
  partitioned_node_coords_one_to_one_distributed->doImport(*partitioned_node_coords_distributed, importer_one_to_one, Tpetra::INSERT);
  node_coords_distributed = partitioned_node_coords_one_to_one_distributed;
  partitioned_map = Teuchos::rcp(new Tpetra::Map<LO,GO,node_type>(*partitioned_map_one_to_one));

  //migrate density vector if this is a restart file read
  if(simparam->restart_file){
    Teuchos::RCP<MV> partitioned_node_densities_distributed = Teuchos::rcp(new MV(partitioned_map, 1));

    //create import object using local node indices map and all indices map
    Tpetra::Import<LO, GO> importer(map, partitioned_map);

    //comms to get ghosts
    partitioned_node_densities_distributed->doImport(*design_node_densities_distributed, importer, Tpetra::INSERT);
    design_node_densities_distributed = partitioned_node_densities_distributed;
  }

  //update nlocal_nodes and node map
  map = partitioned_map;
  nlocal_nodes = map->getLocalNumElements();
  
}
