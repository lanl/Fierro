
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
#include "FEA_Module_Heat_Conduction.h"
#include "Implicit_Solver.h"

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

#define MAX_ELEM_NODES 8
#define FLUX_EPSILON 0.000000001
#define DENSITY_EPSILON 0.0001
#define BC_EPSILON 1.0e-8

using namespace utils;


FEA_Module_Heat_Conduction::FEA_Module_Heat_Conduction(Implicit_Solver *Solver_Pointer) :FEA_Module(Solver_Pointer){
  //create parameter object
  simparam = new Simulation_Parameters();
  // ---- Read input file, define state and boundary conditions ---- //
  simparam->input();
  //create ref element object
  ref_elem = new elements::ref_element();
  //create mesh objects
  //init_mesh = new swage::mesh_t(simparam);
  //mesh = new swage::mesh_t(simparam);
  hessvec_count = update_count = 0;
  linear_solve_time = hessvec_time = hessvec_linear_time = 0;

  //preconditioner construction
  Hierarchy_Constructed = false;

  Matrix_alloc=0;
  gradient_print_sync = 0;

  //boundary condition data
  current_bdy_id = 0;

  //boundary condition flags
  body_term_flag = thermal_flag = electric_flag = false;

  //construct globally distributed temperature and heat flux vectors
  int num_dim = simparam->num_dim;
  local_dof_map = map;
  all_dof_map = all_node_map;
  node_temperatures_distributed = Teuchos::rcp(new MV(local_dof_map, 1));
  node_temperature_gradients_distributed = Teuchos::rcp(new MV(local_dof_map, 1));
  all_node_temperatures_distributed = Teuchos::rcp(new MV(all_dof_map, 1));
  node_heat_fluxes_distributed = Teuchos::rcp(new MV(map, num_dim));
  Global_Nodal_RHS = Teuchos::rcp(new MV(local_dof_map, 1));
  all_node_temperature_gradients_distributed = Teuchos::rcp(new MV(all_node_map, num_dim));
  all_node_heat_fluxes_distributed = Teuchos::rcp(new MV(all_node_map, num_dim));

  //initialize temperatures to 0
  //local variable for host view in the dual view
  host_vec_array all_node_temperatures = all_node_temperatures_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  host_vec_array node_temperatures = node_temperatures_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  for(int init = 0; init < local_dof_map->getNodeNumElements(); init++)
    node_temperatures(init,0) = 0;
  for(int init = 0; init < all_dof_map->getNodeNumElements(); init++)
    all_node_temperatures(init,0) = 0;
}

FEA_Module_Heat_Conduction::~FEA_Module_Heat_Conduction(){
   delete simparam;
}

/* ----------------------------------------------------------------------------
   Initialize sets of element boundary surfaces and arrays for input conditions
------------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::init_boundaries(){
  int num_boundary_sets = simparam->NB;
  int num_surface_flux_sets = simparam->NBSF;
  int num_surface_temp_sets = simparam->NBT;
  int num_dim = simparam->num_dim;
  
  // set the number of boundary sets
  if(myrank == 0)
    std::cout << "building boundary sets " << std::endl;
  
  init_boundary_sets(num_boundary_sets);
  Boundary_Condition_Type_List = CArrayKokkos<int, array_layout, HostSpace, memory_traits>(num_boundary_sets); 
  Boundary_Surface_Heat_Flux = CArrayKokkos<real_t, array_layout, HostSpace, memory_traits>(num_surface_flux_sets,3);
  Boundary_Surface_Temperatures = CArrayKokkos<real_t, array_layout, HostSpace, memory_traits>(num_surface_temp_sets,1);

  //initialize
  for(int ibdy=0; ibdy < num_boundary_sets; ibdy++) Boundary_Condition_Type_List(ibdy) = NONE;

  //allocate nodal data
  Node_DOF_Boundary_Condition_Type = CArrayKokkos<int, array_layout, device_type, memory_traits>(nall_nodes, "Node_DOF_Boundary_Condition_Type");
  Node_Temperature_Boundary_Conditions = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(nall_nodes);
  Node_DOF_Flux_Boundary_Conditions = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(nall_nodes);

  //initialize
  for(int init=0; init < nall_nodes*num_dim; init++)
    Node_DOF_Boundary_Condition_Type(init) = NONE;
}

/* ----------------------------------------------------------------------
   initialize storage for element boundary surfaces corresponding to user BCs
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::init_boundary_sets (int num_sets){
    
  num_boundary_conditions = num_sets;
  if(num_sets == 0){
    std::cout << " Warning: number of boundary conditions = 0";
    return;
  }
  NBoundary_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_sets, "NBoundary_Condition_Patches");
  Boundary_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_sets, nboundary_patches, "Boundary_Condition_Patches");

  //initialize data
  for(int iset = 0; iset < num_sets; iset++) NBoundary_Condition_Patches(iset) = 0;
}


/* ----------------------------------------------------------------------
   Assign sets of element boundary surfaces corresponding to user BCs
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::generate_bcs(){
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
  Boundary_Condition_Type_List(bdy_set_id) = TEMPERATURE_CONDITION;
  Boundary_Surface_Temperatures(surf_disp_set_id,0) = 0;
  Boundary_Surface_Temperatures(surf_disp_set_id,1) = 0;
  Boundary_Surface_Temperatures(surf_disp_set_id,2) = 0;
  surf_disp_set_id++;
    
  *fos << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
  *fos << std::endl;

  //Tag nodes for Boundary conditions such as displacements
  Temperature_Boundary_Conditions();
} // end generate_bcs

/* ----------------------------------------------------------------------
   Assign sets of element boundary surfaces corresponding to user BCs
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::generate_applied_loads(){
  int num_dim = simparam->num_dim;
  int bdy_set_id;
  int surf_flux_set_id = 0;
  int bc_tag;
  real_t value;
  
  //Surface Fluxes Section
  
  *fos << "tagging beam +z heat flux " << std::endl;
  bc_tag = 2;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  //value = 0;
  value = 100;
  bdy_set_id = current_bdy_id++;
  //find boundary patches this BC corresponds to
  tag_boundaries(bc_tag, value, bdy_set_id);
  Boundary_Condition_Type_List(bdy_set_id) = SURFACE_LOADING_CONDITION;
  Boundary_Surface_Heat_Flux(surf_force_set_id,0) = 0.5/simparam->unit_scaling/simparam->unit_scaling;
  Boundary_Surface_Heat_Flux(surf_force_set_id,1) = 0;
  Boundary_Surface_Heat_Flux(surf_force_set_id,2) = 0;
  surf_flux_set_id++;
  *fos << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(bdy_set_id) << std::endl;
  *fos << std::endl;
  
  //Body Term Section

  //apply gravity
  gravity_flag = simparam->gravity_flag;
  gravity_vector = simparam->gravity_vector;

  if(electric_flag||body_force_flag||thermal_flag) body_term_flag = true;

}

/* ----------------------------------------------------------------------
   Initialize global vectors and array maps needed for matrix assembly
------------------------------------------------------------------------- */
void FEA_Module_Heat_Conduction::init_assembly(){
  int num_dim = simparam->num_dim;
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Graph_Fill(nall_nodes, "nall_nodes");
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> current_row_nodes_scanned;
  int current_row_n_nodes_scanned;
  int local_node_index, global_node_index, current_column_index;
  int max_stride = 0;
  size_t nodes_per_element;
  
  //allocate stride arrays
  CArrayKokkos <size_t, array_layout, device_type, memory_traits> Graph_Matrix_Strides_initial(nlocal_nodes, "Graph_Matrix_Strides_initial");
  Graph_Matrix_Strides = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(nlocal_nodes, "Graph_Matrix_Strides");

  //allocate storage for the sparse conductivity matrix map used in the assembly process
  Global_Conductivity_Matrix_Assembly_Map = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(rnum_elem,
                                         max_nodes_per_element,max_nodes_per_element, "Global_Conductivity_Matrix_Assembly_Map");

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
          Global_Conductivity_Matrix_Assembly_Map(ielem,lnode,jnode) = current_column_index;
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
          Global_Conductivity_Matrix_Assembly_Map(ielem,lnode,jnode) = current_column_index;
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
        Global_Conductivity_Matrix_Assembly_Map(current_element_index,element_row_index, element_column_index) 
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

        Global_Conductivity_Matrix_Assembly_Map(current_element_index,element_row_index, element_column_index) 
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
    is used to loop over each element's local conductivity matrix in the assembly process.*/
  
  Conductivity_Matrix_Strides = Graph_Matrix_Strides;

  Conductivity_Matrix = RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout>(Conductivity_Matrix_Strides);
  DOF_Graph_Matrix = Graph_Matrix;
  
}

/* ----------------------------------------------------------------------
   Assemble the Sparse Stiffness Matrix
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::assemble_matrix(){
  int num_dim = simparam->num_dim;
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  int nodes_per_element;
  int current_row_n_nodes_scanned;
  int local_dof_index, global_node_index, current_row, current_column;
  int max_stride = 0;
  
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> Local_Conductivity_Matrix(max_nodes_per_element,max_nodes_per_element);

  //initialize conductivity Matrix entries to 0
  //debug print
    //std::cout << "DOF GRAPH MATRIX ENTRIES ON TASK " << myrank << std::endl;
  for (int idof = 0; idof < nlocal_nodes; idof++){
    for (int istride = 0; istride < Conductivity_Matrix_Strides(idof); istride++){
      Conductivity_Matrix(idof,istride) = 0;
      //debug print
      //std::cout << "{" <<istride + 1 << "," << DOF_Graph_Matrix(idof,istride) << "} ";
    }
    //debug print
    //std::cout << std::endl;
  }

  //assemble the global conductivity matrix
  if(num_dim==2)
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_2Delem_type(Element_Types(ielem), elem2D);
    nodes_per_element = elem2D->num_nodes();
    //construct local conductivity matrix for this element
    local_matrix(ielem, Local_Conductivity_Matrix);
    //assign entries of this local matrix to the sparse global matrix storage;
    for (int inode = 0; inode < nodes_per_element; inode++){
      //see if this node is local
      global_node_index = nodes_in_elem(ielem,inode);
      if(!map->isNodeGlobalElement(global_node_index)) continue;
      //set dof row start index
      current_row = map->getLocalElement(global_node_index);
      for(int jnode = 0; jnode < nodes_per_element; jnode++){
        
        current_column = Global_Conductivity_Matrix_Assembly_Map(ielem,inode,jnode);
        Conductivity_Matrix(current_row, current_column) += Local_Conductivity_Matrix(inode,jnode);
      }
    }
  }

  if(num_dim==3)
  for (int ielem = 0; ielem < rnum_elem; ielem++){
    element_select->choose_3Delem_type(Element_Types(ielem), elem);
    nodes_per_element = elem->num_nodes();
    //construct local conductivity matrix for this element
    local_matrix(ielem, Local_Conductivity_Matrix);
    //assign entries of this local matrix to the sparse global matrix storage;
    for (int inode = 0; inode < nodes_per_element; inode++){
      //see if this node is local
      global_node_index = nodes_in_elem(ielem,inode);
      if(!map->isNodeGlobalElement(global_node_index)) continue;
      //set dof row start index
      current_row = map->getLocalElement(global_node_index);
      for(int jnode = 0; jnode < nodes_per_element; jnode++){
        current_column = Global_Conductivity_Matrix_Assembly_Map(ielem,inode,jnode);
        Conductivity_Matrix(current_row, current_column) += Local_Conductivity_Matrix(inode,jnode);
      }
    }
  }

  //construct distributed conductivity matrix and force vector from local kokkos data
  
  //build column map for the global conductivity matrix
  Teuchos::RCP<const Tpetra::Map<LO,GO,node_type> > colmap;
  const Teuchos::RCP<const Tpetra::Map<LO,GO,node_type> > dommap = local_dof_map;
  
  //debug print
  /*
    std::cout << "DOF GRAPH MATRIX ENTRIES ON TASK " << myrank << std::endl;
  for (int idof = 0; idof < num_dim*nlocal_nodes; idof++){
    for (int istride = 0; istride < Conductivity_Matrix_Strides(idof); istride++){
      //debug print
      std::cout << "{" <<istride + 1 << "," << DOF_Graph_Matrix(idof,istride) << "} ";
    }
    //debug print
    std::cout << std::endl;
  } */

  //debug print of conductivity matrix
  /*
  std::cout << " ------------SPARSE STIFFNESS MATRIX ON TASK"<< myrank << std::endl;
  for (int idof = 0; idof < num_dim*nlocal_nodes; idof++){
      std::cout << "row: " << idof + 1 << " { ";
    for (int istride = 0; istride < Conductivity_Matrix_Strides(idof); istride++){
        std::cout << istride + 1 << " = " << Conductivity_Matrix(idof,istride) << " , " ;
    }
    std::cout << " }"<< std::endl;
  }
  */
  Tpetra::Details::makeColMap<LO,GO,node_type>(colmap,dommap,DOF_Graph_Matrix.get_kokkos_view(), nullptr);

  size_t nnz = DOF_Graph_Matrix.size();

  //debug print
  //std::cout << "DOF GRAPH SIZE ON RANK " << myrank << " IS " << nnz << std::endl;
  
  //local indices in the graph using the constructed column map
  CArrayKokkos<LO, array_layout, device_type, memory_traits> conductivity_local_indices(nnz, "conductivity_local_indices");
  
  //row offsets with compatible template arguments
    row_pointers row_offsets = DOF_Graph_Matrix.start_index_;
    row_pointers row_offsets_pass("row_offsets", nlocal_nodes + 1);
    for(int ipass = 0; ipass < nlocal_nodes + 1; ipass++){
      row_offsets_pass(ipass) = row_offsets(ipass);
    }

  size_t entrycount = 0;
  for(int irow = 0; irow < nlocal_nodes; irow++){
    for(int istride = 0; istride < Conductivity_Matrix_Strides(irow); istride++){
      conductivity_local_indices(entrycount) = colmap->getLocalElement(DOF_Graph_Matrix(irow,istride));
      entrycount++;
    }
  }
  
  if(!Matrix_alloc){
  Global_Conductivity_Matrix = Teuchos::rcp(new MAT(local_dof_map, colmap, row_offsets_pass, conductivity_local_indices.get_kokkos_view(), Conductivity_Matrix.get_kokkos_view()));
  Global_Conductivity_Matrix->fillComplete();
  Matrix_alloc = 1;
  }

  //filter small negative numbers (that should have been 0 from cancellation) from floating point error
  /*
  for (int idof = 0; idof < num_dim*nlocal_nodes; idof++){
    for (int istride = 0; istride < Conductivity_Matrix_Strides(idof); istride++){
      if(Conductivity_Matrix(idof,istride)<0.000000001*simparam->Elastic_Modulus*DENSITY_EPSILON||Conductivity_Matrix(idof,istride)>-0.000000001*simparam->Elastic_Modulus*DENSITY_EPSILON)
      Conductivity_Matrix(idof,istride) = 0;
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
  //Global_Conductivity_Matrix->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;

  //Print solution vector
  //*fos << "Global Nodal Forces :" << std::endl;
  //Global_Nodal_Forces->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
}


/* ----------------------------------------------------------------------
   Construct the global applied force vector
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::assemble_vector(){
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  //local variable for host view in the dual view
  host_vec_array Nodal_RHS = Global_Nodal_RHS->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
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
  int surface_flux_set_id = 0;
  int num_dim = simparam->num_dim;
  int nodes_per_elem = max_nodes_per_element;
  int num_gauss_points = simparam->num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  int current_element_index, local_surface_id, surf_dim1, surf_dim2, surface_sign;
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
  real_t heat_flux[3], wedge_product, Jacobian, current_density, weight_multiply, inner_product, surface_normal[3], surface_norm;
  real_t specific_internal_energy_rate;
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

   //RHS vector initialization
  for(int i=0; i < nlocal_nodes; i++)
    Nodal_RHS(i,0) = 0;

  /*Loop through boundary sets and check if they apply surface fluxes.
  These sets can have overlapping nodes since applied loading conditions
  are assumed to be additive*/
  for(int iboundary = 0; iboundary < num_boundary_sets; iboundary++){
    if(Boundary_Condition_Type_List(iboundary)!=SURFACE_LOADING_CONDITION) continue;
    //std::cout << "I REACHED THE LOADING BOUNDARY CONDITION" <<std::endl;
    num_bdy_patches_in_set = NBoundary_Condition_Patches(iboundary);
    
    heat_flux[0] = Boundary_Surface_Heat_Flux(surface_flux_set_id,0);
    heat_flux[1] = Boundary_Surface_Heat_Flux(surface_flux_set_id,1);
    heat_flux[2] = Boundary_Surface_Heat_Flux(surface_flux_set_id,2);
    surface_flux_set_id++;

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
        if(local_surface_id%2==0){
          quad_coordinate(2) = -1;
          surface_sign = 1;
        }
        else{
          quad_coordinate(2) = 1;
          surface_sign = -1;
        }
      }
      else if(local_surface_id<4){
      surf_dim1 = 0;
      surf_dim2 = 2;
      quad_coordinate(0) = legendre_nodes_1D(x_quad);
      quad_coordinate(2) = legendre_nodes_1D(y_quad);
      //set to -1 or 1 for an isoparametric space
        if(local_surface_id%2==0){
          quad_coordinate(1) = -1;
          surface_sign = 1;
        }
        else{
          quad_coordinate(1) = 1;
          surface_sign = -1;
        }
      }
      else if(local_surface_id<6){
      surf_dim1 = 1;
      surf_dim2 = 2;
      quad_coordinate(1) = legendre_nodes_1D(x_quad);
      quad_coordinate(2) = legendre_nodes_1D(y_quad);
      //set to -1 or 1 for an isoparametric space
        if(local_surface_id%2==0){
          quad_coordinate(0) = -1;
          surface_sign = 1;
        }
        else{
          quad_coordinate(0) = 1;
          surface_sign = -1;
        }
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
      //compute product scaling the reference differential surface element to the deformed differential surface element
      wedge_product = sqrt(pow(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2),2)+
               pow(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2),2)+
               pow(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1),2));
      
      //compute surface normal through cross product of the two tangent vectors
      surface_normal[0] = surface_sign*JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2);
      surface_normal[1] = surface_sign*JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2);
      if(num_dim==3)
        surface_normal[2] = surface_sign*JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1);
      else
        surface_normal[2] = 0;
      
      //normalize
      surface_normal[0] = surface_normal[0]/wedge_product;
      surface_normal[1] = surface_normal[1]/wedge_product;
      surface_normal[2] = surface_normal[2]/wedge_product;
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
        inner_product = surface_normal[0]*heat_flux[0] + surface_normal[1]*heat_flux[1] + surface_normal[2]*heat_flux[2];

        if(inner_product!=0)
          Nodal_RHS(node_id,0) += wedge_product*quad_coordinate_weight(0)*quad_coordinate_weight(1)*inner_product*basis_values(local_nodes[node_count]);
        
      }
      }
    }
  }
  }

    //apply line distribution of forces

    //apply point forces

    //apply contribution from non-zero displacement boundary conditions

    //apply body forces
    if(body_term_flag){
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
      Body_Term(ielem, current_density, specific_internal_energy_rate);
    
      //evaluate contribution to force vector component
      for(int ibasis=0; ibasis < nodes_per_elem; ibasis++){
        if(!map->isNodeGlobalElement(nodes_in_elem(ielem, ibasis))) continue;
        local_node_id = map->getLocalElement(nodes_in_elem(ielem, ibasis));
        Nodal_RHS(local_node_id,0) += Jacobian*weight_multiply*specific_internal_energy_rate*basis_values(ibasis);
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

/* ----------------------------------------------------------------------
   Retrieve body force at a point
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::Body_Term(size_t ielem, real_t density, real_t &specific_internal_energy_rate){
  real_t unit_scaling = simparam->unit_scaling;
  int num_dim = simparam->num_dim;
  
  //init 
  specific_internal_energy_rate = 0;
  if(thermal_flag){
    specific_internal_energy_rate += simparam->specific_internal_energy_rate*density;
  }
  
  /*
  if(body_term_flag){

  }

  if(electric_flag){

  }

  */
}

/* ----------------------------------------------------------------------
   Gradient of body force at a point
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::Gradient_Body_Term(size_t ielem, real_t density, real_t &gradient_specific_internal_energy_rate){
  real_t unit_scaling = simparam->unit_scaling;
  int num_dim = simparam->num_dim;
  
  //init 
  gradient_specific_internal_energy_rate = 0;\
  if(thermal_flag){
    specific_internal_energy_rate += simparam->specific_internal_energy_rate;
  }
  
  /*
  if(thermal_flag){

  }

  if(electric_flag){

  }

  */
}

/* ----------------------------------------------------------------------
   Retrieve material properties associated with a finite element
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::Element_Material_Properties(size_t ielem, real_t &Element_Conductivity, real_t density){
  real_t unit_scaling = simparam->unit_scaling;
  real_t penalty_product = 1;
  if(density < 0) density = 0;
  for(int i = 0; i < simparam->penalty_power; i++)
    penalty_product *= density;
  //relationship between density and conductivity
  Element_Conductivity = (DENSITY_EPSILON + (1 - DENSITY_EPSILON)*penalty_product)*simparam->Thermal_Conductivity/unit_scaling;
}

/* ----------------------------------------------------------------------
   Retrieve derivative of material properties with respect to local density
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::Gradient_Element_Material_Properties(size_t ielem, real_t &Element_Conductivity_Derivative, real_t density){
  real_t unit_scaling = simparam->unit_scaling;
  real_t penalty_product = 1;
  Element_Conductivity_Derivative = 0;
  if(density < 0) density = 0;
  for(int i = 0; i < simparam->penalty_power - 1; i++)
    penalty_product *= density;
  //relationship between density and conductivity
  Element_Conductivity_Derivative = simparam->penalty_power*(1 - DENSITY_EPSILON)*penalty_product*simparam->Thermal_Conductivity/unit_scaling;
}

/* --------------------------------------------------------------------------------
   Retrieve second derivative of material properties with respect to local density
----------------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::Concavity_Element_Material_Properties(size_t ielem, real_t &Element_Conductivity_Derivative, real_t density){
  real_t unit_scaling = simparam->unit_scaling;
  real_t penalty_product = 1;
  ElementConductivity_Derivative = 0;
  if(density < 0) density = 0;
  if(simparam->penalty_power>=2){
    for(int i = 0; i < simparam->penalty_power - 2; i++)
      penalty_product *= density;
    //relationship between density and conductivity
    Element_Conductivity_Derivative = simparam->penalty_power*(simparam->penalty_power-1)*(1 - DENSITY_EPSILON)*penalty_product*simparam->Thermal_Conductivity/unit_scaling;
  }
}

/* ----------------------------------------------------------------------
   Construct the local conductivity matrix
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::local_matrix(int ielem, CArrayKokkos<real_t, array_layout, device_type, memory_traits> &Local_Matrix){
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
  real_t matrix_subterm1, matrix_subterm2, matrix_subterm3, invJacobian, Jacobian, weight_multiply;
  real_t Element_Modulus, Poisson_Ratio;
  //CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
  //CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
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
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix_contribution(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix(Brows,num_dim*elem->num_basis());
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix_contribution(Brows,num_dim*elem->num_basis());
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

  //initialize local conductivity matrix storage
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
    else
    quad_coordinate_weight(2) = 1;
    weight_multiply = quad_coordinate_weight(0)*quad_coordinate_weight(1)*quad_coordinate_weight(2);

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
    invJacobian = 1/Jacobian;
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

    //compute the contributions of this quadrature point to all the local conductivity matrix elements
    for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++)
      for(int jfill=ifill; jfill < num_dim*nodes_per_elem; jfill++){
        matrix_term = 0;
        for(int span = 0; span < Brows; span++){
          matrix_term += B_matrix_contribution(span,ifill)*CB_matrix_contribution(span,jfill);
        }
        Local_Matrix(ifill,jfill) += Elastic_Constant*weight_multiply*matrix_term*invJacobian;
        if(ifill!=jfill)
          Local_Matrix(jfill,ifill) = Local_Matrix(ifill,jfill);
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

    //debug print of local conductivity matrix
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
   Loop through applied boundary conditions and tag node ids to remove 
   necessary rows and columns from the assembled linear system
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::Temperature_Boundary_Conditions(){
  int num_bdy_patches_in_set, patch_id;
  int warning_flag = 0;
  int local_flag;
  int current_node_index, current_node_id;
  int num_boundary_sets = num_boundary_conditions;
  int surface_temp_set_id = 0;
  int num_dim = simparam->num_dim;
  int bc_option;
  int DOF_BC_type;
  real_t temperature;
  CArrayKokkos<int, array_layout, device_type, memory_traits> Temperatures_Conditions(num_dim);
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> first_condition_per_node(nall_nodes);
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Surface_Nodes;
  Number_DOF_BCS = 0;

  //host view of local nodal displacements
  host_vec_array node_temperatures_host = node_temperatures_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);

  //initialize to -1 (DO NOT make -1 an index for bdy sets)
  for(int inode = 0 ; inode < nlocal_nodes; inode++)
    first_condition_per_node(inode) = -1;
  
  //scan for surface method of setting fixed nodal displacements
  for(int iboundary = 0; iboundary < num_boundary_sets; iboundary++){
    
    if(Boundary_Condition_Type_List(iboundary)==TEMPERATURE_CONDITION){
      bc_option=0;
      DOF_BC_type = TEMPERATURE_CONDITION;
    }
    else{
      continue;
    }
      
      //debug print of surface conditions
      //std::cout << " Surface BC types " << Boundary_Condition_Type_List(iboundary) <<std::endl;

      num_bdy_patches_in_set = NBoundary_Condition_Patches(iboundary);
      if(bc_option==0) {
        bc_dim_set[0]=1;
        displacement(0) = Boundary_Surface_Temperatures(surface_temp_set_id,0);
      }
      surface_temp_set_id++;
      
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

          //warning for reapplied a displacement boundary condition (For now there is an output flag being set that triggers output later)
          if(Node_DOF_Boundary_Condition_Type(current_node_id)==DOF_BC_type){
            //if overlap is just due to the loop over patches, a warning is not needed
            if(first_condition_per_node(current_node_id)!=iboundary) warning_flag = 1;
          }
          else{
            first_condition_per_node(current_node_id) = iboundary;
            Node_DOF_Boundary_Condition_Type(current_node_id) = DOF_BC_type;
            Node_Temperature_Boundary_Conditions(current_node_id) = temperature;
            //counts local DOF being constrained
            if(local_flag){
            Number_DOF_BCS++;
            node_temperatures_host(current_node_id,0) = temperature;
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
  //if(warning_flag)
  //std::cout << std::endl << "One or more displacement boundary conditions overlap on a subset of nodes; please revise input" << std::endl << std::endl;

}

/* ----------------------------------------------------------------------
   Compute the gradient of heat capacity potential with respect to nodal densities
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::compute_adjoint_gradients(const_host_vec_array design_variables, host_vec_array design_gradients){
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_vec_array all_node_temperatures = all_node_temperatures_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
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
  real_t Element_Modulus_Gradient, Poisson_Ratio, gradient_force_density[3];
  real_t Elastic_Constant, Shear_Term, Pressure_Term;
  real_t inner_product, matrix_term, Jacobian, invJacobian, weight_multiply;
  //CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
  //CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
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
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix_contribution(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix(Brows,num_dim*elem->num_basis());
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix_contribution(Brows,num_dim*elem->num_basis());
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
      current_nodal_displacements(node_loop*num_dim) = all_node_temperatures(local_dof_idx,0);
      current_nodal_displacements(node_loop*num_dim+1) = all_node_temperatures(local_dof_idy,0);
      current_nodal_displacements(node_loop*num_dim+2) = all_node_temperatures(local_dof_idz,0);
      
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
    invJacobian = 1/Jacobian;

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
    
    //compute the contributions of this quadrature point to all the local conductivity matrix elements
    for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++){
      for(int jfill=ifill; jfill < num_dim*nodes_per_elem; jfill++){
        matrix_term = 0;
        for(int span = 0; span < Brows; span++){
          matrix_term += B_matrix_contribution(span,ifill)*CB_matrix_contribution(span,jfill);
        }
        Local_Matrix_Contribution(ifill,jfill) = matrix_term;
        if(ifill!=jfill)
          Local_Matrix_Contribution(jfill,ifill) = Local_Matrix_Contribution(ifill,jfill);
      }
    }

    //compute inner product for this quadrature point contribution
    inner_product = 0;
    for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++){
      for(int jfill=ifill; jfill < num_dim*nodes_per_elem; jfill++){
        if(ifill==jfill)
          inner_product += Local_Matrix_Contribution(ifill, jfill)*current_nodal_displacements(ifill)*current_nodal_displacements(jfill);
        else
          inner_product += 2*Local_Matrix_Contribution(ifill, jfill)*current_nodal_displacements(ifill)*current_nodal_displacements(jfill);
        //debug
        //if(Local_Matrix_Contribution(ifill, jfill)<0) Local_Matrix_Contribution(ifill, jfill) = - Local_Matrix_Contribution(ifill, jfill);
        //inner_product += Local_Matrix_Contribution(ifill, jfill);
      }
    }

    //evaluate local conductivity matrix gradient with respect to igradient
    for(int igradient=0; igradient < nodes_per_elem; igradient++){
      if(!map->isNodeGlobalElement(nodes_in_elem(ielem, igradient))) continue;
      local_node_id = map->getLocalElement(nodes_in_elem(ielem, igradient));
      
      //debug print
      //std::cout << "contribution for " << igradient + 1 << " is " << inner_product << std::endl;
      design_gradients(local_node_id,0) -= inner_product*Elastic_Constant*basis_values(igradient)*weight_multiply*0.5*invJacobian;
    }

      //evaluate gradient of body force (such as gravity which depends on density) with respect to igradient
    if(body_term_flag){
      //look up element material properties at this point as a function of density
      Gradient_Body_Term(ielem, current_density, gradient_force_density);
      for(int igradient=0; igradient < nodes_per_elem; igradient++){
      if(!map->isNodeGlobalElement(nodes_in_elem(ielem, igradient))) continue;
      local_node_id = map->getLocalElement(nodes_in_elem(ielem, igradient));
      
      //compute inner product for this quadrature point contribution
      inner_product = 0;
      for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++){
        inner_product += gradient_force_density[ifill%num_dim]*current_nodal_displacements(ifill)*basis_values(ifill/num_dim);
      }
      
      //debug print
      //std::cout << "contribution for " << igradient + 1 << " is " << inner_product << std::endl;
      design_gradients(local_node_id,0) += inner_product*basis_values(igradient)*weight_multiply*Jacobian;
      }
    }
    }
  }
  
    
  //debug print

}

/* ----------------------------------------------------------------------
   Compute the gradient of heat capacity potential with respect to nodal densities
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::compute_adjoint_hessian_vec(const_host_vec_array design_densities, host_vec_array hessvec, Teuchos::RCP<const MV> direction_vec_distributed){
  //local variable for host view in the dual view
  real_t current_cpu_time = Solver_Pointer_->CPU_Time();
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_vec_array all_node_temperatures = all_node_temperatures_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array nodes_in_elem = nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  
  const_host_vec_array Element_Densities;
  //local variable for host view of densities from the dual view
  bool nodal_density_flag = simparam->nodal_density_flag;
  const_host_vec_array all_node_densities;
  if(nodal_density_flag)
  all_node_densities = all_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  else
  Element_Densities = Global_Element_Densities->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
  host_vec_array unbalanced_B_view = unbalanced_B->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
  const_host_vec_array direction_vec = direction_vec_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
  Teuchos::RCP<Xpetra::MultiVector<real_t,LO,GO,node_type>> xlambda = xX;
  Teuchos::RCP<MV> lambda = X;
  const_host_vec_array lambda_view = lambda->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);

  int num_dim = simparam->num_dim;
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam->num_gauss_points;
  int strain_max_flag = simparam->strain_max_flag;
  int z_quad,y_quad,x_quad, direct_product_count;
  int solve_flag, zero_strain_flag;
  LO local_node_id, jlocal_node_id, temp_id, local_dof_id, local_reduced_dof_id, local_dof_idx, local_dof_idy, local_dof_idz;
  GO current_global_index, global_dof_id;

  direct_product_count = std::pow(num_gauss_points,num_dim);
  real_t Element_Modulus_Gradient, Element_Modulus_Concavity, Poisson_Ratio, gradient_force_density[3];
  real_t Elastic_Constant, Gradient_Elastic_Constant, Concavity_Elastic_Constant, Shear_Term, Pressure_Term;
  real_t inner_product, matrix_term, Jacobian, invJacobian, weight_multiply;
  real_t direction_vec_reduce, local_direction_vec_reduce;
  //CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
  //CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
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
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_adjoint_displacements(elem->num_basis()*num_dim);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());

  size_t Brows;
  if(num_dim==2) Brows = 3;
  if(num_dim==3) Brows = 6;
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix_contribution(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix(Brows,num_dim*elem->num_basis());
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix_contribution(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix(Brows,num_dim*elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> C_matrix(Brows,Brows);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> Local_Matrix_Contribution(num_dim*nodes_per_elem,num_dim*nodes_per_elem);

  //initialize weights
  elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
  elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);
  
  real_t current_density = 1;
  
  //direction_vec_distributed->describe(*fos,Teuchos::VERB_EXTREME);

  //initialize gradient value to zero
  for(size_t inode = 0; inode < nlocal_nodes; inode++)
    hessvec(inode,0) = 0;
  
  //initialize RHS vector
  for(int i=0; i < local_reduced_dof_map->getNodeNumElements(); i++)
    unbalanced_B_view(i,0) = 0;
  
  //sum components of direction vector
  direction_vec_reduce = local_direction_vec_reduce = 0;
  for(int i = 0; i < nlocal_nodes; i++)
    local_direction_vec_reduce += direction_vec(i,0);
  
  MPI_Allreduce(&local_direction_vec_reduce,&direction_vec_reduce,1,MPI_DOUBLE,MPI_SUM,world);

  //comms to get ghost components of direction vector needed for matrix inner products
  Tpetra::Import<LO, GO> node_importer(map, all_node_map);
  
  Teuchos::RCP<MV> all_direction_vec_distributed = Teuchos::rcp(new MV(all_node_map, 1));
  //comms to get ghosts
  all_direction_vec_distributed->doImport(*direction_vec_distributed, node_importer, Tpetra::INSERT);
  
  const_host_vec_array all_direction_vec = all_direction_vec_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);

  //loop through each element to contribute to the RHS of the hessvec adjoint equation
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
      current_nodal_displacements(node_loop*num_dim) = all_node_temperatures(local_dof_idx,0);
      current_nodal_displacements(node_loop*num_dim+1) = all_node_temperatures(local_dof_idy,0);
      current_nodal_displacements(node_loop*num_dim+2) = all_node_temperatures(local_dof_idz,0);
      
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
      invJacobian = 1/Jacobian;

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
    
    //compute the contributions of this quadrature point to all the local conductivity matrix elements
    for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++){
      for(int jfill=ifill; jfill < num_dim*nodes_per_elem; jfill++){
        matrix_term = 0;
        for(int span = 0; span < Brows; span++){
          matrix_term += B_matrix_contribution(span,ifill)*CB_matrix_contribution(span,jfill);
        }
        Local_Matrix_Contribution(ifill,jfill) = matrix_term;
        if(ifill!=jfill)
          Local_Matrix_Contribution(jfill,ifill) = Local_Matrix_Contribution(ifill,jfill);
      }
    }

    //evaluate local conductivity matrix gradient with respect to igradient
    for(int igradient=0; igradient < nodes_per_elem; igradient++){
      //if(!map->isNodeGlobalElement(nodes_in_elem(ielem, igradient))) continue;
      local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, igradient));
      
      //compute rhs product for this quadrature point contribution  
      for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++){
        local_dof_id = all_dof_map->getLocalElement(nodes_in_elem(ielem, ifill/num_dim)*num_dim);
        local_dof_id += ifill%num_dim;
        global_dof_id = all_dof_map->getGlobalElement(local_dof_id);
        if(Node_DOF_Boundary_Condition_Type(local_dof_id)!=DISPLACEMENT_CONDITION&&local_reduced_dof_original_map->isNodeGlobalElement(global_dof_id)){
          local_reduced_dof_id = local_reduced_dof_original_map->getLocalElement(global_dof_id);
          inner_product = 0;
          for(int jfill=0; jfill < num_dim*nodes_per_elem; jfill++){
            inner_product += Local_Matrix_Contribution(ifill, jfill)*current_nodal_displacements(jfill);
            //debug
            //if(Local_Matrix_Contribution(ifill, jfill)<0) Local_Matrix_Contribution(ifill, jfill) = - Local_Matrix_Contribution(ifill, jfill);
            //inner_product += Local_Matrix_Contribution(ifill, jfill);
          }
          unbalanced_B_view(local_reduced_dof_id,0) += inner_product*Elastic_Constant*basis_values(igradient)*weight_multiply*all_direction_vec(local_node_id,0)*invJacobian;
        }
      }
      } //density gradient loop
    }//quadrature loop
  }//element index loop
  
  //*fos << "Elastic Modulus Gradient" << Element_Modulus_Gradient <<std::endl;
  //*fos << "DISPLACEMENT" << std::endl;
  //all_node_temperatures_distributed->describe(*fos,Teuchos::VERB_EXTREME);

  //*fos << "RHS vector" << std::endl;
  //unbalanced_B->describe(*fos,Teuchos::VERB_EXTREME);
  //balance RHS vector due to missing BC dofs
  //import object to rebalance force vector
  Tpetra::Import<LO, GO> Bvec_importer(local_reduced_dof_map, local_balanced_reduced_dof_map);
  
  //comms to rebalance force vector
  balanced_B->doImport(*unbalanced_B, Bvec_importer, Tpetra::INSERT);
  
  //solve for adjoint vector
  int num_iter = 2000;
  double solve_tol = 1e-05;
  int cacheSize = 0;
  std::string solveType         = "belos";
  std::string belosType         = "cg";
  // =========================================================================
  // Preconditioner construction
  // =========================================================================
  //bool useML   = Linear_Solve_Params->isParameter("use external multigrid package") && (Linear_Solve_Params->get<std::string>("use external multigrid package") == "ml");
  //out<<"*********** MueLu ParameterList ***********"<<std::endl;
  //out<<*Linear_Solve_Params;
  //out<<"*******************************************"<<std::endl;
  
  //H->Write(-1, -1);
  //H->describe(*fos,Teuchos::VERB_EXTREME);
  
  // =========================================================================
  // System solution (Ax = b)
  // =========================================================================
  //since matrix graph and A are the same from the last update solve, the Hierarchy H need not be rebuilt
  //xwrap_balanced_A->describe(*fos,Teuchos::VERB_EXTREME);
  if(simparam->equilibrate_matrix_flag){
    Solver_Pointer_->preScaleRightHandSides(*balanced_B,"diag");
    Solver_Pointer_->preScaleInitialGuesses(*lambda,"diag");
  }
  real_t current_cpu_time2 = Solver_Pointer_->CPU_Time();
  comm->barrier();
  SystemSolve(xwrap_balanced_A,xlambda,xbalanced_B,H,Prec,*fos,solveType,belosType,false,false,false,cacheSize,0,true,true,num_iter,solve_tol);
  comm->barrier();
  hessvec_linear_time += Solver_Pointer_->CPU_Time() - current_cpu_time2;

  if(simparam->equilibrate_matrix_flag){
    Solver_Pointer_->postScaleSolutionVectors(*lambda,"diag");
  }
  //scale by reciprocal ofdirection vector sum
  lambda->scale(1/direction_vec_reduce);
  //*fos << "LAMBDA" << std::endl;
  //lambda->describe(*fos,Teuchos::VERB_EXTREME);
  //communicate adjoint vector to original all dof map for simplicity now (optimize out later)
  Teuchos::RCP<MV> adjoint_distributed = Teuchos::rcp(new MV(local_dof_map, 1));
  Teuchos::RCP<MV> all_adjoint_distributed = Teuchos::rcp(new MV(all_dof_map, 1));
  //communicate solution on reduced map to the all node map vector for post processing of strain etc.
  //intermediate storage on the unbalanced reduced system
  Teuchos::RCP<MV> reduced_adjoint_distributed = Teuchos::rcp(new MV(local_reduced_dof_map, 1));
  //create import object using local node indices map and all indices map
  Tpetra::Import<LO, GO> reduced_adjoint_importer(local_balanced_reduced_dof_map, local_reduced_dof_map);

  //comms to get displacements on reduced unbalanced displacement vector
  reduced_adjoint_distributed->doImport(*lambda, reduced_adjoint_importer, Tpetra::INSERT);

  //populate node displacement multivector on the local dof map
  const_host_vec_array reduced_adjoint_host = reduced_adjoint_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  host_vec_array adjoint_host = adjoint_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);

  for(int init = 0; init < local_dof_map->getNodeNumElements(); init++)
    adjoint_host(init,0) = 0;

  for(LO i=0; i < local_reduced_dof_original_map->getNodeNumElements(); i++){
   local_reduced_dof_id = local_dof_map->getLocalElement(Free_Indices(i));
    adjoint_host(local_reduced_dof_id,0) = reduced_adjoint_host(i,0);
  }
  
  //import for displacement of ghosts
  Tpetra::Import<LO, GO> ghost_displacement_importer(local_dof_map, all_dof_map);

  //comms to get displacements on all node map
  all_adjoint_distributed->doImport(*adjoint_distributed, ghost_displacement_importer, Tpetra::INSERT);
  host_vec_array all_adjoint = all_adjoint_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  //*fos << "ALL ADJOINT" << std::endl;
  //all_adjoint_distributed->describe(*fos,Teuchos::VERB_EXTREME);
//now that adjoint is computed, calculate the hessian vector product
//loop through each element and assign the contribution to Hessian vector product for each of its local nodes

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
      current_nodal_displacements(node_loop*num_dim) = all_node_temperatures(local_dof_idx,0);
      current_nodal_displacements(node_loop*num_dim+1) = all_node_temperatures(local_dof_idy,0);
      current_nodal_displacements(node_loop*num_dim+2) = all_node_temperatures(local_dof_idz,0);
      current_adjoint_displacements(node_loop*num_dim) = all_adjoint(local_dof_idx,0);
      current_adjoint_displacements(node_loop*num_dim+1) = all_adjoint(local_dof_idy,0);
      current_adjoint_displacements(node_loop*num_dim+2) = all_adjoint(local_dof_idz,0);
      
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
    invJacobian = 1/Jacobian;

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

    //look up element material properties at this point as a function of density
    Concavity_Element_Material_Properties(ielem, Element_Modulus_Concavity, Poisson_Ratio, current_density);
    Gradient_Element_Material_Properties(ielem, Element_Modulus_Gradient, Poisson_Ratio, current_density);
    
    Gradient_Elastic_Constant = Element_Modulus_Gradient/((1 + Poisson_Ratio)*(1 - 2*Poisson_Ratio));
    Concavity_Elastic_Constant = Element_Modulus_Concavity/((1 + Poisson_Ratio)*(1 - 2*Poisson_Ratio));
    Shear_Term = 0.5 - Poisson_Ratio;
    Pressure_Term = 1 - Poisson_Ratio;
    //*fos << "Elastic Modulus Concavity" << Element_Modulus_Concavity << " " << Element_Modulus_Gradient << std::endl;
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
    
    //compute the contributions of this quadrature point to all the local conductivity matrix elements
      for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++){
        for(int jfill=ifill; jfill < num_dim*nodes_per_elem; jfill++){
          matrix_term = 0;
          for(int span = 0; span < Brows; span++){
            matrix_term += B_matrix_contribution(span,ifill)*CB_matrix_contribution(span,jfill);
          }
          Local_Matrix_Contribution(ifill,jfill) = matrix_term;
          if(jfill!=ifill)
            Local_Matrix_Contribution(jfill,ifill) = Local_Matrix_Contribution(ifill,jfill);
        }
      }
      
    //compute inner product for this quadrature point contribution
    inner_product = 0;
    for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++){
      for(int jfill=ifill; jfill < num_dim*nodes_per_elem; jfill++){
        if(ifill==jfill)
          inner_product += Local_Matrix_Contribution(ifill, jfill)*current_nodal_displacements(ifill)*current_nodal_displacements(jfill);
        else
          inner_product += 2*Local_Matrix_Contribution(ifill, jfill)*current_nodal_displacements(ifill)*current_nodal_displacements(jfill);
      }
    }

    //evaluate local conductivity matrix concavity with respect to igradient and jgradient
    for(int igradient=0; igradient < nodes_per_elem; igradient++){
      local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, igradient));
      for(int jgradient=igradient; jgradient < nodes_per_elem; jgradient++){
        jlocal_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, jgradient));
        //debug print
        //std::cout << "contribution for " << igradient + 1 << " is " << inner_product << std::endl;
        if(map->isNodeGlobalElement(nodes_in_elem(ielem, igradient))){
        temp_id = map->getLocalElement(nodes_in_elem(ielem, igradient));
          hessvec(temp_id,0) -= inner_product*Concavity_Elastic_Constant*basis_values(igradient)*all_direction_vec(jlocal_node_id,0)*
                                  basis_values(jgradient)*weight_multiply*0.5*invJacobian;
        }
        if(igradient!=jgradient&&map->isNodeGlobalElement(nodes_in_elem(ielem, jgradient))){
          temp_id = map->getLocalElement(nodes_in_elem(ielem, jgradient));
          hessvec(jlocal_node_id,0) -= inner_product*Concavity_Elastic_Constant*basis_values(igradient)*all_direction_vec(local_node_id,0)*
                                      basis_values(jgradient)*weight_multiply*0.5*invJacobian;

        }
      }
    }
    
    //compute inner product for this quadrature point contribution
    inner_product = 0;
    for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++){
      for(int jfill=0; jfill < num_dim*nodes_per_elem; jfill++){
        inner_product += Local_Matrix_Contribution(ifill, jfill)*current_adjoint_displacements(ifill)*current_nodal_displacements(jfill);
        //debug
        //if(Local_Matrix_Contribution(ifill, jfill)<0) Local_Matrix_Contribution(ifill, jfill) = - Local_Matrix_Contribution(ifill, jfill);
        //inner_product += Local_Matrix_Contribution(ifill, jfill);
      }
    }

    //evaluate local conductivity matrix gradient with respect to igradient (augmented term with adjoint vector)
    for(int igradient=0; igradient < nodes_per_elem; igradient++){
      if(!map->isNodeGlobalElement(nodes_in_elem(ielem, igradient))) continue;
      local_node_id = map->getLocalElement(nodes_in_elem(ielem, igradient));
      
      //debug print
      //std::cout << "contribution for " << igradient + 1 << " is " << inner_product << std::endl;
      hessvec(local_node_id,0) += inner_product*direction_vec_reduce*Gradient_Elastic_Constant*basis_values(igradient)*weight_multiply*invJacobian;
      }

      //evaluate gradient of body force (such as gravity which depends on density) with respect to igradient
      if(body_term_flag){
        for(int igradient=0; igradient < nodes_per_elem; igradient++){
        if(!map->isNodeGlobalElement(nodes_in_elem(ielem, igradient))) continue;
        local_node_id = map->getLocalElement(nodes_in_elem(ielem, igradient));
        //look up element material properties at this point as a function of density
        Gradient_Body_Term(ielem, current_density, gradient_force_density);
      
        //compute inner product for this quadrature point contribution
        inner_product = 0;
        for(int ifill=0; ifill < num_dim*nodes_per_elem; ifill++){
          inner_product -= gradient_force_density[ifill%num_dim]*
                           current_adjoint_displacements(ifill)*basis_values(ifill/num_dim);
        }
      
        //debug print
        //std::cout << "contribution for " << igradient + 1 << " is " << inner_product << std::endl;
        hessvec(local_node_id,0) += inner_product*direction_vec_reduce*basis_values(igradient)*weight_multiply*Jacobian;
        }
      }
    }
  }//end element loop for hessian vector product
  hessvec_time += Solver_Pointer_->CPU_Time() - current_cpu_time;
}

/* -------------------------------------------------------------------------------------------
   Compute the maximum nodal heat fluxes resulting from minimizing the L2 error
   between flux (subspace solution) and a nodal interpolation (nodal fluxes defined at each node)
   for each element. Mainly used for output and is approximate.
---------------------------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::compute_nodal_heat_fluxes(){
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_vec_array all_node_temperatures = all_node_temperatures_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
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
  real_t matrix_subterm1, matrix_subterm2, matrix_subterm3, Jacobian, weight_multiply;
  //CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_nodes_1D(num_gauss_points);
  //CArrayKokkos<real_t, array_layout, device_type, memory_traits> legendre_weights_1D(num_gauss_points);
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
      current_nodal_displacements(node_loop*num_dim) = all_node_temperatures(local_dof_idx,0);
      current_nodal_displacements(node_loop*num_dim+1) = all_node_temperatures(local_dof_idy,0);
      current_nodal_displacements(node_loop*num_dim+2) = all_node_temperatures(local_dof_idz,0);

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
        projection_vector(irow,icol) += weight_multiply*quad_strain(irow)*basis_values(icol);
        
      }

    //compute contribution to projection matrix (only upper part is set)
    for(int irow=0; irow < nodes_per_elem; irow++)
      for(int icol=0; icol < nodes_per_elem; icol++){
        //if(irow<=icol)
        projection_matrix(irow,icol) += weight_multiply*basis_values(irow)*basis_values(icol)*Jacobian;
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
        if(fabs(projection_vector(istrain,icol))>FLUX_EPSILON) zero_strain_flag = 0;
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
   Solve the FEA linear system
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::linear_solver_parameters(){
  if(simparam->direct_solver_flag){
    Linear_Solve_Params = Teuchos::rcp(new Teuchos::ParameterList("Amesos2"));
    auto superlu_params = Teuchos::sublist(Teuchos::rcpFromRef(*Linear_Solve_Params), "SuperLU_DIST");
    superlu_params->set("Equil", true);
    //superlu_params.set("Trans","TRANS","Whether to solve with A^T");
    //superlu_params.set("ColPerm","NATURAL","Use 'natural' ordering of columns");
  
  }
  else{
    Linear_Solve_Params = Teuchos::rcp(new Teuchos::ParameterList("MueLu"));
    std::string xmlFileName = "MueLu_Thermal_3D_Params.xml";
    //std::string xmlFileName = "simple_test.xml";
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&(*Linear_Solve_Params)), *comm);
  }
}

/* ----------------------------------------------------------------------
   Solve the FEA linear system
------------------------------------------------------------------------- */

int FEA_Module_Heat_Conduction::solve(){
  //local variable for host view in the dual view
  const_host_vec_array node_coords = node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  //local variable for host view in the dual view
  const_host_vec_array Nodal_RHS = Global_Nodal_RHS->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  int num_dim = simparam->num_dim;
  int nodes_per_elem = max_nodes_per_element;
  int local_node_index, current_row, current_column;
  size_t reduced_index;
  int max_stride = 0;
  size_t access_index, row_access_index, row_counter;
  global_size_t reduced_row_count;
  GO global_index, global_dof_index;
  LO reduced_local_dof_index;

  //*fos << Amesos2::version() << std::endl << std::endl;

  bool printTiming   = true;
  bool verbose       = false;
  std::string filename("arc130.mtx");
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("print-timing","no-print-timing",&printTiming,"Print solver timing statistics");

  const size_t numVectors = 1;
  
  //construct global data for example; normally this would be from file input and distributed
  //according to the row map at that point
  global_size_t nrows = num_nodes;
  
 
  //number of boundary conditions on this mpi rank
  global_size_t local_nboundaries = Number_DOF_BCS;
  size_t local_nrows_reduced = nlocal_nodes - local_nboundaries;

  //obtain total number of boundary conditions on all ranks
  global_size_t global_nboundaries = 0;
  MPI_Allreduce(&local_nboundaries,&global_nboundaries,1,MPI_INT,MPI_SUM,world);
  global_size_t nrows_reduced = nrows - global_nboundaries;
  //std::cout << "GLOBAL BC DATA " << nrows_reduced << " " << nrows <<"\n ";  

  //Rebalance distribution of the global conductivity matrix rows here later since
  //rows and columns are being removed.

  //global_size_t *entries_per_row = new global_size_t[local_nrows];
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Reduced_Conductivity_Matrix_Strides(local_nrows_reduced,"Reduced_Conductivity_Matrix_Strides");
  row_pointers reduced_row_offsets_pass = row_pointers("reduced_row_offsets_pass", local_nrows_reduced+1);

  //init row_offsets
  for(int i=0; i < local_nrows_reduced+1; i++){
    reduced_row_offsets_pass(i) = 0;
  }
  
  //stores global indices belonging to this MPI rank from the non-reduced map corresponding to the reduced system
  Free_Indices = CArrayKokkos<GO, array_layout, device_type, memory_traits>(local_nrows_reduced,"Free_Indices");
  reduced_index = 0;
  for(LO i=0; i < nlocal_nodes; i++)
    if((Node_DOF_Boundary_Condition_Type(i)!=TEMPERATURE_CONDITION)){
        Free_Indices(reduced_index) = local_dof_map->getGlobalElement(i);
        reduced_index++;
      }
    
  
  //compute reduced matrix strides
  for(LO i=0; i < local_nrows_reduced; i++){
    access_index = local_dof_map->getLocalElement(Free_Indices(i));
    reduced_row_count = 0;
    for(LO j=0; j < Conductivity_Matrix_Strides(access_index); j++){
      global_dof_index = DOF_Graph_Matrix(access_index,j);
      row_access_index = all_dof_map->getLocalElement(global_dof_index);
      if((Node_DOF_Boundary_Condition_Type(row_access_index)!=TEMPERATURE_CONDITION)){
        reduced_row_count++;
      }
    }
    Reduced_Conductivity_Matrix_Strides(i) = reduced_row_count;
  }
  
  RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits> Reduced_DOF_Graph_Matrix(Reduced_Conductivity_Matrix_Strides); //stores global dof indices
  RaggedRightArrayKokkos<LO, array_layout, device_type, memory_traits> Reduced_Local_DOF_Graph_Matrix(Reduced_Conductivity_Matrix_Strides); //stores local dof indices
  RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout> Reduced_Conductivity_Matrix(Reduced_Conductivity_Matrix_Strides);

  //template compatible row offsets (may change to shallow copy later if it works on device types etc.)
  row_pointers reduced_row_offsets = Reduced_DOF_Graph_Matrix.start_index_;
    for(int ipass = 0; ipass < local_nrows_reduced+1; ipass++){
      reduced_row_offsets_pass(ipass) = reduced_row_offsets(ipass);
    }
  
  //construct maps to define set of global indices for the reduced node set
  //stores indices that aren't contiguous due to removal of BCS
  local_reduced_dof_original_map =
    Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(nrows_reduced,Free_Indices.get_kokkos_view(),0,comm) );

  //stores contiguous indices with an unbalanced local distribution
  local_reduced_dof_map = 
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
  size_t all_nrows_reduced = local_nrows_reduced + nghost_nodes;
  for(LO i=nlocal_nodes; i < nall_nodes; i++){
      if((Node_DOF_Boundary_Condition_Type(i)==TEMPERATURE_CONDITION))
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
  for(LO i=0; i < nall_nodes; i++)
    if((Node_DOF_Boundary_Condition_Type(i)!=TEMPERATURE_CONDITION)){
        All_Free_Indices(reduced_index) = all_dof_map->getGlobalElement(i);
        reduced_index++;
      }
  
  //construct map to define set of global indices for the reduced node set including ghosts
  //passing invalid forces the map to count the global elements
  all_reduced_dof_original_map =
    Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(Teuchos::OrdinalTraits<GO>::invalid(),All_Free_Indices.get_kokkos_view(),0,comm) );

  //debug print
  /*
  *fos << "All reduced dof original indices :" << std::endl;
  local_reduced_dof_original_map->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;
  std::fflush(stdout);

  //debug print
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
  *fos << "All reduced global indices :" << std::endl;
  all_reduced_global_indices->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;
  std::fflush(stdout);
  */
  
  //store the new global indices for the reduced matrix graph
  for(LO i=0; i < local_nrows_reduced; i++){
    access_index = local_dof_map->getLocalElement(Free_Indices(i));
    row_counter = 0;
    for(LO j=0; j < Conductivity_Matrix_Strides(access_index); j++){
      global_dof_index = DOF_Graph_Matrix(access_index,j);
      row_access_index = all_dof_map->getLocalElement(global_dof_index);
      if((Node_DOF_Boundary_Condition_Type(row_access_index)!=TEMPERATURE_CONDITION)){
        reduced_local_dof_index = all_reduced_dof_original_map->getLocalElement(global_dof_index);
        //std::cout << "REDUCED LOCAL INDEX ON TASK " << myrank << " is " << Reduced_Conductivity_Matrix_Strides(i) << Reduced_DOF_Graph_Matrix(i,row_counter++) << std::endl;
        Reduced_DOF_Graph_Matrix(i,row_counter++) = all_reduced_global_indices_host(reduced_local_dof_index,0);
      }
    }
  }
  
  //store reduced conductivity matrix values
  for(LO i=0; i < local_nrows_reduced; i++){
    access_index = local_dof_map->getLocalElement(Free_Indices(i));
    row_counter = 0;
    for(LO j=0; j < Conductivity_Matrix_Strides(access_index); j++){
      global_dof_index = DOF_Graph_Matrix(access_index,j);
      row_access_index = all_dof_map->getLocalElement(global_dof_index);
      if((Node_DOF_Boundary_Condition_Type(row_access_index)!=TEMPERATURE_CONDITION)){
        Reduced_Conductivity_Matrix(i,row_counter++) = Conductivity_Matrix(access_index,j);
      }
    }
  }
  
  // create a Map for the reduced global conductivity matrix that is evenly distributed amongst mpi ranks
  local_balanced_reduced_dof_map = 
    Teuchos::rcp( new Tpetra::Map<LO,GO,node_type>(nrows_reduced,0,comm));

    //debug print
  /*
  *fos << "local_reduced_dof_map :" << std::endl;
  local_reduced_dof_map->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;
  std::fflush(stdout);

  //debug print
  *fos << "local_balanced_reduced_dof_map :" << std::endl;
  local_balanced_reduced_dof_map->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;
  std::fflush(stdout);
  */

  //build column map
  Teuchos::RCP<const Tpetra::Map<LO,GO,node_type> > colmap;
  const Teuchos::RCP<const Tpetra::Map<LO,GO,node_type> > dommap = local_reduced_dof_map;

  //debug print
  /*
  if(myrank==4){
  std::cout << "Reduced DOF Graph Matrix on Rank " << myrank << std::endl;
  for(LO i=0; i < local_nrows_reduced; i++){
    for(LO j=0; j < Reduced_Conductivity_Matrix_Strides(i); j++){
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
    for(LO j=0; j < Reduced_Conductivity_Matrix_Strides(i); j++)
      Reduced_Local_DOF_Graph_Matrix(i,j) = colmap->getLocalElement(Reduced_DOF_Graph_Matrix(i,j));
  
  Teuchos::RCP<MAT> unbalanced_A = Teuchos::rcp(new MAT(local_reduced_dof_map, colmap, reduced_row_offsets_pass,
                   Reduced_Local_DOF_Graph_Matrix.get_kokkos_view(), Reduced_Conductivity_Matrix.get_kokkos_view()));
  unbalanced_A->fillComplete();
  Teuchos::RCP<const_MAT> const_unbalanced_A(new const_MAT(*unbalanced_A));
  
  //This completes the setup for A matrix of the linear system

  //debug print of A matrix before balancing
  //*fos << "Reduced Stiffness Matrix :" << std::endl;
  //const_unbalanced_A->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;

  //communicate reduced conductivity matrix entries for better load balancing
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
  X = Teuchos::rcp(new MV(local_balanced_reduced_dof_map, 1));
  //return !EXIT_SUCCESS;
  //X->randomize();

  // Create Kokkos view of RHS B vector (Force Vector)  
  vec_array Bview_pass = vec_array("Bview_pass", local_nrows_reduced,1);

  //set bview to force vector data
  for(LO i=0; i < local_nrows_reduced; i++){
    access_index = local_dof_map->getLocalElement(Free_Indices(i));
    Bview_pass(i,0) = Nodal_Forces(access_index,0);
  }
  
  unbalanced_B = Teuchos::rcp(new MV(local_reduced_dof_map, Bview_pass));
  
  //import object to rebalance force vector
  Tpetra::Import<LO, GO> Bvec_importer(local_reduced_dof_map, local_balanced_reduced_dof_map);

  balanced_B = Teuchos::rcp(new MV(local_balanced_reduced_dof_map, 1));
  
  //comms to rebalance force vector
  balanced_B->doImport(*unbalanced_B, Bvec_importer, Tpetra::INSERT);
  
  //if(myrank==0)
  //*fos << "RHS :" << std::endl;
  //balanced_B->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;

  //debug print
  //if(update_count==42){
    //Tpetra::MatrixMarket::Writer<MAT> market_writer();
    //Tpetra::MatrixMarket::Writer<MAT>::writeSparseFile("A_matrix.txt", *balanced_A, "A_matrix", "Stores conductivity matrix values");
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
    using impl_scalar_type =
      typename Kokkos::Details::ArithTraits<real_t>::val_type;
    using mag_type = typename Kokkos::ArithTraits<impl_scalar_type>::mag_type;
    // Instead of checking each time for rank, create a rank 0 stream

    xbalanced_B = Teuchos::rcp(new Xpetra::TpetraMultiVector<real_t,LO,GO,node_type>(balanced_B));
    xX = Teuchos::rcp(new Xpetra::TpetraMultiVector<real_t,LO,GO,node_type>(X));
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
    Teuchos::RCP<MV> tnullspace = Teuchos::rcp(new MV(local_balanced_reduced_dof_map, 1));
    tnullspace->putScalar(1);
    Teuchos::RCP<Xpetra::MultiVector<real_t,LO,GO,node_type>> nullspace = Teuchos::rcp(new Xpetra::TpetraMultiVector<real_t,LO,GO,node_type>(tnullspace));

    Teuchos::RCP<Xpetra::MultiVector<real_t,LO,GO,node_type>> material = Teuchos::null;
    Teuchos::RCP<Xpetra::CrsMatrix<real_t,LO,GO,node_type>> xbalanced_A = Teuchos::rcp(new Xpetra::TpetraCrsMatrix<real_t,LO,GO,node_type>(balanced_A));
    xwrap_balanced_A = Teuchos::rcp(new Xpetra::CrsMatrixWrap<real_t,LO,GO,node_type>(xbalanced_A));
    //xwrap_balanced_A->SetFixedBlockSize(1);
   
    //randomize initial vector
    xX->setSeed(100);
    xX->randomize();
    
     if(simparam->equilibrate_matrix_flag){
      Solver_Pointer_->equilibrateMatrix(xwrap_balanced_A,"diag");
      Solver_Pointer_->preScaleRightHandSides(*balanced_B,"diag");
      Solver_Pointer_->preScaleInitialGuesses(*X,"diag");
     }

    //debug print
    //if(myrank==0)
    //*fos << "Xpetra A matrix :" << std::endl;
    //xX->describe(*fos,Teuchos::VERB_EXTREME);
    //*fos << std::endl;
    //std::fflush(stdout);
    
    int num_iter = 2000;
    double solve_tol = 1e-06;
    int cacheSize = 0;
    std::string solveType         = "belos";
    std::string belosType         = "cg";
    // =========================================================================
    // Preconditioner construction
    // =========================================================================
    //bool useML   = Linear_Solve_Params->isParameter("use external multigrid package") && (Linear_Solve_Params->get<std::string>("use external multigrid package") == "ml");
    //out<<"*********** MueLu ParameterList ***********"<<std::endl;
    //out<<*Linear_Solve_Params;
    //out<<"*******************************************"<<std::endl;
    
    //xwrap_balanced_A->describe(*fos,Teuchos::VERB_EXTREME);
    
    comm->barrier();
    //PreconditionerSetup(A,coordinates,nullspace,material,paramList,false,false,useML,0,H,Prec);
    if(Hierarchy_Constructed){
      ReuseXpetraPreconditioner(xwrap_balanced_A, H);
    }
    else{
      PreconditionerSetup(xwrap_balanced_A,coordinates,nullspace,material,*Linear_Solve_Params,false,false,false,0,H,Prec);
      Hierarchy_Constructed = true;
    }
    comm->barrier();
    
    //H->Write(-1, -1);
    //H->describe(*fos,Teuchos::VERB_EXTREME);
    

    // =========================================================================
    // System solution (Ax = b)
    // =========================================================================

    real_t current_cpu_time = Solver_Pointer_->CPU_Time();
    SystemSolve(xwrap_balanced_A,xX,xbalanced_B,H,Prec,*fos,solveType,belosType,false,false,false,cacheSize,0,true,true,num_iter,solve_tol);
    linear_solve_time += Solver_Pointer_->CPU_Time() - current_cpu_time;
    comm->barrier();

    if(simparam->equilibrate_matrix_flag){
      Solver_Pointer_->postScaleSolutionVectors(*X,"diag");
    }

    if(simparam->multigrid_timers){
    Teuchos::RCP<Teuchos::ParameterList> reportParams = rcp(new Teuchos::ParameterList);
        reportParams->set("How to merge timer sets",   "Union");
        reportParams->set("alwaysWriteLocal",          false);
        reportParams->set("writeGlobalStats",          true);
        reportParams->set("writeZeroTimers",           false);
    std::ios_base::fmtflags ff(fos->flags());
    *fos << std::fixed;
    Teuchos::TimeMonitor::report(comm.ptr(), *fos, "", reportParams);
    *fos << std::setiosflags(ff);
    //xwrap_balanced_A->describe(*fos,Teuchos::VERB_EXTREME);
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
  Teuchos::RCP<MV> reduced_node_temperatures_distributed = Teuchos::rcp(new MV(local_reduced_dof_map, 1));
  //create import object using local node indices map and all indices map
  Tpetra::Import<LO, GO> reduced_temperatures_importer(local_balanced_reduced_dof_map, local_reduced_dof_map);

  //comms to get displacements on reduced unbalanced displacement vector
  reduced_node_temperatures_distributed->doImport(*X, reduced_temperatures_importer, Tpetra::INSERT);

  //populate node displacement multivector on the local dof map
  const_host_vec_array reduced_node_temperatures_host = reduced_node_temperatures_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  host_vec_array node_temperatures_host = node_temperatures_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);

  for(int init = 0; init < local_dof_map->getNodeNumElements(); init++)
    node_temperatures_host(init,0) = 0;

  for(LO i=0; i < local_nrows_reduced; i++){
    access_index = local_dof_map->getLocalElement(Free_Indices(i));
    node_temperatures_host(access_index,0) = reduced_node_temperatures_host(i,0);
  }
  
  //import for displacement of ghosts
  Tpetra::Import<LO, GO> ghost_temperature_importer(local_dof_map, all_dof_map);

  //comms to get displacements on all node map
  all_node_temperatures_distributed->doImport(*node_temperatures_distributed, ghost_displacement_importer, Tpetra::INSERT);

  //if(myrank==0)
  //*fos << "All displacements :" << std::endl;
  //all_node_temperatures_distributed->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  
  return !EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------------------------
   Communicate ghosts using the current optimization design data
---------------------------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::comm_variables(Teuchos::RCP<const MV> zp){
  
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

void FEA_Module_Heat_Conduction::update_linear_solve(Teuchos::RCP<const MV> zp){
  
  //set density vector to the current value chosen by the optimizer
  test_node_densities_distributed = zp;

  assemble_matrix();

  if(body_term_flag)
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
