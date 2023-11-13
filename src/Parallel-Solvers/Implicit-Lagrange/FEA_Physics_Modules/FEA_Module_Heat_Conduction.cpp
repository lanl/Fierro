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
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"
#include "Tpetra_Details_FixedHashTable.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Import_Util2.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include <set>

#include "elements.h"
#include "swage.h"
#include "matar.h"
#include "utilities.h"
#include "node_combination.h"
#include "Simulation_Parameters/FEA_Module/Heat_Conduction_Parameters.h"
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
#define BC_EPSILON 1.0e-6

using namespace utils;


FEA_Module_Heat_Conduction::FEA_Module_Heat_Conduction(
    Heat_Conduction_Parameters& in_params,
    Solver *Solver_Pointer, 
    const int my_fea_module_index
  ) : FEA_Module(Solver_Pointer) {

  //assign interfacing index
  my_fea_module_index_ = my_fea_module_index;
  Module_Type = FEA_MODULE_TYPE::Heat_Conduction;

  //recast solver pointer for non-base class access
  Implicit_Solver_Pointer_ = dynamic_cast<Implicit_Solver*>(Solver_Pointer);

  module_params = in_params;
  simparam = Implicit_Solver_Pointer_->simparam;
  
  //TO parameters
  penalty_power = simparam.optimization_options.simp_penalty_power;
  nodal_density_flag = simparam.nodal_density_flag;

  //create ref element object
  //ref_elem = new elements::ref_element();
  //create mesh objects
  //init_mesh = new swage::mesh_t(simparam);
  //mesh = new swage::mesh_t(simparam);
  hessvec_count = update_count = 0;
  linear_solve_time = hessvec_time = hessvec_linear_time = 0;

  //preconditioner construction
  Hierarchy_Constructed = false;

  gradient_print_sync = 0;
  
  //boundary condition data
  max_boundary_sets = max_temp_boundary_sets = max_load_boundary_sets = num_surface_temp_sets = num_surface_flux_sets = 0;

  //boundary condition flags
  matrix_bc_reduced = body_term_flag = thermal_flag = electric_flag = false;

  //construct globally distributed temperature and heat flux vectors
  int num_dim = simparam.num_dims;
  local_dof_map = map;
  all_dof_map = all_node_map;
  node_temperatures_distributed = Teuchos::rcp(new MV(local_dof_map, 1));
  node_temperature_gradients_distributed = Teuchos::rcp(new MV(local_dof_map, 1));
  all_node_temperatures_distributed = Teuchos::rcp(new MV(all_dof_map, 1));
  node_heat_fluxes_distributed = Teuchos::rcp(new MV(map, num_dim));
  Global_Nodal_RHS = Teuchos::rcp(new MV(local_dof_map, 1));
  Global_Nodal_Heat = Teuchos::rcp(new MV(local_dof_map, 1));
  all_node_temperature_gradients_distributed = Teuchos::rcp(new MV(all_node_map, num_dim));
  all_node_heat_fluxes_distributed = Teuchos::rcp(new MV(all_node_map, num_dim));
  adjoints_allocated = false;

  //initialize temperatures to 0
  //local variable for host view in the dual view
  host_vec_array all_node_temperatures = all_node_temperatures_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  host_vec_array node_temperatures = node_temperatures_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  for(int init = 0; init < local_dof_map->getLocalNumElements(); init++)
    node_temperatures(init,0) = 0;
  for(int init = 0; init < all_dof_map->getLocalNumElements(); init++)
    all_node_temperatures(init,0) = 0;

  //construct specific dof map here if needed (more than one dof per node and not the coordinate dof map already provided) from node map
  //this thermal module just equates the dof map and node map since theres one scalar dof per node

  //output setup
  init_output();
}

FEA_Module_Heat_Conduction::~FEA_Module_Heat_Conduction(){ }

/* ----------------------------------------------------------------------------
   Initialize sets of element boundary surfaces and arrays for input conditions
------------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::init_boundaries(){
  max_load_boundary_sets = module_params.loading_conditions.size();
  max_temp_boundary_sets = module_params.boundary_conditions.size();
  max_boundary_sets = max_load_boundary_sets + max_temp_boundary_sets;
  
  // set the number of boundary sets
  if(myrank == 0)
    std::cout << "building boundary sets " << std::endl;

  //initialize to 1 since there must be at least 1 boundary set anyway; read in may occure later
  if(max_boundary_sets==0) max_boundary_sets = 1;
  //std::cout << "NUM BOUNDARY CONDITIONS ON RANK " << myrank << " FOR INIT " << num_boundary_conditions <<std::endl;
  init_boundary_sets(max_boundary_sets);

  //allocate nodal data
  Node_DOF_Boundary_Condition_Type = CArrayKokkos<int, array_layout, device_type, memory_traits>(nall_nodes, "Node_DOF_Boundary_Condition_Type");
  Node_Temperature_Boundary_Conditions = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(nall_nodes);
  Node_Heat_Flux_Boundary_Conditions = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(nall_nodes);

  //initialize
  for(int init=0; init < nall_nodes; init++)
    Node_DOF_Boundary_Condition_Type(init) = NONE;

  Number_DOF_BCS = 0;
}

/* ----------------------------------------------------------------------
   initialize storage for element boundary surfaces corresponding to user BCs
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::init_boundary_sets (int num_sets){

  if(num_sets == 0){
    std::cout << " Warning: number of boundary conditions = 0";
    return;
  }
  //initialize maximum
  max_boundary_sets = num_sets;
  if(max_load_boundary_sets == 0) max_load_boundary_sets = num_sets;
  if(max_temp_boundary_sets == 0) max_temp_boundary_sets = num_sets;
  //std::cout << " DEBUG PRINT "<<num_sets << " " << nboundary_patches << std::endl;
  Boundary_Condition_Type_List = CArrayKokkos<int, array_layout, HostSpace, memory_traits>(num_sets, "Boundary_Condition_Type_List");
  NBoundary_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_sets, "NBoundary_Condition_Patches");
  //std::cout << "NBOUNDARY PATCHES ON RANK " << myrank << " FOR INIT IS " << nboundary_patches <<std::endl;
  Boundary_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_sets, nboundary_patches, "Boundary_Condition_Patches");
  Boundary_Surface_Heat_Flux = CArrayKokkos<real_t, array_layout, HostSpace, memory_traits>(max_load_boundary_sets, 3, "Boundary_Surface_Heat_Flux");
  Boundary_Surface_Temperatures = CArrayKokkos<real_t, array_layout, HostSpace, memory_traits>(max_temp_boundary_sets, 3, "Boundary_Surface_Temperatures");

  //initialize data
  for(int iset = 0; iset < num_sets; iset++) NBoundary_Condition_Patches(iset) = 0;

   //initialize
  for(int ibdy=0; ibdy < num_sets; ibdy++) Boundary_Condition_Type_List(ibdy) = NONE;
}

/* ----------------------------------------------------------------------------
   Grow boundary conditions sets of element boundary surfaces
------------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::grow_boundary_sets(int num_sets){
  int num_dim = simparam.num_dims;

  if(num_sets == 0){
    std::cout << " Warning: number of boundary conditions being set to 0";
    return;
  }

  //std::cout << " DEBUG PRINT "<<num_sets << " " << nboundary_patches << std::endl;
  if(num_sets>max_boundary_sets){
    //temporary storage for previous data
    CArrayKokkos<int, array_layout, HostSpace, memory_traits> Temp_Boundary_Condition_Type_List = Boundary_Condition_Type_List;
    CArrayKokkos<size_t, array_layout, device_type, memory_traits> Temp_NBoundary_Condition_Patches = NBoundary_Condition_Patches;
    CArrayKokkos<size_t, array_layout, device_type, memory_traits> Temp_Boundary_Condition_Patches = Boundary_Condition_Patches;
    
    max_boundary_sets = num_sets + 5; //5 is an arbitrary buffer
    Boundary_Condition_Type_List = CArrayKokkos<int, array_layout, HostSpace, memory_traits>(max_boundary_sets, "Boundary_Condition_Type_List");
    NBoundary_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(max_boundary_sets, "NBoundary_Condition_Patches");
    //std::cout << "NBOUNDARY PATCHES ON RANK " << myrank << " FOR GROW " << nboundary_patches <<std::endl;
    Boundary_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(max_boundary_sets, nboundary_patches, "Boundary_Condition_Patches");

    //copy previous data back over
    //std::cout << "NUM BOUNDARY CONDITIONS ON RANK " << myrank << " FOR COPY " << max_boundary_sets <<std::endl;
    for(int iset = 0; iset < num_boundary_conditions; iset++){
      Boundary_Condition_Type_List(iset) = Temp_Boundary_Condition_Type_List(iset);
      NBoundary_Condition_Patches(iset) = Temp_NBoundary_Condition_Patches(iset);
      for(int ipatch = 0; ipatch < nboundary_patches; ipatch++){
        Boundary_Condition_Patches(iset, ipatch) = Temp_Boundary_Condition_Patches(iset, ipatch);
      }
    }
    
    //initialize data
    for(int iset = num_boundary_conditions; iset < max_boundary_sets; iset++) NBoundary_Condition_Patches(iset) = 0;

    //initialize
    for(int ibdy = num_boundary_conditions; ibdy < max_boundary_sets; ibdy++) Boundary_Condition_Type_List(ibdy) = NONE;
  }
}

/* ----------------------------------------------------------------------------
   Grow storage for temperature boundary conditions
------------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::grow_temperature_condition_sets(int num_sets){
  int num_dim = simparam.num_dims;
  
  if(num_sets == 0){
    std::cout << " Warning: number of boundary conditions being set to 0";
    return;
  }
  
  //std::cout << " DEBUG PRINT "<<num_sets << " " << nboundary_patches << std::endl;
  if(num_sets>max_temp_boundary_sets){
    //temporary storage for previous data
    CArrayKokkos<real_t, array_layout, HostSpace, memory_traits> Temp_Boundary_Surface_Temperatures = Boundary_Surface_Temperatures;
    
    max_temp_boundary_sets = num_sets + 5; //5 is an arbitrary buffer
    Boundary_Surface_Temperatures = CArrayKokkos<real_t, array_layout, HostSpace, memory_traits>(max_temp_boundary_sets, 3, "Boundary_Surface_Temperatures");

    //copy previous data back over
    for(int iset = 0; iset < num_surface_temp_sets; iset++){
      Boundary_Surface_Temperatures(iset,0) = Temp_Boundary_Surface_Temperatures(iset,0);
      Boundary_Surface_Temperatures(iset,1) = Temp_Boundary_Surface_Temperatures(iset,1);
      Boundary_Surface_Temperatures(iset,2) = Temp_Boundary_Surface_Temperatures(iset,2);
    }
  }
  
}

/* ----------------------------------------------------------------------------
   Grow boundary conditions sets of element boundary surfaces
------------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::grow_loading_condition_sets(int num_sets){
  int num_dim = simparam.num_dims;
  
  if(num_sets == 0){
    std::cout << " Warning: number of boundary conditions being set to 0";
    return;
  }
  
  //std::cout << " DEBUG PRINT "<<num_sets << " " << nboundary_patches << std::endl;
  if(num_sets>max_load_boundary_sets){
    //temporary storage for previous data
    CArrayKokkos<real_t, array_layout, HostSpace, memory_traits> Temp_Boundary_Surface_Heat_Flux = Boundary_Surface_Heat_Flux;
    
    max_load_boundary_sets = num_sets + 5; //5 is an arbitrary buffer
    Boundary_Surface_Heat_Flux = CArrayKokkos<real_t, array_layout, HostSpace, memory_traits>(max_load_boundary_sets, 3, "Boundary_Surface_Heat_Flux");

    //copy previous data back over
    for(int iset = 0; iset < num_surface_flux_sets; iset++){
      Boundary_Surface_Heat_Flux(iset,0) = Temp_Boundary_Surface_Heat_Flux(iset,0);
      Boundary_Surface_Heat_Flux(iset,1) = Temp_Boundary_Surface_Heat_Flux(iset,1);
      Boundary_Surface_Heat_Flux(iset,2) = Temp_Boundary_Surface_Heat_Flux(iset,2);
    }
  }
}

/* ----------------------------------------------------------------------
   Assign sets of element boundary surfaces corresponding to user BCs
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::generate_bcs(){
  int num_dim = simparam.num_dims;
  int bc_tag;
  real_t value, temp_temp;
  real_t surface_limits[4];

  for (auto bc : module_params.boundary_conditions) {
    switch (bc.surface.type) {
      case BOUNDARY_TYPE::x_plane:
        bc_tag = 0;
        break;
      case BOUNDARY_TYPE::y_plane:
        bc_tag = 1;
        break;
      case BOUNDARY_TYPE::z_plane:
        bc_tag = 2;
        break;
      default:
        throw std::runtime_error("Invalid surface type: " + to_string(bc.surface.type));
    }

    value = bc.surface.plane_position * simparam.get_unit_scaling();
    
    if(num_boundary_conditions + 1>max_boundary_sets) grow_boundary_sets(num_boundary_conditions+1);
    if(num_surface_temp_sets + 1>max_load_boundary_sets) grow_loading_condition_sets(num_surface_temp_sets+1);
    //tag_boundaries(bc_tag, value, num_boundary_conditions, surface_limits);
    if(bc.surface.use_limits){
      surface_limits[0] = bc.surface.surface_limits_sl;
      surface_limits[1] = bc.surface.surface_limits_su;
      surface_limits[2] = bc.surface.surface_limits_tl;
      surface_limits[3] = bc.surface.surface_limits_tu;
      tag_boundaries(bc_tag, value, num_boundary_conditions, surface_limits);
    }
    else{
      tag_boundaries(bc_tag, value, num_boundary_conditions);
    }

    if (bc.type == BOUNDARY_CONDITION_TYPE::temperature) {
      Boundary_Condition_Type_List(num_boundary_conditions) = TEMPERATURE_CONDITION;
      Boundary_Surface_Temperatures(num_surface_temp_sets,0) = bc.value;
    }
    
    if(Boundary_Surface_Temperatures(num_surface_temp_sets,0)) nonzero_bc_flag = true;
    *fos << "tagging " << bc_tag << " at " << value <<  std::endl;
    *fos << "tagged a set " << std::endl;
    std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(num_boundary_conditions) << std::endl;
    *fos << std::endl;
    num_boundary_conditions++;
    num_surface_temp_sets++;
    
  }
  
  /*
  // tag the z=0 plane,  (Direction, value, bdy_set)
  *fos << "tagging z = 0 " << std::endl;
  bc_tag = 2;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  value = 0.0 * simparam.get_unit_scaling();
  surface_limits[0] = surface_limits[2] = 4;
  surface_limits[1] = surface_limits[3] = 6;
  if(num_boundary_conditions + 1>max_boundary_sets) grow_boundary_sets(num_boundary_conditions+1);
  if(num_surface_temp_sets + 1>max_load_boundary_sets) grow_loading_condition_sets(num_surface_temp_sets+1);
  //tag_boundaries(bc_tag, value, bdy_set_id, surface_limits);
  tag_boundaries(bc_tag, value, num_boundary_conditions);
  Boundary_Condition_Type_List(num_boundary_conditions) = TEMPERATURE_CONDITION;
  Boundary_Surface_Temperatures(num_surface_temp_sets,0) = 293;
  if(Boundary_Surface_Temperatures(num_surface_temp_sets,0)) nonzero_bc_flag = true;
    
  *fos << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(num_boundary_conditions) << std::endl;
  *fos << std::endl;
  
  num_boundary_conditions++;
  num_surface_temp_sets++;
  */
  //Tag nodes for Boundary conditions such as temperatures
  Temperature_Boundary_Conditions();
} // end generate_bcs

/* ----------------------------------------------------------------------
   Assign sets of element boundary surfaces corresponding to user BCs
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::generate_applied_loads() {
  int num_dim = simparam.num_dims;
  int bc_tag, dim1_other, dim2_other;
  real_t value;
  real_t surface_limits[4];
  
  //Surface Fluxes Section
  //find user flux settings
  std::string fea_module_base = "fea_module_";
  std::string bc_base = ":loading_condition_";
  std::string index, bc_index;
  std::string fea_module_name, bc_name;
  
  index = std::to_string(my_fea_module_index_+1);
  bc_index = std::to_string(num_surface_flux_sets+1);
  fea_module_name = fea_module_base + index;
  bc_name = fea_module_name + bc_base + bc_index;

  double unit_scaling = simparam.get_unit_scaling();
  for (auto lc : module_params.loading_conditions) {
    switch (lc->surface.type) {
      case BOUNDARY_TYPE::x_plane:
        bc_tag = 0;
        break;
      case BOUNDARY_TYPE::y_plane:
        bc_tag = 1;
        break;
      case BOUNDARY_TYPE::z_plane:
        bc_tag = 2;
        break;
      default:
        throw std::runtime_error("Invalid surface type: " + to_string(lc->surface.type));
    }
    
    value = lc->surface.plane_position * unit_scaling;
    
    if(num_boundary_conditions + 1>max_boundary_sets) grow_boundary_sets(num_boundary_conditions+1);
    if(num_surface_flux_sets + 1>max_load_boundary_sets) grow_loading_condition_sets(num_surface_flux_sets+1);
    //tag_boundaries(bc_tag, value, num_boundary_conditions, surface_limits);
    if(lc->surface.use_limits){
      surface_limits[0] = lc->surface.surface_limits_sl;
      surface_limits[1] = lc->surface.surface_limits_su;
      surface_limits[2] = lc->surface.surface_limits_tl;
      surface_limits[3] = lc->surface.surface_limits_tu;
      tag_boundaries(bc_tag, value, num_boundary_conditions, surface_limits);
    }
    else{
      tag_boundaries(bc_tag, value, num_boundary_conditions);
    }
    lc->apply(
      [&](const Surface_Flux_Condition& lc) { 
        Boundary_Condition_Type_List(num_boundary_conditions) = SURFACE_LOADING_CONDITION;
        if(lc.specification == LOADING_SPECIFICATION::normal){
          if(bc_tag==0){
            dim1_other = 1;
            dim2_other = 2;
          }
          else if(bc_tag==1){
            dim1_other = 0;
            dim2_other = 2;
          }
          else if(bc_tag==2){
            dim1_other = 0;
            dim2_other = 1;
          }
          if (lc.flux_value){
            Boundary_Surface_Heat_Flux(num_surface_flux_sets,bc_tag) = lc.flux_value / unit_scaling / unit_scaling;
            Boundary_Surface_Heat_Flux(num_surface_flux_sets,dim1_other) = 0;
            Boundary_Surface_Heat_Flux(num_surface_flux_sets,dim2_other) = 0;
          }
        }
      }
    );

    *fos << "tagging " << bc_tag << " at " << value <<  std::endl;
    
    *fos << "tagged a set " << std::endl;
    std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(num_boundary_conditions) << std::endl;
    *fos << std::endl;
    num_boundary_conditions++;
    num_surface_flux_sets++;
    
    bc_index = std::to_string(num_surface_flux_sets+1);
    bc_name = fea_module_name + bc_base + bc_index;
  }
  /*
  *fos << "tagging beam +z heat flux " << std::endl;
  bc_tag = 2;  // bc_tag = 0 xplane, 1 yplane, 2 zplane, 3 cylinder, 4 is shell
  //value = 0;
  value = 100;
  
  //grow arrays as needed
  if(num_boundary_conditions + 1>max_boundary_sets) grow_boundary_sets(num_boundary_conditions+1);
  if(num_surface_flux_sets + 1>max_load_boundary_sets) grow_loading_condition_sets(num_surface_flux_sets+1);

  //find boundary patches this BC corresponds to
  tag_boundaries(bc_tag, value, num_boundary_conditions);
  Boundary_Condition_Type_List(num_boundary_conditions) = SURFACE_LOADING_CONDITION;
  Boundary_Surface_Heat_Flux(num_surface_flux_sets,0) = 0;
  Boundary_Surface_Heat_Flux(num_surface_flux_sets,1) = 0;
  Boundary_Surface_Heat_Flux(num_surface_flux_sets,2) = -0.1/simparam.get_unit_scaling()/simparam.get_unit_scaling();

  *fos << "tagged a set " << std::endl;
  std::cout << "number of bdy patches in this set = " << NBoundary_Condition_Patches(num_boundary_conditions) << std::endl;
  *fos << std::endl;
  
  num_boundary_conditions++;
  num_surface_flux_sets++;
  */
  //Body Term Section

  //apply body terms
  thermal_flag = module_params.thermal_flag;

  if(electric_flag||thermal_flag) body_term_flag = true;

}

/* ----------------------------------------------------------------------
   Initialize global vectors and array maps needed for matrix assembly
------------------------------------------------------------------------- */
void FEA_Module_Heat_Conduction::init_assembly(){
  int num_dim = simparam.num_dims;
  const_host_elem_conn_array nodes_in_elem = global_nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
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
  Graph_Matrix = RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits>(Graph_Matrix_Strides);
  for(int inode = 0; inode < nlocal_nodes; inode++)
    for(int istride = 0; istride < Graph_Matrix_Strides(inode); istride++)
      Graph_Matrix(inode,istride) = Repeat_Graph_Matrix(inode,istride);

  //deallocate repeat matrix
  
  /*At this stage the sparse graph should have unique global indices on each row.
    The constructed Assembly map (to the global sparse matrix)
    is used to loop over each element's local conductivity matrix in the assembly process.*/
  
  Conductivity_Matrix_Strides = Graph_Matrix_Strides;

  Conductivity_Matrix = RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout>(Conductivity_Matrix_Strides);
  DOF_Graph_Matrix = RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits> (Conductivity_Matrix_Strides);
  for (int idof = 0; idof < nlocal_nodes; idof++){
    for (int istride = 0; istride < Conductivity_Matrix_Strides(idof); istride++){
      DOF_Graph_Matrix(idof,istride) = Graph_Matrix(idof,istride);
      //debug print
      //std::cout << "{" <<istride + 1 << "," << DOF_Graph_Matrix(idof,istride) << "} ";
    }
    //debug print
    //std::cout << std::endl;
  }

  //construct distributed conductivity matrix and force vector from local kokkos data
  
  //build column map for the global conductivity matrix
  Teuchos::RCP<const Tpetra::Map<LO,GO,node_type> > colmap;
  const Teuchos::RCP<const Tpetra::Map<LO,GO,node_type> > dommap = local_dof_map;

  Tpetra::Details::makeColMap<LO,GO,node_type>(colmap,dommap,DOF_Graph_Matrix.get_kokkos_view(), nullptr);

  size_t nnz = DOF_Graph_Matrix.size();

  //debug print
  //std::cout << "DOF GRAPH SIZE ON RANK " << myrank << " IS " << nnz << std::endl;
  
  //local indices in the graph using the constructed column map
  CArrayKokkos<LO, array_layout, device_type, memory_traits> conductivity_local_indices(nnz, "conductivity_local_indices");
  
  //row offsets with compatible template arguments
    Kokkos::View<size_t *,array_layout, device_type, memory_traits> row_offsets = DOF_Graph_Matrix.start_index_;
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
  
  //sort values and indices
  Tpetra::Import_Util::sortCrsEntries<row_pointers, indices_array, values_array>(row_offsets_pass, conductivity_local_indices.get_kokkos_view(), Conductivity_Matrix.get_kokkos_view());

  Global_Conductivity_Matrix = Teuchos::rcp(new MAT(local_dof_map, colmap, row_offsets_pass, conductivity_local_indices.get_kokkos_view(), Conductivity_Matrix.get_kokkos_view()));
  Global_Conductivity_Matrix->fillComplete();
  
}

/* ----------------------------------------------------------------------
   Assemble the Sparse Conductivity Matrix
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::assemble_matrix(){
  int num_dim = simparam.num_dims;
  const_host_elem_conn_array nodes_in_elem = global_nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
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

  //reset unsorted DOF Graph corresponding to assembly mapped values
  for (int idof = 0; idof < nlocal_nodes; idof++){
    for (int istride = 0; istride < Conductivity_Matrix_Strides(idof); istride++){
      DOF_Graph_Matrix(idof,istride) = Graph_Matrix(idof,istride);
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

  matrix_bc_reduced = false;

  
  Teuchos::RCP<const Tpetra::Map<LO,GO,node_type> > colmap = Global_Conductivity_Matrix->getCrsGraph()->getColMap();

  //unsorted local column indices

  size_t nnz = DOF_Graph_Matrix.size();

  //debug print
  //std::cout << "DOF GRAPH SIZE ON RANK " << myrank << " IS " << nnz << std::endl;
  
  //local indices in the graph using the constructed column map
  CArrayKokkos<LO, array_layout, device_type, memory_traits> conductivity_local_indices(nnz, "conductivity_local_indices");

  //row offsets
  Kokkos::View<size_t *,array_layout, device_type, memory_traits> row_offsets = DOF_Graph_Matrix.start_index_;
  row_pointers row_offsets_pass("row_offsets", nlocal_nodes+1);
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

  //sort values and indices
  Tpetra::Import_Util::sortCrsEntries<row_pointers, indices_array, values_array>(row_offsets_pass, conductivity_local_indices.get_kokkos_view(), Conductivity_Matrix.get_kokkos_view());

  //set global indices for DOF graph from sorted local indices
  entrycount = 0;
  for(int irow = 0; irow < nlocal_nodes; irow++){
    for(int istride = 0; istride < Conductivity_Matrix_Strides(irow); istride++){
      DOF_Graph_Matrix(irow,istride) = colmap->getGlobalElement(conductivity_local_indices(entrycount));
      entrycount++;
    }
  }

  //filter small negative numbers (that should have been 0 from cancellation) from floating point error
  /*
  for (int idof = 0; idof < num_dim*nlocal_nodes; idof++){
    for (int istride = 0; istride < Conductivity_Matrix_Strides(idof); istride++){
      if(Conductivity_Matrix(idof,istride)<0.000000001*simparam.Elastic_Modulus*density_epsilon||Conductivity_Matrix(idof,istride)>-0.000000001*simparam.Elastic_Modulus*density_epsilon)
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
  //*fos << "Global Conductivity Matrix :" << std::endl;
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
  const_host_elem_conn_array nodes_in_elem = global_nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  //local variable for host view in the dual view
  host_vec_array Nodal_RHS = Global_Nodal_RHS->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  const_host_vec_array Element_Densities;
  //local variable for host view of densities from the dual view
  //bool nodal_density_flag = simparam.nodal_density_flag;
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
  int num_dim = simparam.num_dims;
  int nodes_per_elem = max_nodes_per_element;
  int num_gauss_points = simparam.num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  int current_element_index, local_surface_id, surf_dim1, surf_dim2, surface_sign, normal_sign;
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
  real_t heat_flux[3], wedge_product, Jacobian, current_density, weight_multiply, inner_product, surface_normal[3], surface_norm, normal_displacement;
  real_t specific_internal_energy_rate;
  CArray<GO> Surface_Nodes;
  
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
          surface_sign = -1;
        }
        else{
          quad_coordinate(2) = 1;
          surface_sign = 1;
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
          surface_sign = -1;
        }
        else{
          quad_coordinate(1) = 1;
          surface_sign = 1;
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
          surface_sign = -1;
        }
        else{
          quad_coordinate(0) = 1;
          surface_sign = 1;
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
      for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
        //current_node_index = Surface_Nodes(node_loop);
        //local_node_id = all_node_map->getLocalElement(current_node_index);
        local_node_id = all_node_map->getLocalElement(nodes_in_elem(current_element_index, node_loop));
        nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
        nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
        nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
      }

      if(local_surface_id<2){
        //compute shape function derivatives
        elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
        elem->partial_eta_basis(basis_derivative_s2,quad_coordinate);
        elem->partial_mu_basis(basis_derivative_s3,quad_coordinate);
      }
      else if(local_surface_id<4){
        //compute shape function derivatives
        elem->partial_xi_basis(basis_derivative_s1,quad_coordinate);
        elem->partial_mu_basis(basis_derivative_s2,quad_coordinate);
        elem->partial_eta_basis(basis_derivative_s3,quad_coordinate);
      }
      else if(local_surface_id<6){
        //compute shape function derivatives
        elem->partial_eta_basis(basis_derivative_s1,quad_coordinate);
        elem->partial_mu_basis(basis_derivative_s2,quad_coordinate);
        elem->partial_xi_basis(basis_derivative_s3,quad_coordinate);
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

      //derivative of x,y,z w.r.t t
      JT_row3(0) = 0;
      JT_row3(1) = 0;
      JT_row3(2) = 0;
      for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
        JT_row3(0) += nodal_positions(node_loop,0)*basis_derivative_s3(node_loop);
        JT_row3(1) += nodal_positions(node_loop,1)*basis_derivative_s3(node_loop);
        JT_row3(2) += nodal_positions(node_loop,2)*basis_derivative_s3(node_loop);
      }
      

      //compute jacobian for this surface
      //compute product scaling the reference differential surface element to the deformed differential surface element
      wedge_product = sqrt(pow(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2),2)+
               pow(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2),2)+
               pow(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1),2));

      //compute the determinant of the Jacobian
      Jacobian = JT_row1(0)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                 JT_row1(1)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                 JT_row1(2)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1));
      
      //compute surface normal through cross product of the two tangent vectors
      surface_normal[0] = (JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2));
      surface_normal[1] = -(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2));
      if(num_dim==3)
        surface_normal[2] = (JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1));
      else
        surface_normal[2] = 0;
      
      //compute sign of the normal by finding s3 displacement corresponding to x,y,z displacement along normal (Cramer's rule)
      normal_displacement = (surface_normal[0]*surface_normal[0]+surface_normal[1]*surface_normal[1]+surface_normal[2]*surface_normal[2])/Jacobian;
      normal_sign = -1;
      if(signbit(normal_displacement)==signbit(surface_sign)){
        normal_sign = 1;
      }
      //normalize
      surface_normal[0] *= normal_sign/wedge_product;
      surface_normal[1] *= normal_sign/wedge_product;
      surface_normal[2] *= normal_sign/wedge_product;
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

    //apply body terms
    if(body_term_flag){
      //initialize quadrature data
      elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
      elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);
      direct_product_count = std::pow(num_gauss_points,num_dim);

      for(size_t ielem = 0; ielem < rnum_elem; ielem++){

      //acquire set of nodes and nodal temperatures for this local element
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

  //apply contribution from non-zero temperature boundary conditions
    if(nonzero_bc_flag){
      for(int irow = 0; irow < nlocal_nodes; irow++){
        for(int istride = 0; istride < Conductivity_Matrix_Strides(irow); istride++){
          node_id = all_node_map->getLocalElement(DOF_Graph_Matrix(irow,istride));
          if(Node_DOF_Boundary_Condition_Type(node_id)==TEMPERATURE_CONDITION&&Node_Temperature_Boundary_Conditions(node_id)){
            Nodal_RHS(irow,0) -= Conductivity_Matrix(irow,istride)*Node_Temperature_Boundary_Conditions(node_id);  
          }
        }//for
      }//for
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
  real_t unit_scaling = simparam.get_unit_scaling();
  int num_dim = simparam.num_dims;
  
  //init 
  specific_internal_energy_rate = 0;
  if(thermal_flag){
    specific_internal_energy_rate += module_params.material.specific_internal_energy_rate*density;
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
  real_t unit_scaling = simparam.get_unit_scaling();
  int num_dim = simparam.num_dims;
  
  //init 
  gradient_specific_internal_energy_rate = 0;\
  if(thermal_flag){
    gradient_specific_internal_energy_rate += module_params.material.specific_internal_energy_rate;
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
  real_t unit_scaling = simparam.get_unit_scaling();
  real_t penalty_product = 1;
  real_t density_epsilon = simparam.optimization_options.density_epsilon;
  if(density < 0) density = 0;
  for(int i = 0; i < penalty_power; i++)
    penalty_product *= density;
  //relationship between density and conductivity
  Element_Conductivity = (density_epsilon + (1 - density_epsilon)*penalty_product)*module_params.material.thermal_conductivity/unit_scaling;
}

/* ----------------------------------------------------------------------
   Retrieve derivative of material properties with respect to local density
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::Gradient_Element_Material_Properties(size_t ielem, real_t &Element_Conductivity_Derivative, real_t density){
  real_t unit_scaling = simparam.get_unit_scaling();
  real_t penalty_product = 1;
  real_t density_epsilon = simparam.optimization_options.density_epsilon;
  Element_Conductivity_Derivative = 0;
  if(density < 0) density = 0;
  for(int i = 0; i < penalty_power - 1; i++)
    penalty_product *= density;
  //relationship between density and conductivity
  Element_Conductivity_Derivative = penalty_power*(1 - density_epsilon)*penalty_product*module_params.material.thermal_conductivity/unit_scaling;
}

/* --------------------------------------------------------------------------------
   Retrieve second derivative of material properties with respect to local density
----------------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::Concavity_Element_Material_Properties(size_t ielem, real_t &Element_Conductivity_Derivative, real_t density){
  real_t unit_scaling = simparam.get_unit_scaling();
  real_t penalty_product = 1;
  real_t density_epsilon = simparam.optimization_options.density_epsilon;
  Element_Conductivity_Derivative = 0;
  if(density < 0) density = 0;
  if(penalty_power>=2){
    for(int i = 0; i < penalty_power - 2; i++)
      penalty_product *= density;
    //relationship between density and conductivity
    Element_Conductivity_Derivative = penalty_power*(penalty_power-1)*(1 - density_epsilon)*penalty_product*module_params.material.thermal_conductivity/unit_scaling;
  }
}

/* ----------------------------------------------------------------------
   Construct the local conductivity matrix
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::local_matrix(int ielem, CArrayKokkos<real_t, array_layout, device_type, memory_traits> &Local_Matrix){
  //local variable for host view in the dual view
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array nodes_in_elem = global_nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_vec_array Element_Densities;
  //local variable for host view of densities from the dual view
  //bool nodal_density_flag = simparam.nodal_density_flag;
  const_host_vec_array all_node_densities;
  if(nodal_density_flag)
  all_node_densities = all_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  else
  Element_Densities = Global_Element_Densities->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
  int num_dim = simparam.num_dims;
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam.num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  size_t local_node_id;

  real_t unit_scaling = simparam.get_unit_scaling();
  bool topology_optimization_on = simparam.topology_optimization_on;
  direct_product_count = std::pow(num_gauss_points,num_dim);
  real_t matrix_term;
  real_t matrix_subterm1, matrix_subterm2, matrix_subterm3, invJacobian, Jacobian, weight_multiply;
  real_t Element_Conductivity;
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

  size_t Brows = num_dim;
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix_contribution(Brows,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix(Brows,elem->num_basis());
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix_contribution(Brows,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix(Brows,elem->num_basis());
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
  for(int ifill=0; ifill < nodes_per_elem; ifill++)
      for(int jfill=0; jfill < nodes_per_elem; jfill++)
      Local_Matrix(ifill,jfill) = 0;

  //B matrix initialization
  for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < nodes_per_elem; icol++){
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
    if(topology_optimization_on){
      Element_Material_Properties((size_t) ielem,Element_Conductivity, current_density);
    }
    else{
      Element_Conductivity = module_params.material.thermal_conductivity/unit_scaling;
    }

    //debug print
    //std::cout << "Element Material Params " << Elastic_Constant << std::endl;

    //compute Isotropic Conductivity (C) matrix (not multiplied by conductivity)
    if(num_dim==2){
      C_matrix(0,0) = -1;
      C_matrix(1,1) = -1;
    }
    if(num_dim==3){
      C_matrix(0,0) = -1;
      C_matrix(1,1) = -1;
      C_matrix(2,2) = -1;
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
    
    for(int ishape=0; ishape < nodes_per_elem; ishape++){
      B_matrix_contribution(0,ishape) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                                         basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
                                         basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(1,ishape) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                                         basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
                                         basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      if(num_dim==3){
        B_matrix_contribution(2,ishape) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
                                           basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
                                           basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      }
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
      for(int icol=0; icol < nodes_per_elem; icol++)
      B_matrix(irow,icol) += B_matrix_contribution(irow,icol);

    //compute the previous multiplied by the Elastic (C) Matrix
    for(int irow=0; irow < Brows; irow++){
      for(int icol=0; icol < nodes_per_elem; icol++){
        CB_matrix_contribution(irow,icol) = 0;
        for(int span=0; span < Brows; span++){
          CB_matrix_contribution(irow,icol) += C_matrix(irow,span)*B_matrix_contribution(span,icol);
        }
      }
    }

    //accumulate CB matrix
    for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < nodes_per_elem; icol++)
      CB_matrix(irow,icol) += CB_matrix_contribution(irow,icol);

    //compute the contributions of this quadrature point to all the local conductivity matrix elements
    for(int ifill=0; ifill < nodes_per_elem; ifill++)
      for(int jfill=ifill; jfill < nodes_per_elem; jfill++){
        matrix_term = 0;
        for(int span = 0; span < Brows; span++){
          matrix_term += B_matrix_contribution(span,ifill)*CB_matrix_contribution(span,jfill);
        }
        Local_Matrix(ifill,jfill) += Element_Conductivity*weight_multiply*matrix_term*invJacobian;
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
      std::cout << " ------------LOCAL Conductivity MATRIX "<< ielem + 1 <<"--------------"<<std::endl;
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
  int num_dim = simparam.num_dims;
  int bc_option;
  int DOF_BC_type;
  real_t temperature;
  CArrayKokkos<int, array_layout, device_type, memory_traits> Temperatures_Conditions(num_dim);
  CArrayKokkos<int, array_layout, device_type, memory_traits> first_condition_per_node(nall_nodes);
  CArray<GO> Surface_Nodes;

  //host view of local nodal temperatures
  host_vec_array node_temperatures_host = node_temperatures_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);

  //initialize to -1 (DO NOT make -1 an index for bdy sets)
  for(int inode = 0 ; inode < nlocal_nodes; inode++)
    first_condition_per_node(inode) = -1;
  
  //scan for surface method of setting fixed nodal temperatures
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
        temperature = Boundary_Surface_Temperatures(surface_temp_set_id,0);
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

  //scan for direct setting of nodal temperatures from input
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
  const_host_elem_conn_array nodes_in_elem = global_nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  
  const_host_vec_array Element_Densities;
  //local variable for host view of densities from the dual view
  //bool nodal_density_flag = simparam.nodal_density_flag;
  const_host_vec_array all_node_densities;
  if(nodal_density_flag)
  all_node_densities = all_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  else
  Element_Densities = Global_Element_Densities->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
  int num_dim = simparam.num_dims;
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam.num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  size_t local_node_id, local_dof_idx, local_dof_idy, local_dof_idz;
  GO current_global_index;

  direct_product_count = std::pow(num_gauss_points,num_dim);
  real_t Element_Conductivity_Gradient, Poisson_Ratio, gradient_specific_internal_energy_rate;
  real_t Element_Conductivity;
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
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_nodal_temperatures(elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());

  size_t Brows = num_dim;
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix_contribution(Brows,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix(Brows,elem->num_basis());
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix_contribution(Brows,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix(Brows,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> C_matrix(Brows,Brows);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> Local_Matrix_Contribution(nodes_per_elem,nodes_per_elem);

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
      for(int icol=0; icol < nodes_per_elem; icol++){
        CB_matrix(irow,icol) = 0;
      }

    //acquire set of nodes and nodal temperatures for this local element
    for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
      local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, node_loop));
      nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
      nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
      nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
      current_nodal_temperatures(node_loop) = all_node_temperatures(local_node_id,0);
      
      if(nodal_density_flag) nodal_density(node_loop) = all_node_densities(local_node_id,0);
      //debug print
      /*
      std::cout << "node index access x "<< local_node_id << std::endl;
      std::cout << "local index access x "<< local_dof_idx << " displacement x " << current_nodal_temperatures(node_loop*num_dim) <<std::endl;
      std::cout << "local index access y "<< local_dof_idy << " displacement y " << current_nodal_temperatures(node_loop*num_dim + 1) << std::endl;
      std::cout << "local index access z "<< local_dof_idz << " displacement z " << current_nodal_temperatures(node_loop*num_dim + 2) << std::endl; 
      */
    }

    //debug print of current_nodal_temperatures
    /*
    std::cout << " ------------nodal temperatures for Element "<< ielem + 1 <<"--------------"<<std::endl;
    std::cout << " { ";
    for (int idof = 0; idof < num_dim*nodes_per_elem; idof++){
      std::cout << idof + 1 << " = " << current_nodal_temperatures(idof) << " , " ;
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
    for(int ishape=0; ishape < nodes_per_elem; ishape++){
      B_matrix_contribution(0,ishape) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                                         basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
                                         basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(1,ishape) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                                         basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
                                         basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      if(num_dim==3){
        B_matrix_contribution(2,ishape) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
                                           basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
                                           basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      }
    }
    
    //look up element material properties at this point as a function of density
    Gradient_Element_Material_Properties(ielem, Element_Conductivity_Gradient, current_density);

    //debug print
    //std::cout << "Element Material Params " << Elastic_Constant << std::endl;

    //compute Isotropic Conductivity (C) matrix; not scaled by k for gradient calculation later
    if(num_dim==2){
      C_matrix(0,0) = -1;
      C_matrix(1,1) = -1;
    }
    if(num_dim==3){
      C_matrix(0,0) = -1;
      C_matrix(1,1) = -1;
      C_matrix(2,2) = -1;
    }

    //compute the previous multiplied by the Elastic (C) Matrix
    for(int irow=0; irow < Brows; irow++){
      for(int icol=0; icol < nodes_per_elem; icol++){
        CB_matrix_contribution(irow,icol) = 0;
        for(int span=0; span < Brows; span++){
          CB_matrix_contribution(irow,icol) += C_matrix(irow,span)*B_matrix_contribution(span,icol);
        }
      }
    }
    
    //compute the contributions of this quadrature point to all the local conductivity matrix elements
    for(int ifill=0; ifill < nodes_per_elem; ifill++){
      for(int jfill=ifill; jfill < nodes_per_elem; jfill++){
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
    for(int ifill=0; ifill < nodes_per_elem; ifill++){
      for(int jfill=ifill; jfill < nodes_per_elem; jfill++){
        if(ifill==jfill)
          inner_product += Local_Matrix_Contribution(ifill, jfill)*current_nodal_temperatures(ifill)*current_nodal_temperatures(jfill);
        else
          inner_product += 2*Local_Matrix_Contribution(ifill, jfill)*current_nodal_temperatures(ifill)*current_nodal_temperatures(jfill);
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
      design_gradients(local_node_id,0) -= inner_product*Element_Conductivity_Gradient*basis_values(igradient)*weight_multiply*invJacobian;
    }

      //evaluate gradient of body force (such as gravity which depends on density) with respect to igradient
    if(body_term_flag){
      //look up element material properties at this point as a function of density
      Gradient_Body_Term(ielem, current_density, gradient_specific_internal_energy_rate);
      for(int igradient=0; igradient < nodes_per_elem; igradient++){
      if(!map->isNodeGlobalElement(nodes_in_elem(ielem, igradient))) continue;
      local_node_id = map->getLocalElement(nodes_in_elem(ielem, igradient));
      
      //compute inner product for this quadrature point contribution
      inner_product = 0;
      for(int ifill=0; ifill < nodes_per_elem; ifill++){
        inner_product += gradient_specific_internal_energy_rate*current_nodal_temperatures(ifill)*basis_values(ifill);
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

/* ----------------------------------------------------------------------------------------------
   Compute the hessian*vector product of heat capacity potential with respect to nodal densities
------------------------------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::compute_adjoint_hessian_vec(const_host_vec_array design_densities, host_vec_array hessvec, Teuchos::RCP<const MV> direction_vec_distributed){
  //local variable for host view in the dual view
  real_t current_cpu_time = Implicit_Solver_Pointer_->CPU_Time();
  const_host_vec_array all_node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_vec_array all_node_temperatures = all_node_temperatures_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_elem_conn_array nodes_in_elem = global_nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  
  const_host_vec_array Element_Densities;
  //local variable for host view of densities from the dual view
  //bool nodal_density_flag = simparam.nodal_density_flag;
  const_host_vec_array all_node_densities;
  if(nodal_density_flag)
  all_node_densities = all_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  else
  Element_Densities = Global_Element_Densities->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
  const_host_vec_array direction_vec = direction_vec_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);

  if(!adjoints_allocated){
    adjoint_temperatures_distributed = Teuchos::rcp(new MV(local_dof_map, 1));
    adjoint_equation_RHS_distributed = Teuchos::rcp(new MV(local_dof_map, 1));
    all_adjoint_temperatures_distributed = Teuchos::rcp(new MV(all_dof_map, 1));
    adjoints_allocated = true;
  }

  Teuchos::RCP<MV> lambda = adjoint_temperatures_distributed;
  Teuchos::RCP<Xpetra::MultiVector<real_t,LO,GO,node_type>> xlambda = Teuchos::rcp(new Xpetra::TpetraMultiVector<real_t,LO,GO,node_type>(lambda));
  Teuchos::RCP<Xpetra::MultiVector<real_t,LO,GO,node_type>> xB = Teuchos::rcp(new Xpetra::TpetraMultiVector<real_t,LO,GO,node_type>(adjoint_equation_RHS_distributed));
  
  host_vec_array adjoint_equation_RHS_view = adjoint_equation_RHS_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
  
  const_host_vec_array lambda_view = lambda->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);

  int num_dim = simparam.num_dims;
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam.num_gauss_points;
  int z_quad,y_quad,x_quad, direct_product_count;
  LO local_node_id, jlocal_node_id, temp_id, local_dof_id, local_reduced_dof_id, local_dof_idx, local_dof_idy, local_dof_idz;
  GO current_global_index, global_dof_id;

  direct_product_count = std::pow(num_gauss_points,num_dim);
  real_t Element_Conductivity_Gradient, Element_Conductivity_Concavity, gradient_specific_internal_energy_rate;
  real_t Element_Conductivity, specific_internal_energy_rate;
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
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_nodal_temperatures(elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_adjoint_temperatures(elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());

  size_t Brows = num_dim;
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix_contribution(Brows,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix(Brows,elem->num_basis());
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix_contribution(Brows,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix(Brows,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> C_matrix(Brows,Brows);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> Local_Matrix_Contribution(nodes_per_elem,nodes_per_elem);

  //initialize weights
  elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
  elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);
  
  real_t current_density = 1;
  
  //direction_vec_distributed->describe(*fos,Teuchos::VERB_EXTREME);

  //initialize gradient value to zero
  for(size_t inode = 0; inode < nlocal_nodes; inode++)
    hessvec(inode,0) = 0;
  
  //initialize RHS vector
  //initialize RHS vector
  for(int i=0; i < local_dof_map->getLocalNumElements(); i++)
    adjoint_equation_RHS_view(i,0) = 0;
  
  //sum components of direction vector
  direction_vec_reduce = local_direction_vec_reduce = 0;
  for(int i = 0; i < nlocal_nodes; i++)
    local_direction_vec_reduce += direction_vec(i,0);
  
  MPI_Allreduce(&local_direction_vec_reduce,&direction_vec_reduce,1,MPI_DOUBLE,MPI_SUM,world);

  //comms to get ghost components of direction vector needed for matrix inner products
  //Tpetra::Import<LO, GO> node_importer(map, all_node_map);
  
  Teuchos::RCP<MV> all_direction_vec_distributed = Teuchos::rcp(new MV(all_node_map, 1));
  //comms to get ghosts
  all_direction_vec_distributed->doImport(*direction_vec_distributed, *importer, Tpetra::INSERT);
  
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
      for(int icol=0; icol < nodes_per_elem; icol++){
        CB_matrix(irow,icol) = 0;
      }

    //acquire set of nodes and nodal temperatures for this local element
    for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
      local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, node_loop));
      nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
      nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
      nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
      current_nodal_temperatures(node_loop) = all_node_temperatures(local_node_id,0);
      
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
      for(int ishape=0; ishape < nodes_per_elem; ishape++){
      B_matrix_contribution(0,ishape) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                                         basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
                                         basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(1,ishape) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                                         basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
                                         basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      if(num_dim==3){
        B_matrix_contribution(2,ishape) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
                                           basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
                                           basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      }
    }
    
    //look up element material properties at this point as a function of density
    Gradient_Element_Material_Properties(ielem, Element_Conductivity_Gradient, current_density);

    //debug print
    //std::cout << "Element Material Params " << Elastic_Constant << std::endl;

    //compute Isotropic Conductivity (C) matrix
    if(num_dim==2){
      C_matrix(0,0) = -1;
      C_matrix(1,1) = -1;
    }
    if(num_dim==3){
      C_matrix(0,0) = -1;
      C_matrix(1,1) = -1;
      C_matrix(2,2) = -1;
    }

    //compute the previous multiplied by the Elastic (C) Matrix
    for(int irow=0; irow < Brows; irow++){
      for(int icol=0; icol < nodes_per_elem; icol++){
        CB_matrix_contribution(irow,icol) = 0;
        for(int span=0; span < Brows; span++){
          CB_matrix_contribution(irow,icol) += C_matrix(irow,span)*B_matrix_contribution(span,icol);
        }
      }
    }
    
    //compute the contributions of this quadrature point to all the local conductivity matrix elements
    for(int ifill=0; ifill < nodes_per_elem; ifill++){
      for(int jfill=ifill; jfill < nodes_per_elem; jfill++){
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
      for(int ifill=0; ifill < nodes_per_elem; ifill++){
        global_dof_id = nodes_in_elem(ielem, ifill);
        local_dof_id = all_node_map->getLocalElement(global_dof_id);
        if(Node_DOF_Boundary_Condition_Type(local_dof_id)!=TEMPERATURE_CONDITION&&local_dof_map->isNodeGlobalElement(global_dof_id)){
          inner_product = 0;
          for(int jfill=0; jfill < nodes_per_elem; jfill++){
            inner_product += Local_Matrix_Contribution(ifill, jfill)*current_nodal_temperatures(jfill);
            //debug
            //if(Local_Matrix_Contribution(ifill, jfill)<0) Local_Matrix_Contribution(ifill, jfill) = - Local_Matrix_Contribution(ifill, jfill);
            //inner_product += Local_Matrix_Contribution(ifill, jfill);
          }
          adjoint_equation_RHS_view(local_dof_id,0) += 2*inner_product*Element_Conductivity_Gradient*basis_values(igradient)*weight_multiply*all_direction_vec(local_node_id,0)*invJacobian;
        }
      }
      } //density gradient loop
    }//quadrature loop
  }//element index loop
  
  //set adjoint equation RHS terms to 0 if they correspond to a boundary constraint DOF index
  for(int i=0; i < local_dof_map->getLocalNumElements(); i++){
    if(Node_DOF_Boundary_Condition_Type(i)==TEMPERATURE_CONDITION)
      adjoint_equation_RHS_view(i,0) = 0;
  }


  //assign old Conductivity matrix entries
  if(!matrix_bc_reduced){
  LO stride_index;
  for(LO i=0; i < nlocal_nodes; i++){
    if((Node_DOF_Boundary_Condition_Type(i)==TEMPERATURE_CONDITION)){
      for(LO j = 0; j < Conductivity_Matrix_Strides(i); j++){
        global_dof_id = DOF_Graph_Matrix(i,j);
        local_dof_id = all_dof_map->getLocalElement(global_dof_id);
        Original_Conductivity_Entries(i,j) = Conductivity_Matrix(i,j);
        Original_Conductivity_Entry_Indices(i,j) = j;
        if(local_dof_id == i){
          Conductivity_Matrix(i,j) = 1;
        }
        else{     
          Conductivity_Matrix(i,j) = 0;
        }
      }//stride for
    }
    else{
      stride_index = 0;
      for(LO j = 0; j < Conductivity_Matrix_Strides(i); j++){
        global_dof_id = DOF_Graph_Matrix(i,j);
        local_dof_id = all_dof_map->getLocalElement(global_dof_id);
        if((Node_DOF_Boundary_Condition_Type(local_dof_id)==TEMPERATURE_CONDITION)){
          Original_Conductivity_Entries(i,stride_index) = Conductivity_Matrix(i,j);
          Original_Conductivity_Entry_Indices(i,stride_index) = j;   
          Conductivity_Matrix(i,j) = 0;
          stride_index++;
        }
      }
    }
  }//row for
  
  matrix_bc_reduced = true;
  }
  
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
  //xA->describe(*fos,Teuchos::VERB_EXTREME);
  if(module_params.equilibrate_matrix_flag){
    Implicit_Solver_Pointer_->preScaleRightHandSides(*Global_Nodal_RHS,"diag");
    Implicit_Solver_Pointer_->preScaleInitialGuesses(*lambda,"diag");
  }
  real_t current_cpu_time2 = Implicit_Solver_Pointer_->CPU_Time();
  comm->barrier();
  SystemSolve(xA,xlambda,xB,H,Prec,*fos,solveType,belosType,false,false,false,cacheSize,0,true,true,num_iter,solve_tol);
  comm->barrier();
  hessvec_linear_time += Implicit_Solver_Pointer_->CPU_Time() - current_cpu_time2;

  if(module_params.equilibrate_matrix_flag){
    Implicit_Solver_Pointer_->postScaleSolutionVectors(*lambda,"diag");
  }
  //scale by reciprocal ofdirection vector sum
  lambda->scale(1/direction_vec_reduce);
  
  //import for displacement of ghosts
  //Tpetra::Import<LO, GO> ghost_displacement_importer(local_dof_map, all_dof_map);

  //comms to get temperatures on all node map
  all_adjoint_temperatures_distributed->doImport(*adjoint_temperatures_distributed, *importer, Tpetra::INSERT);
  host_vec_array all_adjoint = all_adjoint_temperatures_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
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
      for(int icol=0; icol < nodes_per_elem; icol++){
        CB_matrix(irow,icol) = 0;
      }

    //acquire set of nodes and nodal temperatures for this local element
    for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
      local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, node_loop));
      nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
      nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
      nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
      current_nodal_temperatures(node_loop) = all_node_temperatures(local_node_id,0);
      current_adjoint_temperatures(node_loop) = all_adjoint(local_node_id,0);
      
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
    for(int ishape=0; ishape < nodes_per_elem; ishape++){
      B_matrix_contribution(0,ishape) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                                         basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
                                         basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(1,ishape) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                                         basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
                                         basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      if(num_dim==3){
        B_matrix_contribution(2,ishape) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
                                           basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
                                           basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      }
    }

    //look up element material properties at this point as a function of density
    Concavity_Element_Material_Properties(ielem, Element_Conductivity_Concavity, current_density);
    Gradient_Element_Material_Properties(ielem, Element_Conductivity_Gradient, current_density);

    //*fos << "Elastic Modulus Concavity" << Element_Modulus_Concavity << " " << Element_Modulus_Gradient << std::endl;
    //debug print
    //std::cout << "Element Material Params " << Elastic_Constant << std::endl;

    //compute Isotropic Conductivity (C) matrix
    if(num_dim==2){
      C_matrix(0,0) = -1;
      C_matrix(1,1) = -1;
    }
    if(num_dim==3){
      C_matrix(0,0) = -1;
      C_matrix(1,1) = -1;
      C_matrix(2,2) = -1;
    }

    //compute the previous multiplied by the Elastic (C) Matrix
    for(int irow=0; irow < Brows; irow++){
      for(int icol=0; icol < nodes_per_elem; icol++){
        CB_matrix_contribution(irow,icol) = 0;
        for(int span=0; span < Brows; span++){
          CB_matrix_contribution(irow,icol) += C_matrix(irow,span)*B_matrix_contribution(span,icol);
        }
      }
    }
    
    //compute the contributions of this quadrature point to all the local conductivity matrix elements
      for(int ifill=0; ifill < nodes_per_elem; ifill++){
        for(int jfill=ifill; jfill < nodes_per_elem; jfill++){
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
    for(int ifill=0; ifill < nodes_per_elem; ifill++){
      for(int jfill=ifill; jfill < nodes_per_elem; jfill++){
        if(ifill==jfill)
          inner_product += Local_Matrix_Contribution(ifill, jfill)*current_nodal_temperatures(ifill)*current_nodal_temperatures(jfill);
        else
          inner_product += 2*Local_Matrix_Contribution(ifill, jfill)*current_nodal_temperatures(ifill)*current_nodal_temperatures(jfill);
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
          hessvec(temp_id,0) -= inner_product*Element_Conductivity_Concavity*basis_values(igradient)*all_direction_vec(jlocal_node_id,0)*
                                  basis_values(jgradient)*weight_multiply*invJacobian;
        }
        if(igradient!=jgradient&&map->isNodeGlobalElement(nodes_in_elem(ielem, jgradient))){
          //temp_id = map->getLocalElement(nodes_in_elem(ielem, jgradient));
          hessvec(jlocal_node_id,0) -= inner_product*Element_Conductivity_Concavity*basis_values(igradient)*all_direction_vec(local_node_id,0)*
                                      basis_values(jgradient)*weight_multiply*invJacobian;

        }
      }
    }
    
    //compute inner product for this quadrature point contribution
    inner_product = 0;
    for(int ifill=0; ifill < nodes_per_elem; ifill++){
      for(int jfill=0; jfill < nodes_per_elem; jfill++){
        inner_product += Local_Matrix_Contribution(ifill, jfill)*current_adjoint_temperatures(ifill)*current_nodal_temperatures(jfill);
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
      hessvec(local_node_id,0) += inner_product*direction_vec_reduce*Element_Conductivity_Gradient*basis_values(igradient)*weight_multiply*invJacobian;
    }

      //evaluate gradient of body term (such as energy source which depends on density) with respect to igradient
      if(body_term_flag){
        for(int igradient=0; igradient < nodes_per_elem; igradient++){
        if(!map->isNodeGlobalElement(nodes_in_elem(ielem, igradient))) continue;
        local_node_id = map->getLocalElement(nodes_in_elem(ielem, igradient));
        //look up element material properties at this point as a function of density
        Gradient_Body_Term(ielem, current_density, gradient_specific_internal_energy_rate);
      
        //compute inner product for this quadrature point contribution
        inner_product = 0;
        for(int ifill=0; ifill < nodes_per_elem; ifill++){
          inner_product -= gradient_specific_internal_energy_rate*
                           current_adjoint_temperatures(ifill)*basis_values(ifill/num_dim);
        }
      
        //debug print
        //std::cout << "contribution for " << igradient + 1 << " is " << inner_product << std::endl;
        hessvec(local_node_id,0) += inner_product*direction_vec_reduce*basis_values(igradient)*weight_multiply*Jacobian;
        }
      }
    }
  }//end element loop for hessian vector product
  hessvec_time += Implicit_Solver_Pointer_->CPU_Time() - current_cpu_time;
}

/* ----------------------------------------------------------------------------
   Initialize output data structures
------------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::init_output(){
  //check user parameters for output
  bool output_temperature_flag = simparam.output(FIELD::temperature);
  bool output_temperature_gradient_flag = simparam.output(FIELD::temperature_gradient);
  bool output_heat_flux_flag = simparam.output(FIELD::heat_flux);
  int num_dim = simparam.num_dims;
  
  if(output_temperature_flag){
    output_temperature_index = noutput;
    noutput += 1;
    module_outputs.resize(noutput);

    vector_style.resize(noutput);
    vector_style[noutput-1] = NODAL;

    output_vector_sizes.resize(noutput);
    output_vector_sizes[noutput-1] = 1;

    output_dof_names.resize(noutput);
    output_dof_names[noutput-1].resize(1);
    output_dof_names[noutput-1][0] = "Temperature";
  }
  if(output_temperature_gradient_flag){
    output_temperature_gradient_index = noutput;
    noutput += 1;
    module_outputs.resize(noutput);

    vector_style.resize(noutput);
    vector_style[noutput-1] = NODAL;

    output_vector_sizes.resize(noutput);
    output_vector_sizes[noutput-1] = num_dim;

    output_dof_names.resize(noutput);
    output_dof_names[noutput-1].resize(num_dim);
    output_dof_names[noutput-1][0] = "temperature_gradient_x";
    output_dof_names[noutput-1][1] = "temperature_gradient_y";
    output_dof_names[noutput-1][2] = "temperature_gradient_z";
  }
  if(output_heat_flux_flag){
    output_heat_flux_index = noutput;
    noutput += 1;
    module_outputs.resize(noutput);

    vector_style.resize(noutput);
    vector_style[noutput-1] = NODAL;

    output_vector_sizes.resize(noutput);
    output_vector_sizes[noutput-1] = num_dim;

    output_dof_names.resize(noutput);
    output_dof_names[noutput-1].resize(num_dim);
    output_dof_names[noutput-1][0] = "heat_flux_x";
    output_dof_names[noutput-1][1] = "heat_flux_y";
    output_dof_names[noutput-1][2] = "heat_flux_z";
  }
}

/* -------------------------------------------------------------------------------------------
   Prompts sorting of thermal response output data. For now, nodal temperatures and heat fluxes.
---------------------------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::sort_output(Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > sorted_map){
  
  bool output_temperature_flag = simparam.output(FIELD::temperature);
  bool output_temperature_gradient_flag = simparam.output(FIELD::temperature_gradient);
  bool output_heat_flux_flag = simparam.output(FIELD::heat_flux);
  int num_dim = simparam.num_dims;
  
  //reset modules so that host view falls out of scope
  for(int init = 0; init < noutput; init++){
    const_host_vec_array dummy;
    module_outputs[init] = dummy;
  }
  
  //collect nodal temperature information
  if(output_temperature_flag){
    //importer from local node distribution to collected distribution
    Tpetra::Import<LO, GO> node_sorting_importer(map, sorted_map);

    sorted_node_temperatures_distributed = Teuchos::rcp(new MV(sorted_map, 1));

    //comms to collect
    sorted_node_temperatures_distributed->doImport(*(node_temperatures_distributed), node_sorting_importer, Tpetra::INSERT);

    //set host views of the collected data to print out from
    module_outputs[output_temperature_index] = sorted_node_temperatures_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  
  }

  //collect strain data
  if(output_heat_flux_flag){

    //importer for strains, all nodes to global node set on rank 0
    //Tpetra::Import<LO, GO> strain_collection_importer(all_node_map, global_reduce_map);

    //collected nodal density information
    Teuchos::RCP<MV> sorted_node_heat_fluxes_distributed = Teuchos::rcp(new MV(sorted_map, num_dim));

    //importer from local node distribution to collected distribution
    Tpetra::Import<LO, GO> node_sorting_importer(map, sorted_map);

    //comms to collect
    sorted_node_heat_fluxes_distributed->doImport(*(node_heat_fluxes_distributed), node_sorting_importer, Tpetra::INSERT);

    //debug print
    //std::ostream &out = std::cout;
    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    //if(myrank==0)
    //*fos << "Collected nodal temperatures :" << std::endl;
    //collected_node_strains_distributed->describe(*fos,Teuchos::VERB_EXTREME);
    //*fos << std::endl;
    //std::fflush(stdout);

    //host view to print from
    module_outputs[output_heat_flux_index] = sorted_node_heat_fluxes_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  }
}

/* ------------------------------------------------------------------------------------------------
   Prompts computation of elastic response output data. For now, nodal temperatures and heat fluxes.
--------------------------------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::collect_output(Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > global_reduce_map){
  
  bool output_temperature_flag = simparam.output(FIELD::temperature);
  bool output_temperature_gradient_flag = simparam.output(FIELD::temperature_gradient);
  bool output_heat_flux_flag = simparam.output(FIELD::heat_flux);
  int num_dim = simparam.num_dims;
  
  //reset modules so that host view falls out of scope
  for(int init = 0; init < noutput; init++){
    const_host_vec_array dummy;
    module_outputs[init] = dummy;
  }

  //collect nodal temperature information
  if(output_temperature_flag){
  //importer from local node distribution to collected distribution
  Tpetra::Import<LO, GO> node_collection_importer(map, global_reduce_map);

  collected_node_temperatures_distributed = Teuchos::rcp(new MV(global_reduce_map, 1));

  //comms to collect
  collected_node_temperatures_distributed->doImport(*(node_temperatures_distributed), node_collection_importer, Tpetra::INSERT);

  //set host views of the collected data to print out from
  if(myrank==0){
   module_outputs[output_temperature_index] = collected_node_temperatures_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  }
  }

  //collect strain data
  if(output_heat_flux_flag){

    //importer for strains, all nodes to global node set on rank 0
    //Tpetra::Import<LO, GO> strain_collection_importer(all_node_map, global_reduce_map);

    //collected nodal density information
    collected_node_heat_fluxes_distributed = Teuchos::rcp(new MV(global_reduce_map, num_dim));

    //importer from local node distribution to collected distribution
    Tpetra::Import<LO, GO> node_collection_importer(map, global_reduce_map);

    //comms to collect
    collected_node_heat_fluxes_distributed->doImport(*(node_heat_fluxes_distributed), node_collection_importer, Tpetra::INSERT);

    //debug print
    //std::ostream &out = std::cout;
    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    //if(myrank==0)
    //*fos << "Collected nodal temperatures :" << std::endl;
    //collected_node_strains_distributed->describe(*fos,Teuchos::VERB_EXTREME);
    //*fos << std::endl;
    //std::fflush(stdout);

    //host view to print from
    if(myrank==0)
      module_outputs[output_heat_flux_index] = collected_node_heat_fluxes_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  }
}

/* -------------------------------------------------------------------------------------------
   Prompts computation of thermal response output data. For now, nodal heat fluxes.
---------------------------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::compute_output(){
  compute_nodal_heat_fluxes();
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
  const_host_elem_conn_array nodes_in_elem = global_nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  host_vec_array all_node_heat_fluxes = all_node_heat_fluxes_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  host_vec_array node_heat_fluxes = node_heat_fluxes_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  const_host_elem_conn_array node_nconn = node_nconn_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  const_host_vec_array Element_Densities;
  //local variable for host view of densities from the dual view
  //bool nodal_density_flag = simparam.nodal_density_flag;
  const_host_vec_array all_node_densities;
  if(nodal_density_flag)
  all_node_densities = all_node_densities_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  else
  Element_Densities = Global_Element_Densities->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);

  int num_dim = simparam.num_dims;
  int nodes_per_elem = elem->num_basis();
  int num_gauss_points = simparam.num_gauss_points;
  int flux_max_flag = module_params.flux_max_flag;
  int z_quad,y_quad,x_quad, direct_product_count;
  int solve_flag, zero_flux_flag;
  size_t local_node_id, local_dof_idx, local_dof_idy, local_dof_idz;
  //real_t J_min = std::numeric_limits<real_t>::max();
  GO current_global_index;
  real_t Element_Conductivity;

  direct_product_count = std::pow(num_gauss_points,num_dim);
  real_t matrix_term, current_heat_flux;
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
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> current_nodal_temperatures(elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> nodal_density(elem->num_basis());

  size_t Brows = num_dim;
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix_contribution(Brows,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> B_matrix(Brows,elem->num_basis());
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix_contribution(Brows,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> CB_matrix(Brows,elem->num_basis());
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> C_matrix(Brows,Brows);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> quad_temperature_gradient(Brows);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> quad_heat_flux(Brows);
  FArrayKokkos<real_t, array_layout, device_type, memory_traits> projection_matrix(max_nodes_per_element,max_nodes_per_element);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> projection_vector(Brows,max_nodes_per_element);
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> heat_flux_vector(max_nodes_per_element);
  //Teuchos::SerialSymDenseMatrix<LO,real_t> projection_matrix_pass;
  Teuchos::RCP<Teuchos::SerialDenseMatrix<LO,real_t>> projection_matrix_pass;
  //Teuchos::SerialDenseVector<LO,real_t> projection_vector_pass;
  Teuchos::RCP<Teuchos::SerialDenseVector<LO,real_t>> projection_vector_pass;
  Teuchos::RCP<Teuchos::SerialDenseVector<LO,real_t>> heat_flux_vector_pass;
  //Teuchos::SerialSpdDenseSolver<LO,real_t> projection_solver;
  Teuchos::SerialDenseSolver<LO,real_t> projection_solver;

  //initialize weights
  elements::legendre_nodes_1D(legendre_nodes_1D,num_gauss_points);
  elements::legendre_weights_1D(legendre_weights_1D,num_gauss_points);

  //initialize fluxes to 0
  //local variable for host view in the dual view
  for(int init = 0; init < map->getLocalNumElements(); init++)
    for(int iflux = 0; iflux < Brows; iflux++)
      node_heat_fluxes(init,iflux) = 0;

  for(int init = 0; init < all_node_map->getLocalNumElements(); init++)
    for(int iflux = 0; iflux < Brows; iflux++)
      all_node_heat_fluxes(init,iflux) = 0;
  

  real_t current_density = 1;

  //compute nodal fluxes as the result of minimizing the L2 error between the flux field and the nodal interpolant over each element
  //assign the maximum or average nodal flux computed from connected elements to nodes sharing elements (safer to know the largest positive or negative values)
  
  for(int ielem = 0; ielem < rnum_elem; ielem++){
    nodes_per_elem = elem->num_basis();
    //B matrix initialization
    for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < nodes_per_elem; icol++){
        B_matrix(irow,icol) = 0;
      }

    //initialize C matrix
    for(int irow = 0; irow < Brows; irow++)
      for(int icol = 0; icol < Brows; icol++)
        C_matrix(irow,icol) = 0;

    //B matrix initialization
    for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < nodes_per_elem; icol++){
        CB_matrix(irow,icol) = 0;
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

    //acquire set of nodes and nodal temperatures for this local element
    for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
      local_node_id = all_node_map->getLocalElement(nodes_in_elem(ielem, node_loop));
      nodal_positions(node_loop,0) = all_node_coords(local_node_id,0);
      nodal_positions(node_loop,1) = all_node_coords(local_node_id,1);
      nodal_positions(node_loop,2) = all_node_coords(local_node_id,2);
      current_nodal_temperatures(node_loop) = all_node_temperatures(local_node_id,0);
      if(nodal_density_flag) nodal_density(node_loop) = all_node_densities(local_node_id,0);
    }

    //debug print of current_nodal_temperatures
    /*
    std::cout << " ------------nodal temperatures for Element "<< ielem + 1 <<"--------------"<<std::endl;
    std::cout << " { ";
    for (int idof = 0; idof < num_dim*nodes_per_elem; idof++){
      std::cout << idof + 1 << " = " << current_nodal_temperatures(idof) << " , " ;
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
    for(int ishape=0; ishape < nodes_per_elem; ishape++){
      B_matrix_contribution(0,ishape) = (basis_derivative_s1(ishape)*(JT_row2(1)*JT_row3(2)-JT_row3(1)*JT_row2(2))-
                                         basis_derivative_s2(ishape)*(JT_row1(1)*JT_row3(2)-JT_row3(1)*JT_row1(2))+
                                         basis_derivative_s3(ishape)*(JT_row1(1)*JT_row2(2)-JT_row2(1)*JT_row1(2)));
      B_matrix_contribution(1,ishape) = (-basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(2)-JT_row3(0)*JT_row2(2))+
                                         basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(2)-JT_row3(0)*JT_row1(2))-
                                         basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(2)-JT_row2(0)*JT_row1(2)));
      if(num_dim==3){
        B_matrix_contribution(2,ishape) = (basis_derivative_s1(ishape)*(JT_row2(0)*JT_row3(1)-JT_row3(0)*JT_row2(1))-
                                           basis_derivative_s2(ishape)*(JT_row1(0)*JT_row3(1)-JT_row3(0)*JT_row1(1))+
                                           basis_derivative_s3(ishape)*(JT_row1(0)*JT_row2(1)-JT_row2(0)*JT_row1(1)));
      }
    }

    //compute Isotropic Conductivity (C) matrix
    if(num_dim==2){
      C_matrix(0,0) = -1;
      C_matrix(1,1) = -1;
    }
    if(num_dim==3){
      C_matrix(0,0) = -1;
      C_matrix(1,1) = -1;
      C_matrix(2,2) = -1;
    }
    
    //multiply by temperature vector to get temperature gradient at this quadrature point
    //division by J ommited since it is multiplied out later
    for(int irow=0; irow < Brows; irow++){
      quad_temperature_gradient(irow) = 0;
      for(int icol=0; icol < nodes_per_elem; icol++){
        quad_temperature_gradient(irow) += B_matrix_contribution(irow,icol)*current_nodal_temperatures(icol);
      }
    }

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

    //look up element material properties at this point as a function of density
    Element_Material_Properties(ielem, Element_Conductivity, current_density);
    
    //compute heat flux with Fourier's law
    for(int irow=0; irow < Brows; irow++){
      quad_heat_flux(irow) = 0;
      for(int icol=0; icol < Brows; icol++){
        quad_heat_flux(irow) += Element_Conductivity*C_matrix(irow,icol)*quad_temperature_gradient(icol);
      }
    }

    //compute contribution to RHS projection vector
    for(int irow=0; irow < Brows; irow++)
      for(int icol=0; icol < nodes_per_elem; icol++){
        projection_vector(irow,icol) += weight_multiply*quad_heat_flux(irow)*basis_values(icol);
      }

    //compute contribution to projection matrix (only upper part is set)
    for(int irow=0; irow < nodes_per_elem; irow++)
      for(int icol=0; icol < nodes_per_elem; icol++){
        //if(irow<=icol)
        projection_matrix(irow,icol) += weight_multiply*basis_values(irow)*basis_values(icol)*Jacobian;
      }
    
    }
    
    //construct matrix and vector wrappers for dense solver
    //projection_matrix_pass = Teuchos::rcp( new Teuchos::SerialSymDenseMatrix<LO,real_t>(Teuchos::View, true, projection_matrix.pointer(), nodes_per_elem, nodes_per_elem));

    //debug print of matrix
    //projection_matrix_pass->print(std::cout);

    heat_flux_vector_pass = Teuchos::rcp( new Teuchos::SerialDenseVector<LO,real_t>(Teuchos::View, heat_flux_vector.pointer(), nodes_per_elem));
    //loop through flux components and solve for nodal values of that component
    for(int iflux = 0; iflux < Brows; iflux++){
      //check if projection vector is zero due to zero fluxes
      zero_flux_flag = 1;
      for(int icol=0; icol < nodes_per_elem; icol++){
        if(fabs(projection_vector(iflux,icol))>FLUX_EPSILON) zero_flux_flag = 0;
      }
      if(!zero_flux_flag){
        projection_vector_pass = Teuchos::rcp( new Teuchos::SerialDenseVector<LO,real_t>(Teuchos::View, &projection_vector(iflux,0), nodes_per_elem));
        projection_matrix_pass = Teuchos::rcp( new Teuchos::SerialDenseMatrix<LO,real_t>(Teuchos::Copy, projection_matrix.pointer(), nodes_per_elem, nodes_per_elem, nodes_per_elem));
        //debug print of vectors
        //projection_vector_pass->print(std::cout);
        //projection_matrix_pass->print(std::cout);
        projection_solver.setMatrix(projection_matrix_pass);
        projection_solver.setVectors(heat_flux_vector_pass, projection_vector_pass);
        projection_solver.factorWithEquilibration(true);
        solve_flag = projection_solver.solve();
        
        //debug print
        //std::cout << "HEAT FLUX COMPONENT ON NODES " << iflux + 1 << std::endl;
        //heat_flux_vector_pass->print(std::cout);
        if(solve_flag) std::cout << "Projection Solve Failed With: " << solve_flag << std::endl;

        //contribute nodal flux for this element to corresponding global nodes
        //replace if abs greater than abs of currently assigned flux; accumulate average if flagged for node average
        for(int node_loop=0; node_loop < nodes_per_elem; node_loop++){
          current_global_index = nodes_in_elem(ielem, node_loop);
          local_node_id = all_node_map->getLocalElement(current_global_index);
          current_heat_flux = (*heat_flux_vector_pass)(node_loop);
          if(flux_max_flag){
            if(fabs(current_heat_flux) > all_node_heat_fluxes(local_node_id, iflux)){
              all_node_heat_fluxes(local_node_id, iflux) = current_heat_flux;

              if(map->isNodeGlobalElement(current_global_index)){
                local_node_id = map->getLocalElement(current_global_index);
                node_heat_fluxes(local_node_id, iflux) = current_heat_flux;
              }
            }
          }
          else{
            //debug print
            if(map->isNodeGlobalElement(current_global_index)){
              local_node_id = map->getLocalElement(current_global_index);
              node_heat_fluxes(local_node_id, iflux) += current_heat_flux/(double)node_nconn(local_node_id,0);
            }
          }
        }
      }
    }
    
  }
    
}

/* ----------------------------------------------------------------------
   Solve the FEA linear system
------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::linear_solver_parameters(){
  if(module_params.direct_solver_flag){
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
  int num_dim = simparam.num_dims;
  int nodes_per_elem = max_nodes_per_element;
  int local_node_index, current_row, current_column;
  int max_stride = 0;
  size_t access_index, row_access_index, row_counter;
  GO global_index, global_dof_index;
  LO local_dof_index;

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
  Original_RHS_Entries = CArrayKokkos<real_t, array_layout, device_type, memory_traits>(local_nboundaries);
  real_t diagonal_bc_scaling = Conductivity_Matrix(0,0);

  //alter rows of RHS to be the boundary condition value on that node
  //first pass counts strides for storage
  row_counter = 0;
  {//view scope
    host_vec_array Nodal_RHS = Global_Nodal_RHS->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    for(LO i=0; i < nlocal_nodes; i++){
      if((Node_DOF_Boundary_Condition_Type(i)==TEMPERATURE_CONDITION)){
        Original_RHS_Entries(row_counter) = Nodal_RHS(i,0);
        row_counter++;
        Nodal_RHS(i,0) = Node_Temperature_Boundary_Conditions(i)*diagonal_bc_scaling;
      }
    }//row for
  }//end view scope
  //change entries of Conductivity matrix corresponding to BCs to 0s (off diagonal elements) and 1 (diagonal elements)
  //storage for original Conductivity matrix values
  Original_Conductivity_Entries_Strides = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(nlocal_nodes);

  //debug print of A matrix before applying BCS
  //*fos << "Reduced Conductivity Matrix :" << std::endl;
  //Global_Conductivity_Matrix->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  //Tpetra::MatrixMarket::Writer<MAT> market_writer();
  //Tpetra::MatrixMarket::Writer<MAT>::writeSparseFile("A_matrix.txt", *Global_Conductivity_Matrix, "A_matrix", "Stores conductivity matrix values");

  //first pass counts strides for storage
  if(!matrix_bc_reduced){
  for(LO i=0; i < nlocal_nodes; i++){
    Original_Conductivity_Entries_Strides(i) = 0;
    if((Node_DOF_Boundary_Condition_Type(i)==TEMPERATURE_CONDITION)){
      Original_Conductivity_Entries_Strides(i) = Conductivity_Matrix_Strides(i);
    }
    else{
      for(LO j = 0; j < Conductivity_Matrix_Strides(i); j++){
        global_dof_index = DOF_Graph_Matrix(i,j);
        local_dof_index = all_dof_map->getLocalElement(global_dof_index);
        if((Node_DOF_Boundary_Condition_Type(local_dof_index)==TEMPERATURE_CONDITION)){
          Original_Conductivity_Entries_Strides(i)++;
        }
      }//stride for
    }
  }//row for
  
  //assign old Conductivity matrix entries
  LO stride_index;
  Original_Conductivity_Entries = RaggedRightArrayKokkos<real_t, array_layout, device_type, memory_traits>(Original_Conductivity_Entries_Strides);
  Original_Conductivity_Entry_Indices = RaggedRightArrayKokkos<LO, array_layout, device_type, memory_traits>(Original_Conductivity_Entries_Strides);
  for(LO i=0; i < nlocal_nodes; i++){
    if((Node_DOF_Boundary_Condition_Type(i)==TEMPERATURE_CONDITION)){
      for(LO j = 0; j < Conductivity_Matrix_Strides(i); j++){
        global_dof_index = DOF_Graph_Matrix(i,j);
        local_dof_index = all_dof_map->getLocalElement(global_dof_index);
        Original_Conductivity_Entries(i,j) = Conductivity_Matrix(i,j);
        Original_Conductivity_Entry_Indices(i,j) = j;
        if(local_dof_index == i){
          Conductivity_Matrix(i,j) = diagonal_bc_scaling;
        }
        else{     
          Conductivity_Matrix(i,j) = 0;
        }
      }//stride for
    }
    else{
      stride_index = 0;
      for(LO j = 0; j < Conductivity_Matrix_Strides(i); j++){
        global_dof_index = DOF_Graph_Matrix(i,j);
        local_dof_index = all_dof_map->getLocalElement(global_dof_index);
        if((Node_DOF_Boundary_Condition_Type(local_dof_index)==TEMPERATURE_CONDITION)){
          Original_Conductivity_Entries(i,stride_index) = Conductivity_Matrix(i,j);
          Original_Conductivity_Entry_Indices(i,stride_index) = j;   
          Conductivity_Matrix(i,j) = 0;
          stride_index++;
        }
      }
    }
  }//row for

  matrix_bc_reduced = true;
  }
  //This completes the setup for A matrix of the linear system

  //debug print of A matrix
  /*
  *fos << "Reduced Conductivity Matrix :" << std::endl;
  Global_Conductivity_Matrix->describe(*fos,Teuchos::VERB_EXTREME);
  *fos << std::endl;
  */

  //debug print
  //if(update_count==42){
    //Tpetra::MatrixMarket::Writer<MAT> market_writer();
    //Tpetra::MatrixMarket::Writer<MAT>::writeSparseFile("A_matrix.txt", *Global_Conductivity_Matrix, "A_matrix", "Stores conductivity matrix values");
  //}
  using impl_scalar_type =
    typename Kokkos::Details::ArithTraits<real_t>::val_type;
  using mag_type = typename Kokkos::ArithTraits<impl_scalar_type>::mag_type;
    // Instead of checking each time for rank, create a rank 0 stream

  xB = Teuchos::rcp(new Xpetra::TpetraMultiVector<real_t,LO,GO,node_type>(Global_Nodal_RHS));
  X = node_temperatures_distributed;
  xX = Teuchos::rcp(new Xpetra::TpetraMultiVector<real_t,LO,GO,node_type>(X));

  Teuchos::RCP<MV> tcoordinates;
  tcoordinates = node_coords_distributed;
  Teuchos::RCP<Xpetra::MultiVector<real_t,LO,GO,node_type>> coordinates = Teuchos::rcp(new Xpetra::TpetraMultiVector<real_t,LO,GO,node_type>(tcoordinates));
    
  //nullspace vector
  Teuchos::RCP<MV> tnullspace = Teuchos::rcp(new MV(local_dof_map, 1));
  tnullspace->putScalar(1);
  Teuchos::RCP<Xpetra::MultiVector<real_t,LO,GO,node_type>> nullspace = Teuchos::rcp(new Xpetra::TpetraMultiVector<real_t,LO,GO,node_type>(tnullspace));

  Teuchos::RCP<Xpetra::MultiVector<real_t,LO,GO,node_type>> material = Teuchos::null;
  Teuchos::RCP<Xpetra::CrsMatrix<real_t,LO,GO,node_type>> xcrs_A = Teuchos::rcp(new Xpetra::TpetraCrsMatrix<real_t,LO,GO,node_type>(Global_Conductivity_Matrix));
  xA = Teuchos::rcp(new Xpetra::CrsMatrixWrap<real_t,LO,GO,node_type>(xcrs_A));
  //xA->SetFixedBlockSize(num_dim);
   
  //randomize initial vector
  xX->setSeed(100);
  xX->randomize();

  //initialize BC components
  /*
  host_vec_array X_view = X->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
  for(LO i=0; i < nlocal_nodes; i++){
    if((Node_DOF_Boundary_Condition_Type(i)==TEMPERATURE_CONDITION)){
      X_view(i,0) = Node_DOF_Displacement_Boundary_Conditions(i);
    }
  }//row for
  */
  if(module_params.equilibrate_matrix_flag){
    Implicit_Solver_Pointer_->equilibrateMatrix(xA,"diag");
    Implicit_Solver_Pointer_->preScaleRightHandSides(*Global_Nodal_RHS,"diag");
    Implicit_Solver_Pointer_->preScaleInitialGuesses(*X,"diag");
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
    
  //xA->describe(*fos,Teuchos::VERB_EXTREME);
  //debug print
  //Tpetra::MatrixMarket::Writer<MAT> market_writer();
  //Tpetra::MatrixMarket::Writer<MAT>::writeSparseFile("A_matrix_reduced.txt", *Global_Conductivity_Matrix, "A_matrix_reduced", "Stores Conductivity matrix values");  
  //Xpetra::IO<real_t,LO,GO,node_type>WriteLocal("A_matrixlocal.txt", *xA);
  comm->barrier();
  //PreconditionerSetup(A,coordinates,nullspace,material,paramList,false,false,useML,0,H,Prec);
  //xA->describe(*fos,Teuchos::VERB_EXTREME);
  Teuchos::RCP<Tpetra::Vector<real_t,LO,GO,node_type>> tdiagonal = Teuchos::rcp(new Tpetra::Vector<real_t,LO,GO,node_type>(local_dof_map));
  //Teuchos::RCP<Xpetra::Vector<real_t,LO,GO,node_type>> diagonal = Teuchos::rcp(new Xpetra::Vector<real_t,LO,GO,node_type>(tdiagonal));
  //Global_Conductivity_Matrix->getLocalDiagCopy(*tdiagonal);
  //tdiagonal->describe(*fos,Teuchos::VERB_EXTREME);
  if(Hierarchy_Constructed){
    ReuseXpetraPreconditioner(xA, H);
  }
  else{
    PreconditionerSetup(xA,coordinates,nullspace,material,*Linear_Solve_Params,false,false,false,0,H,Prec);
    //Hierarchy_Constructed = true;
  }
  comm->barrier();
    
  //H->Write(-1, -1);
  //H->describe(*fos,Teuchos::VERB_EXTREME);
    
  // =========================================================================
  // System solution (Ax = b)
  // =========================================================================

  real_t current_cpu_time = Implicit_Solver_Pointer_->CPU_Time();
  SystemSolve(xA,xX,xB,H,Prec,*fos,solveType,belosType,false,false,false,cacheSize,0,true,true,num_iter,solve_tol);
  linear_solve_time += Implicit_Solver_Pointer_->CPU_Time() - current_cpu_time;
  comm->barrier();

  if(module_params.equilibrate_matrix_flag){
    Implicit_Solver_Pointer_->postScaleSolutionVectors(*X,"diag");
    Implicit_Solver_Pointer_->postScaleSolutionVectors(*Global_Nodal_RHS,"diag");
  }

  if(module_params.multigrid_timers){
    Teuchos::RCP<Teuchos::ParameterList> reportParams = rcp(new Teuchos::ParameterList);
    reportParams->set("How to merge timer sets",   "Union");
    reportParams->set("alwaysWriteLocal",          false);
    reportParams->set("writeGlobalStats",          true);
    reportParams->set("writeZeroTimers",           false);
    std::ios_base::fmtflags ff(fos->flags());
    *fos << std::fixed;
    Teuchos::TimeMonitor::report(comm.ptr(), *fos, "", reportParams);
    *fos << std::setiosflags(ff);
    //xA->describe(*fos,Teuchos::VERB_EXTREME);
  }
  //return !EXIT_SUCCESS;
  //timing statistics for LU solver
  //solver->printTiming(*fos);
  
  //import for displacement of ghosts
  Tpetra::Import<LO, GO> ghost_temperature_importer(local_dof_map, all_dof_map);

  //comms to get temperatures on all node map
  all_node_temperatures_distributed->doImport(*node_temperatures_distributed, *importer, Tpetra::INSERT);

  //reinsert global conductivity values corresponding to BC indices to facilitate heat potential calculation
  if(matrix_bc_reduced){
    for(LO i = 0; i < nlocal_nodes; i++){
      for(LO j = 0; j < Original_Conductivity_Entries_Strides(i); j++){
        access_index = Original_Conductivity_Entry_Indices(i,j);
        Conductivity_Matrix(i,access_index) = Original_Conductivity_Entries(i,j);
      }
    }//row for
    matrix_bc_reduced = false;
  }

  //compute nodal heat vector (used by other functions such as TO) due to inputs and constraints
  Global_Conductivity_Matrix->apply(*node_temperatures_distributed,*Global_Nodal_Heat);

  //if(myrank==0)
  //*fos << "All temperatures :" << std::endl;
  //all_node_temperatures_distributed->describe(*fos,Teuchos::VERB_EXTREME);
  //*fos << std::endl;
  
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------------------------
   Communicate ghosts using the current optimization design data
---------------------------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::comm_variables(Teuchos::RCP<const MV> zp){
  
  if(simparam.topology_optimization_on)
    comm_densities(zp);
}

/* -------------------------------------------------------------------------------------------
   update nodal displacement information in accordance with current optimization vector
---------------------------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::update_linear_solve(Teuchos::RCP<const MV> zp, int compute_step){
  
  if(compute_step!=last_compute_step){
    //set density vector to the current value chosen by the optimizer
    test_node_densities_distributed = zp;

    assemble_matrix();

    if(body_term_flag||nonzero_bc_flag)
      assemble_vector();
  
    //solve for new nodal temperatures
    int solver_exit = solve();
    if(solver_exit != EXIT_SUCCESS){
      std::cout << "Linear Solver Error" << std::endl <<std::flush;
      return;
    }
  
    update_count++;
    last_compute_step = compute_step;
  }
}

/* -------------------------------------------------------------------------------------------
   enforce constraints on nodes due to BCS
---------------------------------------------------------------------------------------------- */

void FEA_Module_Heat_Conduction::node_density_constraints(host_vec_array node_densities_lower_bound){
  LO local_node_index;
  int num_dim = simparam.num_dims;

  if(simparam.optimization_options.thick_condition_boundary){
    for(int i = 0; i < nlocal_nodes; i++){
      if(Node_DOF_Boundary_Condition_Type(i) == TEMPERATURE_CONDITION){
        for(int j = 0; j < Graph_Matrix_Strides(i); j++){
          if(map->isNodeGlobalElement(Graph_Matrix(i,j))){
            local_node_index = map->getLocalElement(Graph_Matrix(i,j));
            node_densities_lower_bound(local_node_index,0) = 1;
          }
        }
      }
    }
  }
  else{
    for(int i = 0; i < nlocal_nodes; i++){
      if(Node_DOF_Boundary_Condition_Type(i) == TEMPERATURE_CONDITION){
        node_densities_lower_bound(i,0) = 1;
      }
    }
  }
}
