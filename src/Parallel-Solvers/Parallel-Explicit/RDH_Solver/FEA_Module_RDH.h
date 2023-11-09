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


#ifndef FEA_MODULE_RDH_H
#define FEA_MODULE_RDH_H

#include "mesh.h"
#include "ref_elem.h"
#include "state.h"
#include "matar.h"
#include "elements.h"
#include "node_combination.h"
#include "Solver.h"
#include "FEA_Module.h"
#include "Simulation_Parameters/FEA_Module/RDH_Parameters.h"
#include "Simulation_Parameters/Simulation_Parameters_Explicit.h"
//#include "Simulation_Parameters_Elasticity.h"
//#include "Simulation_Parameters_Dynamic_Optimization.h"
#include "material_models.h"

class Explicit_Solver;

class FEA_Module_RDH: public FEA_Module{

public:
  
  FEA_Module_RDH( RDH_Parameters& params, Solver *Solver_Pointer, std::shared_ptr<mesh_t> mesh_in, const int my_fea_module_index = 0);
  ~FEA_Module_RDH();
  
  //initialize data for boundaries of the model and storage for boundary conditions and applied loads
  void rdh_interface_setup(node_t &node, elem_t &elem, mesh_t &mesh, corner_t &corner);

  void setup(mesh_t &mesh);

  void cleanup_material_models();

  int solve();

  void module_cleanup();

  void rdh_solve();

  Simulation_Parameters_Explicit simparam;
  RDH_Parameters module_params;
  Explicit_Solver *Explicit_Solver_Pointer_;

  //mesh_t mesh;
  fe_ref_elem_t ref_elem;
  //std::shared_ptr<fe_ref_elem_t> ref_elem;
  std::shared_ptr<mesh_t> mesh;

  //shallow copies of mesh class views
  size_t num_nodes_in_elem;
  // corner ids in node
  RaggedRightArrayKokkos <size_t> corners_in_node;
  CArrayKokkos <size_t> num_corners_in_node;
    
  // elem ids in node
  RaggedRightArrayKokkos <size_t> elems_in_node;
    
  // node ids in node
  RaggedRightArrayKokkos <size_t> nodes_in_node;
  CArrayKokkos <size_t> num_nodes_in_node;
    
  // node ids in elem
  DCArrayKokkos <size_t> nodes_in_elem;
    
  // corner ids in elem
  CArrayKokkos <size_t> corners_in_elem;
    
  // elem ids in elem
  RaggedRightArrayKokkos <size_t> elems_in_elem;
  CArrayKokkos <size_t> num_elems_in_elem;
    
  // patch ids in elem
  CArrayKokkos <size_t> patches_in_elem;
    
  // node ids in a patch
  CArrayKokkos <size_t> nodes_in_patch;
    
  // element ids in a patch
  CArrayKokkos <size_t> elems_in_patch;
  
  // bdy nodes
  CArrayKokkos <size_t> bdy_nodes;
  
  //Topology optimization filter variable
  DCArrayKokkos<double> relative_element_densities;
  
  //Local FEA data
  host_elem_conn_array interface_nodes_in_elem; //host view of element connectivity to nodes

  //Global FEA data
  Teuchos::RCP<MV> node_velocities_distributed;
  Teuchos::RCP<MV> initial_node_coords_distributed;
  Teuchos::RCP<MV> all_initial_node_coords_distributed;
  Teuchos::RCP<MV> initial_node_velocities_distributed;
  Teuchos::RCP<MV> all_node_velocities_distributed;
  Teuchos::RCP<MV> all_cached_node_velocities_distributed;
  Teuchos::RCP<MV> node_masses_distributed;
  Teuchos::RCP<MV> ghost_node_masses_distributed;
  Teuchos::RCP<MV> adjoint_vector_distributed;
  Teuchos::RCP<MV> phi_adjoint_vector_distributed;
  Teuchos::RCP<MV> psi_adjoint_vector_distributed;
  Teuchos::RCP<std::vector<Teuchos::RCP<MV>>> forward_solve_velocity_data;
  Teuchos::RCP<std::vector<Teuchos::RCP<MV>>> forward_solve_coordinate_data;
  Teuchos::RCP<std::vector<Teuchos::RCP<MV>>> adjoint_vector_data;
  Teuchos::RCP<std::vector<Teuchos::RCP<MV>>> phi_adjoint_vector_data;
  Teuchos::RCP<std::vector<Teuchos::RCP<MV>>> psi_adjoint_vector_data;
  Teuchos::RCP<MV> force_gradient_design;
  Teuchos::RCP<MV> force_gradient_position;
  Teuchos::RCP<MV> force_gradient_velocity;

  //Local FEA data
  DCArrayKokkos<size_t, array_layout, device_type, memory_traits> Global_Gradient_Matrix_Assembly_Map;
  RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits> Graph_Matrix; //stores global indices
  RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits> DOF_Graph_Matrix; //stores global indices
  RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout> Force_Gradient_Positions;
  RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout> Force_Gradient_Velocities;
  RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout> Force_Gradient_Energies;
  RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout> Power_Gradient_Positions;
  RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout> Power_Gradient_Velocities;
  CArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits> Power_Gradient_Energies;
  DCArrayKokkos<size_t, array_layout, device_type, memory_traits> Gradient_Matrix_Strides;
  DCArrayKokkos<size_t, array_layout, device_type, memory_traits> DOF_to_Elem_Matrix_Strides;
  DCArrayKokkos<size_t, array_layout, device_type, memory_traits> Elem_to_Elem_Matrix_Strides;
  DCArrayKokkos<size_t, array_layout, device_type, memory_traits> Graph_Matrix_Strides;
  RaggedRightArrayKokkos<real_t, array_layout, device_type, memory_traits> Original_Gradient_Entries;
  RaggedRightArrayKokkos<LO, array_layout, device_type, memory_traits> Original_Gradient_Entry_Indices;
  DCArrayKokkos<size_t, array_layout, device_type, memory_traits> Original_Gradient_Entries_Strides;

  //distributed matrices
  Teuchos::RCP<MAT> distributed_force_gradient_positions;
  Teuchos::RCP<MAT> distributed_force_gradient_velocities;

  std::vector<real_t> time_data;
  int max_time_steps, last_time_step;

  // ---------------------------------------------------------------------
  //    state data type declarations (must stay in scope for output after run)
  // ---------------------------------------------------------------------
  mesh_t mesh_interface;
  node_t  node_interface;
  elem_t  elem_interface;
  corner_t  corner_interface;

  //Dual View wrappers
  // Dual Views of the individual node struct variables
  DViewCArrayKokkos <double> node_coords;
  DViewCArrayKokkos <double> node_vel;
  DViewCArrayKokkos <double> node_mass;
             
  // Dual Views of the individual elem struct variables
  DViewCArrayKokkos <double> elem_den;
  DViewCArrayKokkos <double> elem_pres;
  DViewCArrayKokkos <double> elem_stress; // always 3D even in 2D-RZ
  DViewCArrayKokkos <double> elem_sspd;
  DViewCArrayKokkos <double> elem_sie;
  DViewCArrayKokkos <double> elem_vol;
  DViewCArrayKokkos <double> elem_div;    
  DViewCArrayKokkos <double> elem_mass;
  DViewCArrayKokkos <size_t> elem_mat_id;

  // Element velocity gradient 
  DCArrayKokkos <double> elem_vel_grad;

  // for storing global variables used in user material model
  DCArrayKokkos <double> global_vars;

  //elem_user_output_vars allow users to output variables of interest per element
  DCArrayKokkos <double> elem_user_output_vars;

  //material models
  DCArrayKokkos <eos_t> elem_eos;
  DCArrayKokkos <strength_t> elem_strength;

  // Dual Views of the corner struct variables
  DViewCArrayKokkos <double> corner_force;
  DViewCArrayKokkos <double> corner_mass;
  
  //Boundary Conditions Data
  DCArrayKokkos<size_t> Local_Index_Boundary_Patches;
  //CArray <Nodal_Combination> Patch_Nodes;
  enum bc_type {NONE, POINT_LOADING_CONDITION, LINE_LOADING_CONDITION, SURFACE_LOADING_CONDITION};
  
  //Boundary Conditions Data
  int max_boundary_sets;

  //output dof data
  //Global arrays with collected data used to print
  int output_velocity_index, output_strain_index, output_stress_index;
  
  //parameters
  double time_value, time_final, dt, dt_max, dt_min, dt_cfl, graphics_time, graphics_dt_ival;
  size_t graphics_cyc_ival, cycle_stop, rk_num_stages, graphics_id;
  double fuzz, tiny, small;
  CArray <double> graphics_times;

  //optimization flags
  bool kinetic_energy_objective;
  
  bool have_loading_conditions;
};

#endif
