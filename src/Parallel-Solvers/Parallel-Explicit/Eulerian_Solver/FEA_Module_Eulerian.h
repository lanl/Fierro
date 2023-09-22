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

#ifndef FEA_MODULE_EULERIAN_H
#define FEA_MODULE_EULERIAN_H

#include "mesh.h"
#include "state.h"
#include "matar.h"
#include "elements.h"
#include "node_combination.h"
#include "Solver.h"
#include "FEA_Module.h"
#include "Simulation_Parameters_Eulerian.h"
#include "Simulation_Parameters_Elasticity.h"
#include "Simulation_Parameters_Dynamic_Optimization.h"

class Explicit_Solver_Eulerian;

class FEA_Module_Eulerian: public FEA_Module{

public:
  
  FEA_Module_Eulerian(Solver *Solver_Pointer, mesh_t& mesh, const int my_fea_module_index = 0);
  ~FEA_Module_Eulerian();

  void setup();

  void euler_solve();

  void comm_node_masses();

  void comm_adjoint_vectors(int cycle);

  void comm_variables(Teuchos::RCP<const MV> zp);

  void read_conditions_ansys_dat(std::ifstream *in, std::streampos before_condition_header);

  //interfaces between user input and creating data structures for bcs
  void generate_bcs();

  void Displacement_Boundary_Conditions();

  void init_output();

  void compute_output();
  
  void sort_output(Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > sorted_map);

  void collect_output(Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > global_reduce_map);

  void init_boundaries();

  //initializes memory for arrays used in the global stiffness matrix assembly
  void init_boundary_sets(int num_boundary_sets);

  void grow_boundary_sets(int num_boundary_sets);

  void rk_init(DViewCArrayKokkos <double> &node_coords,
               DViewCArrayKokkos <double> &node_vel,
               DViewCArrayKokkos <double> &elem_sie,
               DViewCArrayKokkos <double> &elem_stress,
               const size_t num_elems,
               const size_t num_nodes);

  void get_timestep(mesh_t &mesh,
                    DViewCArrayKokkos <double> &node_coords,
                    DViewCArrayKokkos <double> &node_vel,
                    DViewCArrayKokkos <double> &elem_sspd,
                    DViewCArrayKokkos <double> &elem_vol);


  void get_timestep2D(mesh_t &mesh,
                      DViewCArrayKokkos <double> &node_coords,
                      DViewCArrayKokkos <double> &node_vel,
                      DViewCArrayKokkos <double> &elem_sspd,
                      DViewCArrayKokkos <double> &elem_vol);

  void write_outputs (const mesh_t &mesh,
                      DViewCArrayKokkos <double> &node_coords,
                      DViewCArrayKokkos <double> &node_vel,
                      DViewCArrayKokkos <double> &node_mass,
                      DViewCArrayKokkos <double> &elem_den,
                      DViewCArrayKokkos <double> &elem_pres,
                      DViewCArrayKokkos <double> &elem_stress,
                      DViewCArrayKokkos <double> &elem_sspd,
                      DViewCArrayKokkos <double> &elem_sie,
                      DViewCArrayKokkos <double> &elem_vol,
                      DViewCArrayKokkos <double> &elem_mass,
                      DViewCArrayKokkos <size_t> &elem_mat_id);

  void node_density_constraints(host_vec_array node_densities_lower_bound);

  void example_function(double rk_alpha,
                         const size_t num_nodes,
                         DViewCArrayKokkos <double> &node_coords,
                         const DViewCArrayKokkos <double> &node_vel);

  KOKKOS_FUNCTION
  void example_device_function(const ViewCArrayKokkos <double> &B_matrix,
                 const size_t elem_gid,
                 const DViewCArrayKokkos <double> &node_coords,
                 const ViewCArrayKokkos <size_t>  &elem_node_gids,
                 const size_t rk_level) const;
  
  
  bool nodal_density_flag;
  real_t penalty_power;
  Teuchos::RCP<MAT> Global_Stiffness_Matrix;
  RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout> Stiffness_Matrix;
  DCArrayKokkos<size_t, array_layout, device_type, memory_traits> Stiffness_Matrix_Strides;
  DCArrayKokkos<size_t, array_layout, device_type, memory_traits> Global_Stiffness_Matrix_Assembly_Map;
  //end elastic TO data
  
  Simulation_Parameters_Eulerian simparam;
  Simulation_Parameters_Elasticity simparam_elasticity;
  Simulation_Parameters_Dynamic_Optimization simparam_dynamic_opt;
  Explicit_Solver_Eulerian *Explicit_Solver_Pointer_;

  elements::ref_element  *ref_elem;
  
  mesh_t& mesh;
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

  // patch ids in bdy set
  size_t num_bdy_sets;
  DynamicRaggedRightArrayKokkos <size_t> bdy_patches_in_set;
  
  // bdy nodes
  CArrayKokkos <size_t> bdy_nodes;

  // node ids in bdy_patch set
  RaggedRightArrayKokkos <size_t> bdy_nodes_in_set;
  DCArrayKokkos <size_t> num_bdy_nodes_in_set;
  
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
  Teuchos::RCP<std::vector<Teuchos::RCP<MV>>> forward_solve_velocity_data;
  Teuchos::RCP<std::vector<Teuchos::RCP<MV>>> forward_solve_coordinate_data;
  Teuchos::RCP<std::vector<Teuchos::RCP<MV>>> adjoint_vector_data;
  Teuchos::RCP<std::vector<Teuchos::RCP<MV>>> phi_adjoint_vector_data;
  Teuchos::RCP<MV> force_gradient_design;
  Teuchos::RCP<MV> force_gradient_position;
  Teuchos::RCP<MV> force_gradient_velocity;

  //Local FEA data
  DCArrayKokkos<size_t, array_layout, device_type, memory_traits> Global_Gradient_Matrix_Assembly_Map;
  RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits> Graph_Matrix; //stores global indices
  RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits> DOF_Graph_Matrix; //stores global indices
  RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout> Force_Gradient_Positions;
  RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout> Force_Gradient_Velocities;
  DCArrayKokkos<size_t, array_layout, device_type, memory_traits> Gradient_Matrix_Strides;
  DCArrayKokkos<size_t, array_layout, device_type, memory_traits> Graph_Matrix_Strides;
  RaggedRightArrayKokkos<real_t, array_layout, device_type, memory_traits> Original_Gradient_Entries;
  RaggedRightArrayKokkos<LO, array_layout, device_type, memory_traits> Original_Gradient_Entry_Indices;
  DCArrayKokkos<size_t, array_layout, device_type, memory_traits> Original_Gradient_Entries_Strides;

  //distributed matrices
  Teuchos::RCP<MAT> distributed_force_gradient_positions;
  Teuchos::RCP<MAT> distributed_force_gradient_velocities;

  std::vector<real_t> time_data;
  int max_time_steps, last_time_step;

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
  
  //file parameters
  DCArrayKokkos <size_t>read_from_file;

  //parameters
  double time_value, time_final, dt, dt_max, dt_min, dt_cfl, graphics_time, graphics_dt_ival;
  size_t graphics_cyc_ival, cycle_stop, rk_num_stages, graphics_id;
  double fuzz, tiny, small;
  CArray <double> graphics_times;

  //optimization flags
  bool kinetic_energy_objective;
};

#endif // end HEADER_H
