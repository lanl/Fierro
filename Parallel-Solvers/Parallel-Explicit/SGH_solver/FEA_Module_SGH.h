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

#ifndef FEA_MODULE_SGH_H
#define FEA_MODULE_SGH_H

#include "mesh.h"
#include "state.h"
#include "matar.h"
#include "elements.h"
#include "node_combination.h"
#include "Solver.h"
#include "FEA_Module.h"

class Explicit_Solver_SGH;
class Simulation_Parameters_SGH;
class Simulation_Parameters_Topology_Optimization;

class FEA_Module_SGH: public FEA_Module{

public:
  
  FEA_Module_SGH(Solver *Solver_Pointer);
  ~FEA_Module_SGH();
  
  //initialize data for boundaries of the model and storage for boundary conditions and applied loads
  void sgh_interface_setup(mesh_t &mesh, node_t &node, elem_t &elem, corner_t &corner);

  void setup(mesh_t &mesh,
             const DViewCArrayKokkos <double> &node_coords,
             DViewCArrayKokkos <double> &node_vel,
             DViewCArrayKokkos <double> &node_mass,
             const DViewCArrayKokkos <double> &elem_den,
             const DViewCArrayKokkos <double> &elem_pres,
             const DViewCArrayKokkos <double> &elem_stress,
             const DViewCArrayKokkos <double> &elem_sspd,
             const DViewCArrayKokkos <double> &elem_sie,
             const DViewCArrayKokkos <double> &elem_vol,
             const DViewCArrayKokkos <double> &elem_mass,
             const DViewCArrayKokkos <size_t> &elem_mat_id,
             const DViewCArrayKokkos <double> &elem_statev,
             const DViewCArrayKokkos <double> &corner_mass);

  void sgh_solve(mesh_t &mesh,
                 DViewCArrayKokkos <double> &node_coords,
                 DViewCArrayKokkos <double> &node_vel,
                 DViewCArrayKokkos <double> &node_mass,
                 DViewCArrayKokkos <double> &elem_den,
                 DViewCArrayKokkos <double> &elem_pres,
                 DViewCArrayKokkos <double> &elem_stress,
                 DViewCArrayKokkos <double> &elem_sspd,
                 DViewCArrayKokkos <double> &elem_sie,
                 DViewCArrayKokkos <double> &elem_vol,
                 DViewCArrayKokkos <double> &elem_div,
                 DViewCArrayKokkos <double> &elem_mass,
                 DViewCArrayKokkos <size_t> &elem_mat_id,
                 DViewCArrayKokkos <double> &elem_statev,
                 DViewCArrayKokkos <double> &corner_force,
                 DViewCArrayKokkos <double> &corner_mass);

  void get_force_sgh(const CArrayKokkos <material_t> &material,
                     const mesh_t &mesh,
                     const DViewCArrayKokkos <double> &node_coords,
                     const DViewCArrayKokkos <double> &node_vel,
                     const DViewCArrayKokkos <double> &elem_den,
                     const DViewCArrayKokkos <double> &elem_sie,
                     const DViewCArrayKokkos <double> &elem_pres,
                     const DViewCArrayKokkos <double> &elem_stress,
                     const DViewCArrayKokkos <double> &elem_sspd,
                     const DViewCArrayKokkos <double> &elem_vol,
                     const DViewCArrayKokkos <double> &elem_div,
                     const DViewCArrayKokkos <size_t> &elem_mat_id,
                     DViewCArrayKokkos <double> &corner_force,
                     const double fuzz,
                     const double small,
                     const DViewCArrayKokkos <double> &elem_statev,
                     const double dt,
                     const double rk_alpha);

  void get_force_sgh2D(const CArrayKokkos <material_t> &material,
                       const mesh_t &mesh,
                       const DViewCArrayKokkos <double> &node_coords,
                       const DViewCArrayKokkos <double> &node_vel,
                       const DViewCArrayKokkos <double> &elem_den,
                       const DViewCArrayKokkos <double> &elem_sie,
                       const DViewCArrayKokkos <double> &elem_pres,
                       const DViewCArrayKokkos <double> &elem_stress,
                       const DViewCArrayKokkos <double> &elem_sspd,
                       const DViewCArrayKokkos <double> &elem_vol,
                       const DViewCArrayKokkos <double> &elem_div,
                       const DViewCArrayKokkos <size_t> &elem_mat_id,
                       DViewCArrayKokkos <double> &corner_force,
                       const double fuzz,
                       const double small,
                       const DViewCArrayKokkos <double> &elem_statev,
                       const double dt,
                       const double rk_alpha);
  
  void update_position_sgh(double rk_alpha,
                           double dt,
                           const size_t num_dims,
                           const size_t num_nodes,
                           DViewCArrayKokkos <double> &node_coords,
                           const DViewCArrayKokkos <double> &node_vel);

  void get_vol(const DViewCArrayKokkos <double> &elem_vol,
               const DViewCArrayKokkos <double> &node_coords,
               const mesh_t &mesh);

  KOKKOS_INLINE_FUNCTION
  void get_vol_hex(const DViewCArrayKokkos <double> &elem_vol,
                   const size_t elem_gid,
                   const DViewCArrayKokkos <double> &node_coords,
                   const ViewCArrayKokkos <size_t>  &elem_node_gids) const;


  KOKKOS_INLINE_FUNCTION
  void get_vol_quad(const DViewCArrayKokkos <double> &elem_vol,
                    const size_t elem_gid,
                    const DViewCArrayKokkos <double> &node_coords,
                    const ViewCArrayKokkos <size_t>  &elem_node_gids) const;


  KOKKOS_FUNCTION
  double get_area_quad(const size_t elem_gid,
                       const DViewCArrayKokkos <double> &node_coords,
                       const ViewCArrayKokkos <size_t>  &elem_node_gids) const;


  KOKKOS_FUNCTION
  void get_bmatrix(const ViewCArrayKokkos <double> &B_matrix,
                   const size_t elem_gid,
                   const DViewCArrayKokkos <double> &node_coords,
                   const ViewCArrayKokkos <size_t>  &elem_node_gids) const;


  KOKKOS_FUNCTION
  void get_bmatrix2D(const ViewCArrayKokkos <double> &B_matrix,
                     const size_t elem_gid,
                     const DViewCArrayKokkos <double> &node_coords,
                     const ViewCArrayKokkos <size_t>  &elem_node_gids) const;

  KOKKOS_FUNCTION
  void get_area_weights2D(const ViewCArrayKokkos <double> &corner_areas,
                          const size_t elem_gid,
                          const DViewCArrayKokkos <double> &node_coords,
                          const ViewCArrayKokkos <size_t>  &elem_node_gids) const;


  KOKKOS_INLINE_FUNCTION
  double heron(const double x1,
               const double y1,
               const double x2,
               const double y2,
               const double x3,
               const double y3) const;

  void get_divergence(DViewCArrayKokkos <double> &elem_div,
                      const mesh_t mesh,
                      const DViewCArrayKokkos <double> &node_coords,
                      const DViewCArrayKokkos <double> &node_vel,
                      const DViewCArrayKokkos <double> &elem_vol);


  void get_divergence2D(DViewCArrayKokkos <double> &elem_div,
                        const mesh_t mesh,
                        const DViewCArrayKokkos <double> &node_coords,
                        const DViewCArrayKokkos <double> &node_vel,
                        const DViewCArrayKokkos <double> &elem_vol);


  KOKKOS_FUNCTION
  void get_velgrad(ViewCArrayKokkos <double> &vel_grad,
                   const ViewCArrayKokkos <size_t>  &elem_node_gids,
                   const DViewCArrayKokkos <double> &node_vel,
                   const ViewCArrayKokkos <double> &b_matrix,
                   const double elem_vol,
                   const size_t elem_gid) const;


  KOKKOS_FUNCTION
  void get_velgrad2D(ViewCArrayKokkos <double> &vel_grad,
                     const ViewCArrayKokkos <size_t>  &elem_node_gids,
                     const DViewCArrayKokkos <double> &node_vel,
                     const ViewCArrayKokkos <double> &b_matrix,
                     const double elem_vol,
                     const double elem_area,
                     const size_t elem_gid) const;

  KOKKOS_INLINE_FUNCTION
  void decompose_vel_grad(ViewCArrayKokkos <double> &D_tensor,
                          ViewCArrayKokkos <double> &W_tensor,
                          const ViewCArrayKokkos <double> &vel_grad,
                          const ViewCArrayKokkos <size_t>  &elem_node_gids,
                          const size_t elem_gid,
                          const DViewCArrayKokkos <double> &node_coords,
                          const DViewCArrayKokkos <double> &node_vel,
                          const double vol) const;


  void update_velocity_sgh(double rk_alpha,
                           double dt,
                           const mesh_t &mesh,
                           DViewCArrayKokkos <double> &node_vel,
                           const DViewCArrayKokkos <double> &node_mass,
                           const DViewCArrayKokkos <double> &corner_force);
  
  void tag_bdys(const CArrayKokkos <boundary_t> &boundary,
                mesh_t &mesh,
                const DViewCArrayKokkos <double> &node_coords);

  void boundary_velocity(const mesh_t &mesh,
                         const CArrayKokkos <boundary_t> &boundary,
                         DViewCArrayKokkos <double> &node_vel);
  
  KOKKOS_INLINE_FUNCTION 
  size_t check_bdy(const size_t patch_gid,
                   const int this_bc_tag,
                   const double val,
                   const mesh_t &mesh,
                   const DViewCArrayKokkos <double> &node_coords) const;

  void rk_init(DViewCArrayKokkos <double> &node_coords,
               DViewCArrayKokkos <double> &node_vel,
               DViewCArrayKokkos <double> &elem_sie,
               DViewCArrayKokkos <double> &elem_stress,
               const size_t num_dims,
               const size_t num_elems,
               const size_t num_nodes);


  void get_timestep(mesh_t &mesh,
                    DViewCArrayKokkos <double> &node_coords,
                    DViewCArrayKokkos <double> &node_vel,
                    DViewCArrayKokkos <double> &elem_sspd,
                    DViewCArrayKokkos <double> &elem_vol,
                    double time_value,
                    const double graphics_time,
                    const double time_final,
                    const double dt_max,
                    const double dt_min,
                    const double dt_cfl,
                    double &dt,
                    const double fuzz);


  void get_timestep2D(mesh_t &mesh,
                      DViewCArrayKokkos <double> &node_coords,
                      DViewCArrayKokkos <double> &node_vel,
                      DViewCArrayKokkos <double> &elem_sspd,
                      DViewCArrayKokkos <double> &elem_vol,
                      double time_value,
                      const double graphics_time,
                      const double time_final,
                      const double dt_max,
                      const double dt_min,
                      const double dt_cfl,
                      double &dt,
                      const double fuzz);
  
  void update_energy_sgh(double rk_alpha,
                         double dt,
                         const mesh_t &mesh,
                         const DViewCArrayKokkos <double> &node_vel,
                         const DViewCArrayKokkos <double> &node_coords,
                         DViewCArrayKokkos <double> &elem_sie,
                         const DViewCArrayKokkos <double> &elem_mass,
                         const DViewCArrayKokkos <double> &corner_force);
                   
  void update_state(const CArrayKokkos <material_t> &material,
                    const mesh_t &mesh,
                    const DViewCArrayKokkos <double> &node_coords,
                    const DViewCArrayKokkos <double> &node_vel,
                    DViewCArrayKokkos <double> &elem_den,
                    DViewCArrayKokkos <double> &elem_pres,
                    DViewCArrayKokkos <double> &elem_stress,
                    DViewCArrayKokkos <double> &elem_sspd,
                    const DViewCArrayKokkos <double> &elem_sie,
                    const DViewCArrayKokkos <double> &elem_vol,
                    const DViewCArrayKokkos <double> &elem_mass,
                    const DViewCArrayKokkos <size_t> &elem_mat_id,
                    const DViewCArrayKokkos <double> &elem_statev,
                    const double dt,
                    const double rk_alpha);


  void update_state2D(const CArrayKokkos <material_t> &material,
                      const mesh_t &mesh,
                      const DViewCArrayKokkos <double> &node_coords,
                      const DViewCArrayKokkos <double> &node_vel,
                      DViewCArrayKokkos <double> &elem_den,
                      DViewCArrayKokkos <double> &elem_pres,
                      DViewCArrayKokkos <double> &elem_stress,
                      DViewCArrayKokkos <double> &elem_sspd,
                      const DViewCArrayKokkos <double> &elem_sie,
                      const DViewCArrayKokkos <double> &elem_vol,
                      const DViewCArrayKokkos <double> &elem_mass,
                      const DViewCArrayKokkos <size_t> &elem_mat_id,
                      const DViewCArrayKokkos <double> &elem_statev,
                      const double dt,
                      const double rk_alpha);

  KOKKOS_INLINE_FUNCTION
  void user_eos_model(const DViewCArrayKokkos <double> &elem_pres,
                      const DViewCArrayKokkos <double> &elem_stress,
                      const size_t elem_gid,
                      const size_t mat_id,
                      const DViewCArrayKokkos <double> &elem_state_vars,
                      const DViewCArrayKokkos <double> &elem_sspd,
                      const double den,
                      const double sie) const;


  KOKKOS_INLINE_FUNCTION
  void user_strength_model(const DViewCArrayKokkos <double> &elem_pres,
                           const DViewCArrayKokkos <double> &elem_stress,
                           const size_t elem_gid,
                           const size_t mat_id,
                           const DViewCArrayKokkos <double> &elem_state_vars,
                           const DViewCArrayKokkos <double> &elem_sspd,
                           const double den,
                           const double sie,
                           const ViewCArrayKokkos <double> &vel_grad,
                           const ViewCArrayKokkos <size_t> &elem_node_gids,
                           const DViewCArrayKokkos <double> &node_coords,
                           const DViewCArrayKokkos <double> &node_vel,
                           const double vol,
                           const double dt,
                           const double rk_alpha) const;


  KOKKOS_INLINE_FUNCTION
  void user_strength_model_vpsc(const DViewCArrayKokkos <double> &elem_pres,
                              const DViewCArrayKokkos <double> &elem_stress,
                              const size_t elem_gid,
                              const size_t mat_id,
                              const DViewCArrayKokkos <double> &elem_state_vars,
                              const DViewCArrayKokkos <double> &elem_sspd,
                              const double den,
                              const double sie,
                              const ViewCArrayKokkos <double> &vel_grad,
                              const ViewCArrayKokkos <size_t> &elem_node_gids,
                              const DViewCArrayKokkos <double> &node_coords,
                              const DViewCArrayKokkos <double> &node_vel,
                              const double vol,
                              const double dt,
                              const double rk_alpha) const;


void user_model_init(const DCArrayKokkos <double> &file_state_vars,
                     const size_t num_state_vars,
                     const size_t mat_id,
                     const size_t num_elems);

  void build_boundry_node_sets(const CArrayKokkos <boundary_t> &boundary, mesh_t &mesh);
  
  void init_boundaries();

  //initializes memory for arrays used in the global stiffness matrix assembly
  void init_boundary_sets(int num_boundary_sets);

  void grow_boundary_sets(int num_boundary_sets);

  int solve();

  void comm_variables(Teuchos::RCP<const MV> zp);

  void read_conditions_ansys_dat(std::ifstream *in, std::streampos before_condition_header);

  //interfaces between user input and creating data structures for bcs
  void generate_bcs();

  void Displacement_Boundary_Conditions();

  void init_output();

  void compute_output();
  
  void sort_output(Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > sorted_map);

  void collect_output(Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > global_reduce_map);

  
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
                      DViewCArrayKokkos <size_t> &elem_mat_id,
                      CArray <double> &graphics_times,
                      size_t &graphics_id,
                      const double time_value);


  void ensight(const mesh_t &mesh,
               const DViewCArrayKokkos <double> &node_coords,
               const DViewCArrayKokkos <double> &node_vel,
               const DViewCArrayKokkos <double> &node_mass,
               const DViewCArrayKokkos <double> &elem_den,
               const DViewCArrayKokkos <double> &elem_pres,
               const DViewCArrayKokkos <double> &elem_stress,
               const DViewCArrayKokkos <double> &elem_sspd,
               const DViewCArrayKokkos <double> &elem_sie,
               const DViewCArrayKokkos <double> &elem_vol,
               const DViewCArrayKokkos <double> &elem_mass,
               const DViewCArrayKokkos <size_t> &elem_mat_id,
               CArray <double> &graphics_times,
               size_t &graphics_id,
               const double time_value);


  void state_file(const mesh_t &mesh,
                  const DViewCArrayKokkos <double> &node_coords,
                  const DViewCArrayKokkos <double> &node_vel,
                  const DViewCArrayKokkos <double> &node_mass,
                  const DViewCArrayKokkos <double> &elem_den,
                  const DViewCArrayKokkos <double> &elem_pres,
                  const DViewCArrayKokkos <double> &elem_stress,
                  const DViewCArrayKokkos <double> &elem_sspd,
                  const DViewCArrayKokkos <double> &elem_sie,
                  const DViewCArrayKokkos <double> &elem_vol,
                  const DViewCArrayKokkos <double> &elem_mass,
                  const DViewCArrayKokkos <size_t> &elem_mat_id,
                  const double time_value );

  void node_density_constraints(host_vec_array node_densities_lower_bound);
  
  Simulation_Parameters_SGH *simparam;
  Simulation_Parameters_Topology_Optimization *simparam_TO;
  Explicit_Solver_SGH *Explicit_Solver_Pointer_;
  
  //output stream
  Teuchos::RCP<Teuchos::FancyOStream> fos;
  
  elements::element_selector *element_select;
  elements::Element3D *elem;
  elements::Element2D *elem2D;
  elements::ref_element  *ref_elem;
  
  Explicit_Solver_SGH *explicit_solver_pointer;
  
  //Local FEA data
  size_t nlocal_nodes;
  dual_vec_array dual_node_coords; //coordinates of the nodes
  dual_vec_array dual_node_densities; //topology optimization design variable
  dual_elem_conn_array dual_nodes_in_elem; //dual view of element connectivity to nodes
  host_elem_conn_array nodes_in_elem; //host view of element connectivity to nodes
  CArrayKokkos<elements::elem_types::elem_type, array_layout, HostSpace, memory_traits> Element_Types;
  CArrayKokkos<size_t, array_layout, HostSpace, memory_traits> Nodes_Per_Element_Type;

  //Ghost data on this MPI rank
  size_t nghost_nodes;
  CArrayKokkos<GO, Kokkos::LayoutLeft, node_type::device_type> ghost_nodes;
  CArrayKokkos<int, array_layout, device_type, memory_traits> ghost_node_ranks;

  //Local FEA data including ghosts
  size_t nall_nodes;
  size_t rnum_elem;

  //Global FEA data
  long long int num_nodes, num_elem;
  Teuchos::RCP<const Teuchos::Comm<int> > comm;
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > map; //map of node indices
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > ghost_node_map; //map of node indices with ghosts on each rank
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > all_node_map; //map of node indices with ghosts on each rank
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > element_map; //non overlapping map of elements owned by each rank used in reduction ops
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > all_element_map; //overlapping map of elements connected to the local nodes in each rank
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > local_dof_map; //map of local dofs (typically num_node_local*num_dim)
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > all_dof_map; //map of local and ghost dofs (typically num_node_all*num_dim)
  Teuchos::RCP<MCONN> nodes_in_elem_distributed; //element to node connectivity table
  Teuchos::RCP<MCONN> node_nconn_distributed; //how many elements a node is connected to
  Teuchos::RCP<MV> node_coords_distributed;
  Teuchos::RCP<MV> all_node_coords_distributed;
  Teuchos::RCP<MV> design_node_densities_distributed;
  Teuchos::RCP<const MV> test_node_densities_distributed;
  Teuchos::RCP<MV> all_node_densities_distributed;
  Teuchos::RCP<MV> Global_Element_Densities;
  
  //Boundary Conditions Data
  //CArray <Nodal_Combination> Patch_Nodes;
  size_t nboundary_patches;
  size_t num_boundary_conditions;
  int current_bdy_id;
  CArrayKokkos<Node_Combination, array_layout, device_type, memory_traits> Boundary_Patches;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Boundary_Condition_Patches; //set of patches corresponding to each boundary condition
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> NBoundary_Condition_Patches;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Boundary_Condition_Patches_strides;
  enum bc_type {NONE, POINT_LOADING_CONDITION, LINE_LOADING_CONDITION, SURFACE_LOADING_CONDITION};

  //element selection parameters and data
  size_t max_nodes_per_element;

  //determines if rhs gets a contribution from bcs
  bool nonzero_bc_flag;

  //lists what kind of boundary condition the nodal DOF is subjected to if any
  CArrayKokkos<int, array_layout, device_type, memory_traits> Node_DOF_Boundary_Condition_Type;
  //lists what kind of boundary condition each boundary set is assigned to
  CArrayKokkos<int, array_layout, HostSpace, memory_traits> Boundary_Condition_Type_List;
  
  //number of displacement boundary conditions acting on nodes; used to size the reduced global stiffness map
  size_t Number_DOF_BCS;

  //MPI data
  int myrank; //index of this mpi rank in the world communicator
  int nranks; //number of mpi ranks in the world communicator
  MPI_Comm world; //stores the default communicator object (MPI_COMM_WORLD)

  //! mapping used to get local ghost index from the global ID.
  //typedef ::Tpetra::Details::FixedHashTable<GO, LO, Kokkos::HostSpace::device_type>
    //global_to_local_table_host_type;

  //global_to_local_table_host_type global2local_map;
  //CArrayKokkos<int, Kokkos::LayoutLeft, Kokkos::HostSpace::device_type> active_ranks;

  //Pertains to local mesh information being stored as prescribed by the row map
  global_size_t min_gid;
  global_size_t max_gid;
  global_size_t index_base;

  //allocation flags to avoid repeat MV and global matrix construction
  int Matrix_alloc;

  //debug flags
  int gradient_print_sync;

  //Topology Optimization parameter
  int penalty_power;
  bool nodal_density_flag;

  //runtime and counters for performance output
  double linear_solve_time, hessvec_time, hessvec_linear_time;
  int update_count, hessvec_count;

  //nodal DOF output data
  enum vector_styles {NODAL, DOF}; //multivector can store as ndof by 1 or nnode by vector_size
  int noutput;
  bool displaced_mesh_flag;
  int displacement_index;
  std::vector<std::vector<std::string>> output_dof_names;
  std::vector<const_host_vec_array> module_outputs;
  std::vector<vector_styles> vector_style;
  std::vector<int> output_vector_sizes;
  
  //Boundary Conditions Data
  int max_boundary_sets;

  //output dof data
  //Global arrays with collected data used to print
  int output_velocity_index, output_strain_index, output_stress_index;
};

#endif // end HEADER_H
