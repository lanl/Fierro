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
#ifndef AO_CONTACT_H
#define AO_CONTACT_H

#include "matar.h"
#include "simulation_parameters.hpp"
#include "state.hpp"

using namespace mtr;

struct AO_contact_state_t
{
    // ---------------------------------------------------------------
    // boundary node variables
    // ---------------------------------------------------------------
    DCArrayKokkos <double> bdy_node_coords;     // subset of coords only including boundary nodes
    CArrayKokkos <double> bdy_node_vels;        // boundary node velocities (TLQS: delta_u/dt)
    CArrayKokkos <double> bdy_node_accels;      // boundary node accelerations (TLQS: zero)
    CArrayKokkos <double> bdy_node_mass;        // (num_bdy_nodes) lumped nodal mass

    // ---------------------------------------------------------------
    // sort variables (rebuilt once per load step)
    // ---------------------------------------------------------------
    CArrayKokkos <double> bounding_boxes;                       // coords of bounding box of each boundary surface
    RaggedRightArrayKokkos <size_t> nodes_in_bounding_boxes;    // nodes that lie in the bounding box of each boundary surface
    DCArrayKokkos <size_t> num_nodes_in_bounding_boxes;         // stride array for nodes_in_bounding_boxes
    const size_t num_bins = 10; // TODO: THIS NEEDS TO BE EITHER CALCULATED BASED ON MESH OR SET AS AN INPUT FROM THE YAML
    double lebesgue_overshoot;
    swage::PointCloud_t bdy_node_point_cloud;
    DCArrayKokkos<size_t>          num_surfs_in_node_boxes;     // (num_bdy_nodes) number of surface bounding boxes containing each node
    RaggedRightArrayKokkos<size_t> surfs_in_node_boxes;         // (bdy node, slot) -> bdy_surf_lid
    RaggedRightArrayKokkos<size_t> cols_in_node_boxes;          // (bdy node, slot) -> node_lid column in pairing_check_vars

    // ---------------------------------------------------------------
    // detection variables (reset every sweep)
    // ---------------------------------------------------------------
    DRaggedRightArrayKokkos <double> pairing_check_vars;        // (surf, box node, 3) gap, xi, eta; sentinel if not penetrating
    CArrayKokkos<double> bdy_surf_node_normals;                 // outward normal at each node in each surface
    const double filter_tol = 0.0; // TODO: THIS NEEDS TO BE EITHER CALCULATED BASED ON MESH OR SET AS AN INPUT FROM THE YAML
    double max_gap;
    CArrayKokkos<size_t> num_pen_surfs_in_node;                 // (num_bdy_nodes) surfaces each node penetrates this sweep
    CArrayKokkos<size_t> pen_nodes;                             // (num_bdy_nodes) penetrating bdy node lids this sweep, valid to contact_counts(0)
    DCArrayKokkos <double> bdy_node_pred_coords;                // (num_bdy_nodes, 3) predicted coords for repeated sweeps

    // ---------------------------------------------------------------
    // counters, host copy current after penetration_sweep
    //   [0] penetrating nodes this sweep
    //   [1] pairing errors this sweep (should be 0)
    //   [2] active nodes this load step (valid length of active_nodes)
    //   [3] pair-set changes this sweep (pairs added + pairs deleted by refresh)
    // ---------------------------------------------------------------
    DCArrayKokkos<size_t> contact_counts;

    // ---------------------------------------------------------------
    // persistent pair variables
    // reset once per load step in AO_contact_sort; persist across Picard
    // iterations, refreshed every sweep, new pairs appended by the sweep
    // ---------------------------------------------------------------
    CArrayKokkos<size_t>            is_active;                  // (num_bdy_nodes) 1 once the node is in active_nodes this load step
    CArrayKokkos<size_t>            active_nodes;               // (num_bdy_nodes) bdy node lids with pairs, valid to contact_counts(2)
    CArrayKokkos<size_t>            num_pairs_in_node;          // (num_bdy_nodes) current pairs per node
    RaggedRightArrayKokkos<size_t>  contact_pair_surfs;         // (bdy node, pair) -> target bdy_surf_lid
    DRaggedRightArrayKokkos<double> contact_pair_vars;          // (bdy node, pair, 7) -> xi, eta, nx, ny, nz, F_tot, f_inc
    CArrayKokkos<size_t>            num_pairs_touching_node;    // (num_bdy_nodes) pairs involving each node (relaxation weight), recounted every sweep

    // ---------------------------------------------------------------
    // contact forces
    // ---------------------------------------------------------------
    CArrayKokkos<double> contact_forces;                        // (num_bdy_nodes, 3) assembled from the pairs
    CArrayKokkos<double> contact_forces_used;                   // (num_bdy_nodes, 3) copy handed to the last Picard solve

};

// ********************************************************
// STARTING FUNCTIONS FOR INITIALIZATION OF CONTACT STATE
// ********************************************************

// gets the max of the lebesgue function
double get_lebesgue_constant_1d(const elements::ReferenceElement_t& ref_elem,
                                const size_t num_samples = 100001);

// getting the overshoot factor for sizing bounding boxes
void get_surf_overshoot_factor(const elements::ReferenceElement_t& ref_elem, double& lebesgue_overshoot);

// row-sum lumped mass of the boundary nodes from the material point densities
void get_bdy_node_mass(const elements::ReferenceElement_t& ref_elem,
                       const elements::Quadrature_t& Quad,
                       const swage::Mesh_t& mesh,
                       const State_t& State,
                       const MPICArrayKokkos<double>& node_coords,
                       CArrayKokkos<double>& bdy_node_mass);

void AO_contact_initialize(AO_contact_state_t& AO_contact_state,
                           const elements::ReferenceElement_t& ref_elem,
                           const elements::Quadrature_t& Quad,
                           const swage::Mesh_t& mesh,
                           const State_t& State,
                           const MPICArrayKokkos<double>& node_coords,
                           const size_t num_bdy_nodes,
                           const size_t num_bdy_surfs,
                           const size_t num_nodes_in_surf);

// ********************************************************
// ENDING FUNCTIONS FOR INITIALIZATION OF CONTACT STATE
// ********************************************************



// ********************************************************
// STARTING FUNCTIONS FOR SORTING NODES FOR PAIRING
// ********************************************************

// uses a kokkos parallel reduce to get the max values in one kernel launch
void get_max_vel_and_accel(double& vx_max, double& vy_max, double& vz_max,
                           double& ax_max, double& ay_max, double& az_max,
                           const CArrayKokkos <double>& bdy_node_vels,
                           const CArrayKokkos <double>& bdy_node_accels);

// gets the bounding box for all boundary surfaces
void get_bounding_boxes(const DCArrayKokkos <double>& bdy_node_coords,
                        const CArrayKokkos <double>& bdy_node_vels,
                        const CArrayKokkos <double>& bdy_node_accels,
                        const size_t num_bdy_surfs,
                        const size_t num_nodes_in_surf,
                        const double dt,
                        const double lebesgue_overshoot,
                        const CArrayKokkos <size_t>& bdy_nodes_in_bdy_surf,
                        CArrayKokkos <double>& bounding_boxes);

double get_min_bdy_elem_edge_length(const size_t num_bdy_surfs,
                                    const size_t num_dofs_1d,
                                    const CArrayKokkos<size_t>& bdy_surfs,
                                    const CArrayKokkos<int>& elems_in_surf,
                                    const DCArrayKokkos<size_t>& nodes_in_elem,
                                    const MPICArrayKokkos<double>& node_coords);

void invert_bounding_boxes(const size_t num_bdy_surfs,
                           const size_t num_bdy_nodes,
                           const RaggedRightArrayKokkos<size_t>& nodes_in_bounding_boxes,
                           DCArrayKokkos<size_t>& num_surfs_in_node_boxes,
                           RaggedRightArrayKokkos<size_t>& surfs_in_node_boxes,
                           RaggedRightArrayKokkos<size_t>& cols_in_node_boxes);

// updates bdy_node_coords and nodes_in_bounding_boxes
void AO_contact_sort(DCArrayKokkos <double>& bdy_node_coords,
                     const CArrayKokkos <double>& bdy_node_vels,
                     const CArrayKokkos <double>& bdy_node_accels,
                     const size_t num_bdy_nodes,
                     const MPICArrayKokkos <double>& node_coords,
                     const CArrayKokkos <size_t>& bdy_nodes,
                     swage::PointCloud_t& bdy_node_point_cloud,
                     const size_t num_bins,
                     const size_t num_bdy_surfs,
                     const size_t num_nodes_in_surf,
                     const double dt,
                     const CArrayKokkos <size_t>& bdy_nodes_in_bdy_surf,
                     CArrayKokkos <double>& bounding_boxes,
                     const double lebesgue_overshoot,
                     DCArrayKokkos <size_t>& num_nodes_in_bounding_boxes,
                     RaggedRightArrayKokkos <size_t>& nodes_in_bounding_boxes,
                     const size_t num_dofs_1d,
                     const CArrayKokkos<size_t>& bdy_surfs,
                     const CArrayKokkos<int>& elems_in_surf,
                     const DCArrayKokkos<size_t>& nodes_in_elem,
                     double& max_gap,
                     DRaggedRightArrayKokkos<double>& pairing_check_vars,
                     DCArrayKokkos<size_t>& num_surfs_in_node_boxes,
                     RaggedRightArrayKokkos<size_t>& surfs_in_node_boxes,
                     RaggedRightArrayKokkos<size_t>& cols_in_node_boxes,
                     RaggedRightArrayKokkos<size_t>& contact_pair_surfs,
                     DRaggedRightArrayKokkos<double>& contact_pair_vars,
                     CArrayKokkos<size_t>& num_pairs_in_node,
                     CArrayKokkos<size_t>& is_active,
                     CArrayKokkos<double>& contact_forces,
                     CArrayKokkos<double>& contact_forces_used,
                     DCArrayKokkos<size_t>& contact_counts);

// ********************************************************
// ENDING FUNCTIONS FOR SORTING NODES FOR PAIRING
// ********************************************************



// ********************************************************
// STARTING FUNCTIONS FOR CHECKING PENETRATION
// ********************************************************

KOKKOS_FUNCTION
void lagrange_1D(const CArrayKokkos<double>& dof_positions_1d,
                 const size_t num_dofs_1d,
                 const size_t a,
                 const double x,
                 double& val,
                 double& dval);

// 1D Lagrange basis values and first derivatives at x
KOKKOS_FUNCTION
void lagrange_val_and_deriv_1D(double* val,
                               double* dval,
                               const CArrayKokkos<double>& dof_positions_1d,
                               const size_t num_dofs_1d,
                               const double x);

// build the cross product to get the normal direction
KOKKOS_FUNCTION
void get_normal(const CArrayKokkos<double>& dof_positions_1d,
                const ViewCArrayKokkos<size_t>& bdy_nodes_in_the_surf,
                const DCArrayKokkos<double>& bdy_node_coords,
                const size_t face_lid,
                const double xi,
                const double eta,
                double* normal);

// outward unit normal at each GLL node of each boundary surface
void get_bdy_surf_node_normals(const size_t num_bdy_surfs,
                               const size_t num_nodes_in_surf,
                               const CArrayKokkos<size_t>& bdy_surfs,
                               const CArrayKokkos<int>& faces_in_surf,
                               const CArrayKokkos<size_t>& bdy_nodes_in_bdy_surf,
                               const CArrayKokkos<double>& dof_positions_1d,
                               const DCArrayKokkos<double>& bdy_node_coords,
                               CArrayKokkos<double>& bdy_surf_node_normals);

// returns true if the node should proceed to the Newton solve
// filter 1: node is part of this surface's connectivity         -> reject
// filter 2: node is outside the tangent plane at every
//           surface node by more than filter_tol                -> reject
KOKKOS_FUNCTION
bool check_filters(const size_t bdy_node_lid,
                   const ViewCArrayKokkos<size_t>& bdy_nodes_in_the_surf,
                   const ViewCArrayKokkos<double>& surf_node_normals,
                   const DCArrayKokkos<double>& bdy_node_coords,
                   const size_t num_nodes_in_surf,
                   const double filter_tol);

KOKKOS_FUNCTION
void get_face_dims(const size_t face_lid,
                   const size_t num_dofs_1d,
                   size_t& fixed_dim,
                   size_t& xi_dim,
                   size_t& eta_dim,
                   size_t& face_layer);

KOKKOS_FUNCTION
void lagrange_1D_d2(const CArrayKokkos<double>& dof_positions_1d,
                    const size_t num_dofs_1d,
                    const size_t a,
                    const double x,
                    double& val,
                    double& dval,
                    double& ddval);

KOKKOS_FUNCTION
void get_surf_point_and_derivs(const CArrayKokkos<double>& dof_positions_1d,
                               const ViewCArrayKokkos<size_t>& bdy_nodes_in_the_surf,
                               const DCArrayKokkos<double>& bdy_node_coords,
                               const double xi,
                               const double eta,
                               double* x_s,
                               double* x_xi,
                               double* x_eta,
                               double* x_xixi,
                               double* x_xieta,
                               double* x_etaeta);

KOKKOS_FUNCTION
bool closest_point_projected_newton(const double* x_I,
                                    const ViewCArrayKokkos<size_t>& bdy_nodes_in_the_surf,
                                    const DCArrayKokkos<double>& bdy_node_coords,
                                    const CArrayKokkos<double>& dof_positions_1d,
                                    double& xi,
                                    double& eta,
                                    double& f);

double seed_vel_alignment(const double* x_I,
                          const double* x_c,
                          const double* v,
                          const double  v_mag);

// is the node penetrating the surface
KOKKOS_FUNCTION
bool penetration_check(const size_t bdy_node_lid,
                       const size_t bdy_surf_lid,
                       const size_t node_lid,
                       const ViewCArrayKokkos<size_t>& bdy_nodes_in_the_surf,
                       const DCArrayKokkos<double>& bdy_node_coords,
                       const CArrayKokkos<double>& bdy_node_vels,
                       const CArrayKokkos<double>& dof_positions_1d,
                       const CArrayKokkos<double>& surf_qpt_basis,
                       const CArrayKokkos<double>& surf_qpt_positions,
                       const size_t face_lid,
                       const double max_gap,
                       const DRaggedRightArrayKokkos<double>& pairing_check_vars);

KOKKOS_FUNCTION
bool refresh_closest_point(const size_t bdy_node_lid,
                           const ViewCArrayKokkos<size_t>& bdy_nodes_in_the_surf,
                           const DCArrayKokkos<double>& bdy_node_coords,
                           const CArrayKokkos<double>& dof_positions_1d,
                           const double slide_tol,
                           double& xi,
                           double& eta);

// find contact pairs from nodes_in_bounding_boxes
void penetration_sweep(const size_t num_bdy_surfs,
                       const size_t num_bdy_nodes,
                       const size_t num_nodes_in_surf,
                       const size_t num_nodes_in_elem,
                       const CArrayKokkos<size_t>& bdy_surfs,
                       const CArrayKokkos<int>& faces_in_surf,
                       const CArrayKokkos<int>& elems_in_surf,
                       const CArrayKokkos<size_t>& bdy_nodes,
                       const RaggedRightArrayKokkos<size_t>& elems_in_node,
                       const DCArrayKokkos<size_t>& nodes_in_elem,
                       const CArrayKokkos<size_t>& bdy_nodes_in_bdy_surf,
                       const CArrayKokkos<double>& dof_positions_1d,
                       const CArrayKokkos<double>& surf_qpt_basis,
                       const CArrayKokkos<double>& surf_qpt_positions,
                       const MPICArrayKokkos<double>& node_coords,
                       const DCArrayKokkos<double>& bdy_node_coords,
                       const CArrayKokkos<double>& bdy_node_vels,
                       const CArrayKokkos<double>& bdy_surf_node_normals,
                       const RaggedRightArrayKokkos<size_t>& nodes_in_bounding_boxes,
                       const RaggedRightArrayKokkos<size_t>& surfs_in_node_boxes,
                       const RaggedRightArrayKokkos<size_t>& cols_in_node_boxes,
                       const double filter_tol,
                       const double max_gap,
                       DRaggedRightArrayKokkos<double>& pairing_check_vars,
                       CArrayKokkos<size_t>& num_pen_surfs_in_node,
                       CArrayKokkos<size_t>& pen_nodes,
                       DCArrayKokkos<size_t>& contact_counts,
                       CArrayKokkos<size_t>& is_active,
                       CArrayKokkos<size_t>& active_nodes,
                       CArrayKokkos<size_t>& num_pairs_in_node,
                       RaggedRightArrayKokkos<size_t>& contact_pair_surfs,
                       DRaggedRightArrayKokkos<double>& contact_pair_vars,
                       CArrayKokkos<size_t>& num_pairs_touching_node);

// ********************************************************
// ENDING FUNCTIONS FOR CHECKING PENETRATION
// ********************************************************



// ********************************************************
// STARTING FUNCTIONS FOR GETTING CONTACT FORCES
// ********************************************************

void AO_contact_update_bdy_coords(const size_t num_bdy_nodes,
                                  const MPICArrayKokkos<double>& node_coords,
                                  const CArrayKokkos<size_t>& bdy_nodes,
                                  DCArrayKokkos<double>& bdy_node_coords);

void AO_contact_snapshot_forces(const size_t num_bdy_nodes,
                                const CArrayKokkos<double>& contact_forces,
                                CArrayKokkos<double>& contact_forces_used);

void predict_bdy_node_coords(const size_t num_bdy_nodes,
                             const double dt,
                             const DCArrayKokkos<double>& bdy_node_coords,
                             const CArrayKokkos<double>& bdy_node_mass,
                             const CArrayKokkos<double>& contact_forces,
                             const CArrayKokkos<double>& contact_forces_used,
                             DCArrayKokkos<double>& bdy_node_pred_coords);

KOKKOS_FUNCTION
void scatter_pair_force(const size_t bdy_node_lid,
                        const ViewCArrayKokkos<size_t>& surf_nodes,
                        const CArrayKokkos<double>& dof_positions_1d,
                        const double xi,
                        const double eta,
                        const double* n,
                        const double F,
                        const CArrayKokkos<double>& contact_forces);

void assemble_contact_forces(const size_t num_active,
                             const size_t num_nodes_in_surf,
                             const CArrayKokkos<size_t>& active_nodes,
                             const CArrayKokkos<size_t>& num_pairs_in_node,
                             const RaggedRightArrayKokkos<size_t>& contact_pair_surfs,
                             const DRaggedRightArrayKokkos<double>& contact_pair_vars,
                             const CArrayKokkos<size_t>& bdy_nodes_in_bdy_surf,
                             const CArrayKokkos<double>& dof_positions_1d,
                             CArrayKokkos<double>& contact_forces);

KOKKOS_FUNCTION
int frictionless_increment(const size_t bdy_node_lid,
                           const ViewCArrayKokkos<size_t>& surf_nodes,
                           const DCArrayKokkos<double>& bdy_node_coords,
                           const CArrayKokkos<double>& bdy_node_mass,
                           const CArrayKokkos<double>& contact_forces,
                           const CArrayKokkos<double>& contact_forces_used,
                           const CArrayKokkos<double>& dof_positions_1d,
                           const double dt,
                           const double slide_tol,
                           const double res_tol,
                           double& xi,
                           double& eta,
                           double* n,
                           double& f_inc);

bool force_resolution(const size_t num_active,
                      const size_t num_nodes_in_surf,
                      const double dt,
                      const double slide_tol,
                      const double res_tol,
                      const double rel_tol,
                      const size_t max_global_iter,
                      const CArrayKokkos<size_t>& active_nodes,
                      CArrayKokkos<size_t>& num_pairs_in_node,
                      RaggedRightArrayKokkos<size_t>& contact_pair_surfs,
                      DRaggedRightArrayKokkos<double>& contact_pair_vars,
                      const CArrayKokkos<size_t>& num_pairs_touching_node,
                      const CArrayKokkos<size_t>& bdy_nodes_in_bdy_surf,
                      const CArrayKokkos<double>& dof_positions_1d,
                      const DCArrayKokkos<double>& bdy_node_coords,
                      const CArrayKokkos<double>& bdy_node_mass,
                      const CArrayKokkos<double>& contact_forces_used,
                      CArrayKokkos<double>& contact_forces,
                      size_t& num_deleted,
                      size_t& num_failed);

void AO_contact_resolve(AO_contact_state_t& AO_contact_state,
                        const size_t num_bdy_surfs,
                        const size_t num_bdy_nodes,
                        const size_t num_nodes_in_surf,
                        const size_t num_nodes_in_elem,
                        const CArrayKokkos<size_t>& bdy_surfs,
                        const CArrayKokkos<int>& faces_in_surf,
                        const CArrayKokkos<int>& elems_in_surf,
                        const CArrayKokkos<size_t>& bdy_nodes,
                        const RaggedRightArrayKokkos<size_t>& elems_in_node,
                        const DCArrayKokkos<size_t>& nodes_in_elem,
                        const CArrayKokkos<size_t>& bdy_nodes_in_bdy_surf,
                        const CArrayKokkos<double>& dof_positions_1d,
                        const CArrayKokkos<double>& surf_qpt_basis,
                        const CArrayKokkos<double>& surf_qpt_positions,
                        const MPICArrayKokkos<double>& node_coords,
                        const double dt);

// ********************************************************
// ENDING FUNCTIONS FOR GETTING CONTACT FORCES
// ********************************************************

#endif  // AO_CONTACT_H