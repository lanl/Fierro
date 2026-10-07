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

using namespace mtr;

struct AO_contact_state_t
{
    // bounding node variables
    DCArrayKokkos <double> bdy_node_coords;                     // subset of coords only including boundary nodes
    CArrayKokkos <double> bdy_node_vels;                        // subset of coords only including boundary nodes
    CArrayKokkos <double> bdy_node_accels;                      // subset of coords only including boundary nodes

    // sort variables
    CArrayKokkos <double> bounding_boxes;                       // coords of bounding box of each boundary surface
    RaggedRightArrayKokkos <size_t> nodes_in_bounding_boxes;    // nodes that lie in the bounding box of each boundary surface
    DCArrayKokkos <size_t> num_nodes_in_bounding_boxes;         // stride array for nodes_in_bounding_boxes
    const size_t num_bins = 10; // TODO: THIS NEEDS TO BE EITHER CALCULATED BASED ON MESH OR SET AS AN INPUT FROM THE YAML
    double lebesgue_overshoot;
    swage::PointCloud_t bdy_node_point_cloud;

    // pairing variables
    DRaggedRightArrayKokkos <double> pairing_check_vars;        // stores min(distance_to_surf), xi, and eta for a node compared to a surface
    CArrayKokkos<double> bdy_surf_node_normals;                 // stores outward normal at each node in surface to avoid redundant calculations
    const double filter_tol = 0.0;

};

// ********************************************************
// STARTING FUNCTIONS FOR INITIALIZATION OF CONTACT STATE
// ********************************************************

// gets the max of the lebesgue function
double get_lebesgue_constant_1d(const elements::ReferenceElement_t& ref_elem,
                                const size_t num_samples = 100001);

// getting the overshoot factor for sizing bounding boxes
void get_surf_overshoot_factor(const elements::ReferenceElement_t& ref_elem, double& lebesgue_overshoot);

// sizes necessary arrays and gets static variables
void AO_contact_initialize(DCArrayKokkos <double>& bdy_node_coords,
                           CArrayKokkos <double>& bdy_node_vels,
                           CArrayKokkos <double>& bdy_node_accels,
                           const size_t num_bdy_nodes,
                           DCArrayKokkos <size_t>& num_nodes_in_bounding_boxes,
                           const size_t num_bdy_surfs,
                           const elements::ReferenceElement_t& ref_elem,
                           double& lebesgue_overshoot,
                           CArrayKokkos <double>& bounding_boxes,
                           CArrayKokkos<double>& bdy_surf_node_normals,
                           const size_t num_nodes_in_elem);

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
                     RaggedRightArrayKokkos <size_t>& nodes_in_bounding_boxes);

// ********************************************************
// ENDING FUNCTIONS FOR SORTING NODES FOR PAIRING
// ********************************************************



// ********************************************************
// STARTING FUNCTIONS FOR CHECKING PENETRATION
// ********************************************************

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
void get_bdy_surf_node_normals(const swage::Mesh_t& mesh,
                               const CArrayKokkos<double>& dof_positions_1d,
                               const DCArrayKokkos<double>& bdy_node_coords,
                               const CArrayKokkos<size_t>& bdy_nodes_in_bdy_surf,
                               const size_t num_bdy_surfs,
                               const size_t num_nodes_in_surf,
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

// is the node penetrating the surface
void penetration_check();

// find contact pairs from nodes_in_bounding_boxes
void penetration_sweep(const swage::Mesh_t& mesh,
                       const CArrayKokkos <double>& dof_positions_1d,
                       const RaggedRightArrayKokkos <size_t>& nodes_in_bounding_boxes,
                       DCArrayKokkos <size_t>& num_nodes_in_bounding_boxes,
                       DRaggedRightArrayKokkos <double>& pairing_check_vars,
                       const DCArrayKokkos <double>& bdy_node_coords,
                       const CArrayKokkos <size_t>& bdy_nodes_in_bdy_surf,
                       CArrayKokkos <double>& bdy_surf_node_normals,
                       const size_t num_bdy_surfs,
                       const size_t num_nodes_in_surf,
                       const double filter_tol);

// ********************************************************
// ENDING FUNCTIONS FOR CHECKING PENETRATION
// ********************************************************



// ********************************************************
// STARTING FUNCTIONS FOR GETTING CONTACT FORCES
// ********************************************************

// ********************************************************
// ENDING FUNCTIONS FOR GETTING CONTACT FORCES
// ********************************************************

#endif  // AO_CONTACT_H