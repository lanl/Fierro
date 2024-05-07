/**********************************************************************************************
 � 2020. Triad National Security, LLC. All rights reserved.
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

#include "sgh_solver.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn boundary_velocity
///
/// \brief Evolves the boundary according to a give velocity
///
/// \param The simulation mesh
/// \param An array of boundary_condition_t that contain information about BCs
/// \param A view into the nodal velocity array
/// \param The current simulation time
///
/////////////////////////////////////////////////////////////////////////////
void SGH::boundary_velocity(const mesh_t&     mesh,
    const CArrayKokkos<boundary_condition_t>& boundary,
    DCArrayKokkos<double>& node_vel,
    const double time_value)
{
    // Loop over boundary sets
    for (size_t bdy_set = 0; bdy_set < mesh.num_bdy_sets; bdy_set++) {
        // Loop over boundary nodes in a boundary set
        FOR_ALL(bdy_node_lid, 0, mesh.num_bdy_nodes_in_set.host(bdy_set), {
            // reflected (boundary array is on the device)
            if (boundary(bdy_set).type == boundary_conds::reflected) {
                // directions with type:
                // x_plane  = 0,
                // y_plane  = 1,
                // z_plane  = 2,
                size_t direction = boundary(bdy_set).direction;

                size_t bdy_node_gid = mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid);

                // Set velocity to zero in that directdion
                node_vel(1, bdy_node_gid, direction) = 0.0;
            }
            else if (boundary(bdy_set).type == boundary_conds::fixed) {
                size_t bdy_node_gid = mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid);

                for (size_t dim = 0; dim < mesh.num_dims; dim++) {
                    // Set velocity to zero
                    node_vel(1, bdy_node_gid, dim) = 0.0;
                }
            } // end if
            else if (boundary(bdy_set).type == boundary_conds::velocity) {
                size_t bdy_node_gid = mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid);

                // directions with type:
                // x_plane  = 0,
                // y_plane  = 1,
                // z_plane  = 2,
                size_t direction = boundary(bdy_set).direction;

                // Set velocity to that directdion to specified value
                // if t_end > time > t_start
                // v(t) = v0 exp(-v1*(time - time_start) )
                if (time_value >= boundary(bdy_set).hydro_bc_vel_t_start
                    && time_value <= boundary(bdy_set).hydro_bc_vel_t_end) {
                    double time_delta = time_value - boundary(bdy_set).hydro_bc_vel_t_start;

                    node_vel(1, bdy_node_gid, direction) =
                        boundary(bdy_set).hydro_bc_vel_0 *
                        exp(-boundary(bdy_set).hydro_bc_vel_1 * time_delta);
                } // end if on time
            } // end if
        }); // end for bdy_node_lid
    } // end for bdy_set

    return;
} // end boundary_velocity function


/////////////////////////////////////////////////////////////////////////////
///
/// \fn boundary_velocity
///
/// \brief Evolves the boundary according to a give velocity
///
/// \param The simulation mesh
/// \param An array of boundary_condition_t that contain information about BCs
/// \param A view into the nodal velocity array
/// \param The current simulation time
///
/////////////////////////////////////////////////////////////////////////////
void SGH::boundary_contact(const mesh_t &mesh, const node_t &nodes, const corner_t &corner)
{
    contact_bank.sort(mesh, nodes, corner);
} // end boundary_contact function