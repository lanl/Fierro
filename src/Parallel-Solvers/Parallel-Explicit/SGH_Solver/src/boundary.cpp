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
#include "mesh.h"
#include "state.h"
#include "FEA_Module_SGH.h"
#include "Simulation_Parameters/FEA_Module/Boundary_Conditions.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn boundary_velocity
///
/// \brief Evolves the boundary according to a give velocity
///
/// \param The simulation mesh
/// \param An array of boundary_t that contain information about BCs
/// \param A view into the nodal velocity array
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::boundary_velocity(const mesh_t& mesh,
    const DCArrayKokkos<boundary_t>& boundary,
    DViewCArrayKokkos<double>& node_vel)
{
    const size_t rk_level = rk_num_bins - 1;
    int num_dims = num_dim;
    // Loop over boundary sets
    for (size_t bdy_set = 0; bdy_set < num_bdy_sets; bdy_set++)
    {
        // Loop over boundary nodes in a boundary set
        FOR_ALL_CLASS(bdy_node_lid, 0, num_bdy_nodes_in_set.host(bdy_set), {
            // reflected (boundary array is on the device)
            if (boundary(bdy_set).type == BOUNDARY_CONDITION_TYPE::reflected)
            {
                // directions with hydro_bc:
                // x_plane  = 0,
                // y_plane  = 1,
                // z_plane  = 2,
                size_t direction = boundary(bdy_set).surface.planar_surface_index();

                size_t bdy_node_gid = bdy_nodes_in_set(bdy_set, bdy_node_lid);

                // Set velocity to zero in that directdion
                node_vel(rk_level, bdy_node_gid, direction) = 0.0;
            }
            else if (boundary(bdy_set).type == BOUNDARY_CONDITION_TYPE::fixed_position)
            {
                size_t bdy_node_gid = bdy_nodes_in_set(bdy_set, bdy_node_lid);

                // debug clause
                // if(bdy_node_gid==549412) print_flag(0) = true;

                for (size_t dim = 0; dim < num_dims; dim++)
                {
                    // Set velocity to zero
                    node_vel(rk_level, bdy_node_gid, dim) = 0.0;
                }
            }
            else if (boundary(bdy_set).type == BOUNDARY_CONDITION_TYPE::velocity)
            {
                size_t bdy_node_gid = mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid);

                node_vel(rk_level, bdy_node_gid, 0) = boundary(bdy_set).u;
                node_vel(rk_level, bdy_node_gid, 1) = boundary(bdy_set).v;
                if (mesh.num_dims == 3)
                {
                    node_vel(rk_level, bdy_node_gid, 2) = boundary(bdy_set).w;
                }
            } // end if
        }); // end for bdy_node_lid
    } // end for bdy_set

    return;
} // end boundary_velocity function