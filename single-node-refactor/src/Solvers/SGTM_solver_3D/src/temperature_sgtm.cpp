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

#include "sgtm_solver_3D.h"
#include "mesh.h"
#include "state.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn update_temperature
///
/// \brief Evolves the specific internal energy
///
/// \param The current Runge Kutta alpha value
/// \param Time step size
/// \param The simulation mesh
/// \param A view into the nodal velocity data
/// \param A view into the nodal position data
/// \param A view into the element specific internal energy
/// \param A view into the element mass
/// \param A view into the corner force data
///
/////////////////////////////////////////////////////////////////////////////
void SGTM3D::update_temperature(
    const Mesh_t& mesh,
    const DCArrayKokkos<double>& corner_flux,
    const DCArrayKokkos<double>& node_temp,
    const DCArrayKokkos<double>& node_mass,
    const DCArrayKokkos<double>& node_flux,
    const double rk_alpha,
    const double dt) const
{
    // loop over all the nodes in the mesh
    FOR_ALL(node_gid, 0, mesh.num_nodes, {
        
        double flux = 0.0;

        // node_flux(1, node_gid) = 0.0;

        // loop over all corners around the node and calculate the nodal gradient
        for (size_t corner_lid = 0; corner_lid < mesh.num_corners_in_node(node_gid); corner_lid++) {
            
            // Get corner gid
            size_t corner_gid = mesh.corners_in_node(node_gid, corner_lid);

            node_flux(1, node_gid) += corner_flux(1, corner_gid);

            flux += corner_flux(1, corner_gid);

        } // end for corner_lid

        // update the temperature
        node_temp(1, node_gid) = node_temp(0, node_gid) + rk_alpha * dt * node_flux(1, node_gid) / (node_mass(node_gid)*903.0);

    }); // end for parallel for over nodes

    return;
} // end subroutine
