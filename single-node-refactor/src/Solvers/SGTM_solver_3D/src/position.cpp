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

/////////////////////////////////////////////////////////////////////////////
///
/// \fn update_position
///
/// \brief Updates the nodal positions based on the nodal velocity
///
/// \param Runge Kutta time integration alpha value
/// \param Time step size
/// \param Number of dimensions in the mesh (REMOVE)
/// \param Number of nodes in the mesh
/// \param View of nodal position data
/// \param View of nodal velocity data
///
/////////////////////////////////////////////////////////////////////////////
void SGTM3D::update_position(
    double rk_alpha,
    double dt,
    const size_t num_dims,
    const size_t num_nodes,
    DCArrayKokkos<double>& node_coords,
    const DCArrayKokkos<double>& node_coords_n0,
    const DCArrayKokkos<double>& node_vel,
    const DCArrayKokkos<double>& node_vel_n0) const
{
    // loop over all the nodes in the mesh
    FOR_ALL(node_gid, 0, num_nodes, {
        for (int dim = 0; dim < num_dims; dim++) {
            double half_vel = (node_vel(node_gid, dim) + node_vel_n0(node_gid, dim)) * 0.5;
            node_coords(node_gid, dim) = node_coords_n0(node_gid, dim) + rk_alpha * dt * half_vel;
        }
    }); // end parallel for over nodes
} // end subroutine