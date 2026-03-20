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

    #include "tlqs_solver_3D.hpp"
    //#include "mesh.hpp""

/////////////////////////////////////////////////////////////////////////////
///
/// \fn timestep_init
///
/// \brief This function saves the variables at rk_stage = 0, which is t_n
///
/// \param View of nodal position data
/// \param View of nodal velocity data
/// \param View of element specific internal energy data
/// \param View of element stress
/// \param Number of dimension (REMOVE)
/// \param Number of elements
/// \param Number of nodes
///
/////////////////////////////////////////////////////////////////////////////
void TLQS3D::timestep_init(
    DCArrayKokkos<double>& node_coords,
    DCArrayKokkos<double>& node_coords_n0,
    DCArrayKokkos<double>& node_vel,
    DCArrayKokkos<double>& node_vel_n0,
    DRaggedRightArrayKokkos<double>& MaterialPoints_sie,
    DRaggedRightArrayKokkos<double>& MaterialPoints_sie_n0,
    DRaggedRightArrayKokkos<double>& MaterialPoints_stress,
    DRaggedRightArrayKokkos<double>& MaterialPoints_stress_n0,
    const size_t num_dims,
    const size_t num_elems,
    const size_t num_nodes,
    const size_t mat_id) const
{
    // // save elem quantities
    // FOR_ALL(matpt_lid, 0, num_mat_points, {
    //     // stress is always 3D even with 2D-RZ
    //     for (size_t i = 0; i < 3; i++) {
    //         for (size_t j = 0; j < 3; j++) {
    //             MaterialPoints_stress_n0(mat_id, matpt_lid, i, j) = MaterialPoints_stress(mat_id, matpt_lid, i, j);
    //         }
    //     }  // end for

    //     MaterialPoints_sie_n0(mat_id, matpt_lid) = MaterialPoints_sie(mat_id, matpt_lid);
    // }); // end parallel for

    // // save nodal quantities
    // FOR_ALL(node_gid, 0, num_nodes, {
    //     for (size_t i = 0; i < num_dims; i++) {
    //         node_coords_n0(node_gid, i) = node_coords(node_gid, i);
    //         node_vel_n0(node_gid, i)    = node_vel(node_gid, i);
    //     }
    // }); // end parallel for
    // Kokkos::fence();

    return;
} // end rk_init