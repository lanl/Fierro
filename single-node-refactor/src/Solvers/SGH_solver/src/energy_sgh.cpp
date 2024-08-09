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
#include "mesh.h"
#include "state.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn update_energy
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
void SGH::update_energy(const double rk_alpha,
    const double dt,
    const mesh_t& mesh,
    const DCArrayKokkos<double>& node_vel,
    const DCArrayKokkos<double>& node_coords,
    const DCArrayKokkos<double>& MaterialPoints_sie,
    const DCArrayKokkos<double>& MaterialPoints_mass,
    const DCArrayKokkos<double>& MaterialCorners_force,
    const corners_in_mat_t corners_in_mat_elem,
    const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
    const size_t num_mat_elems
    ) const
{

    // loop over all the elements in the mesh
    FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

        // get elem gid
        size_t elem_gid = MaterialToMeshMaps_elem(mat_elem_lid); 

        // the material point index = the material elem index for a 1-point element
        size_t mat_point_lid = mat_elem_lid;

        
        double MaterialPoints_power = 0.0;


        // --- tally the contribution from each corner to the element ---

        // Loop over the nodes in the element
        for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {

            // corner lid and node lid
            size_t corner_lid = node_lid;

            // Get node global id for the local node id
            size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

            // Get the corner global id for the local corner id
            //size_t corner_gid = mesh.corners_in_elem(elem_gid, corner_lid);

            // Get the material corner lid
            size_t mat_corner_lid = corners_in_mat_elem(mat_elem_lid, corner_lid);


            double node_radius = 1;
            if (mesh.num_dims == 2) {
                node_radius = node_coords(1, node_gid, 1);
            }

            // calculate the Power=F dot V for this corner
            for (size_t dim = 0; dim < mesh.num_dims; dim++) {

                double half_vel = (node_vel(1, node_gid, dim) + node_vel(0, node_gid, dim)) * 0.5;
                MaterialPoints_power += MaterialCorners_force(mat_corner_lid, dim) * node_radius * half_vel;

            } // end for dim
        } // end for node_lid

        // update the specific energy
        MaterialPoints_sie(1, mat_point_lid) = MaterialPoints_sie(0, mat_point_lid) -
                                rk_alpha * dt / MaterialPoints_mass(mat_point_lid) * MaterialPoints_power;

       
    }); // end parallel loop over the elements

    return;
} // end subroutine
