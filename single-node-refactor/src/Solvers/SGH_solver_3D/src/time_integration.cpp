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

#include "sgh_solver_3D.h"
#include "mesh.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn rk_init
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
void SGH3D::rk_init(
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
    const size_t num_mat_points,
    const size_t mat_id) const
{
    // save elem quantities
    FOR_ALL(matpt_lid, 0, num_mat_points, {
        // stress is always 3D even with 2D-RZ
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 3; j++) {
                MaterialPoints_stress_n0(mat_id, matpt_lid, i, j) = MaterialPoints_stress(mat_id, matpt_lid, i, j);
            }
        }  // end for

        MaterialPoints_sie_n0(mat_id, matpt_lid) = MaterialPoints_sie(mat_id, matpt_lid);
    }); // end parallel for

    // save nodal quantities
    FOR_ALL(node_gid, 0, num_nodes, {
        for (size_t i = 0; i < num_dims; i++) {
            node_coords_n0(node_gid, i) = node_coords(node_gid, i);
            node_vel_n0(node_gid, i)    = node_vel(node_gid, i);
        }
    }); // end parallel for
    Kokkos::fence();

    return;
} // end rk_init

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_timestep
///
/// \brief This function calculates the time step by finding the shortest distance
///        between any two nodes in the mesh.
///
/// WARNING WARNING :  Only works for 3D, 8 node elements
///
/// \param Simulation mesh
/// \param View of nodal position data
/// \param View of nodal velocity data
/// \param View of element sound speed
/// \param View of element volume
///
/// REMOVE EXCESS TIME RELATED VARIABLES
///
/////////////////////////////////////////////////////////////////////////////
void SGH3D::get_timestep(Mesh_t& mesh,
                       DCArrayKokkos<double>& node_coords,
                       DCArrayKokkos<double>& node_vel,
                       DCArrayKokkos<double>& GaussPoints_vol,
                       DRaggedRightArrayKokkos<double>& MaterialPoints_sspd,
                       DRaggedRightArrayKokkos<bool>&   MaterialPoints_eroded,
                       DRaggedRightArrayKokkos<size_t>& MaterialToMeshMaps_elem,
                       size_t num_mat_elems,
                       double time_value,
                       const double graphics_time,
                       const double time_final,
                       const double dt_max,
                       const double dt_min,
                       const double dt_cfl,
                       double&      dt,
                       const double fuzz,
                       const double tiny,
                       const size_t mat_id) const
{
    // increase dt by 10%, that is the largest dt value
    dt = dt * 1.1;

    double dt_lcl;
    double min_dt_calc;
    FOR_REDUCE_MIN(mat_elem_lid, 0, num_mat_elems, dt_lcl, {
        size_t elem_gid = MaterialToMeshMaps_elem(mat_id, mat_elem_lid);

        double coords0[24];  // element coords
        ViewCArrayKokkos<double> coords(coords0, 8, 3);

        double distance0[28];  // array for holding distances between each node
        ViewCArrayKokkos<double> dist(distance0, 28);

        // Getting the coordinates of the element
        for (size_t node_lid = 0; node_lid < 8; node_lid++) {
            for (size_t dim = 0; dim < mesh.num_dims; dim++) {
                coords(node_lid, dim) = node_coords(mesh.nodes_in_elem(elem_gid, node_lid), dim);
            } // end for dim
        } // end for loop over node_lid

        // loop conditions needed for distance calculation
        size_t countA = 0;
        size_t countB = 1;
        size_t a;
        size_t b;
        size_t loop = 0;

        // Only works for 3D
        // Solving for the magnitude of distance between each node
        for (size_t i = 0; i < 28; i++) {
            a = countA;
            b = countB;

            // returns magnitude of distance between each node, 28 total options
            dist(i) = fabs(sqrt((pow((coords(b, 0) - coords(a, 0)), 2.0)
            + pow((coords(b, 1) - coords(a, 1)), 2.0)
            + pow((coords(b, 2) - coords(a, 2)), 2.0))));

            countB++;
            countA++;

            // tricky indexing
            if (countB > 7) {
                loop++;
                countB = 1 + loop;
                countA = 0;
            }
        } // end for i

        double dist_min = dist(0);

        for (int i = 0; i < 28; ++i) {
            dist_min = fmin(dist(i), dist_min);
        }

        // local dt calc based on CFL
        double dt_lcl_ = dt_cfl * dist_min / (MaterialPoints_sspd(mat_id, mat_elem_lid) + fuzz);

        if (MaterialPoints_eroded(mat_id, mat_elem_lid) == true) {
            dt_lcl_ = 1.0e32;  // a huge time step as this element doesn't exist
        }

        // make dt be in bounds
        dt_lcl_ = fmin(dt_lcl_, dt_max);    // make dt small than dt_max
        dt_lcl_ = fmax(dt_lcl_, dt_min);    // make dt larger than dt_min

        if (dt_lcl_ < dt_lcl) {
            dt_lcl = dt_lcl_;
        }
    }, min_dt_calc);  // end parallel reduction
    Kokkos::fence();

    // save the min dt
    if (min_dt_calc < dt) {
        dt = min_dt_calc;
    }

    // ensure time step hits the graphics time intervals, adding tiny to ensure dt passes graphics time
    dt = fmin(dt, fmax(fuzz,(graphics_time - time_value)) + tiny);

    // make dt hit final time, adding tiny to ensure dt passes it by a little bit
    dt = fmin(dt, time_final - time_value + tiny);
    return;
} // end get_timestep
