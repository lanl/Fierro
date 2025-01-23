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

#include "sgh_solver_rz.h"
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
void SGHRZ::rk_init_rz(DCArrayKokkos<double>& node_coords,
    DCArrayKokkos<double>& node_vel,
    DCArrayKokkos<double>& MaterialPoints_sie,
    DCArrayKokkos<double>& MaterialPoints_stress,
    const size_t num_dims,
    const size_t num_elems,
    const size_t num_nodes,
    const size_t num_mat_points) const
{

    // save elem quantities
    FOR_ALL(matpt_lid, 0, num_mat_points, {

        // stress is always 3D even with 2D-RZ
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 3; j++) {
                MaterialPoints_stress(0, matpt_lid, i, j) = MaterialPoints_stress(1, matpt_lid, i, j);
            }
        }  // end for

        MaterialPoints_sie(0, matpt_lid) = MaterialPoints_sie(1, matpt_lid);
    }); // end parallel for

    // save nodal quantities
    FOR_ALL(node_gid, 0, num_nodes, {
        for (size_t i = 0; i < num_dims; i++) {
            node_coords(0, node_gid, i) = node_coords(1, node_gid, i);
            node_vel(0, node_gid, i)    = node_vel(1, node_gid, i);
        }
    }); // end parallel for
    Kokkos::fence();

    return;
} // end rk_init



/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_timestep_rz
///
/// \brief This function calculates the time step by finding the shortest distance
///        between any two nodes in the mesh.
///
/// WARNING WARNING :  Only works for 2D, 4 node elements
///
/// \param Simulation mesh
/// \param View of nodal position data
/// \param View of nodal velocity data
/// \param View of element sound speed
/// \param View of element volume
///
/////////////////////////////////////////////////////////////////////////////
void SGHRZ::get_timestep_rz(Mesh_t& mesh,
                            DCArrayKokkos<double>& node_coords,
                            DCArrayKokkos<double>& node_vel,
                            DCArrayKokkos<double>& GaussPoints_vol,
                            DCArrayKokkos<double>& MaterialPoints_sspd,
                            DCArrayKokkos<bool>&   MaterialPoints_eroded,
                            DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
                            size_t num_mat_elems,
                            double time_value,
                            const double graphics_time,
                            const double time_final,
                            const double dt_max,
                            const double dt_min,
                            const double dt_cfl,
                            double&      dt,
                            const double fuzz) const
{
    // increase dt by 10%, that is the largest dt value
    dt = dt * 1.1;

    double dt_lcl;
    double min_dt_calc;
    FOR_REDUCE_MIN(mat_elem_lid, 0, num_mat_elems, dt_lcl, {

        size_t elem_gid = MaterialToMeshMaps_elem(mat_elem_lid); 

        double coords0[8];  // element coords
        ViewCArrayKokkos<double> coords(coords0, 4, 2);

        double distance0[6];  // array for holding distances between each node
        ViewCArrayKokkos<double> dist(distance0, 6);

        // Getting the coordinates of the nodes of the element
        for (size_t node_lid = 0; node_lid < 4; node_lid++) {
            for (size_t dim = 0; dim < mesh.num_dims; dim++) {
                coords(node_lid, dim) = node_coords(1, mesh.nodes_in_elem(elem_gid, node_lid), dim);
            } // end for dim
        } // end for loop over node_lid

        // Only works for 2D
        // Solving for the magnitude of distance between each node
        size_t count = 0;
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = i + 1; j <= 3; j++) {
                // returns magnitude of distance between each node, 6 total options
                dist(count) = fabs(
                                 sqrt(pow((coords(i, 0) - coords(j, 0)), 2.0)
                        + pow((coords(i, 1) - coords(j, 1)), 2.0) )
                    );
                count++;
            } // end for j
        } // end for i

        double dist_min = dist(0);

        for (int i = 0; i < 6; ++i) {
            dist_min = fmin(dist(i), dist_min);
        }

        // local dt calc based on CFL
        double dt_lcl_ = dt_cfl * dist_min / (MaterialPoints_sspd(mat_elem_lid) + fuzz);


        if (MaterialPoints_eroded(mat_elem_lid) == true){
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

    // ensure time step hits the graphics time intervals
    dt = fmin(dt, (graphics_time - time_value) + fuzz);

    // make dt be exact for final time
    dt = fmin(dt, time_final - time_value);

    return;
} // end get_timestep_rz
