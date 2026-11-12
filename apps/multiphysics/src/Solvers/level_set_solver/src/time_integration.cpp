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

#include "level_set_solver.hpp"
//#include "mesh.hpp""

/////////////////////////////////////////////////////////////////////////////
///
/// \fn rk_init
///
/// \brief This function saves the variables at rk_stage = 0, which is t_n
///

/// \param Levelset field
/// \param Levelset field
/// \param MaterialtoMeshMaps elem
/// \param Number of dimension
/// \param Number of material elems
///
/////////////////////////////////////////////////////////////////////////////
void LevelSet::rk_init(
    MPICArrayKokkos<real_t>& GaussPoints_level_set,
    MPICArrayKokkos<real_t>& GaussPoints_level_set_n0,
    DRaggedRightArrayKokkos<size_t>& elem_mat_elem,
    const size_t num_dims,
    const size_t num_mat_elems,
    const size_t mat_id) const
{
    // save elem quantities
    FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

        // get the element
        size_t elem_gid = elem_mat_elem(mat_id, mat_elem_sid);

        // Note:
        // loop gauss point, this method is for a single quadrature point element
        const size_t gauss_pt = elem_gid;

        GaussPoints_level_set_n0(gauss_pt) = GaussPoints_level_set(gauss_pt);
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
/// \param coordinates nodal data
/// \param volume velocity data
/// \param MaterialtoMeshMaps elem
///
/// REMOVE EXCESS TIME RELATED VARIABLES
///
/////////////////////////////////////////////////////////////////////////////
void LevelSet::get_timestep(
    const swage::Mesh_t& mesh,
    const Material_t& Materials,
    const MPICArrayKokkos<real_t>& node_coords,
    const DCArrayKokkos<real_t>& GaussPoints_vol,
    const DRaggedRightArrayKokkos<size_t>& elem_mat_elem,
    const size_t num_mat_elems,
    const size_t mat_id,
    const real_t time_value,
    const real_t graphics_time,
    const real_t time_final,
    const real_t dt_max,
    const real_t dt_min,
    const real_t dt_cfl,
    real_t&      dt,
    const real_t fuzz,
    const real_t tiny) const
{

    // increase dt by 10%, that is the largest dt value
    dt = dt * real_t(1.1);

    double dt_lcl;
    double min_dt_calc;
    FOR_REDUCE_MIN(mat_elem_sid, 0, num_mat_elems, dt_lcl, {
        size_t elem_gid = elem_mat_elem(mat_id, mat_elem_sid);

        real_t coords0[24];  // element coords
        ViewCArrayKokkos<real_t> coords(coords0, 8, 3);

        real_t distance0[28];  // array for holding distances between each node
        ViewCArrayKokkos<real_t> dist(distance0, 28);

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
            dist(i) = Kokkos::fabs(Kokkos::sqrt((Kokkos::pow((coords(b, 0) - coords(a, 0)), real_t(2.0))
            + Kokkos::pow((coords(b, 1) - coords(a, 1)), real_t(2.0))
            + Kokkos::pow((coords(b, 2) - coords(a, 2)), real_t(2.0)))));

            countB++;
            countA++;

            // tricky indexing
            if (countB > 7) {
                loop++;
                countB = 1 + loop;
                countA = 0;
            }
        } // end for i

        real_t dist_min = dist(0);

        for (int i = 0; i < 28; ++i) {
            dist_min = Kokkos::fmin(dist(i), dist_min);
        }

        // get the velocity
        const real_t normal_velocity = Materials.MaterialFunctions(mat_id).normal_velocity;
        const real_t curvature_velocity = Materials.MaterialFunctions(mat_id).curvature_velocity;
        
        real_t kappa = real_t(1.0)/(real_t(3.0)*dist_min); // max curvature is around 3 DeltaX, curvature = 1/R
        real_t front_velocity = real_t(1.5)+normal_velocity + curvature_velocity*kappa;
        // the 1.5 is for redistancing velocity that has a value around 1.0, so using 1.5

        // local dt calc based on CFL
        real_t dt_lcl_ = dt_cfl * dist_min / (front_velocity + fuzz);


        // make dt be in bounds
        dt_lcl_ = Kokkos::fmin(dt_lcl_, dt_max);    // make dt small than dt_max
        dt_lcl_ = Kokkos::fmax(dt_lcl_, dt_min);    // make dt larger than dt_min

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
    dt = Kokkos::fmin(dt, Kokkos::fmax(fuzz,(graphics_time - time_value)) + tiny);

    // make dt hit final time, adding tiny to ensure dt passes it by a little bit
    dt = Kokkos::fmin(dt, time_final - time_value + tiny);

    return;
} // end get_timestep

/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_timestep_2D
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
void LevelSet::get_timestep_2D(
    const swage::Mesh_t& mesh,
    const Material_t& Materials,
    const MPICArrayKokkos<real_t>& node_coords,
    const DCArrayKokkos<real_t>& GaussPoints_vol,
    const DRaggedRightArrayKokkos<size_t>& elem_mat_elem,
    const size_t num_mat_elems,
    const size_t mat_id,
    const real_t time_value,
    const real_t graphics_time,
    const real_t time_final,
    const real_t dt_max,
    const real_t dt_min,
    const real_t dt_cfl,
    real_t&      dt,
    const real_t fuzz,
    const real_t tiny) const
{
    // increase dt by 10%, that is the largest dt value
    dt = dt * real_t(1.1);

    double dt_lcl;
    double min_dt_calc;
    FOR_REDUCE_MIN(mat_elem_sid, 0, num_mat_elems, dt_lcl, {

        size_t elem_gid = elem_mat_elem(mat_id, mat_elem_sid); 

        real_t coords0[8];  // element coords
        ViewCArrayKokkos<real_t> coords(coords0, 4, 2);

        real_t distance0[6];  // array for holding distances between each node
        ViewCArrayKokkos<real_t> dist(distance0, 6);

        // Getting the coordinates of the nodes of the element
        for (size_t node_lid = 0; node_lid < 4; node_lid++) {
            for (size_t dim = 0; dim < mesh.num_dims; dim++) {
                coords(node_lid, dim) = node_coords(mesh.nodes_in_elem(elem_gid, node_lid), dim);
            } // end for dim
        } // end for loop over node_lid

        // Only works for 2D
        // Solving for the magnitude of distance between each node
        size_t count = 0;
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = i + 1; j <= 3; j++) {
                // returns magnitude of distance between each node, 6 total options
                dist(count) = Kokkos::fabs(
                        Kokkos::sqrt(Kokkos::pow((coords(i, 0) - coords(j, 0)), real_t(2.0))
                + Kokkos::pow((coords(i, 1) - coords(j, 1)), real_t(2.0)) )
                );
                count++;
            } // end for j
        } // end for i

        real_t dist_min = dist(0);

        for (int i = 0; i < 6; ++i) {
        dist_min = Kokkos::fmin(dist(i), dist_min);
        }

        // get the velocity

        const real_t normal_velocity = Materials.MaterialFunctions(mat_id).normal_velocity;
        const real_t curvature_velocity = Materials.MaterialFunctions(mat_id).curvature_velocity;

        real_t kappa = real_t(1.0)/(real_t(3.0)*dist_min); // max curvature is around 3 DeltaX, curvature = 1/R
        real_t front_velocity = real_t(1.5)+normal_velocity + curvature_velocity*kappa;
        // the 1.5 is for redistancing velocity that has a value around 1.0, so using 1.5

        // local dt calc based on CFL
        real_t dt_lcl_ = dt_cfl * dist_min / (front_velocity + fuzz);


        // make dt be in bounds
        dt_lcl_ = Kokkos::fmin(dt_lcl_, dt_max);    // make dt small than dt_max
        dt_lcl_ = Kokkos::fmax(dt_lcl_, dt_min);    // make dt larger than dt_min

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
    dt = Kokkos::fmin(dt, Kokkos::fmax(fuzz,(graphics_time - time_value)) + tiny);

    // make dt hit final time, adding tiny to ensure dt passes it by a little bit
    dt = Kokkos::fmin(dt, time_final - time_value + tiny);

    return;
} // end get_timestep_2D