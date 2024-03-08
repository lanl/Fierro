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
#include "state.h"
#include "mesh.h"
#include "FEA_Module_SGH.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn update_state
///
/// \brief This calls the models to update state
///
/// \param An array of material_t that contains material specific data
/// \param The simulation mesh
/// \param A view into the nodal position array
/// \param A view into the nodal velocity array
/// \param A view into the element density array
/// \param A view into the element specific internal energy array
/// \param A view into the element pressure array
/// \param A view into the element stress array
/// \param A view into the element sound speed array
/// \param A view into the element volume array
/// \param A view into the element divergence of velocity array
/// \param A view into the element material identifier array
/// \param The current Runge Kutta integration alpha value
/// \param The current cycle index
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::update_state(const DCArrayKokkos<material_t>& material,
    const mesh_t& mesh,
    const DViewCArrayKokkos<double>& node_coords,
    const DViewCArrayKokkos<double>& node_vel,
    DViewCArrayKokkos<double>& elem_den,
    DViewCArrayKokkos<double>& elem_pres,
    DViewCArrayKokkos<double>& elem_stress,
    DViewCArrayKokkos<double>& elem_sspd,
    const DViewCArrayKokkos<double>& elem_sie,
    const DViewCArrayKokkos<double>& elem_vol,
    const DViewCArrayKokkos<double>& elem_mass,
    const DViewCArrayKokkos<size_t>& elem_mat_id,
    const double rk_alpha,
    const size_t cycle
    )
{
    const size_t rk_level = rk_num_bins - 1;
    int num_dims = num_dim;

    // loop over all the elements in the mesh
    FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
        const size_t num_dims = num_dims;
        const size_t num_nodes_in_elem = num_nodes_in_elem;

        // cut out the node_gids for this element
        ViewCArrayKokkos<size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), num_nodes_in_elem);

        // --- Density ---
        elem_den(elem_gid) = elem_mass(elem_gid) / elem_vol(elem_gid);

        size_t mat_id = elem_mat_id(elem_gid);

        // initialize elem pressure
        elem_pres(elem_gid) = 0;

        // --- Stress ---
        // hyper elastic plastic model
        if (material(mat_id).strength_type == STRENGTH_TYPE::hyper) {
            // cut out the node_gids for this element
            ViewCArrayKokkos<size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), num_nodes_in_elem);

            // --- Density ---
            elem_den(elem_gid) = elem_mass(elem_gid) / elem_vol(elem_gid);

            // corner area normals
            double area_array[24];
            ViewCArrayKokkos<double> area(area_array, num_nodes_in_elem, num_dims);

            // velocity gradient
            double vel_grad_array[9];
            ViewCArrayKokkos<double> vel_grad(vel_grad_array, num_dims, num_dims);

            // get the B matrix which are the OUTWARD corner area normals
            get_bmatrix(area,
                        elem_gid,
                        node_coords,
                        elem_node_gids,
                        rk_level);

            // --- Calculate the velocity gradient ---
            get_velgrad(vel_grad,
                        elem_node_gids,
                        node_vel,
                        area,
                        elem_vol(elem_gid),
                        elem_gid,
                        rk_level);

            // --- call strength model ---
            elem_strength(elem_gid).calc_stress(elem_pres,
                                        elem_stress,
                                        elem_gid,
                                        mat_id,
                                        eos_state_vars,
                                        strength_state_vars,
                                        eos_global_vars,
                                        strength_global_vars,
                                        elem_user_output_vars,
                                        elem_sspd,
                                        elem_den(elem_gid),
                                        elem_sie(elem_gid),
                                        vel_grad,
                                        elem_node_gids,
                                        node_coords,
                                        node_vel,
                                        elem_vol(elem_gid),
                                        dt,
                                        rk_alpha,
                                        cycle,
                                        rk_level,
                                        time_value);
        } // end logical on hyper strength model

        // --- Pressure ---
        elem_eos(elem_gid).calc_pressure(elem_pres,
                                 elem_stress,
                                 elem_gid,
                                 elem_mat_id(elem_gid),
                                 eos_state_vars,
                                 strength_state_vars,
                                 eos_global_vars,
                                 strength_global_vars,
                                 elem_user_output_vars,
                                 elem_sspd,
                                 elem_den(elem_gid),
                                 elem_sie(rk_level, elem_gid));

        // --- Sound speed ---
        elem_eos(elem_gid).calc_sound_speed(elem_pres,
                                    elem_stress,
                                    elem_gid,
                                    elem_mat_id(elem_gid),
                                    eos_state_vars,
                                    strength_state_vars,
                                    eos_global_vars,
                                    strength_global_vars,
                                    elem_user_output_vars,
                                    elem_sspd,
                                    elem_den(elem_gid),
                                    elem_sie(rk_level, elem_gid));
    }); // end parallel for
    Kokkos::fence();

    return;
} // end method to update state

/////////////////////////////////////////////////////////////////////////////
///
/// \fn update_state2D
///
/// \brief Updates the state for 2D elements
///
/// \param An array of material_t that contains material specific data
/// \param The simulation mesh
/// \param A view into the nodal position array
/// \param A view into the nodal velocity array
/// \param A view into the element density array
/// \param A view into the element specific internal energy array
/// \param A view into the element pressure array
/// \param A view into the element stress array
/// \param A view into the element sound speed array
/// \param A view into the element volume array
/// \param A view into the element divergence of velocity array
/// \param A view into the element material identifier array
/// \param The current Runge Kutta integration alpha value
/// \param The current cycle index
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::update_state2D(const DCArrayKokkos<material_t>& material,
    const mesh_t& mesh,
    const DViewCArrayKokkos<double>& node_coords,
    const DViewCArrayKokkos<double>& node_vel,
    DViewCArrayKokkos<double>& elem_den,
    DViewCArrayKokkos<double>& elem_pres,
    DViewCArrayKokkos<double>& elem_stress,
    DViewCArrayKokkos<double>& elem_sspd,
    const DViewCArrayKokkos<double>& elem_sie,
    const DViewCArrayKokkos<double>& elem_vol,
    const DViewCArrayKokkos<double>& elem_mass,
    const DViewCArrayKokkos<size_t>& elem_mat_id,
    const double rk_alpha,
    const size_t cycle
    )
{
    const size_t rk_level = rk_num_bins - 1;
    int num_dims = num_dim;

    // loop over all the elements in the mesh
    FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
        const size_t num_dims = num_dims;
        const size_t num_nodes_in_elem = num_nodes_in_elem;

        // cut out the node_gids for this element
        ViewCArrayKokkos<size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), num_nodes_in_elem);

        // --- Density ---
        elem_den(elem_gid) = elem_mass(elem_gid) / elem_vol(elem_gid);

        size_t mat_id = elem_mat_id(elem_gid);

        // initialize elem pressure
        elem_pres(elem_gid) = 0;

        // --- Stress ---
        // hyper elastic plastic model
        if (material(mat_id).strength_type == STRENGTH_TYPE::hyper) {
            // cut out the node_gids for this element
            ViewCArrayKokkos<size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), num_nodes_in_elem);

            // --- Density ---
            elem_den(elem_gid) = elem_mass(elem_gid) / elem_vol(elem_gid);

            // corner area normals
            double area_array[8];
            ViewCArrayKokkos<double> area(area_array, num_nodes_in_elem, num_dims);

            // velocity gradient
            double vel_grad_array[4];
            ViewCArrayKokkos<double> vel_grad(vel_grad_array, num_dims, num_dims);

            // get the B matrix which are the OUTWARD corner area normals
            get_bmatrix(area,
                        elem_gid,
                        node_coords,
                        elem_node_gids,
                        rk_level);

            // --- Calculate the velocity gradient ---
            get_velgrad(vel_grad,
                        elem_node_gids,
                        node_vel,
                        area,
                        elem_vol(elem_gid),
                        elem_gid,
                        rk_level);

            // --- call strength model ---
            elem_strength(elem_gid).calc_stress(elem_pres,
                                        elem_stress,
                                        elem_gid,
                                        mat_id,
                                        eos_state_vars,
                                        strength_state_vars,
                                        eos_global_vars,
                                        strength_global_vars,
                                        elem_user_output_vars,
                                        elem_sspd,
                                        elem_den(elem_gid),
                                        elem_sie(elem_gid),
                                        vel_grad,
                                        elem_node_gids,
                                        node_coords,
                                        node_vel,
                                        elem_vol(elem_gid),
                                        dt,
                                        rk_alpha,
                                        cycle,
                                        rk_level,
                                        time_value);
        } // end logical on hyper strength model

        // --- Pressure ---
        elem_eos(elem_gid).calc_pressure(elem_pres,
                                 elem_stress,
                                 elem_gid,
                                 elem_mat_id(elem_gid),
                                 eos_state_vars,
                                 strength_state_vars,
                                 eos_global_vars,
                                 strength_global_vars,
                                 elem_user_output_vars,
                                 elem_sspd,
                                 elem_den(elem_gid),
                                 elem_sie(rk_level, elem_gid));

        // --- Sound speed ---
        elem_eos(elem_gid).calc_sound_speed(elem_pres,
                                    elem_stress,
                                    elem_gid,
                                    elem_mat_id(elem_gid),
                                    eos_state_vars,
                                    strength_state_vars,
                                    eos_global_vars,
                                    strength_global_vars,
                                    elem_user_output_vars,
                                    elem_sspd,
                                    elem_den(elem_gid),
                                    elem_sie(rk_level, elem_gid));
    }); // end parallel for
    Kokkos::fence();

    return;
} // end method to update state
