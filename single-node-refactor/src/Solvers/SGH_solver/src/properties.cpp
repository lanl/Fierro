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
/// \fn update_state
///
/// \brief This calls the models to update state
///
/// \param Material that contains material specific data
/// \param The simulation mesh
/// \param A view into the nodal position array
/// \param A view into the nodal velocity array
/// \param A view into the element density array
/// \param A view into the element pressure array
/// \param A view into the element stress array
/// \param A view into the element sound speed array
/// \param A view into the element specific internal energy array
/// \param A view into the element volume array
/// \param A view into the element mass
/// \param A view into the element material identifier array
/// \param A view into the element state variables
/// \param Time step size
/// \param The current Runge Kutta integration alpha value
///
/////////////////////////////////////////////////////////////////////////////
void SGH::update_state(const Material_t& Materials,
    const mesh_t& mesh,
    const DCArrayKokkos<double>& node_coords,
    const MPIArrayKokkos<double>& node_vel,
    DCArrayKokkos<double>& MaterialPoints_den,
    DCArrayKokkos<double>& MaterialPoints_pres,
    DCArrayKokkos<double>& MaterialPoints_stress,
    DCArrayKokkos<double>& MaterialPoints_sspd,
    const DCArrayKokkos<double>& MaterialPoints_sie,
    const DCArrayKokkos<double>& GaussPoints_vol,
    const DCArrayKokkos<double>& MaterialPoints_mass,
    const DCArrayKokkos<size_t>& GaussPoints_mat_id,
    const DCArrayKokkos<double>& MaterialPoints_statev,
    const DCArrayKokkos<bool>&   GaussPoints_eroded,
    const double dt,
    const double rk_alpha) const
{
    // loop over all the elements in the mesh
    FOR_ALL(elem_gid, 0, mesh.num_elems, {
        const size_t num_dims = mesh.num_dims;
        const size_t num_nodes_in_elem = mesh.num_nodes_in_elem;

        // cut out the node_gids for this element
        ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), num_nodes_in_elem);

        // --- Density ---
        MaterialPoints_den(elem_gid) = MaterialPoints_mass(elem_gid) / GaussPoints_vol(elem_gid);

        size_t mat_id = GaussPoints_mat_id(elem_gid);

        // --- Stress ---
        // state_based elastic plastic model
        if (Materials.MaterialEnums(mat_id).StrengthType == model::stateBased) {
            // cut out the node_gids for this element
            ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), num_nodes_in_elem);

            // --- Density ---
            MaterialPoints_den(elem_gid) = MaterialPoints_mass(elem_gid) / GaussPoints_vol(elem_gid);

            // corner area normals
            double area_array[24];
            ViewCArrayKokkos<double> area(area_array, num_nodes_in_elem, num_dims);

            // velocity gradient
            double vel_grad_array[9];
            ViewCArrayKokkos<double> vel_grad(vel_grad_array, num_dims, num_dims);

            // get the B matrix which are the OUTWARD corner area normals
            geometry::get_bmatrix(area, elem_gid, node_coords, elem_node_gids);

            // --- Calculate the velocity gradient ---
            get_velgrad(vel_grad,
                        elem_node_gids,
                        node_vel,
                        area,
                        GaussPoints_vol(elem_gid),
                        elem_gid);

            // --- call strength model ---
            Materials.MaterialFunctions(mat_id).calc_stress(
                                         MaterialPoints_pres,
                                         MaterialPoints_stress,
                                         elem_gid,
                                         mat_id,
                                         MaterialPoints_statev,
                                         MaterialPoints_sspd,
                                         MaterialPoints_den(elem_gid),
                                         MaterialPoints_sie(elem_gid),
                                         vel_grad,
                                         elem_node_gids,
                                         node_coords,
                                         node_vel,
                                         GaussPoints_vol(elem_gid),
                                         dt,
                                         rk_alpha);
        } // end logical on state_based strength model

        // apply the element erosion model
        //if (material(mat_id).erosion_type == model::erosion) {
        //    // starting simple, but in the future call an erosion model
        //    if (MaterialPoints_pres(elem_gid) <= material(mat_id).erode_tension_val
        //        || MaterialPoints_den(elem_gid) <= material(mat_id).erode_density_val) {
        //        GaussPoints_mat_id(elem_gid) = material(mat_id).void_mat_id;
        //
        //        GaussPoints_eroded(elem_gid) = true;
        //    } // end if
        //} // end if
        if (Materials.MaterialFunctions(mat_id).erode != NULL) {
            

            // --- Element erosion model ---
            Materials.MaterialFunctions(mat_id).erode(
                                   MaterialPoints_pres,
                                   MaterialPoints_stress,
                                   GaussPoints_eroded,
                                   GaussPoints_mat_id,
                                   elem_gid,
                                   Materials.MaterialFunctions(mat_id).void_mat_id,
                                   Materials.MaterialFunctions(mat_id).erode_tension_val,
                                   Materials.MaterialFunctions(mat_id).erode_density_val,
                                   MaterialPoints_sspd,
                                   MaterialPoints_den,
                                   MaterialPoints_sie(1, elem_gid));

        } // end if

        if (Materials.MaterialEnums(mat_id).EOSType == model::decoupledEOSType) {

            // --- Pressure ---
            Materials.MaterialFunctions(mat_id).calc_pressure(
                                           MaterialPoints_pres,
                                           MaterialPoints_stress,
                                           elem_gid,
                                           GaussPoints_mat_id(elem_gid),
                                           MaterialPoints_statev,
                                           MaterialPoints_sspd,
                                           MaterialPoints_den(elem_gid),
                                           MaterialPoints_sie(1, elem_gid),
                                           Materials.eos_global_vars);   
            // --- Sound Speed ---                               
            Materials.MaterialFunctions(mat_id).calc_sound_speed(
                                              MaterialPoints_pres,
                                              MaterialPoints_stress,
                                              elem_gid,
                                              GaussPoints_mat_id(elem_gid),
                                              MaterialPoints_statev,
                                              MaterialPoints_sspd,
                                              MaterialPoints_den(elem_gid),
                                              MaterialPoints_sie(1, elem_gid),
                                              Materials.eos_global_vars);
        }
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
/// \param A view into the element pressure array
/// \param A view into the element stress array
/// \param A view into the element sound speed array
/// \param A view into the element specific internal energy array
/// \param A view into the element volume array
/// \param A view into the element mass
/// \param A view into the element material identifier array
/// \param A view into the element state variables
/// \param Time step size
/// \param The current Runge Kutta integration alpha value
///
/////////////////////////////////////////////////////////////////////////////
void SGH::update_state2D(const Material_t& Materials,
    const mesh_t& mesh,
    const DCArrayKokkos<double>& node_coords,
    const MPIArrayKokkos<double>& node_vel,
    DCArrayKokkos<double>& MaterialPoints_den,
    DCArrayKokkos<double>& MaterialPoints_pres,
    DCArrayKokkos<double>& MaterialPoints_stress,
    DCArrayKokkos<double>& MaterialPoints_sspd,
    const DCArrayKokkos<double>& MaterialPoints_sie,
    const DCArrayKokkos<double>& GaussPoints_vol,
    const DCArrayKokkos<double>& MaterialPoints_mass,
    const DCArrayKokkos<size_t>& GaussPoints_mat_id,
    const DCArrayKokkos<double>& MaterialPoints_statev,
    const double dt,
    const double rk_alpha) const
{
    // loop over all the elements in the mesh
    FOR_ALL(elem_gid, 0, mesh.num_elems, {
        const size_t num_dims = mesh.num_dims;
        const size_t num_nodes_in_elem = mesh.num_nodes_in_elem;

        // cut out the node_gids for this element
        ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), num_nodes_in_elem);

        // --- Density ---
        MaterialPoints_den(elem_gid) = MaterialPoints_mass(elem_gid) / GaussPoints_vol(elem_gid);

        size_t mat_id = GaussPoints_mat_id(elem_gid);

        // --- Stress ---
        // state_based elastic plastic model
        if (Materials.MaterialEnums(mat_id).StrengthType == model::stateBased) {
            // cut out the node_gids for this element
            ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), num_nodes_in_elem);

            // --- Density ---
            MaterialPoints_den(elem_gid) = MaterialPoints_mass(elem_gid) / GaussPoints_vol(elem_gid);

            // corner area normals
            double area_array[8];
            ViewCArrayKokkos<double> area(area_array, num_nodes_in_elem, num_dims);

            // velocity gradient
            double vel_grad_array[4];
            ViewCArrayKokkos<double> vel_grad(vel_grad_array, num_dims, num_dims);

            // get the B matrix which are the OUTWARD corner area normals
            geometry::get_bmatrix(area, elem_gid, node_coords, elem_node_gids);

            // --- Calculate the velocity gradient ---
            get_velgrad(vel_grad,
                        elem_node_gids,
                        node_vel,
                        area,
                        GaussPoints_vol(elem_gid),
                        elem_gid);

            // --- call strength model ---
            // Material.MaterialFunctions(mat_id).strength_model(MaterialPoints_pres,
            //                                 MaterialPoints_stress,
            //                                 elem_gid,
            //                                 mat_id,
            //                                 MaterialPoints_statev,
            //                                 MaterialPoints_sspd,
            //                                 MaterialPoints_den(elem_gid),
            //                                 MaterialPoints_sie(elem_gid),
            //                                 vel_grad,
            //                                 elem_node_gids,
            //                                 node_coords,
            //                                 node_vel,
            //                                 GaussPoints_vol(elem_gid),
            //                                 dt,
            //                                 rk_alpha);
        } // end logical on state_based strength model

        // --- Erosion ---
        // apply the element erosion model
        //if (Materials.MaterialFunctions(mat_id).erode != NULL) {
        //
        //    // --- Element erosion model ---
        //    material.MaterialFunctions(mat_id).erode(MaterialPoints_pres,
        //                           MaterialPoints_stress,
        //                           GaussPoints_eroded,
        //                           GaussPoints_mat_id,
        //                           elem_gid,
        //                           material(mat_id).void_mat_id,
        //                           material(mat_id).erode_tension_val,
        //                           material(mat_id).erode_density_val,
        //                           MaterialPoints_sspd,
        //                           MaterialPoints_den,
        //                           MaterialPoints_sie(1, elem_gid));
        //} // end if

        // --- Pressure ---
        if (Materials.MaterialEnums(mat_id).EOSType == model::decoupledEOSType) {

            // --- Pressure ---
            Materials.MaterialFunctions(mat_id).calc_pressure(
                                           MaterialPoints_pres,
                                           MaterialPoints_stress,
                                           elem_gid,
                                           GaussPoints_mat_id(elem_gid),
                                           MaterialPoints_statev,
                                           MaterialPoints_sspd,
                                           MaterialPoints_den(elem_gid),
                                           MaterialPoints_sie(1, elem_gid),
                                           Materials.eos_global_vars);   
            // --- Sound Speed ---                               
            Materials.MaterialFunctions(mat_id).calc_sound_speed(
                                              MaterialPoints_pres,
                                              MaterialPoints_stress,
                                              elem_gid,
                                              GaussPoints_mat_id(elem_gid),
                                              MaterialPoints_statev,
                                              MaterialPoints_sspd,
                                              MaterialPoints_den(elem_gid),
                                              MaterialPoints_sie(1, elem_gid),
                                              Materials.eos_global_vars);
        }
    }); // end parallel for
    Kokkos::fence();

    return;
} // end method to update state
