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
#include "material.h"
#include "mesh.h"
#include "geometry_new.h"

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
void SGTM3D::update_state(
    const Material_t& Materials,
    const Mesh_t&     mesh,
    const DCArrayKokkos<double>& node_coords,
    const DCArrayKokkos<double>& node_vel,
    const DCArrayKokkos<double>& MaterialPoints_den,
    const DCArrayKokkos<double>& MaterialPoints_pres,
    const DCArrayKokkos<double>& MaterialPoints_stress,
    const DCArrayKokkos<double>& MaterialPoints_sspd,
    const DCArrayKokkos<double>& MaterialPoints_sie,
    const DCArrayKokkos<double>& GaussPoints_vol,
    const DCArrayKokkos<double>& MaterialPoints_mass,
    const DCArrayKokkos<double>& MaterialPoints_statev,
    const DCArrayKokkos<bool>&   MaterialPoints_eroded,
    const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
    const double dt,
    const double rk_alpha,
    const size_t num_material_elems,
    const size_t mat_id) const
{
    const size_t num_dims = mesh.num_dims;

    // --- pressure ---
    if (Materials.MaterialEnums.host(mat_id).EOSType == model::decoupledEOSType) {
        // loop over all the elements the material lives in
        FOR_ALL(mat_elem_lid, 0, num_material_elems, {
            // get elem gid
            size_t elem_gid = MaterialToMeshMaps_elem(mat_elem_lid);

            // get the material points for this material
            // Note, with the SGH method, they are equal
            size_t mat_point_lid = mat_elem_lid;

            // for this method, gauss point is equal to elem_gid
            size_t gauss_gid = elem_gid;

            // --- Density ---
            MaterialPoints_den(mat_point_lid) = MaterialPoints_mass(mat_point_lid) / GaussPoints_vol(gauss_gid);

            // --- Pressure ---
            Materials.MaterialFunctions(mat_id).calc_pressure(
                                        MaterialPoints_pres,
                                        MaterialPoints_stress,
                                        mat_point_lid,
                                        mat_id,
                                        MaterialPoints_statev,
                                        MaterialPoints_sspd,
                                        MaterialPoints_den(mat_point_lid),
                                        MaterialPoints_sie(1, mat_point_lid),
                                        Materials.eos_global_vars);

            // --- Sound Speed ---
            Materials.MaterialFunctions(mat_id).calc_sound_speed(
                                        MaterialPoints_pres,
                                        MaterialPoints_stress,
                                        mat_point_lid,
                                        mat_id,
                                        MaterialPoints_statev,
                                        MaterialPoints_sspd,
                                        MaterialPoints_den(mat_point_lid),
                                        MaterialPoints_sie(1, mat_point_lid),
                                        Materials.eos_global_vars);
        }); // end parallel for over mat elem lid
    } // if decoupled EOS
    else {
        // only calculate density as pressure and sound speed come from the coupled strength model

        // --- Density ---
        // loop over all the elements the material lives in
        FOR_ALL(mat_elem_lid, 0, num_material_elems, {
            // get elem gid
            size_t elem_gid = MaterialToMeshMaps_elem(mat_elem_lid);

            // get the material points for this material
            // Note, with the SGH method, they are equal
            size_t mat_point_lid = mat_elem_lid;

            // for this method, gauss point is equal to elem_gid
            size_t gauss_gid = elem_gid;

            // --- Density ---
            MaterialPoints_den(mat_point_lid) = MaterialPoints_mass(mat_point_lid) / GaussPoints_vol(gauss_gid);
        }); // end parallel for over mat elem lid
        Kokkos::fence();
    } // end if

    // --- Stress ---

    // state_based elastic plastic model
    if (Materials.MaterialEnums.host(mat_id).StrengthType == model::stateBased) {
        const size_t num_nodes_in_elem = mesh.num_nodes_in_elem;

        // loop over all the elements the material lives in
        FOR_ALL(mat_elem_lid, 0, num_material_elems, {
            // get elem gid
            size_t elem_gid = MaterialToMeshMaps_elem(mat_elem_lid);

            // get the material points for this material
            // Note, with the SGH method, they are equal
            size_t mat_point_lid = mat_elem_lid;

            // for this method, gauss point is equal to elem_gid
            size_t gauss_gid = elem_gid;

            // cut out the node_gids for this element
            ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), num_nodes_in_elem);

            // --- Density ---
            MaterialPoints_den(mat_point_lid) = MaterialPoints_mass(mat_point_lid) / GaussPoints_vol(gauss_gid);

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
            // Materials.MaterialFunctions(mat_id).calc_stress(
            //                                 MaterialPoints_pres,
            //                                 MaterialPoints_stress,
            //                                 mat_point_lid,
            //                                 mat_id,
            //                                 MaterialPoints_statev,
            //                                 MaterialPoints_sspd,
            //                                 MaterialPoints_den(mat_point_lid),
            //                                 MaterialPoints_sie(1, mat_point_lid),
            //                                 vel_grad,
            //                                 elem_node_gids,
            //                                 node_coords,
            //                                 node_vel,
            //                                 GaussPoints_vol(gauss_gid),
            //                                 dt,
            //                                 rk_alpha);
        }); // end parallel for over mat elem lid
    } // end if state_based strength model

    // --- mat point erosion ---
    if (Materials.MaterialEnums.host(mat_id).ErosionModels != model::noErosion) {
        // loop over all the elements the material lives in
        FOR_ALL(mat_elem_lid, 0, num_material_elems, {
            // get elem gid
            size_t elem_gid = MaterialToMeshMaps_elem(mat_elem_lid);

            // get the material points for this material
            // Note, with the SGH method, they are equal
            size_t mat_point_lid = mat_elem_lid;

            // for this method, gauss point is equal to elem_gid
            size_t gauss_gid = elem_gid;

            // --- Element erosion model ---
            // Materials.MaterialFunctions(mat_id).erode(
            //                        MaterialPoints_eroded,
            //                        MaterialPoints_stress,
            //                        MaterialPoints_pres(mat_point_lid),
            //                        MaterialPoints_den(mat_point_lid),
            //                        MaterialPoints_sie(1, mat_point_lid),
            //                        MaterialPoints_sspd(mat_point_lid),
            //                        Materials.MaterialFunctions(mat_id).erode_tension_val,
            //                        Materials.MaterialFunctions(mat_id).erode_density_val,
            //                        mat_point_lid);

            // apply a void eos if mat_point is eroded
            if (MaterialPoints_eroded(mat_point_lid)) {
                MaterialPoints_pres(mat_point_lid) = 0.0;
                MaterialPoints_sspd(mat_point_lid) = 1.0e-32;
                MaterialPoints_den(mat_point_lid) = 1.0e-32;

                for (size_t i = 0; i < 3; i++) {
                    for (size_t j = 0; j < 3; j++) {
                        MaterialPoints_stress(1, mat_point_lid, i, j) = 0.0;
                    }
                }  // end for i,j
            } // end if on eroded
        }); // end parallel for
    } // end if elem erosion

    return;
} // end method to update state
