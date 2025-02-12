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
/// \param DualArrays for the nodal position 
/// \param DualArrays for the nodal velocity 
/// \param DualArrays for the material point density 
/// \param DualArrays for the material point pressure 
/// \param DualArrays for the material point stress 
/// \param DualArrays for the material point sound speed 
/// \param DualArrays for the material point specific internal energy 
/// \param DualArrays for the gauss point volume 
/// \param DualArrays for the material point mass
/// \param DualArrays for the material point eos state vars
/// \param DualArrays for the material point strength state vars
/// \param DualArrays for the material point identifier for erosion
/// \param DualArrays for the element that the material lives inside
/// \param Time step size
/// \param The current Runge Kutta integration alpha value
/// \param The number of material elems
/// \param The material id
///
/////////////////////////////////////////////////////////////////////////////
void SGH3D::update_state(
    const Material_t& Materials,
    const Mesh_t&     mesh,
    const DCArrayKokkos<double>& node_coords,
    const DCArrayKokkos<double>& node_vel,
    const DCArrayKokkos<double>& GaussPoints_vel_grad,
    const DCArrayKokkos<double>& MaterialPoints_den,
    const DCArrayKokkos<double>& MaterialPoints_pres,
    const DCArrayKokkos<double>& MaterialPoints_stress,
    const DCArrayKokkos<double>& MaterialPoints_sspd,
    const DCArrayKokkos<double>& MaterialPoints_sie,
    const DCArrayKokkos<double>& GaussPoints_vol,
    const DCArrayKokkos<double>& MaterialPoints_mass,
    const DCArrayKokkos<double>& MaterialPoints_eos_state_vars,
    const DCArrayKokkos<double>& MaterialPoints_strength_state_vars,
    const DCArrayKokkos<bool>&   MaterialPoints_eroded,
    const DCArrayKokkos<double>& MaterialPoints_shear_modulii,
    const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
    const double time_value,
    const double dt,
    const double rk_alpha,
    const size_t cycle,
    const size_t num_material_elems,
    const size_t mat_id) const
{
    const size_t num_dims = mesh.num_dims;
    const size_t num_nodes_in_elem = mesh.num_nodes_in_elem;

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
                                        MaterialPoints_eos_state_vars,
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
                                        MaterialPoints_eos_state_vars,
                                        MaterialPoints_sspd,
                                        MaterialPoints_den(mat_point_lid),
                                        MaterialPoints_sie(1, mat_point_lid),
                                        MaterialPoints_shear_modulii,
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

            // --- call strength model ---
            Materials.MaterialFunctions(mat_id).calc_stress(
                                        GaussPoints_vel_grad,
                                        node_coords,
                                        node_vel,
                                        mesh.nodes_in_elem,
                                        MaterialPoints_pres,
                                        MaterialPoints_stress,
                                        MaterialPoints_sspd,
                                        MaterialPoints_eos_state_vars,
                                        MaterialPoints_strength_state_vars,
                                        MaterialPoints_den(mat_point_lid),
                                        MaterialPoints_sie(1,mat_point_lid),
                                        MaterialPoints_shear_modulii,
                                        MaterialToMeshMaps_elem,
                                        Materials.eos_global_vars,
                                        Materials.strength_global_vars,
                                        GaussPoints_vol(elem_gid),
                                        dt,
                                        rk_alpha,
                                        time_value,
                                        cycle,
                                        mat_point_lid,
                                        mat_id,
                                        gauss_gid,
                                        elem_gid);
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
            Materials.MaterialFunctions(mat_id).erode(
                                   MaterialPoints_eroded,
                                   MaterialPoints_stress,
                                   MaterialPoints_pres(mat_point_lid),
                                   MaterialPoints_den(mat_point_lid),
                                   MaterialPoints_sie(1, mat_point_lid),
                                   MaterialPoints_sspd(mat_point_lid),
                                   Materials.MaterialFunctions(mat_id).erode_tension_val,
                                   Materials.MaterialFunctions(mat_id).erode_density_val,
                                   mat_point_lid);

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
    } // end if elem errosion

    return;
} // end method to update state



/////////////////////////////////////////////////////////////////////////////
///
/// \fn update_stress
///
/// \brief This function calculates the corner forces and the evolves stress
///
/// \param Material that contains material specific data
/// \param The simulation mesh
/// \param DualArray for gauss point vol
/// \param DualArray for nodal node coords
/// \param DualArray for nodal velocity
/// \param DualArray for mat point density
/// \param DualArray for mat point specific internal energy 
/// \param DualArray for mat point pressure 
/// \param DualArray for mat point stress 
/// \param DualArray for mat point sound speed 
/// \param DualArray for mat point eos state vars
/// \param DualArray for mat point strength state vars
/// \param DualArray for the mapping from mat lid to elem
/// \param num_mat_elems
/// \param material id
/// \param fuzz
/// \param small
/// \param time_value
/// \param Time step size
/// \param The current Runge Kutta integration alpha value
/// \param Cycle in the calculation
///
/////////////////////////////////////////////////////////////////////////////
void SGH3D::update_stress(
    const Material_t& Materials,
    const Mesh_t& mesh,
    const DCArrayKokkos<double>& GaussPoints_vol,
    const DCArrayKokkos<double>& node_coords,
    const DCArrayKokkos<double>& node_vel,
    const DCArrayKokkos<double>& GaussPoints_vel_grad,
    const DCArrayKokkos<double>& MaterialPoints_den,
    const DCArrayKokkos<double>& MaterialPoints_sie,
    const DCArrayKokkos<double>& MaterialPoints_pres,
    const DCArrayKokkos<double>& MaterialPoints_stress,
    const DCArrayKokkos<double>& MaterialPoints_sspd,
    const DCArrayKokkos<double>& MaterialPoints_eos_state_vars,
    const DCArrayKokkos<double>& MaterialPoints_strength_state_vars,
    const DCArrayKokkos<double>& MaterialPoints_shear_modulii,
    const DCArrayKokkos<size_t>& MaterialToMeshMaps_elem,
    const size_t num_mat_elems,
    const size_t mat_id,
    const double fuzz,
    const double small,
    const double time_value,
    const double dt,
    const double rk_alpha,
    const size_t cycle) const
{
    // --- Update Stress ---
    // calculate the new stress at the next rk level, if it is an increment_based model
    // increment_based strength model

    const size_t num_dims = 3;
    const size_t num_nodes_in_elem = 8;


    // ==================================================
    // launching another solver, which then calls the material model interface
    // ==================================================

    
    if (Materials.MaterialEnums.host(mat_id).StrengthRunLocation == model::host ||
        Materials.MaterialEnums.host(mat_id).StrengthRunLocation == model::dual){

        CArrayKokkos<size_t> elem_node_gids(num_nodes_in_elem); 
        for (size_t mat_elem_lid=0; mat_elem_lid<num_mat_elems; mat_elem_lid++){

            // get elem gid
            size_t elem_gid = MaterialToMeshMaps_elem.host(mat_elem_lid); 


            size_t gauss_gid = elem_gid;

            // the material point index = the material elem index for a 1-point element
            size_t mat_point_lid = mat_elem_lid;


            // --- call strength model from the host side ---
            Materials.MaterialFunctions.host(mat_id).calc_stress(
                                            GaussPoints_vel_grad,
                                            node_coords,
                                            node_vel,
                                            mesh.nodes_in_elem,
                                            MaterialPoints_pres,
                                            MaterialPoints_stress,
                                            MaterialPoints_sspd,
                                            MaterialPoints_eos_state_vars,
                                            MaterialPoints_strength_state_vars,
                                            MaterialPoints_den(mat_point_lid),
                                            MaterialPoints_sie(1,mat_point_lid),
                                            MaterialPoints_shear_modulii,
                                            MaterialToMeshMaps_elem,
                                            Materials.eos_global_vars,
                                            Materials.strength_global_vars,
                                            GaussPoints_vol(elem_gid),
                                            dt,
                                            rk_alpha,
                                            time_value,
                                            cycle,
                                            mat_point_lid,
                                            mat_id,
                                            gauss_gid,
                                            elem_gid);

        } // end serial loop over the material_elem_lids


    } // call another solver on the host, which then calls strength

    // ============================================
    // --- Device launched model is here
    // ============================================
    else {

        // --- calculate the forces acting on the nodes from the element ---
        FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

            // get elem gid
            size_t elem_gid = MaterialToMeshMaps_elem(mat_elem_lid); 

            size_t gauss_gid = elem_gid;

            // the material point index = the material elem index for a 1-point element
            size_t mat_point_lid = mat_elem_lid;


            // --- call strength model ---
            Materials.MaterialFunctions(mat_id).calc_stress(
                                            GaussPoints_vel_grad,
                                            node_coords,
                                            node_vel,
                                            mesh.nodes_in_elem,
                                            MaterialPoints_pres,
                                            MaterialPoints_stress,
                                            MaterialPoints_sspd,
                                            MaterialPoints_eos_state_vars,
                                            MaterialPoints_strength_state_vars,
                                            MaterialPoints_den(mat_point_lid),
                                            MaterialPoints_sie(1,mat_point_lid),
                                            MaterialPoints_shear_modulii,
                                            MaterialToMeshMaps_elem,
                                            Materials.eos_global_vars,
                                            Materials.strength_global_vars,
                                            GaussPoints_vol(elem_gid),
                                            dt,
                                            rk_alpha,
                                            time_value,
                                            cycle,
                                            mat_point_lid,
                                            mat_id,
                                            gauss_gid,
                                            elem_gid);

        });  // end parallel for over elems that have the materials

    } // end if run on device

}; // end function to increment stress tensor
