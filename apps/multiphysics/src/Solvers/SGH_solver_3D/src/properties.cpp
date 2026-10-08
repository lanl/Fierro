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

#include "sgh_solver_3D.hpp"
#include "material.hpp"
//#include "mesh.hpp""
#include "geometry_new.hpp"

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
    const swage::Mesh_t&     mesh,
    const MPICArrayKokkos<double> & node_coords,
    const MPICArrayKokkos<double> & node_coords_t0,
    const MPICArrayKokkos<double> & node_vel,
    const DCArrayKokkos<double>   & GaussPoints_vel_grad,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_den,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_pres,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_stress,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_stress_n0,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_sspd,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_sie,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_volfrac,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_geo_volfrac,
    const DCArrayKokkos<double>& GaussPoints_vol,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_mass,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_eos_state_vars,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_strength_state_vars,
    const DRaggedRightArrayKokkos<bool>&   MaterialPoints_eroded,
    DRaggedRightArrayKokkos<double>& MaterialPoints_deformation_grad,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_deformation_grad_t0,
    const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
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
        FOR_ALL(mat_elem_sid, 0, num_material_elems, {
            // get elem gid
            size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);

            // get the material points for this material
            // Note, with the SGH method, they are equal
            size_t mat_point_sid = mat_elem_sid;

            // for this method, gauss point is equal to elem_gid
            size_t gauss_gid = elem_gid;

            // --- Density ---
            MaterialPoints_den(mat_id, mat_point_sid) = MaterialPoints_mass(mat_id, mat_point_sid) / 
                  (GaussPoints_vol(gauss_gid)*MaterialPoints_volfrac(mat_id, mat_point_sid)*MaterialPoints_geo_volfrac(mat_id, mat_point_sid) + 1.0e-20);

            // --- Pressure ---
            Materials.MaterialFunctions(mat_id).calc_pressure(
                                        MaterialPoints_pres,
                                        MaterialPoints_stress,
                                        mat_point_sid,
                                        mat_id,
                                        MaterialPoints_eos_state_vars,
                                        MaterialPoints_sspd,
                                        MaterialPoints_den(mat_id, mat_point_sid),
                                        MaterialPoints_sie(mat_id, mat_point_sid),
                                        Materials.eos_global_vars);

                                        
            // --- Sound Speed ---
            Materials.MaterialFunctions(mat_id).calc_sound_speed(
                                        MaterialPoints_pres,
                                        MaterialPoints_stress,
                                        mat_point_sid,
                                        mat_id,
                                        MaterialPoints_eos_state_vars,
                                        MaterialPoints_sspd,
                                        MaterialPoints_den(mat_id, mat_point_sid),
                                        MaterialPoints_sie(mat_id, mat_point_sid),
                                        MaterialPoints_deformation_grad,
                                        Materials.eos_global_vars);

        }); // end parallel for over mat elem lid
    } // if decoupled EOS
    else {
        // only calculate density as pressure and sound speed come from the coupled strength model

        // --- Density ---
        // loop over all the elements the material lives in
        FOR_ALL(mat_elem_sid, 0, num_material_elems, {
            // get elem gid
            size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);

            // get the material points for this material
            // Note, with the SGH method, they are equal
            size_t mat_point_sid = mat_elem_sid;

            // for this method, gauss point is equal to elem_gid
            size_t gauss_gid = elem_gid;

            // --- Density ---
            MaterialPoints_den(mat_id, mat_point_sid) = MaterialPoints_mass(mat_id, mat_point_sid) / (GaussPoints_vol(gauss_gid) + 1.0e-20);
        }); // end parallel for over mat elem lid
        Kokkos::fence();
    } // end if

    // --- Stress ---

    // state_based elastic plastic model
    if (Materials.MaterialEnums.host(mat_id).StrengthType == model::stateBased) {

        // ---------------------------------------
        // calculate deformation gradient for these models
        // remember: Fmodel(t) = F0*F(t), where F(t) = grad(displacment)
        get_deformation_grad(MaterialPoints_deformation_grad,
                             mesh,
                             node_coords,
                             node_coords_t0,
                             GaussPoints_vol, // remember: GaussPoint = Elem with this solver
                             elem_in_mat_elem,
                             num_material_elems, 
                             mat_id);


        // loop over all the elements the material lives in
        FOR_ALL(mat_elem_sid, 0, num_material_elems, {
            // get elem gid
            size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);

            // get the material points for this material
            // Note, with the SGH method, they are equal
            size_t mat_point_sid = mat_elem_sid;

            // for this method, gauss point is equal to elem_gid
            size_t gauss_gid = elem_gid;

            // acocunt for the reference deformation 
            double F_total[3][3];
            for (size_t i=0; i<3; i++)
            for (size_t j=0; j<3; j++)
            for (size_t k=0; k<3; k++){
                F_total[i][j] += MaterialPoints_deformation_grad_t0(i,k)*MaterialPoints_deformation_grad(k,j); 
            }

            // save the total elastic deformation gradient
            for (size_t i=0; i<3; i++)
            for (size_t j=0; j<3; j++){
                MaterialPoints_deformation_grad(i,j) = F_total[i][j];
            }

            // --- call strength model ---
            Materials.MaterialFunctions(mat_id).calc_stress(
                                        GaussPoints_vel_grad,
                                        node_coords,
                                        node_coords_t0,
                                        node_vel,
                                        mesh.nodes_in_elem,
                                        MaterialPoints_pres,
                                        MaterialPoints_stress,
                                        MaterialPoints_stress_n0,
                                        MaterialPoints_sspd,
                                        MaterialPoints_eos_state_vars,
                                        MaterialPoints_strength_state_vars,
                                        MaterialPoints_den(mat_id, mat_point_sid),
                                        MaterialPoints_sie(mat_id, mat_point_sid),
                                        MaterialPoints_deformation_grad,
                                        elem_in_mat_elem,
                                        Materials.eos_global_vars,
                                        Materials.strength_global_vars,
                                        GaussPoints_vol(elem_gid),
                                        dt,
                                        rk_alpha,
                                        time_value,
                                        cycle,
                                        mat_point_sid,
                                        mat_id,
                                        gauss_gid,
                                        elem_gid);

        }); // end parallel for over mat elem lid
    } // end if state_based strength model

    // --- mat point erosion ---
    if (Materials.MaterialEnums.host(mat_id).ErosionModels != model::noErosion) {
        // loop over all the elements the material lives in
        FOR_ALL(mat_elem_sid, 0, num_material_elems, {
            // get elem gid
            size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);

            // get the material points for this material
            // Note, with the SGH method, they are equal
            size_t mat_point_sid = mat_elem_sid;

            // for this method, gauss point is equal to elem_gid
            size_t gauss_gid = elem_gid;

            // --- Element erosion model ---
            Materials.MaterialFunctions(mat_id).erode(
                                   MaterialPoints_eroded,
                                   MaterialPoints_stress,
                                   MaterialPoints_pres(mat_id, mat_point_sid),
                                   MaterialPoints_den(mat_id, mat_point_sid),
                                   MaterialPoints_sie(mat_id, mat_point_sid),
                                   MaterialPoints_sspd(mat_id, mat_point_sid),
                                   Materials.MaterialFunctions(mat_id).erode_tension_val,
                                   Materials.MaterialFunctions(mat_id).erode_density_val,
                                   mat_point_sid,
                                   mat_id);

            // apply a void eos if mat_point is eroded
            if (MaterialPoints_eroded(mat_id, mat_point_sid)) {
                MaterialPoints_pres(mat_id, mat_point_sid) = 0.0;
                MaterialPoints_sspd(mat_id, mat_point_sid) = 1.0e-32;
                MaterialPoints_den(mat_id, mat_point_sid) = 1.0e-32;

                for (size_t i = 0; i < 3; i++) {
                    for (size_t j = 0; j < 3; j++) {
                        MaterialPoints_stress(mat_id, mat_point_sid, i, j) = 0.0;
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
    const swage::Mesh_t& mesh,
    const DCArrayKokkos<double>& GaussPoints_vol,
    const MPICArrayKokkos<double>& node_coords,
    const MPICArrayKokkos<double>& node_coords_t0,
    const MPICArrayKokkos<double>& node_vel,
    const DCArrayKokkos<double>& GaussPoints_vel_grad,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_den,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_sie,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_pres,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_stress,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_stress_n0,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_sspd,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_eos_state_vars,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_strength_state_vars,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_deformation_grad,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_deformation_grad_t0,
    const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
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
        for (size_t mat_elem_sid=0; mat_elem_sid<num_mat_elems; mat_elem_sid++){

            // get elem gid
            size_t elem_gid = elem_in_mat_elem.host(mat_id, mat_elem_sid); 


            size_t gauss_gid = elem_gid;

            // the material point index = the material elem index for a 1-point element
            size_t mat_point_sid = mat_elem_sid;


            // --- call strength model from the host side ---
            Materials.MaterialFunctions.host(mat_id).calc_stress(
                                            GaussPoints_vel_grad,
                                            node_coords,
                                            node_coords_t0,
                                            node_vel,
                                            mesh.nodes_in_elem,
                                            MaterialPoints_pres,
                                            MaterialPoints_stress,
                                            MaterialPoints_stress_n0,
                                            MaterialPoints_sspd,
                                            MaterialPoints_eos_state_vars,
                                            MaterialPoints_strength_state_vars,
                                            MaterialPoints_den(mat_id, mat_point_sid),
                                            MaterialPoints_sie(mat_id, mat_point_sid),
                                            MaterialPoints_deformation_grad,
                                            elem_in_mat_elem,
                                            Materials.eos_global_vars,
                                            Materials.strength_global_vars,
                                            GaussPoints_vol(elem_gid),
                                            dt,
                                            rk_alpha,
                                            time_value,
                                            cycle,
                                            mat_point_sid,
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
        FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

            // get elem gid
            size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid); 

            size_t gauss_gid = elem_gid;

            // the material point index = the material elem index for a 1-point element
            size_t mat_point_sid = mat_elem_sid;


            // --- call strength model ---
            Materials.MaterialFunctions(mat_id).calc_stress(
                                            GaussPoints_vel_grad,
                                            node_coords,
                                            node_coords_t0,
                                            node_vel,
                                            mesh.nodes_in_elem,
                                            MaterialPoints_pres,
                                            MaterialPoints_stress,
                                            MaterialPoints_stress_n0,
                                            MaterialPoints_sspd,
                                            MaterialPoints_eos_state_vars,
                                            MaterialPoints_strength_state_vars,
                                            MaterialPoints_den(mat_id, mat_point_sid),
                                            MaterialPoints_sie(mat_id, mat_point_sid),
                                            MaterialPoints_deformation_grad,
                                            elem_in_mat_elem,
                                            Materials.eos_global_vars,
                                            Materials.strength_global_vars,
                                            GaussPoints_vol(elem_gid),
                                            dt,
                                            rk_alpha,
                                            time_value,
                                            cycle,
                                            mat_point_sid,
                                            mat_id,
                                            gauss_gid,
                                            elem_gid);

        });  // end parallel for over elems that have the materials

    } // end if run on device

}; // end function to increment stress tensor



/////////////////////////////////////////////////////////////////////////////
///
/// \fn get_deformation_grad
///
/// \brief This function calculates the element average deformation gradient
///
/// \param deformation gradient
/// \param mesh object
/// \param node coordinates
/// \param node_t0 coordinates at initial, reference configuration
/// \param The volume of the particular element
///
/////////////////////////////////////////////////////////////////////////////
void SGH3D::get_deformation_grad(
    DRaggedRightArrayKokkos<double>& elem_deformation_grad,
    const swage::Mesh_t& mesh,
    const MPICArrayKokkos<double>& node_coords,
    const MPICArrayKokkos<double>& node_coords_t0,
    const DCArrayKokkos<double>& elem_vol,
    const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
    const size_t num_mat_elems,
    const size_t mat_id) const
{
    const size_t num_nodes_in_elem = 8;
    const size_t num_dims = 3;

    // --- loop over material elems ---
    FOR_ALL(mat_elem_sid, 0, num_mat_elems, {

        // get elem gid
        size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid); 

        // the material point storage index = the material elem index for a 1-point element
        size_t mat_point_sid = mat_elem_sid;

        // displacements in x, y, z directions at the nodes
        double u_array[num_nodes_in_elem];
        double v_array[num_nodes_in_elem];
        double w_array[num_nodes_in_elem];

        ViewCArrayKokkos<double> u(u_array, num_nodes_in_elem); // x-dir displacement component
        ViewCArrayKokkos<double> v(v_array, num_nodes_in_elem); // y-dir displacement component
        ViewCArrayKokkos<double> w(w_array, num_nodes_in_elem); // z-dir displacement component

        // cut out the node_gids for this element
        ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), num_nodes_in_elem);

        // The b_matrix are the outward corner area normals
        double b_matrix_array[24];
        ViewCArrayKokkos<double> b_matrix(b_matrix_array, num_nodes_in_elem, num_dims);
        geometry::get_bmatrix(b_matrix, elem_gid, node_coords, elem_node_gids);

        // get the vertex displacments for the elem
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
            // Get node gid
            size_t node_gid = elem_node_gids(node_lid);

            u(node_lid) = node_coords(node_gid, 0) - node_coords_t0(node_gid, 0);
            v(node_lid) = node_coords(node_gid, 1) - node_coords_t0(node_gid, 1);
            w(node_lid) = node_coords(node_gid, 2) - node_coords_t0(node_gid, 2);
        } // end for

        // --- calculate the velocity gradient terms ---
        double inverse_vol = 1.0 / elem_vol(elem_gid);
        // x-dir
        elem_deformation_grad(mat_point_sid, 0, 0) = (u(0) * b_matrix(0, 0) + u(1) * b_matrix(1, 0)
            + u(2) * b_matrix(2, 0) + u(3) * b_matrix(3, 0)
            + u(4) * b_matrix(4, 0) + u(5) * b_matrix(5, 0)
            + u(6) * b_matrix(6, 0) + u(7) * b_matrix(7, 0)) * inverse_vol;

        elem_deformation_grad(mat_point_sid, 0, 1) = (u(0) * b_matrix(0, 1) + u(1) * b_matrix(1, 1)
            + u(2) * b_matrix(2, 1) + u(3) * b_matrix(3, 1)
            + u(4) * b_matrix(4, 1) + u(5) * b_matrix(5, 1)
            + u(6) * b_matrix(6, 1) + u(7) * b_matrix(7, 1)) * inverse_vol;

        elem_deformation_grad(mat_point_sid, 0, 2) = (u(0) * b_matrix(0, 2) + u(1) * b_matrix(1, 2)
            + u(2) * b_matrix(2, 2) + u(3) * b_matrix(3, 2)
            + u(4) * b_matrix(4, 2) + u(5) * b_matrix(5, 2)
            + u(6) * b_matrix(6, 2) + u(7) * b_matrix(7, 2)) * inverse_vol;

        // y-dir
        elem_deformation_grad(mat_point_sid, 1, 0) = (v(0) * b_matrix(0, 0) + v(1) * b_matrix(1, 0)
            + v(2) * b_matrix(2, 0) + v(3) * b_matrix(3, 0)
            + v(4) * b_matrix(4, 0) + v(5) * b_matrix(5, 0)
            + v(6) * b_matrix(6, 0) + v(7) * b_matrix(7, 0)) * inverse_vol;

        elem_deformation_grad(mat_point_sid, 1, 1) = (v(0) * b_matrix(0, 1) + v(1) * b_matrix(1, 1)
            + v(2) * b_matrix(2, 1) + v(3) * b_matrix(3, 1)
            + v(4) * b_matrix(4, 1) + v(5) * b_matrix(5, 1)
            + v(6) * b_matrix(6, 1) + v(7) * b_matrix(7, 1)) * inverse_vol;
        elem_deformation_grad(mat_point_sid, 1, 2) = (v(0) * b_matrix(0, 2) + v(1) * b_matrix(1, 2)
            + v(2) * b_matrix(2, 2) + v(3) * b_matrix(3, 2)
            + v(4) * b_matrix(4, 2) + v(5) * b_matrix(5, 2)
            + v(6) * b_matrix(6, 2) + v(7) * b_matrix(7, 2)) * inverse_vol;

        // z-dir
        elem_deformation_grad(mat_point_sid, 2, 0) = (w(0) * b_matrix(0, 0) + w(1) * b_matrix(1, 0)
            + w(2) * b_matrix(2, 0) + w(3) * b_matrix(3, 0)
            + w(4) * b_matrix(4, 0) + w(5) * b_matrix(5, 0)
            + w(6) * b_matrix(6, 0) + w(7) * b_matrix(7, 0)) * inverse_vol;

        elem_deformation_grad(mat_point_sid, 2, 1) = (w(0) * b_matrix(0, 1) + w(1) * b_matrix(1, 1)
            + w(2) * b_matrix(2, 1) + w(3) * b_matrix(3, 1)
            + w(4) * b_matrix(4, 1) + w(5) * b_matrix(5, 1)
            + w(6) * b_matrix(6, 1) + w(7) * b_matrix(7, 1)) * inverse_vol;

        elem_deformation_grad(mat_point_sid, 2, 2) = (w(0) * b_matrix(0, 2) + w(1) * b_matrix(1, 2)
            + w(2) * b_matrix(2, 2) + w(3) * b_matrix(3, 2)
            + w(4) * b_matrix(4, 2) + w(5) * b_matrix(5, 2)
            + w(6) * b_matrix(6, 2) + w(7) * b_matrix(7, 2)) * inverse_vol;

    });  // end parallel for over mat elems
    Kokkos::fence();


    return;
} // end subroutine