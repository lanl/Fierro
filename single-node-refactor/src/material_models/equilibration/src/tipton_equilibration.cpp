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

#include "tipton_equilibration.hpp"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn TiptonEquilibrationModel
///
/// \brief Using Tipton's model to relax multiple materials in an element
///        to be in equilibrium. The model has functions to equilibrate
///        material volume fractions and geometric volume fractions. The
///        material volfrac is only applied to materials in an element with 
///        volfrac < 1.  The geometric volfrac is only applied to materials 
///        in an element with geo_volfrac < 1.
///
/// \param Material object
/// \param Mesh object
/// \param State object
/// \param GuassPoint average pressure array
/// \param GuassPoint denominater used to calculate average pressure 
/// \param GuassPoint smallest volume fraction of all materials at that point 
/// \param GuassPoint limiter on changing volume fraction 
/// \param dt is time step
/// \param rk_alpha is coefficient used in the RungeKutta time integration
/// \param fuzz is a very small number, like 1e-16
/// \param small is a small number, like 1e-10
///
/////////////////////////////////////////////////////////////////////////////

namespace TiptonEquilibrationModel {
    


    // ----------------------------------------------------------
    // This is the Tipton material pt equilibration model
    // ----------------------------------------------------------
    void mat_equilibration(const Material_t& Materials, 
                      const Mesh_t& mesh, 
                      State_t& State,
                      CArrayKokkos <double>& GaussPoint_pres,
                      CArrayKokkos <double>& GaussPoint_pres_denominator,
                      CArrayKokkos <double>& GaussPoint_volfrac_min,
                      CArrayKokkos <double>& GaussPoint_volfrac_limiter,
                      const double dt,
                      const double rk_alpha,
                      const double fuzz,
                      const double small)
    {

        const double my_fuzz = 1.e-13;

        const size_t num_mats = Materials.num_mats;

        GaussPoint_pres.set_values(0.0);
        GaussPoint_pres_denominator.set_values(0.0);
        GaussPoint_volfrac_min.set_values(1.0);
        GaussPoint_volfrac_limiter.set_values(1.0);

        // calculate weigted average pressure at gauss points
        for(size_t mat_id = 0; mat_id < num_mats; mat_id++){

            // build numerator and denominator
            build_gauss_point_averages( 
                mesh,
                GaussPoint_pres,
                GaussPoint_pres_denominator,
                GaussPoint_volfrac_min,
                State.MaterialPoints.volfrac,  // material volfrac
                State.MaterialPoints.pres,
                State.MaterialPoints.den,
                State.MaterialPoints.sspd,
                State.MaterialToMeshMaps.elem,
                State.points_in_mat_elem,
                dt,
                rk_alpha,
                my_fuzz,
                State.MaterialToMeshMaps.num_material_elems.host(mat_id),
                mat_id);

        } // end for mat_id

        // calculate the average pressure
        calc_gauss_point_averages( 
            mesh,
            GaussPoint_pres,
            GaussPoint_pres_denominator,
            my_fuzz);



        // calculate volfrac change and limiter on the change
        for(size_t mat_id = 0; mat_id < num_mats; mat_id++){

            calc_volfrac_change (
                mesh,
                GaussPoint_pres,
                GaussPoint_pres_denominator,
                GaussPoint_volfrac_min,
                GaussPoint_volfrac_limiter,
                State.GaussPoints.vel_grad,
                State.GaussPoints.vol,
                State.MaterialPoints.volfrac,       // material volfrac
                State.MaterialPoints.delta_volfrac, // the change in the material volfrac
                State.MaterialPoints.pres,
                State.MaterialPoints.den,
                State.MaterialPoints.sspd,
                State.MaterialPoints.mass,
                State.MaterialToMeshMaps.elem,
                State.points_in_mat_elem,
                dt,
                rk_alpha,
                my_fuzz,
                Materials.equilibration_global_vars,
                Materials.num_equilibration_global_vars,
                State.MaterialToMeshMaps.num_material_elems.host(mat_id),
                mat_id);

        } // end for mat_id


        // calculate volfrac and energy change
        for(size_t mat_id = 0; mat_id < num_mats; mat_id++){

            update_volfrac_sie (
                mesh,
                GaussPoint_pres,
                GaussPoint_volfrac_limiter,
                State.GaussPoints.vel_grad,
                State.GaussPoints.vol,
                State.MaterialPoints.volfrac,       // material volfrac
                State.MaterialPoints.delta_volfrac, // the change in material volfrac
                State.MaterialPoints.sie,
                State.MaterialPoints.mass,
                State.MaterialToMeshMaps.elem,
                State.points_in_mat_elem,
                dt,
                rk_alpha,
                my_fuzz,
                State.MaterialToMeshMaps.num_material_elems.host(mat_id),
                mat_id);

        } // end for mat_id

    } // end equilibration function


    // ----------------------------------------------------------
    // This is the Tipton geometric material pt equilibration model
    // ----------------------------------------------------------
    void geo_equilibration(const Material_t& Materials, 
        const Mesh_t& mesh, 
        State_t& State,
        CArrayKokkos <double>& GaussPoint_pres,
        CArrayKokkos <double>& GaussPoint_pres_denominator,
        CArrayKokkos <double>& GaussPoint_volfrac_min,
        CArrayKokkos <double>& GaussPoint_volfrac_limiter,
        const double dt,
        const double rk_alpha,
        const double fuzz,
        const double small)
    {

        const size_t num_mats = Materials.num_mats;

        GaussPoint_pres.set_values(0.0);
        GaussPoint_pres_denominator.set_values(0.0);
        GaussPoint_volfrac_min.set_values(1.0);
        GaussPoint_volfrac_limiter.set_values(1.0);

        // calculate weigted average pressure at gauss points
        for(size_t mat_id = 0; mat_id < num_mats; mat_id++){

            // build numerator and denominator
            build_gauss_point_averages( 
                mesh,
                GaussPoint_pres,
                GaussPoint_pres_denominator,
                GaussPoint_volfrac_min,
                State.MaterialPoints.geo_volfrac,  // geo_volfrac
                State.MaterialPoints.pres,
                State.MaterialPoints.den,
                State.MaterialPoints.sspd,
                State.MaterialToMeshMaps.elem,
                State.points_in_mat_elem,
                dt,
                rk_alpha,
                fuzz,
                State.MaterialToMeshMaps.num_material_elems.host(mat_id),
                mat_id);

        } // end for mat_id

        // calculate the average pressure
        calc_gauss_point_averages( 
            mesh,
            GaussPoint_pres,
            GaussPoint_pres_denominator,
            fuzz);



        // calculate geo_volfrac change and limiter on the change
        for(size_t mat_id = 0; mat_id < num_mats; mat_id++){

            calc_volfrac_change (
                mesh,
                GaussPoint_pres,
                GaussPoint_pres_denominator,
                GaussPoint_volfrac_min,
                GaussPoint_volfrac_limiter,
                State.GaussPoints.vel_grad,
                State.GaussPoints.vol,
                State.MaterialPoints.geo_volfrac,       // geo_volfrac 
                State.MaterialPoints.delta_geo_volfrac, // change in the geo_volfrac 
                State.MaterialPoints.pres,
                State.MaterialPoints.den,
                State.MaterialPoints.sspd,
                State.MaterialPoints.mass,
                State.MaterialToMeshMaps.elem,
                State.points_in_mat_elem,
                dt,
                rk_alpha,
                fuzz,
                Materials.equilibration_global_vars,
                Materials.num_equilibration_global_vars,
                State.MaterialToMeshMaps.num_material_elems.host(mat_id),
                mat_id);

        } // end for mat_id


        // calculate volfrac and energy change
        for(size_t mat_id = 0; mat_id < num_mats; mat_id++){

            update_volfrac_sie (
                mesh,
                GaussPoint_pres,
                GaussPoint_volfrac_limiter,
                State.GaussPoints.vel_grad,
                State.GaussPoints.vol,
                State.MaterialPoints.geo_volfrac,       // geo_volfrac 
                State.MaterialPoints.delta_geo_volfrac, // change in the geo_volfrac
                State.MaterialPoints.sie,
                State.MaterialPoints.mass,
                State.MaterialToMeshMaps.elem,
                State.points_in_mat_elem,
                dt,
                rk_alpha,
                fuzz,
                State.MaterialToMeshMaps.num_material_elems.host(mat_id),
                mat_id);

        } // end for mat_id

    } // end geo_equilibration function


    // ------------------------
    //     Helper functions   
    // ------------------------

    void build_gauss_point_averages (
        const Mesh_t& mesh,
        const CArrayKokkos<double>& GaussPoint_pres,
        const CArrayKokkos<double>& GaussPoint_pres_denominator,
        const CArrayKokkos<double>& GaussPoint_volfrac_min,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_volfrac,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_pres,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_den,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_sspd,
        const DRaggedRightArrayKokkos<size_t>& MaterialToMeshMaps_elem,
        const points_in_mat_t& points_in_mat_elem,
        const double dt,
        const double rk_alpha,
        const double fuzz,
        const size_t num_mat_elems,
        const size_t mat_id)
    {

        // loop over all ellements the material lives in
        FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

            // get elem gid for this material at this lid
            size_t elem_gid = MaterialToMeshMaps_elem(mat_id, mat_elem_lid);

            // loop over gauss points in this element
            for (size_t gauss_pt_lid = 0; gauss_pt_lid < mesh.num_leg_gauss_in_elem; gauss_pt_lid++){

                // get the gauss gid for this point in the element
                size_t gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_pt_lid);

                // get the mat_gauss_pt_storage_lid
                size_t mat_point_storage_lid = points_in_mat_elem(mat_elem_lid, gauss_pt_lid);

                // only do pressure relaxation on materials that have volfrac<1
                if (MaterialPoints_volfrac(mat_id, mat_point_storage_lid )<1.0-fuzz){

                    // the smallest volume fraction is used for limiting volfrac changes
                    GaussPoint_volfrac_min(gauss_gid) = fmin(MaterialPoints_volfrac(mat_id, mat_point_storage_lid), GaussPoint_volfrac_min(gauss_gid)); // min volfrac

                    // calculate average pressure
                    // GaussPoint_avg_press = sum_i(volfrac_i/K_i P_i)/ sum_i(volfrac_i/K_i)
                    // K_i = rho*c^2
                    const double bulk_mod = MaterialPoints_den(mat_id, mat_point_storage_lid)*MaterialPoints_sspd(mat_id, mat_point_storage_lid)*MaterialPoints_sspd(mat_id, mat_point_storage_lid);
                    const double R = MaterialPoints_volfrac(mat_id, mat_point_storage_lid)/(bulk_mod+fuzz); // ratio
                    GaussPoint_pres(gauss_gid) += R * MaterialPoints_pres(mat_id, mat_point_storage_lid);
                    GaussPoint_pres_denominator(gauss_gid) += R; // defined as R_bar

                } // end if volfrac<1
                else {
                    // single material
                    GaussPoint_pres(gauss_gid) = MaterialPoints_pres(mat_id, mat_point_storage_lid);
                    GaussPoint_pres_denominator(gauss_gid) = 1.0;
                } // end if

            } // end for gauss point loop 

        }); // end parallel loop over all material elems in the mesh
        Kokkos::fence();

        return;

    } // end build average fields function



    void  calc_gauss_point_averages( 
        const Mesh_t& mesh,
        const CArrayKokkos<double>&  GaussPoint_pres,
        const CArrayKokkos<double>&  GaussPoint_pres_denominator,
        const double fuzz){

        // loop over all ellements in the mesh
        FOR_ALL(elem_gid, 0, mesh.num_elems, {

            // loop over gauss points in this element
            for (size_t gauss_pt_lid = 0; gauss_pt_lid < mesh.num_leg_gauss_in_elem; gauss_pt_lid++){

                // get the gauss gid for this point in the element
                size_t gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_pt_lid);

                GaussPoint_pres(gauss_gid) /= (GaussPoint_pres_denominator(gauss_gid)+fuzz); 

            } // end for

         });
         Kokkos::fence();

         return;

    } // end function




    void calc_volfrac_change(
        const Mesh_t& mesh,
        const CArrayKokkos<double>& GaussPoint_pres,
        const CArrayKokkos<double>& GaussPoint_pres_denominator,
        const CArrayKokkos<double>& GaussPoint_volfrac_min,
        const CArrayKokkos<double>& GaussPoint_volfrac_limiter,
        const DCArrayKokkos<double>& GaussPoint_vel_grad,
        const DCArrayKokkos<double>& GaussPoint_vol,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_volfrac,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_delta_volfrac,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_pres,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_den,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_sspd,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_mass,
        const DRaggedRightArrayKokkos<size_t>& MaterialToMeshMaps_elem,
        const points_in_mat_t& points_in_mat_elem,
        const double dt,
        const double rk_alpha,
        const double fuzz,
        const CArrayKokkos<double> &equilibration_global_vars,
        const size_t num_global_vars,
        const size_t num_mat_elems,
        const size_t mat_id)
    {

        // loop over all ellements the material lives in
        FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

            // get elem gid for this material at this lid
            size_t elem_gid = MaterialToMeshMaps_elem(mat_id, mat_elem_lid);

            // loop over gauss points in this element
            for (size_t gauss_pt_lid = 0; gauss_pt_lid < mesh.num_leg_gauss_in_elem; gauss_pt_lid++){

                // get the gauss gid for this point in the element
                size_t gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_pt_lid);

                // get the mat_gauss_pt_storage_lid
                size_t mat_point_storage_lid = points_in_mat_elem(mat_elem_lid, gauss_pt_lid);

                // divergence at the Gauss point
                double div = 0.0;
                for (size_t dim=0; dim<mesh.num_dims; dim++){
                    div += GaussPoint_vel_grad(gauss_gid, dim, dim);
                }


                // only do pressure relaxation on materials that have volfrac<1
                if (MaterialPoints_volfrac(mat_id, mat_point_storage_lid)<1.0-fuzz){

                    // volume fraction change, unlimited
                    const double bulk_mod = MaterialPoints_den(mat_id, mat_point_storage_lid)*MaterialPoints_sspd(mat_id, mat_point_storage_lid)*MaterialPoints_sspd(mat_id, mat_point_storage_lid);
                    const double R = MaterialPoints_volfrac(mat_id, mat_point_storage_lid)/(bulk_mod+fuzz);
                    const double R_bar = GaussPoint_pres_denominator(gauss_gid);
                    const double term_1 = R*(MaterialPoints_pres(mat_id, mat_point_storage_lid) - GaussPoint_pres(gauss_gid))/(equilibration_global_vars(1)+fuzz);  // 0.25 is used in V. Chiravalle's IASSD paper
                    const double term_2 = (R/(R_bar+fuzz) - MaterialPoints_volfrac(mat_id, mat_point_storage_lid))*div*rk_alpha*dt;

                    MaterialPoints_delta_volfrac(mat_id, mat_point_storage_lid) = term_1 + term_2;

                    // limit volume fraction change, ensuring positive volumes after relaxation
                    // coef*delta_vol <= param*smallest_volume
                    const double param = fmin(1.0, fmax(0.0, equilibration_global_vars(0)));  // ensuring it is bounded between 0:1

                    GaussPoint_volfrac_limiter(gauss_gid) = fmin(GaussPoint_volfrac_limiter(gauss_gid), 
                                                                 param*GaussPoint_volfrac_min(gauss_gid)/(fabs(MaterialPoints_delta_volfrac(mat_id, mat_point_storage_lid)) + fuzz));
                    

                } // end if volfrac<1
                else {
                    // single material 
                    MaterialPoints_delta_volfrac(mat_id, mat_point_storage_lid) = 0.0; 
                    GaussPoint_volfrac_limiter(mat_id, mat_point_storage_lid) = 1.0;  
                } // end check if multiple materials are present

            } // end for gauss point loop 

        }); // end parallel loop over all material elems in the mesh
        Kokkos::fence();

        return;

    } // end function


    void update_volfrac_sie(
        const Mesh_t& mesh,
        const CArrayKokkos<double>& GaussPoint_pres,
        const CArrayKokkos <double>& GaussPoint_volfrac_limiter,
        const DCArrayKokkos<double>& GaussPoint_vel_grad,
        const DCArrayKokkos<double>& GaussPoint_vol,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_volfrac,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_delta_volfrac,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_sie,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_mass,
        const DRaggedRightArrayKokkos<size_t>& MaterialToMeshMaps_elem,
        const points_in_mat_t& points_in_mat_elem,
        const double dt,
        const double rk_alpha,
        const double fuzz,
        const size_t num_mat_elems,
        const size_t mat_id)
    {

        // loop over all ellements the material lives in
        FOR_ALL(mat_elem_lid, 0, num_mat_elems, {

            // get elem gid for this material at this lid
            size_t elem_gid = MaterialToMeshMaps_elem(mat_id, mat_elem_lid);

            // loop over gauss points in this element
            for (size_t gauss_pt_lid = 0; gauss_pt_lid < mesh.num_leg_gauss_in_elem; gauss_pt_lid++){

                // get the gauss gid for this point in the element
                size_t gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_pt_lid);

                // get the mat_gauss_pt_storage_lid
                size_t mat_point_storage_lid = points_in_mat_elem(mat_elem_lid, gauss_pt_lid);

                // divergence at the Gauss point
                double div = 0.0;
                for (size_t dim=0; dim<mesh.num_dims; dim++){
                    div += GaussPoint_vel_grad(gauss_gid, dim, dim);
                }

                const double delta_volfrac = GaussPoint_volfrac_limiter(gauss_gid)*MaterialPoints_delta_volfrac(mat_id, mat_point_storage_lid);

                // calculating volume fraction change 
                double volfrac_new = MaterialPoints_volfrac(mat_id, mat_point_storage_lid) +  delta_volfrac;  // note: change in volfrac was limited above here
                MaterialPoints_volfrac(mat_id, mat_point_storage_lid) = fmin(1.0, fmax(0.0, volfrac_new));                

                // update internal energy
                // dVol/dt = Vol*div
                double GaussPoint_deltaVol = div*rk_alpha*dt*GaussPoint_vol(gauss_gid);
                MaterialPoints_sie(mat_id, mat_point_storage_lid) -= GaussPoint_pres(gauss_gid)*delta_volfrac*GaussPoint_deltaVol/MaterialPoints_mass(mat_id, mat_point_storage_lid);

            } // end for gauss point loop 

        }); // end parallel loop over all material elems in the mesh
        Kokkos::fence();

        return;

    } // end function

    

} // end namespace