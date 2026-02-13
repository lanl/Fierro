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

#ifndef TIPTON_EQUILIBRATION_H
#define TIPTON_EQUILIBRATION_H
#include "ELEMENTS.h"
#include "material.hpp"
#include "state.hpp"

// -----------------------------------------------------------------------------
// This is the Tipton material pt equilibration model
// ------------------------------------------------------------------------------
namespace TiptonEquilibrationModel {

    
    void mat_equilibration(
        const Material_t& Materials, 
        const swage::Mesh& mesh, 
        State_t& State,
        CArrayKokkos <double>& GaussPoint_pres,
        CArrayKokkos <double>& GaussPoint_pres_denominator,
        CArrayKokkos <double>& GaussPoint_volfrac_min,
        CArrayKokkos <double>& GaussPoint_volfrac_limiter,
        const double dt,
        const double rk_alpha,
        const double fuzz,
        const double small);

    void geo_equilibration(const Material_t& Materials, 
        const swage::Mesh& mesh, 
        State_t& State,
        CArrayKokkos <double>& GaussPoint_pres,
        CArrayKokkos <double>& GaussPoint_pres_denominator,
        CArrayKokkos <double>& GaussPoint_volfrac_min,
        CArrayKokkos <double>& GaussPoint_volfrac_limiter,
        const double dt,
        const double rk_alpha,
        const double fuzz,
        const double small);

    void build_gauss_point_averages (
        const swage::Mesh& mesh,
        const CArrayKokkos<double>& GaussPoint_pres,
        const CArrayKokkos<double>& GaussPoint_pres_denominator,
        const CArrayKokkos<double>& GaussPoint_volfrac_min,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_volfrac,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_pres,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_den,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_sspd,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const points_in_mat_t& points_in_mat_elem,
        const double dt,
        const double rk_alpha,
        const double fuzz,
        const size_t num_mat_elems,
        const size_t mat_id);


    void calc_gauss_point_averages( 
            const swage::Mesh& mesh,
            const CArrayKokkos<double>&  GaussPoint_pres,
            const CArrayKokkos<double>&  GaussPoint_pres_denominator,
            const double fuzz);

    void calc_volfrac_change (
        const swage::Mesh& mesh,
        const CArrayKokkos<double>& GaussPoint_pres,
        const CArrayKokkos<double>& GaussPoint_pres_denominator,
        const CArrayKokkos <double>& GaussPoint_volfrac_min,
        const CArrayKokkos <double>& GaussPoint_volfrac_limiter,
        const DCArrayKokkos<double>& GaussPoint_vel_grad,
        const DCArrayKokkos<double>& GaussPoint_vol,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_volfrac,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_delta_volfrac,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_pres,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_den,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_sspd,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_mass,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const points_in_mat_t& points_in_mat_elem,
        const double dt,
        const double rk_alpha,
        const double fuzz,
        const CArrayKokkos<double> &equilibration_global_vars,
        const size_t num_global_vars,
        const size_t num_mat_elems,
        const size_t mat_id);        


        void update_state_equilibration (
            const swage::Mesh& mesh,
            const Material_t& Materials,
            const CArrayKokkos<double>& GaussPoint_pres,
            const CArrayKokkos <double>& GaussPoint_volfrac_limiter,
            const DCArrayKokkos<double>& GaussPoint_vel_grad,
            const DCArrayKokkos<double>& GaussPoint_vol,
            const DRaggedRightArrayKokkos<double>& MaterialPoints_volfrac_inout,       // the modified value by equilibration
            const DRaggedRightArrayKokkos<double>& MaterialPoints_delta_volfrac_inout, // the modified value by equilibration
            const DRaggedRightArrayKokkos<double>& MaterialPoints_volfrac_in,          // unmodified value
            const DRaggedRightArrayKokkos<double>& MaterialPoints_delta_volfrac_in,    // unmodified value
            const DRaggedRightArrayKokkos<double>& MaterialPoints_pres,
            const DRaggedRightArrayKokkos<double>& MaterialPoints_den,
            const DRaggedRightArrayKokkos<double>& MaterialPoints_sie,
            const DRaggedRightArrayKokkos<double>& MaterialPoints_sspd,
            const DRaggedRightArrayKokkos<double>& MaterialPoints_mass,
            const DRaggedRightArrayKokkos<double>& MaterialPoints_stress,
            const DRaggedRightArrayKokkos<double>& MaterialPoints_shear_modulii,
            const DRaggedRightArrayKokkos<double>& MaterialPoints_eos_state_vars,
            const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
            const points_in_mat_t& points_in_mat_elem,
            const double dt,
            const double rk_alpha,
            const double fuzz,
            const size_t num_mat_elems,
            const size_t mat_id);    

} // end namespace




#endif // end Header Guard