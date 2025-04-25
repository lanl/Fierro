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


// -----------------------------------------------------------------------------
// This is the Tipton material pt equilibration model
// ------------------------------------------------------------------------------
namespace TiptonEquilibrationModel {
    
    KOKKOS_FUNCTION
    static void equilibrate (
        const DCArrayKokkos<bool>& MaterialPoints_volfrac,
        const DCArrayKokkos<double>& MaterialPoints_stress,
        const DCArrayKokkos<double>& MaterialPoint_pres,
        const DCArrayKokkos<double>& MaterialPoint_den,
        const DCArrayKokkos<double>& MaterialPoint_sie,
        const double MaterialPoint_sspd,
        const double GaussPoint_avg_press,
        const double GaussPoint_vol,
        double& de,
        const double dt,
        const double rk_alpha,
        const double length,
        const double fuzz,
        const RaggedRightArrayKokkos<double> &equilibration_global_vars,
        const size_t num_vars,
        const size_t mat_point_lid,
        const size_t mat_id)
    {

        // GaussPoint_avg_press = sum_i(volfrac_i/K_i P_i)/ sum_i(volfrac_i/K_i)
        // K_i = rho*c^2

        const double tau = length/MaterialPoint_sspd;  // relaxation time scale, length is element length scale
        const double dP = rk_alpha*dt/tau*(GaussPoint_avg_press - MaterialPoint_pres(mat_point_lid));  // change in pressure
        
        //MaterialPoint_pres(mat_point_lid) += dP;  // we will update density and energy, which will make update the pressure


        // using r=relaxed

        // find the change in density
        // (P^r - P) = (rho^r - rho) c^2
        // rho^r - rho = dP/c^2 
        const double dden = dP/(MaterialPoint_sspd*MaterialPoint_sspd + fuzz);

        // find the corresponding volume change for the material
        // vol^r - vol = (m/rho^r - m/rho)
        const double den0 = MaterialPoint_den(mat_point_lid) + fuzz;
        const double vol0 = MaterialPoints_volfrac(mat_point_lid)*GaussPoint_vol;
        const double MaterialPoints_mass = vol0*den0;
        double dvol = MaterialPoints_mass*( 1.0/(den0+dden) - vol0 );

        // limit dvol based on max volume fraction change
        // val > dvol/gauss_pt_vol  or
        // val*gauss_pt_vol > dvol
        const double val = fmin(1.0, equilibration_global_vars(mat_id, 0));
        dvol = fmin(dvol, val*GaussPoint_vol);
         

        // limit dvol
        double beta = 1.0;
        if (dvol > fuzz){
            // vol0 + beta*dvol <= gauss_pt_vol
            beta = fmin(1.0, (GaussPoint_vol-vol0)/dvol);
        }
        else if (dvol < -fuzz){
            // vol0 + beta*dvol >= 0
            beta = fmax(0.0, -vol0/dvol);
        } // end if

        // make dvol bounded, ensuring vol is between 0 and gauss pt vol
        dvol *= beta;
        
        // calculate the new volume fraction and new density
        const double volr = vol0 + dvol;
        MaterialPoints_volfrac(mat_point_lid) = volr/GaussPoint_vol;
        MaterialPoint_den(mat_point_lid) = MaterialPoints_mass*volr;

        // calculate specific internal energy change due to equilibration
        de = -MaterialPoint_pres(mat_point_lid)*dvol;  // work = -PdV
        
        // Outside this function:
        //   1) The volume fractions are scaled to ensure they sum to < 1
        //   2) The de is scaled to conserve total energy


        return;
    }

} // end namespace




#endif // end Header Guard