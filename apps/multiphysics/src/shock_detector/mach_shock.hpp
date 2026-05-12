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

#ifndef MACH_SHOCK_H
#define MACH_SHOCK_H

#include "matar.h"
#include "ELEMENTS.h"

namespace MachShockDetector {

    static void detect_shock(const swage::Mesh& mesh,
                             const DCArrayKokkos<double>& GaussPoints_vel_grad,
                             const MPICArrayKokkos<double>& GaussPoints_shock_detector,
                             const DCArrayKokkos<double>& GaussPoints_vol,
                             const DRaggedRightArrayKokkos<double>& MaterialPoints_sspd,
                             const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
                             const DCArrayKokkos <size_t>& num_mat_elems,
                             const double fuzz,
                             const double small,
                             const size_t num_materials)
    {

        // ---- Calculate the shock detector for the Riemann-solver for each gauss point ----
        // ---- choosing the largest for all materials at each gauss point
        //
        // The dissipation from the Riemann problem is limited by phi
        //    phi = (1. - max( 0., min( 1. , r_face ) ))^n
        //  where
        //      r_face = (C* div(u_+)/div(u_z))
        //  The plus denotes the cell center divergence of a neighbor.
        //  The solution will be first order when phi=1 and have
        //  zero dissipation when phi=0.
        //      phi = 0 highest-order solution
        //      phi = 1 first order solution
        //

        for(size_t mat_id = 0; mat_id < num_materials; mat_id++) {
            
            FOR_ALL(mat_elem_sid, 0, num_mat_elems.host(mat_id), {

                size_t elem_gid = elem_in_mat_elem(mat_id, mat_elem_sid);

                for(size_t gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem; gauss_lid++){
                    
                    size_t gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);

                    size_t mat_point_sid = mat_elem_sid; // Note: hard coded to 1 per element

                    double vol = GaussPoints_vol(gauss_gid);
                    
                    double phi    = 0.0;  // the shock detector
                    double r_face = 1.0;  // the ratio on the face
                    double r_min  = 1.0;  // the min ratio for the cell
                    double r_coef = 0.9;  // 0.9; the coefficient on the ratio
                                        //   (1=minmod and 2=superbee)
                    double n_coef = 1.0;  // the power on the limiting coefficient
                                        //   (1=nominal, and n_coeff > 1 oscillatory)

                    double div = GaussPoints_vel_grad(gauss_gid, 0, 0) + 
                                 GaussPoints_vel_grad(gauss_gid, 1, 1) + 
                                 GaussPoints_vel_grad(gauss_gid, 2, 2);

                    // loop over the neighboring cells
                    for (size_t elem_lid = 0; elem_lid < mesh.num_elems_in_elem(elem_gid); elem_lid++) {
                        // Get global index for neighboring cell
                        size_t neighbor_gid = mesh.elems_in_elem(elem_gid, elem_lid);

                        // calculate the velocity divergence in neighbor
                        double div_neighbor = GaussPoints_vel_grad(neighbor_gid, 0, 0) + 
                                              GaussPoints_vel_grad(neighbor_gid, 1, 1) + 
                                              GaussPoints_vel_grad(neighbor_gid, 2, 2);

                        r_face = r_coef * (div_neighbor + small) / (div + small);

                        // store the smallest face ratio
                        r_min = fmin(r_face, r_min);
                    } // end for elem_lid

                    // calculate standard shock detector
                    phi = 1.0 - fmax(0.0, r_min);
                    phi = pow(phi, n_coef);

                    //  Mach number shock detector
                    double omega    = 20.0; // 20.0;    // weighting factor on Mach number
                    double third    = 1.0 / 3.0;
                    double c_length = pow(vol, third); // characteristic length
                    double alpha    = fmin(1.0, omega * (c_length * fabs(div)) / (MaterialPoints_sspd(mat_id, mat_point_sid) + fuzz) );

                    // use Mach based detector with standard shock detector

                    // turn off dissipation in expansion
                    // alpha = fmax(-fabs(div0)/div0 * alpha, 0.0);  // this should be if(div0<0) alpha=alpha else alpha=0

                    phi = alpha * phi;

                    phi = fmax(0.0, fmin(1.0, phi));

                    double prev_phi = GaussPoints_shock_detector(gauss_gid); // Phi for previous material in this element
               

                    GaussPoints_shock_detector(gauss_gid) = fmin(prev_phi, phi);
                }
            });
            MATAR_FENCE();
        }

        return;
    }
} // end MachShockDetector namespace


#endif // end MACH_SHOCK_H Header Guard