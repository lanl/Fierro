#ifndef QR_H
#define QR_H
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

 //////////////////////////

#include <iostream>
#include <stdio.h>
#include <cmath>

#include "tlqs_solver_3D.hpp"

#include "simulation_parameters.hpp"
#include "material.hpp"
#include "boundary_conditions.hpp"
#include "state.hpp"
#include "geometry_new.hpp"
#include "mesh_io.hpp"

#include "matar.h"
using namespace mtr;


// Back substitution to solve Rx = y
void TLQS3D::QR_backsub(const CArrayKokkos <double> &R, 
                     const CArrayKokkos <double> &y,
                     DCArrayKokkos <double> &x,
                     const size_t n_active) {
    
    for (int i = (int)n_active - 1; i >= 0; --i) {
        
        RUN({
            x(i) = y(i);
        });

        double sum = 0.0;
        double sum_lcl = 0.0;

        FOR_REDUCE_SUM(j, i + 1, n_active, 
                       sum_lcl, {
             sum_lcl -= R(i,j) * x(j);
        }, sum);

        RUN({
            x(i) += sum;
            x(i) /= R(i,i);
        });

    } // end for i

    return;

} // end function

// QR Decomposition using Modified Gram-Schmidt
void TLQS3D::QR_decompose(const CArrayKokkos <double> &A, 
                       FArrayKokkos <double> &Q, 
                       CArrayKokkos <double> &R,
                       DCArrayKokkos <double> &v,
                       const size_t n_active) {


    const size_t m = A.dims(0);

    Q.set_values(0.0);
    R.set_values(0.0);

    // Copy columns of A to v, and taking transpose
    FOR_ALL(i, 0, m,
            j, 0, n_active, {
        v(j, i) = A(i,j);  
    });
    Kokkos::fence();


    for (size_t i = 0; i < n_active; ++i) {

        // find the norm of a column in matrix v for row i
        double tally = 0.0;
        double tally_lcl = 0.0;


        FOR_REDUCE_SUM(j, 0, m, 
                       tally_lcl, {
            tally_lcl += v(i,j) * v(i,j);
        }, tally);

        Kokkos::fence();

        RUN({
            R(i,i) = sqrt(tally); // row i norm
        });
        // done with norm calc

        FOR_ALL(j, 0, m, {
            Q(j,i) = v(i,j)/R(i,i);
        });

        Kokkos::fence();

// single parallelism
/*
        FOR_ALL(jj, i+1, n, {

            R(i,jj) = 0.0;

            double sum=0;

            for(size_t k=0; k<m; k++){
                sum += Q(k,i) * A(k,jj);
            }
            R(i,jj) = sum;

            for(size_t k=0; k<m; k++){
                v(jj,k) -= sum * Q(k,i);
            }

        }); // end parallel jj
*/
// nested parallelism

        FOR_FIRST(jj, i+1, n_active, {

            R(i,jj) = 0.0;

            double sum=0;
            double sum_lcl = 0;

            FOR_REDUCE_SUM_SECOND(k, 0, m,
                                  sum_lcl, {
                sum_lcl += Q(k,i) * A(k,jj);
            }, sum);
            R(i,jj) = sum;

           teamMember.team_barrier();


            FOR_SECOND(k, 0, m,{
                v(jj,k) -= sum * Q(k,i);
            });

            teamMember.team_barrier();

        }); // end parallel j

        Kokkos::fence();

    } // end for i

} // end function


double TLQS3D::QR_determinant(const FArrayKokkos <double> &Q,
                           const CArrayKokkos <double> &R)
{
    //const size_t m = Q.dims(0);
    const size_t n = Q.dims(1);
    // Q(m,n,"Q");
    // R(n,n,"R");

    double detR = 1.0;
    int signQ = 1;

    // Accumulate det(R) and track sign adjustments
    // for (size_t i = 0; i < n; ++i) {
    //     if (R(i,i) < 0.0) {
    //         signQ *= -1;             // negative diagonal flips Q's sign
    //         detR *= -R(i,i);         // store magnitude in detR
    //     } else {
    //         detR *= R(i,i);
    //     }
    // }

    double prod_tally;
    double prod_lcl = 1.0;
    FOR_REDUCE_PRODUCT(i, 0, n, 
                    prod_lcl, {
        
        if (R(i,i) < 0.0) {
            prod_lcl *= -R(i,i);         // store magnitude in detR
        } else {
            prod_lcl *= R(i,i);
        }
    }, prod_tally); // end j
    Kokkos::fence();


    detR = prod_tally;



    prod_lcl = 1.0;
    FOR_REDUCE_PRODUCT(i, 0, n, 
                    prod_lcl, {
        
        if (R(i,i) < 0.0) {
            prod_lcl *= -1;              // negative diagonal flips Q's sign
        }
    }, prod_tally); // end j

    signQ = prod_tally;     
    
    Kokkos::fence();

    // Compute sign of det(Q) directly from Q if needed
    // Here we infer sign from R adjustments only.
    double detA = signQ * detR;
    return detA;
}


// Solve for x in Ax = b using QR
// A[m,n]
// x[n]
// b[m]
void TLQS3D::QR_solver(const CArrayKokkos <double> &A, 
                    const CArrayKokkos <double> &b,
                    DCArrayKokkos <double> &x,
                    FArrayKokkos <double> &Q,
                    CArrayKokkos <double> &R,
                    CArrayKokkos <double> &y,
                    DCArrayKokkos <double> &v,
                    const size_t n_active) {
    
    const size_t m = A.dims(0);

    QR_decompose(A, Q, R, v, n_active);

    // Compute Q^t * b
    FOR_FIRST(i, 0, n_active, {

        double sum = 0.0;
        double sum_lcl = 0.0;

        // Q[m,n] so Q^t[n,m] b[m]
        FOR_REDUCE_SUM_SECOND(j, 0, m, 
                              sum_lcl, {
            sum_lcl += Q(j,i) * b(j);
        }, sum);
        y(i) = sum;

    }); // end parallel i

    Kokkos::fence();

    // Solve R x = y
    QR_backsub(R, y, x, n_active);
}

//////////////////////////


#endif // QR