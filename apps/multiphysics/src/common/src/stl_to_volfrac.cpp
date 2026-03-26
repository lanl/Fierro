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

// -----------------------------------------------
// pointcloud reproducing kernels in C++
//  Nathaniel Morgan
// -----------------------------------------------

#include <chrono>   // for timing

#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cmath>

#include <cstdlib> // For rand() and srand()


#include "matar.h"

#include "lu_solver.hpp"


#include <set> // for unorded map testing 



#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

using namespace mtr;

const double PI = 3.14159265358979323846;

// -----------------------------------------------
// inputs:


#define P2 // P1 or P2
#define CUBIC_SPLINE // CUBIC_SPLINE or GAUSS kernel

bool RAND_CLOUD = true; // RAND_CLOUD or uniform

const size_t num_1d_x = 5;
const size_t num_1d_y = 5;
const size_t num_1d_z = 5;

const double h_kernel = 2.0/5.; // it is always 2/num_1d_x
const double num_points_fit = 27; // minimum to P2 fit on structured mesh is 3x3x3

const size_t num_points = num_1d_x*num_1d_y*num_1d_z;

// the bin sizes for finding neighboring points
const double bin_dx = 0.05; // bins in x
const double bin_dy = 0.05; // bins in y
const double bin_dz = 0.05; // bins in z

const double X0 = 0.0;   // origin
const double Y0 = 0.0;
const double Z0 = 0.0;

// length of the domain 
const double LX = 1.0;   // length in x-dir
const double LY = 1.0;
const double LZ = 1.0;


bool check_maps = false; // CPU only!!!!

//
// -----------------------------------------------



struct bin_keys_t{
    size_t i,j,k;
};

KOKKOS_INLINE_FUNCTION
size_t get_gid(size_t i, size_t j, size_t k, size_t num_x, size_t num_y){
    return i + (j + k*num_y)*num_x;
}

KOKKOS_INLINE_FUNCTION
bin_keys_t get_bin_keys(const double x_pt, 
                        const double y_pt, 
                        const double z_pt){
            

    double i_dbl = fmax(0, round((x_pt - X0 - bin_dx*0.5)/bin_dx - 1.0e-10)); // x = ih + X0 + dx_bin*0.5
    double j_dbl = fmax(0, round((y_pt - Y0 - bin_dy*0.5)/bin_dy - 1.0e-10));
    double k_dbl = fmax(0, round((z_pt - Z0 - bin_dz*0.5)/bin_dz - 1.0e-10));

    bin_keys_t bin_keys; // save i,j,k to the bin keys

    // get the integer for the bins
    bin_keys.i = (size_t)i_dbl;
    bin_keys.j = (size_t)j_dbl;
    bin_keys.k = (size_t)k_dbl;

    return bin_keys;

} // end function

KOKKOS_INLINE_FUNCTION
size_t get_bin_gid(const double x_pt, 
                   const double y_pt, 
                   const double z_pt, 
                   const size_t num_bins_x,
                   const size_t num_bins_y,
                   const size_t num_bins_z){
            

    double i_dbl = fmin(num_bins_x-1, fmax(0.0, round((x_pt - X0)/bin_dx - 1.0e-8))); // x = ih + X0
    double j_dbl = fmin(num_bins_y-1, fmax(0.0, round((y_pt - Y0)/bin_dy - 1.0e-8)));
    double k_dbl = fmin(num_bins_z-1, fmax(0.0, round((z_pt - Z0)/bin_dz - 1.0e-8)));

    // get the integers for the bins
    size_t i = (size_t)i_dbl;
    size_t j = (size_t)j_dbl;
    size_t k = (size_t)k_dbl;
    
    // get the 1D index for this bin                               
    return get_gid(i, j, k, num_bins_x, num_bins_y);

} // end function




#if defined(CUBIC_SPLINE)

KOKKOS_FUNCTION
double kernel(const double r[3], double h) {
    
    double diff_sqrd = 0.0;
    for(size_t dim=0; dim<3; dim++){
        diff_sqrd += r[dim]*r[dim];
    } // dim

    const double radius = sqrt(diff_sqrd);
    const double q = radius/h;
    const double alpha = 2.0/(3.0*h);
    if (q < 0.0) return 0.0; // defensive
    if (q < 1.0) return (alpha * (1.0 - 1.5*q*q + 0.75*q*q*q));
    if (q < 2.0) return (alpha * 0.25 * pow(2.0 - q, 3));

    return 0.0;
}


KOKKOS_FUNCTION
// derivative dW/dx_i = - dW/dr where r = xj-xi
void grad_kernel(double *grad_W, const double r[3], const double h, const bool derivative_wrt_i) {

    double diff_sqrd = 0.0;
    for(size_t dim=0; dim<3; dim++){
        diff_sqrd += r[dim]*r[dim];
    } // dim

    const double radius = sqrt(diff_sqrd);
    const double q = radius/h;

    double df_dq = 0.0; // derivative of the dimensionless kernel shape function f(q)
    if (q < 1.0) {
        // f(q) = 1 - 1.5 q^2 + 0.75 q^3
        // f'(q) = -3 q + 2.25 q^2
        df_dq = -3.0 * q + 2.25 * q * q;
    } else if (q < 2.0) {
        // f(q) = 0.25 (2 - q)^3
        // f'(q) = 0.25 * 3 (2 - q)^2 * (-1) = -0.75 (2 - q)^2
        const double two_minus_q = 2.0 - q;
        df_dq = -0.75 * two_minus_q * two_minus_q;
    } else {
        df_dq = 0.0;
    }

    const double alpha = 2.0/(3.0*h);
    const double dW_dr = alpha*(df_dq / h);
    // grad W = dW/dr * (rij / radius)
    const double invr = 1.0 /(radius + 1e-16);

    double drdx = -1.0;  // default is derivative with respect to i
    if(derivative_wrt_i == false){
        drdx = 1.0;  // derivative with respect to j
    } // end if

    for (size_t dim=0; dim<3; ++dim) {
        grad_W[dim] = dW_dr * r[dim] * drdx * invr;
    }

    return;
}

#else 

// Gaussian function part of the RBF
// rbf = exp(-(xj - x)*(xj - x)/h)
KOKKOS_FUNCTION
double kernel(const double r[3], const double h){

    double diff_sqrd = 0.0;

    for(size_t dim=0; dim<3; dim++){
        diff_sqrd += r[dim]*r[dim];
    } // dim

    double norm = 1.0; // / (h * h * h * pow(PI, 1.5));
    return norm * exp(-diff_sqrd / (h * h));
} // end of function



// Gradient Gaussian function
// d/dx rbf = d/dx (exp(-(xj - xi)*(xj - x)/hi^2) 
KOKKOS_FUNCTION
void grad_kernel(double *grad_W, const double r[3], const double h, const bool derivative_wrt_i){

    double diff_sqrd = 0.0;

    for(size_t dim=0; dim<3; dim++){
        diff_sqrd += r[dim]*r[dim];
    } // dim

    double drdx = -1;
    if(derivative_wrt_i == false){
        drdx = 1.0;  // derivative with respect to j
    } // end if

    const double rbf = kernel(r, h);

    // gradient
    for (size_t dim=0; dim<3; ++dim) {
        grad_W[dim] = -2.0/(h*h)*r[dim]*rbf*drdx; 
    }

    return;
} // end of function

#endif

// Gaussian function part of the RBF, symmeterized
// rbf = 0.5*(exp(-(xj - xi)*(xj - xi)/hi^2) + exp(-(xi - xj)*(xi - xj)/hj^2))
KOKKOS_FUNCTION
double kernel_syn(const double r[3], const double hi, const double hj){

    double diff_sqrd = 0.0;

    for(size_t dim=0; dim<3; dim++){
        diff_sqrd += r[dim]*r[dim];
    } // dim

    const double Wi = exp(-diff_sqrd/(hi*hi)); // use kernel func call
    const double Wj = exp(-diff_sqrd/(hj*hj));

    return 0.5*(Wi + Wj);
} // end of function


// d/dx rbf = d/dx ( 0.5(exp(-(xj - xi)*(xj - x)/hi^2) + exp(-(xi - xj)*(xi - xj)/hj^2)) ) 
KOKKOS_FUNCTION
void grad_kernel_sym(double *gradW, const double r[3], const double hi, const double hj) {
    double diff_sqrd = 0.0;
    for (size_t dim=0; dim<3; ++dim){
        diff_sqrd += r[dim]*r[dim];
    }

    const double drdxi = -1;

    double Wi = exp(-diff_sqrd/(hi*hi));
    double Wj = exp(-diff_sqrd/(hj*hj));

    double dWi = -2.0/ (hi*hi) * Wi*drdxi; // it uses xj - xi so a minus one
    double dWj = -2.0/ (hj*hj) * Wj; // it uses xi - xj so it has a +1 for drdxi

    for (size_t dim=0; dim<3; ++dim) {
        gradW[dim] = 0.5 * (dWi * r[dim] - dWj * r[dim]); // second term using -r
    }
}

#if defined(P2)
    // Polynomial basis up to quadratic in 3D (10 terms)
    const size_t num_poly_basis = 10;
    KOKKOS_INLINE_FUNCTION
    void poly_basis(const double r[3], double *p) {

        p[0] = 1.0;
        p[1] = r[0];
        p[2] = r[1];
        p[3] = r[2];
        p[4] = r[0] * r[0];
        p[5] = r[0] * r[1];
        p[6] = r[0] * r[2];
        p[7] = r[1] * r[1];
        p[8] = r[1] * r[2];
        p[9] = r[2] * r[2];

        // for high-order will use (x^a y^b z^c)

        return;
    } // end function


    KOKKOS_INLINE_FUNCTION
    void grad_poly_basis(const double r[3], double (*grad_p)[3], bool derivative_wrt_i) {
        
        // definition, r = r_j - r_i

        double drdx = -1.0;  // default is derivative with respect to i
        double drdy = -1.0;  // default is derivative with respect to i
        double drdz = -1.0;  // default is derivative with respect to i
        if(derivative_wrt_i == false){
            drdx = 1.0;  // derivative with respect to j
            drdy = 1.0;  // derivative with respect to j
            drdz = 1.0;  // derivative with respect to j
        } // end if

        grad_p[0][0] = 0.0;
        grad_p[1][0] = drdx;
        grad_p[2][0] = 0.0;
        grad_p[3][0] = 0.0;
        grad_p[4][0] = 2.0*r[0]*drdx;
        grad_p[5][0] = r[1]*drdx;
        grad_p[6][0] = r[2]*drdx;
        grad_p[7][0] = 0.0;
        grad_p[8][0] = 0.0;
        grad_p[9][0] = 0.0;

        // for high-order will use (x^a y^b z^c)

        grad_p[0][1] = 0.0;
        grad_p[1][1] = 0.0;
        grad_p[2][1] = drdy;
        grad_p[3][1] = 0.0;
        grad_p[4][1] = 0.0;
        grad_p[5][1] = r[0]*drdy;
        grad_p[6][1] = 0.0;
        grad_p[7][1] = 2.0*r[1]*drdy;
        grad_p[8][1] = r[2]*drdy;
        grad_p[9][1] = 0.0;

        // for high-order will use (x^a y^b z^c)

        grad_p[0][2] = 0.0;
        grad_p[1][2] = 0.0;
        grad_p[2][2] = 0.0;
        grad_p[3][2] = drdz;
        grad_p[4][2] = 0.0;
        grad_p[5][2] = 0.0;
        grad_p[6][2] = r[0]*drdz;
        grad_p[7][2] = 0.0;
        grad_p[8][2] = r[1]*drdz;
        grad_p[9][2] = 2.0*r[2]*drdz;

        // for high-order will use (x^a y^b z^c)

        return;
    } // end function
#else
    // Polynomial basis up to quadratic in 3D (10 terms)
    const size_t num_poly_basis = 4;
    KOKKOS_INLINE_FUNCTION
    void poly_basis(const double r[3], double *p) {

        p[0] = 1.0;
        p[1] = r[0];
        p[2] = r[1];
        p[3] = r[2];

        return;
    } // end function


    KOKKOS_INLINE_FUNCTION
    void grad_poly_basis(const double r[3], double (*grad_p)[3], bool derivative_wrt_i) {
        
        double drdx = -1.0;  // default is derivative with respect to i
        double drdy = -1.0;  // default is derivative with respect to i
        double drdz = -1.0;  // default is derivative with respect to i
        if(derivative_wrt_i == false){
            drdx = 1.0;  // derivative with respect to j
            drdy = 1.0;  // derivative with respect to j
            drdz = 1.0;  // derivative with respect to j
        } // end if

        grad_p[0][0] = 0.0;
        grad_p[1][0] = drdx;
        grad_p[2][0] = 0.0;
        grad_p[3][0] = 0.0;

        grad_p[0][1] = 0.0;
        grad_p[1][1] = 0.0;
        grad_p[2][1] = drdy;
        grad_p[3][1] = 0.0;


        grad_p[0][2] = 0.0;
        grad_p[1][2] = 0.0;
        grad_p[2][2] = 0.0;
        grad_p[3][2] = drdz;

        return;
    } // end function
#endif


void calc_basis_functions(
    size_t point_gid,
    const DCArrayKokkos <double>& x,
    const DCArrayKokkos <size_t> points_num_neighbors, 
    const DRaggedRightArrayKokkos <size_t> points_in_point,
    const DCArrayKokkos <double>& vol,
    const CArrayKokkos <double>& p_coeffs,
    const DRaggedRightArrayKokkos <double>& basis,
    const double h)
{

    //---------------------------------------------
    // walk over the neighboring points 
    //---------------------------------------------

    FOR_ALL(neighbor_point_lid, 0, points_num_neighbors(point_gid), {

        size_t neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);

        double p[num_poly_basis];    // array holding polynomial basis [x, y, z, x^2, y^2, ... , yz]
        double r[3];    // vecx_j - vecx_i
        r[0] = x(neighbor_point_gid,0) - x(point_gid,0); // x_j-x_i
        r[1] = x(neighbor_point_gid,1) - x(point_gid,1); // y_j-y_i
        r[2] = x(neighbor_point_gid,2) - x(point_gid,2); // z_j-z_i

        double W = kernel(r, h);
        poly_basis(r,p);

        double correction = 0.0;
        for (size_t a = 0; a < num_poly_basis; ++a){
            correction += p_coeffs(point_gid,a) * p[a];
        } // end for a

        basis(point_gid, neighbor_point_lid) = W * correction;

    }); // neighbor_point_lid

    return;
    
} // end function


// grad_C = -M^-1 * grad_M * C
// grad_M[a][b][i] = V (P[a] * grad_p[i][b] W + grad_p[i][a] P[b] W + P[a] P[b] grad_W[i])
// -M^-1[d][k]*grad_M[k][b][i]*C[b] = [d][i]
// p[d] grad_C[d][i]
void calc_basis_and_grad_basis_functions(
    const DCArrayKokkos <double>& x,
    const DCArrayKokkos <size_t> points_num_neighbors, 
    const DRaggedRightArrayKokkos <size_t> points_in_point,
    const DCArrayKokkos <double>& vol,
    const CArrayKokkos <double>& p_coeffs,
    const CArrayKokkos <double>& M_inv,
    const DRaggedRightArrayKokkos <double>& basis,
    const DRaggedRightArrayKokkos <double>& grad_basis,
    const double h,
    const bool derivative_wrt_i)
{

    // dir 

    // actual number of points
    size_t num_points = x.dims(0);
    
    // loop over all nodes in the problem
    FOR_ALL(point_gid, 0, num_points, {

        // --------------
        // Step 1: assemble grad M, the gradient of the moment matrix

        double grad_m[num_poly_basis][num_poly_basis][3];

        for(size_t a=0;a<num_poly_basis;++a){
            for(size_t b=0;b<num_poly_basis;++b){
                for(size_t dim=0;dim<3;++dim){
                    grad_m[a][b][dim] = 0.0;
                }
            } // end for b
        } // end for a

        // walk over the neighboring points 
        for(size_t neighbor_point_lid=0; neighbor_point_lid<points_num_neighbors(point_gid); neighbor_point_lid++){

            size_t neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);

            double r[3];    // vecx_j - vecx_i
            r[0] = x(neighbor_point_gid,0) - x(point_gid,0); // x_j-x_i
            r[1] = x(neighbor_point_gid,1) - x(point_gid,1); // y_j-y_i
            r[2] = x(neighbor_point_gid,2) - x(point_gid,2); // z_j-z_i



            double W = kernel(r,h);
            //printf("kernel = %f \n", W);
            double grad_W[3]; 
            grad_kernel(grad_W,r,h,derivative_wrt_i);
            //printf("grad kernel = %f, %f, %f \n", grad_W[0], grad_W[1], grad_W[2]);

            double p[num_poly_basis]; 
            poly_basis(r,p);

            double grad_p[num_poly_basis][3]; 
            grad_poly_basis(r, grad_p, derivative_wrt_i);

            double Vj = vol(neighbor_point_gid);

            // grad_M = V (P grad_p W + grad_p P W + P P grad_W)
            for(size_t a=0;a<num_poly_basis;++a){
                for(size_t b=0;b<num_poly_basis;++b){
                    for(size_t dim=0;dim<3;++dim){
                        grad_m[a][b][dim] += Vj * ( p[a]*grad_p[b][dim]*W + grad_p[a][dim]*p[b]*W + p[a]*p[b]*grad_W[dim] );
                    }
                } // end for b
            } // end for a

        } // end for loop over neighboring points


        // -----------
        // Step 2: calculate basis and grad basis

        // walk over the neighboring points 
        for(size_t neighbor_point_lid=0; neighbor_point_lid<points_num_neighbors(point_gid); neighbor_point_lid++){

            size_t neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);

            double r[3];    // vecx_j - vecx_i
            r[0] = x(neighbor_point_gid,0) - x(point_gid,0); // x_j-x_i
            r[1] = x(neighbor_point_gid,1) - x(point_gid,1); // y_j-y_i
            r[2] = x(neighbor_point_gid,2) - x(point_gid,2); // z_j-z_i

            double W = kernel(r, h);

            double p[num_poly_basis];    // array holding polynomial basis [x, y, z, x^2, y^2, ... , yz]
            poly_basis(r,p);
            
            double grad_W[3];
            grad_kernel(grad_W, r, h, derivative_wrt_i);

            double grad_p[num_poly_basis][3]; // matrix holding grad polynomial basis
            grad_poly_basis(r, grad_p, derivative_wrt_i);

            // 
            double correction = 0.0;
            for (size_t a = 0; a < num_poly_basis; ++a){
                correction += p_coeffs(point_gid,a) * p[a];
            } // end for a

            basis(point_gid, neighbor_point_lid) = W * correction;


            // --- gradient contributions ---

            // term from grad poly
            double term[3]; 
            for(size_t dim=0;dim<3;++dim){
                term[dim] = 0.0;
            }

            // calc the grad of poly term
            for(size_t dim=0;dim<3;++dim){
                for (size_t a = 0; a < num_poly_basis; ++a){
                    term[dim] += grad_p[a][dim] * p_coeffs(point_gid,a);
                } // end for a
            } // end for dim

            // saving the grad poly term plus the grad kernel term to the grad basis
            for(size_t dim=0;dim<3;++dim){
                grad_basis(point_gid,neighbor_point_lid,dim) = term[dim]*W + correction*grad_W[dim];
            }

            // --- gradient of correction coefficients (grad_C) ---

            // the last contirubtion to grad_basis is from grad_C
            // sum_a p[a]*W*grad_C[i][a]
            double grad_C[num_poly_basis][3];

            
            for (size_t a = 0; a < num_poly_basis; ++a){
                for(size_t dim=0;dim<3;++dim){
                    grad_C[a][dim]= 0.0;
                }
            } // end for 1
            

            // grad_C = -M^-1 * grad_M * C
            for (size_t d = 0; d < num_poly_basis; ++d) {
                for (size_t b = 0; b < num_poly_basis; ++b) {
                    double Cb = p_coeffs(point_gid,b);
                    for (size_t k = 0; k < num_poly_basis; ++k) {
                        for(size_t dim=0;dim<3;++dim){
                            grad_C[d][dim] -= M_inv(point_gid,d,k) * grad_m[k][b][dim] * Cb;
                        } // end dim
                    } // end for k
                } // end for b
            } // end for d


            // adding grad_C[i][d]* p[d]*W to the grad basis function
            for (size_t d = 0; d < num_poly_basis; ++d){
                for(size_t dim=0;dim<3;++dim){
                    grad_basis(point_gid,neighbor_point_lid,dim) += grad_C[d][dim]*p[d]*W;
                } // end dim
            } // end for d


        } // end for over neighbor_point_lid

    }); // end parallel loop over all points 

    return;
    
} // end function



// Build reproducing kernel poly coefficients for all particles in the domain
void calc_p_coefficients(
    const DCArrayKokkos <double>& x,
    const DCArrayKokkos <size_t> points_num_neighbors, 
    const DRaggedRightArrayKokkos <size_t> points_in_point,
    const DCArrayKokkos <double>& vol,
    const CArrayKokkos <double>& p_coeffs,
    const CArrayKokkos <double>& M_inv,
    double h)
{

    // actual number of points
    size_t num_points = x.dims(0);

    
    // loop over all nodes in the problem
    FOR_ALL(point_gid, 0, num_points, {

        double M_1D[num_poly_basis*num_poly_basis]; 
        ViewCArrayKokkos <double> M(&M_1D[0], num_poly_basis, num_poly_basis);
        M.set_values(0.0);

        // values in rhs after this function will be accessed as p_coeffs(i,0:N)
        ViewCArrayKokkos <double> rhs (&p_coeffs(point_gid,0), num_poly_basis);
        rhs.set_values(0.0);
        rhs(0) = 1.0;   // enforce reproduction of constant 1, everything else is = 0

        double p[num_poly_basis];    // array holding polynomial basis [x, y, z, x^2, y^2, ... , yz]
        double r[3];    // vecx_j - vecx_i

        //---------------------------------------------
        // walk over the neighboring points
        //---------------------------------------------

        for (size_t neighbor_point_lid=0; neighbor_point_lid<points_num_neighbors(point_gid); neighbor_point_lid++){

            size_t neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);

            r[0] = x(neighbor_point_gid,0) - x(point_gid,0); // x_j-x_i
            r[1] = x(neighbor_point_gid,1) - x(point_gid,1); // y_j-y_i
            r[2] = x(neighbor_point_gid,2) - x(point_gid,2); // z_j-z_i

            double W = kernel(r, h);
            poly_basis(r,p);

            // assemble matrix

            for (size_t a = 0; a < num_poly_basis; ++a) {
                for (size_t b = 0; b < num_poly_basis; ++b) {
                    M(a,b) += vol(neighbor_point_gid) * W * p[a] * p[b]; 
                } // end for b
            } // for a

        } // neighbor_point_lid

    
        // -------------
        // solve Ax=B
        // -------------

        size_t perm_1D[num_poly_basis];
        ViewCArrayKokkos <size_t> perm (&perm_1D[0], num_poly_basis);
        for (size_t a = 0; a < num_poly_basis; ++a) {
            perm(a)= 0;
        } // end a

        double vv_1D[num_poly_basis];
        ViewCArrayKokkos <double> vv(&vv_1D[0], num_poly_basis);
        
        // used for LU problem
        int singular = 0;
        int parity = 0;
        singular = LU_decompose(M, perm, vv, parity);  // M is returned as the LU matrix  
        if(singular==0){
            printf("ERROR: matrix is singluar \n");
        }


        // --------------------------------------------------
        // things needed for gradient of the basis function
        double col_1D[num_poly_basis];
        ViewCArrayKokkos <double> col(&col_1D[0], num_poly_basis);

        // making a view, inverting only the matrix at point i
        ViewCArrayKokkos <double> M_inv_point(&M_inv(point_gid,0,0), num_poly_basis,num_poly_basis);

        LU_invert(M,        // input matrix
                  perm,     // permutations
                  M_inv_point, // inverse matrix at point gid
                  col);     // tmp array
        // -------------------------------------------------
        
        // solve for p_coefs
        LU_backsub(M, perm, rhs);  // note: answer is sent back in rhs

    }); // end parallel loop


    return; 
} // end function





int main(int argc, char *argv[])
{
    Kokkos::initialize(argc, argv);
    {  

        printf("Pointcloud Reproducing Kernels \n\n");


        // define a point cloud
        DCArrayKokkos <double> point_positions(num_points, 3, "point_positions");
        DCArrayKokkos <double> point_values(num_points, "point_values"); 

        DCArrayKokkos <double> vol(num_points);
        vol.set_values(0.0);

        // point locations
        if(RAND_CLOUD){
            srand(static_cast<unsigned int>(time(0))); // Seed the random number generator
                for(size_t i=0; i<num_points; i++){
                    point_positions.host(i, 0) = X0 + LX*static_cast<double>(rand())/static_cast<double>(RAND_MAX);
                    point_positions.host(i, 1) = Y0 + LY*static_cast<double>(rand())/static_cast<double>(RAND_MAX);
                    point_positions.host(i, 2) = Z0 + LZ*static_cast<double>(rand())/static_cast<double>(RAND_MAX);
                }
                vol.set_values(1.0);
            }
        else {

            double dx = (LX-X0)/((double)num_1d_x - 1);
            double dy = (LY-Y0)/((double)num_1d_y - 1);
            double dz = (LZ-Z0)/((double)num_1d_z - 1);

            size_t point_gid = 0;  
            for(size_t k=0; k<num_1d_z; k++){
                for(size_t j=0; j<num_1d_y; j++){
                    for(size_t i=0; i<num_1d_x; i++){
                        point_positions.host(point_gid, 0) = X0 + static_cast<double>(i)*dx;
                        point_positions.host(point_gid, 1) = Y0 + static_cast<double>(j)*dy;
                        point_positions.host(point_gid, 2) = Z0 + static_cast<double>(k)*dz;
                        point_gid++;
                    } // end i
                } // end j
            } // end k



            const double elem_dx = LX/((double)num_1d_x);
            const double elem_dy = LY/((double)num_1d_y);
            const double elem_dz = LZ/((double)num_1d_z); 
            const double elem_vol = elem_dx*elem_dy*elem_dz;

            const size_t num_cells_1d_x = num_1d_x-1;
            const size_t num_cells_1d_y = num_1d_y-1;
            const size_t num_cells_1d_z = num_1d_z-1;

            FOR_ALL(k,0,num_cells_1d_z,
                    j,0,num_cells_1d_y,
                    i,0,num_cells_1d_x,{

                for (int kcount=k; kcount<=k+1; kcount++){
                    for (int jcount=j; jcount<=j+1; jcount++){
                        for (int icount=i; icount<=i+1; icount++){
                            size_t point_gid = get_gid(icount, jcount, kcount, num_1d_x, num_1d_y);
                            Kokkos::atomic_add(&vol(point_gid), elem_vol*0.25);
                        } // end i
                    } // end j
                } // end k
                        
            }); // end parallel over k,j,i 
            

        } // end if
        vol.update_host();
        point_positions.update_device();
        Kokkos::fence();

        // point values
        FOR_ALL(i, 0, num_points, {

            //printf("point location at i=%d is (%f, %f, %f) \n", i, point_positions(i, 0), point_positions(i, 1), point_positions(i, 2));
            point_values(i) = sqrt(point_positions(i, 0)*point_positions(i, 0) + 
                                   point_positions(i, 1)*point_positions(i, 1) +
                                   point_positions(i, 2)*point_positions(i, 2));

        }); // end parallel for tri's in the file
        point_values.update_host();
        Kokkos::fence();
        //printf("\n");


        // ----------------------------
        // Make bins here
        // ----------------------------

        printf("making bins \n");
        
        // the number of nodes in the mesh
        size_t num_bins_x = (size_t)( round( (LX - X0)/bin_dx) + 1 );  
        size_t num_bins_y = (size_t)( round( (LY - Y0)/bin_dy) + 1 );  
        size_t num_bins_z = (size_t)( round( (LZ - Z0)/bin_dz) + 1 );  
        //  bin_dx = (LX-X0)/(num_bins_x - 1);
        //  bin_dy = (LY-Y0)/(num_bins_y - 1);
        //  bin_dz = (LZ-Z0)/(num_bins_z - 1);


        size_t num_bins = num_bins_x*num_bins_y*num_bins_z;
        printf("num bins x=%zu, y=%zu, z=%zu \n", num_bins_x, num_bins_y, num_bins_z);

        // bins and their connectivity to each other and points
        DCArrayKokkos <bin_keys_t> keys_in_bin(num_bins, "keys_in_bin"); // mapping from gid to (i,j,k)
        DCArrayKokkos <size_t> num_points_in_bin(num_bins, "num_bins");
        num_points_in_bin.set_values(0);
        DRaggedRightArrayKokkos <size_t> points_in_bin; // allocated later
        

        // connectivity from points to bins
        DCArrayKokkos <size_t> points_bin_gid(num_points, "points_in_gid");
        CArrayKokkos <size_t>  points_bin_lid_storage(num_points, "bin_lid_storage");  // only used to create storage
        DCArrayKokkos <int> points_bin_stencil(num_points, 6, "bin_stencil");   // how imin,imax,jmin,jmax,kmin,kmax range for bins in stencil
        DCArrayKokkos <size_t> points_num_neighbors(num_points, "num_neighbors");
        
        printf("Starting timers \n\n");


        // start timer
        auto time_1 = std::chrono::high_resolution_clock::now();

        // build reverse mapping between gid and i,j,k
        FOR_ALL(i, 0, num_bins_x,
                j, 0, num_bins_y,
                k, 0, num_bins_z, {
            

            // get bin gid for this i,j,k
            size_t bin_gid = get_gid(i, j, k, num_bins_x, num_bins_y);

            // the i,j,k for this bin
            bin_keys_t bin_keys;
            bin_keys.i = i;
            bin_keys.j = j;
            bin_keys.k = k;

            // save mapping from bin_gid to bin_keys i,j,k
            keys_in_bin(bin_gid) = bin_keys;

        });
        Kokkos::fence();
        keys_in_bin.update_host();

        // end timer
        auto time_2 = std::chrono::high_resolution_clock::now();


        // -------------------------------------------------------------------
        // below here, these routine must be called every time particles move
        // -------------------------------------------------------------------

        // start timer
        auto time_3 = std::chrono::high_resolution_clock::now();

        printf("building neighbor point list \n");

        // save bin id to points
        FOR_ALL(point_gid, 0, num_points, {

            // get the 1D index for this bin
            size_t bin_gid = get_bin_gid(point_positions(point_gid,0), 
                                         point_positions(point_gid,1), 
                                         point_positions(point_gid,2),
                                         num_bins_x, 
                                         num_bins_y,
                                         num_bins_z);

            size_t storage_lid = Kokkos::atomic_fetch_add(&num_points_in_bin(bin_gid), 1);
            points_bin_gid(point_gid) = bin_gid; // the id of the bin
            points_bin_lid_storage(point_gid) = storage_lid; // the storage place in the bin

        }); // end for all
        Kokkos::fence();
        points_bin_gid.update_host();
        num_points_in_bin.update_host();


        // allocate points in bin connectivity
        points_in_bin = DRaggedRightArrayKokkos <size_t> (num_points_in_bin, "num_points_in_bin");

        // save points in bin
        FOR_ALL(point_gid, 0, num_points, {

            // get bin gid
            size_t bin_gid = points_bin_gid(point_gid);

            // get it's storage location in the ragged right compressed storage
            size_t storage_lid = points_bin_lid_storage(point_gid);

            // save the point to this bin
            points_in_bin(bin_gid, storage_lid) = point_gid;

        }); // end for all


        // for(size_t bin_gid=0; bin_gid<num_bins; bin_gid++){
        //     if(num_points_in_bin.host(bin_gid) > 0){
        //         size_t i = keys_in_bin(bin_gid).i;
        //         size_t j = keys_in_bin(bin_gid).j; 
        //         size_t k = keys_in_bin(bin_gid).k;

        //         double bin_x = ((double)i)*bin_dx;
        //         double bin_y = ((double)j)*bin_dy;
        //         double bin_z = ((double)k)*bin_dz;
        //         //printf("num points in bin = %zu, bin keys = (%zu, %zu, %zu), bin x = (%f, %f, %f) \n", 
        //         //    num_points_in_bin.host(bin_gid), i, j, k, bin_x, bin_y, bin_z);
        //     }
        // } // end for



        // ------------------------------------------------
        // Find the neighbors around each point using bins
        // ------------------------------------------------
        
        FOR_ALL(point_gid, 0, num_points, {

            // get bin gid
            size_t bin_gid = points_bin_gid(point_gid);
            
            // get i,j,k for this bin
            bin_keys_t bin_keys = keys_in_bin(bin_gid);
            // printf(" keys = %zu, %zu, %zu, bin size = %zu, %zu, %zu \n", 
            //     bin_keys.i, bin_keys.j, bin_keys.k,
            //     num_bins_x, num_bins_y, num_bins_z);

            // loop over neighboring bins
            size_t num_points_found;

            // establish the stencil size to get enough particles
            for(int stencil=1; stencil<100000; stencil++){

                num_points_found = 0;

                const int i = bin_keys.i;
                const int j = bin_keys.j;
                const int k = bin_keys.k;

                const int imin = MAX(0, i-stencil);
                const int imax = MIN(num_bins_x-1, i+stencil);

                const int jmin = MAX(0, j-stencil);
                const int jmax = MIN(num_bins_y-1, j+stencil);

                const int kmin = MAX(0, k-stencil);
                const int kmax = MIN(num_bins_z-1, k+stencil);
                    
                for (int kcount=kmin; kcount<=kmax; kcount++){
                    for (int jcount=jmin; jcount<=jmax; jcount++) {
                        for (int icount=imin; icount<=imax; icount++){

                            // get bin neighbor gid 
                            size_t neighbor_bin_gid = get_gid(icount, jcount, kcount, num_bins_x, num_bins_y);
                            num_points_found += num_points_in_bin(neighbor_bin_gid);

                        } // end for kcount
                    } // end for jcount
                } // end for icount

                // the min number of points required to solve the system is num_poly_basis+1, was 2*num_poly_basis
                if (num_points_found >= num_points_fit  || num_points_found==num_points){

                    const double x_pt_middle = bin_dx*((double)i) + X0; 
                    const double y_pt_middle = bin_dy*((double)j) + Y0; 
                    const double z_pt_middle = bin_dz*((double)k) + Z0; 

                    const double x_pt_minus = bin_dx*((double)imin) + X0; 
                    const double y_pt_minus = bin_dy*((double)jmin) + Y0; 
                    const double z_pt_minus = bin_dz*((double)kmin) + Z0; 
                    
                    const double x_pt_plus = bin_dx*((double)imax) + X0; 
                    const double y_pt_plus = bin_dy*((double)jmax) + Y0; 
                    const double z_pt_plus = bin_dz*((double)kmax) + Z0; 

                    const double dist_minus = sqrt( (x_pt_minus - x_pt_middle)*(x_pt_minus - x_pt_middle) +
                                                    (y_pt_minus - y_pt_middle)*(y_pt_minus - y_pt_middle) +
                                                    (z_pt_minus - z_pt_middle)*(z_pt_minus - z_pt_middle) );

                    const double dist_plus = sqrt( (x_pt_plus - x_pt_middle)*(x_pt_plus - x_pt_middle) +
                                                   (y_pt_plus - y_pt_middle)*(y_pt_plus - y_pt_middle) +
                                                   (z_pt_plus - z_pt_middle)*(z_pt_plus - z_pt_middle) );

                    //printf("h = %f, dist_m = %f, dist_p = %f, num_points=%zu, imin = %d, imax = %d, jmin = %d,  jmax = %d, kmin = %d,  kmax = %d, \n", 
                    //         h_kernel, dist_minus, dist_plus, num_points_found, imin, imax, jmin, jmax, kmin, kmax);

                    // only exit when we exceed kernel distance
                    if (dist_minus >= h_kernel || dist_plus >= h_kernel || num_points_found==num_points){

                        //printf("exiting \n\n");

                        points_bin_stencil(point_gid,0) = imin;
                        points_bin_stencil(point_gid,1) = imax;
                        points_bin_stencil(point_gid,2) = jmin;
                        points_bin_stencil(point_gid,3) = jmax;
                        points_bin_stencil(point_gid,4) = kmin;
                        points_bin_stencil(point_gid,5) = kmax;

                        points_num_neighbors(point_gid) = num_points_found; // including node_i in the list of neighbors
                        points_num_neighbors(point_gid) = num_points_found - 1; // the -1 is because we counted point i as a neighbor

                        break;
                    }
                    // else increase stencil size


                } // end of check
                
            } // end for stencil


        }); // end for all
        Kokkos::fence();
        points_bin_stencil.update_host();



        // account for stencels not overlapping, fixing assymetry in points connectivity
        FOR_ALL(point_gid, 0, num_points, {

            // get bin gid for this point
            size_t bin_gid = points_bin_gid(point_gid);
                    
            // get i,j,k for this bin
            bin_keys_t bin_keys = keys_in_bin(bin_gid);

            const int i = bin_keys.i;
            const int j = bin_keys.j;
            const int k = bin_keys.k;

            // walk over the stencil to get neighbors of this bin
            const int imin = points_bin_stencil(point_gid,0);
            const int imax = points_bin_stencil(point_gid,1);
            const int jmin = points_bin_stencil(point_gid,2);
            const int jmax = points_bin_stencil(point_gid,3);
            const int kmin = points_bin_stencil(point_gid,4);
            const int kmax = points_bin_stencil(point_gid,5);

            // loop over my bin stencil
            for (int kcount=kmin; kcount<=kmax; kcount++){
                for (int jcount=jmin; jcount<=jmax; jcount++) {
                    for (int icount=imin; icount<=imax; icount++){

                        // get bin neighbor gid 
                        size_t neighbor_bin_gid = get_gid(icount, jcount, kcount, num_bins_x, num_bins_y);

                        // save the points in this bin
                        for(size_t neighbor_pt_lid=0; neighbor_pt_lid<num_points_in_bin(neighbor_bin_gid); neighbor_pt_lid++){

                            size_t neighbor_point_gid = points_in_bin(neighbor_bin_gid, neighbor_pt_lid);

                            // check if the point-point pairs have identical, overlapping stencils, if not, increment the number of neighbors
                            const int neighbor_imin = points_bin_stencil(neighbor_point_gid,0);
                            const int neighbor_imax = points_bin_stencil(neighbor_point_gid,1);
                            const int neighbor_jmin = points_bin_stencil(neighbor_point_gid,2);
                            const int neighbor_jmax = points_bin_stencil(neighbor_point_gid,3);
                            const int neighbor_kmin = points_bin_stencil(neighbor_point_gid,4);
                            const int neighbor_kmax = points_bin_stencil(neighbor_point_gid,5);
                            
                            // i,j,k is the bin where point_gid lives
                            bool inside =
                                (i >= neighbor_imin && i <= neighbor_imax) &&
                                (j >= neighbor_jmin && j <= neighbor_jmax) &&
                                (k >= neighbor_kmin && k <= neighbor_kmax);

                            if(!inside){
                                Kokkos::atomic_increment(&points_num_neighbors(neighbor_point_gid)); 
                                // the other stencil didn't see my point because it was smaller, now it does see it
                            }

                        } // neighbor_point_lid

                    } // end for kcount
                } // end for jcount
            } // end for icount        

        }); // end for all
        Kokkos::fence();
        points_num_neighbors.update_host();
        
        // allocate memory for points in point
        DRaggedRightArrayKokkos <size_t> points_in_point(points_num_neighbors, "points_in_point");
        points_num_neighbors.set_values(0);  // this is a num saved counter now

        // ---------------------
        // Save the neighbors
        // ---------------------

        // find neighbors using bins
        FOR_ALL(point_gid, 0, num_points, {

            // get bin gid for this point
            size_t bin_gid = points_bin_gid(point_gid);
                    
            // get i,j,k for this bin
            bin_keys_t bin_keys = keys_in_bin(bin_gid);

            const int i = bin_keys.i;
            const int j = bin_keys.j;
            const int k = bin_keys.k;

            // walk over the stencil to get neighbors
            int imin = points_bin_stencil(point_gid,0);
            int imax = points_bin_stencil(point_gid,1);
            int jmin = points_bin_stencil(point_gid,2);
            int jmax = points_bin_stencil(point_gid,3);
            int kmin = points_bin_stencil(point_gid,4);
            int kmax = points_bin_stencil(point_gid,5);


            for (int kcount=kmin; kcount<=kmax; kcount++){
                for (int jcount=jmin; jcount<=jmax; jcount++) {
                    for (int icount=imin; icount<=imax; icount++){

                        // get bin neighbor gid 
                        size_t neighbor_bin_gid = get_gid(icount, jcount, kcount, num_bins_x, num_bins_y);

                        // save the points in this bin
                        for(size_t neighbor_pt_lid=0; neighbor_pt_lid<num_points_in_bin(neighbor_bin_gid); neighbor_pt_lid++){

                            size_t neighbor_point_gid = points_in_bin(neighbor_bin_gid, neighbor_pt_lid);
                            
                            // make sure its a neighbor
                            if(neighbor_point_gid != point_gid){

                                // save the neighbor
                                size_t num_saved = Kokkos::atomic_fetch_add(&points_num_neighbors(point_gid), 1);
                                points_in_point(point_gid, num_saved) = neighbor_point_gid;
                                
                                
                                // if point j's stencil did not see point i, then save i to j's list
                                const int neighbor_imin = points_bin_stencil(neighbor_point_gid,0);
                                const int neighbor_imax = points_bin_stencil(neighbor_point_gid,1);
                                const int neighbor_jmin = points_bin_stencil(neighbor_point_gid,2);
                                const int neighbor_jmax = points_bin_stencil(neighbor_point_gid,3);
                                const int neighbor_kmin = points_bin_stencil(neighbor_point_gid,4);
                                const int neighbor_kmax = points_bin_stencil(neighbor_point_gid,5);

                                // i,j,k is the bin where point_gid lives
                                bool inside =
                                    (i >= neighbor_imin && i <= neighbor_imax) &&
                                    (j >= neighbor_jmin && j <= neighbor_jmax) &&
                                    (k >= neighbor_kmin && k <= neighbor_kmax);

                                if(!inside){

                                    size_t num_saved_neighbor = Kokkos::atomic_fetch_add(&points_num_neighbors(neighbor_point_gid), 1);
                                    points_in_point(neighbor_point_gid, num_saved_neighbor) = point_gid;
                                    // the other stencil didn't see my point because it was smaller, now it does see it

                                } // end if

                            } // end if neighbor != point_gid

                        } // neighbor_point_lid

                    } // end for kcount
                } // end for jcount
            } // end for icount        

        }); // end for all
        Kokkos::fence();
        points_in_point.update_host();


        // build the reverse map
        DRaggedRightArrayKokkos <size_t> reverse_neighbor_lid(points_num_neighbors); 

        FOR_ALL(point_gid, 0, num_points, {
                
            for(int neighbor_point_lid = 0; neighbor_point_lid<points_num_neighbors(point_gid); neighbor_point_lid++){
                
                // get the point gid for this neighbor
                int neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);
                
                // loop over the neighbors of my neighbor
                size_t found = 0;
                for(int j_lid = 0; j_lid<points_num_neighbors(neighbor_point_gid); j_lid++){

                    // get the neighboring point gid of my neighbor
                    int j_point_gid = points_in_point(neighbor_point_gid, j_lid);
                    if (point_gid == j_point_gid){
                        reverse_neighbor_lid(point_gid, neighbor_point_lid) = j_lid;
                        found = 1;
                        //printf("found \n");
                        break;
                    }
                } // end loop over j's neighboring points
                if(found==0)printf("reverse map for i=%d and j=%d pair not found \n", point_gid, neighbor_point_gid);
            } // end loop over i's neighboring points
                
        });


        // end timer
        auto time_4 = std::chrono::high_resolution_clock::now();

        printf("done building neighbor point list \n");

        // ----------------------------------------
        // Find basis that reconstructs polynomial 
        // ----------------------------------------

        printf("Reconstructing basis using point cloud data \n\n");

        auto time_5 = std::chrono::high_resolution_clock::now();

        CArrayKokkos <double> p_coeffs(num_points, num_poly_basis); // reproducing kernel coefficients at each point
        


        CArrayKokkos <double> M_inv(num_points, num_poly_basis, num_poly_basis);
        CArrayKokkos <double> grad_M(num_points, num_poly_basis, num_poly_basis);
        
        DRaggedRightArrayKokkos <double> basis(points_num_neighbors);        // reproducing kernel basis (num_points, num_neighbors)
        DRaggedRightArrayKokkos <double> grad_basis_i(points_num_neighbors,3); // grad kernel basis j with respect to i (num_points, num_neighbors)
        DRaggedRightArrayKokkos <double> grad_basis_j(points_num_neighbors,3); // grad kernel basis i with respect to j (num_points, num_neighbors)
        


        double h = h_kernel; // kernel width


        printf("building reproducing kernel coefficients \n");

        // build coefficients on basis functions
        calc_p_coefficients(point_positions, 
                            points_num_neighbors, 
                            points_in_point, 
                            vol, 
                            p_coeffs, 
                            M_inv,
                            h);
        
        
        calc_basis_and_grad_basis_functions(
                                    point_positions,
                                    points_num_neighbors, 
                                    points_in_point,
                                    vol,
                                    p_coeffs,
                                    M_inv,
                                    basis,
                                    grad_basis_i,
                                    h,
                                    true);
        
        calc_basis_and_grad_basis_functions(
                                    point_positions,
                                    points_num_neighbors, 
                                    points_in_point,
                                    vol,
                                    p_coeffs,
                                    M_inv,
                                    basis,
                                    grad_basis_j,
                                    h,
                                    false);

        // end timer
        auto time_6 = std::chrono::high_resolution_clock::now();


        // -----------------
        //  Timers
        // -----------------
        printf("\n");
        std::chrono::duration <double, std::milli> ms = time_2 - time_1;
        std::cout << "runtime to create bins = " << ms.count() << "ms\n\n";

        ms = time_4 - time_3;
        std::cout << "runtime to find and save neighbors = " << ms.count() << "ms\n\n";

        ms = time_6 - time_5;
        std::cout << "runtime to calculate basis and grad basis = " << ms.count() << "ms\n\n";



        printf("Checking gradients at points \n\n");

        // performing checks on p_coeffs, basis, and grad_basis
        double partion_unity;
        double partion_unity_lcl;

        double linear_preserving;
        double linear_preserving_lcl;

        double quadratic_preserving;
        double quadratic_preserving_lcl;

        double grad_x_p0; 
        double grad_x_p0_lcl; 
        double grad_y_p0; 
        double grad_y_p0_lcl; 
        double grad_z_p0; 
        double grad_z_p0_lcl; 

        double grad_x_p1; 
        double grad_x_p1_lcl; 
        double grad_y_p1; 
        double grad_y_p1_lcl; 
        double grad_z_p1; 
        double grad_z_p1_lcl; 

        double grad_x_p2; 
        double grad_x_p2_lcl; 
        double grad_y_p2; 
        double grad_y_p2_lcl; 
        double grad_z_p2; 
        double grad_z_p2_lcl; 

        // loop over the particles in the domain
        for(size_t point_gid=0; point_gid<num_points; point_gid++){

            // partition of unity
            FOR_REDUCE_SUM(neighbor_point_lid, 0, points_num_neighbors.host(point_gid), partion_unity_lcl, {
                // get the point gid for this neighboring
                size_t neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);
                partion_unity_lcl += basis(point_gid,neighbor_point_lid)*vol(neighbor_point_gid);
            }, partion_unity);
            

            // linear reproducing
            FOR_REDUCE_SUM(neighbor_point_lid, 0, points_num_neighbors.host(point_gid), linear_preserving_lcl, {
                // get the point gid for this neighboring
                size_t neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);
                linear_preserving_lcl += basis(point_gid,neighbor_point_lid)*vol(neighbor_point_gid)*point_positions(neighbor_point_gid,0);
            }, linear_preserving);


            // quadratic reproducing
            FOR_REDUCE_SUM(neighbor_point_lid, 0, points_num_neighbors.host(point_gid), quadratic_preserving_lcl, {
                // get the point gid for this neighboring
                size_t neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);
                quadratic_preserving_lcl += basis(point_gid,neighbor_point_lid)*vol(neighbor_point_gid)*point_positions(neighbor_point_gid,0)*point_positions(neighbor_point_gid,0);
            }, quadratic_preserving);

            if(fabs(partion_unity-1.0)>1e-13)
                printf("partition unity = %f, ", partion_unity);

            if(fabs(linear_preserving-point_positions(point_gid,0))>1e-13)
                printf("linear fcn error = %f, ", fabs(linear_preserving-point_positions(point_gid,0)));

        #if defined(P2)
            if(fabs(quadratic_preserving-point_positions(point_gid,0)*point_positions(point_gid,0))>1e-13)
                printf("quadratic fcn error = %f at i=%zu \n", fabs(quadratic_preserving-point_positions(point_gid,0)*point_positions(point_gid,0)), point_gid);
        #endif

            // -----------------
            // gradient checks
            // -----------------

            // Sum(grad) = [0]; 
            FOR_REDUCE_SUM(neighbor_point_lid, 0, points_num_neighbors.host(point_gid), grad_x_p0_lcl, {
                // get the point gid for this neighboring
                size_t neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);
                grad_x_p0_lcl += grad_basis_i(point_gid,neighbor_point_lid,0)*vol(neighbor_point_gid);
            }, grad_x_p0);

            FOR_REDUCE_SUM(neighbor_point_lid, 0, points_num_neighbors.host(point_gid), grad_y_p0_lcl, {
                // get the point gid for this neighboring
                size_t neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);
                grad_y_p0_lcl += grad_basis_i(point_gid,neighbor_point_lid,1)*vol(neighbor_point_gid);
            }, grad_y_p0);

            FOR_REDUCE_SUM(neighbor_point_lid, 0, points_num_neighbors.host(point_gid), grad_z_p0_lcl, {
                // get the point gid for this neighboring
                size_t neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);
                grad_z_p0_lcl += grad_basis_i(point_gid,neighbor_point_lid,2)*vol(neighbor_point_gid);
            }, grad_z_p0);

            const double grad_check_P0 = fabs(grad_x_p0)+fabs(grad_y_p0)+fabs(grad_z_p0);
            if(0.333*fabs(grad_check_P0)>1e-8)
                printf("error in grad(P0) = %f, %f, %f, \n", grad_x_p0, grad_y_p0, grad_z_p0);


            // Sum(grad(P1)) = [1]; 
            FOR_REDUCE_SUM(neighbor_point_lid, 0, points_num_neighbors.host(point_gid), grad_x_p1_lcl, {
                // get the point gid for this neighboring
                size_t neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);
                grad_x_p1_lcl += grad_basis_i(point_gid,neighbor_point_lid,0)*vol(neighbor_point_gid)*point_positions(neighbor_point_gid,0);
            }, grad_x_p1);
            FOR_REDUCE_SUM(neighbor_point_lid, 0, points_num_neighbors.host(point_gid), grad_y_p1_lcl, {
                // get the point gid for this neighboring
                size_t neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);
                grad_y_p1_lcl += grad_basis_i(point_gid,neighbor_point_lid,1)*vol(neighbor_point_gid)*point_positions(neighbor_point_gid,1);
            }, grad_y_p1);
            FOR_REDUCE_SUM(neighbor_point_lid, 0, points_num_neighbors.host(point_gid), grad_z_p1_lcl, {
                // get the point gid for this neighboring
                size_t neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);
                grad_z_p1_lcl += grad_basis_i(point_gid,neighbor_point_lid,2)*vol(neighbor_point_gid)*point_positions(neighbor_point_gid,2);
            }, grad_z_p1);

            const double grad_check_P1 = fabs(grad_x_p1 - 1.0)+fabs(grad_y_p1 - 1.0)+fabs(grad_z_p1 - 1.0);
            if(0.333*fabs(grad_check_P1)>1e-8){
                printf("error in grad(P1) = %f, %f, %f, \n", 
                    fabs(grad_x_p1 - 1.0), 
                    fabs(grad_y_p1 - 1.0), 
                    fabs(grad_z_p1 - 1.0));
            }


        #if defined(P2)
            // Sum(grad(P2)) = [2]; 
            FOR_REDUCE_SUM(neighbor_point_lid, 0, points_num_neighbors.host(point_gid), grad_x_p2_lcl, {
                // get the point gid for this neighboring
                size_t neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);
                grad_x_p2_lcl += grad_basis_i(point_gid,neighbor_point_lid,0)*vol(neighbor_point_gid)*point_positions(neighbor_point_gid,0)*point_positions(neighbor_point_gid,0);
            }, grad_x_p2);
            FOR_REDUCE_SUM(neighbor_point_lid, 0, points_num_neighbors.host(point_gid), grad_y_p2_lcl, {
                // get the point gid for this neighboring
                size_t neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);
                grad_y_p2_lcl += grad_basis_i(point_gid,neighbor_point_lid,1)*vol(neighbor_point_gid)*point_positions(neighbor_point_gid,1)*point_positions(neighbor_point_gid,1);
            }, grad_y_p2);
            FOR_REDUCE_SUM(neighbor_point_lid, 0, points_num_neighbors.host(point_gid), grad_z_p2_lcl, {
                // get the point gid for this neighboring
                size_t neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);
                grad_z_p2_lcl += grad_basis_i(point_gid,neighbor_point_lid,2)*vol(neighbor_point_gid)*point_positions(neighbor_point_gid,2)*point_positions(neighbor_point_gid,2);
            }, grad_z_p2);

            const double grad_check_P2 = fabs(grad_x_p2-2.0*point_positions(point_gid,0)) + fabs(grad_y_p2-2.0*point_positions(point_gid,1)) + fabs(grad_z_p2-2.0*point_positions(point_gid,2));
            if(0.333*fabs(grad_check_P2)>1e-8){
                printf("error in grad(P2) = %f, %f, %f, \n", 
                        fabs(grad_x_p2-2.0*point_positions(point_gid,0)), 
                        fabs(grad_y_p2-2.0*point_positions(point_gid,1)), 
                        fabs(grad_z_p2-2.0*point_positions(point_gid,2)));
            }
        #endif

        } // end for point gid


        if(check_maps){
            size_t bad = 0;
            for (size_t i=0; i<num_points; ++i) {
                for (size_t lid=0; lid<points_num_neighbors(i); ++lid) {
                    size_t j = points_in_point(i, lid);
                    size_t rev = reverse_neighbor_lid(i, lid);
                    if (rev == SIZE_MAX) {
                    printf("MISSING reverse for pair (%zu,%zu)\n", i, j); ++bad;
                    if (bad>50) break;
                    } else {
                    size_t check = points_in_point(j, rev);
                    if (check != i) {
                        printf("WRONG reverse: i=%zu j=%zu rev=%zu check=%zu\n", i, j, rev, check);
                        ++bad;
                        if (bad>50) break;
                    }
                    }
                }
                if (bad>50) break;
            }
            if (bad==0) printf("reverse map OK\n");


            for (size_t i=0; i<num_points; ++i) {
                std::set<size_t> seen;
                for (size_t lid=0; lid<points_num_neighbors(i); ++lid) {
                    size_t j = points_in_point(i, lid);
                    //if (j == i) printf("SELF neighbor found at i=%zu lid=%zu\n", i, lid);
                    if (!seen.insert(j).second) printf("DUPLICATE neighbor %zu in list of %zu\n", j, i);
                }
            }
        }


        printf("Building anti-symmetric gradient \n\n");

        // -------------------
        // Anti-sym gradient
        // -------------------
        printf("Testing divergence of vector field u = (x, y, z) \n\n");

        DCArrayKokkos <double> u(num_points, 3);
        FOR_ALL(i, 0, num_points, {
            u(i, 0) = point_positions(i, 0);        
            u(i, 1) = point_positions(i, 1);
            u(i, 2) = point_positions(i, 2);
        });
        u.update_device();

        DCArrayKokkos <double> div(num_points);
        div.set_values(0.0);

        DCArrayKokkos <double> div_fd(num_points);
        div_fd.set_values(0.0);

        FOR_ALL(i_gid, 0, num_points, {

            for(size_t j_lid = 0; j_lid<points_num_neighbors(i_gid); j_lid++){     

                size_t j_gid = points_in_point(i_gid, j_lid);
                size_t i_lid = reverse_neighbor_lid(i_gid, j_lid);

                if(i_gid != points_in_point(j_gid, i_lid)){
                    printf("CHECK: point i = %d, reverse map point i = %zu for j = %zu \n", i_gid, points_in_point(j_gid, i_lid), j_gid);
                }
                //printf("map check: edge points (%d,%zu), rev from j = %zu using i_lid = %zu \n", i_gid, j_gid,  points_in_point(j_gid, i_lid), i_lid);

                double g_ij[3];
                double g_ji[3];
                double sum[3];
                for (int dim=0; dim<3; ++dim) {
                    g_ij[dim] = grad_basis_i(i_gid,j_lid,dim);
                    g_ji[dim] = grad_basis_j(i_gid,j_lid,dim);
                    sum[dim] = g_ij[dim] + g_ji[dim];
                }
                double norm = sqrt(sum[0]*sum[0] + sum[1]*sum[1] + sum[2]*sum[2]);
                //printf("errors: norm=%f, g_ij = (%f, %f, %f), g_ji = (%f, %f, %f) \n", norm, g_ij[0], g_ij[1], g_ij[2], g_ji[0], g_ji[1], g_ji[2]);

                // conservative mesh-free FE
                double contrib = 0.0;
                for (int dim=0; dim<3; ++dim) {
                    contrib += 0.5*(g_ij[dim] - g_ji[dim]) * (u(j_gid, dim) + u(i_gid, dim));
                }
                div(i_gid) += vol(i_gid) * vol(j_gid) * contrib;

                // finite difference
                contrib = 0.0;
                for (int dim=0; dim<3; ++dim) {
                    contrib += g_ij[dim]*u(j_gid, dim)*vol(j_gid);
                }
                div_fd(i_gid) += contrib;

            } // end loop over neighbors

            div(i_gid) /= vol(i_gid);
            // remember: finite difference doesn't have the V_i on the right side, so no division
        });
        div.update_host();
        div_fd.update_host();


        for(size_t point_gid=0; point_gid<num_points; point_gid++){
            double error = fabs(div.host(point_gid) - 3.0);
            if(error > 1e-8){
                printf("div(u) = %f at point %zu, error = %g, vol = %f\n", div.host(point_gid), point_gid, error, vol.host(point_gid));
            }
        } // end for point_gid

        for(size_t point_gid=0; point_gid<num_points; point_gid++){
            double error = fabs(div_fd.host(point_gid) - 3.0);
            if(error > 1e-8){
                printf("div_fd(u) = %f at point %zu, error = %g\n", div_fd.host(point_gid), point_gid, error);
            }
        } // end for point_gid


        double conserve_check;
        double conserve_check_lcl;
        FOR_REDUCE_SUM(point_gid, 0, num_points, 
                       conserve_check_lcl, {

            for(size_t neighbor_point_lid = 0; neighbor_point_lid<points_num_neighbors(point_gid); neighbor_point_lid++){
                
                // get the point gid for this neighbor
                size_t neighbor_point_gid = points_in_point(point_gid, neighbor_point_lid);

                // get the local id of my neighbor that matches my point_gid
                size_t neighbor_lid = reverse_neighbor_lid(point_gid, neighbor_point_lid);

                for (size_t dim=0; dim<3; dim++){
                    // from point i
                    conserve_check_lcl += 0.5*(grad_basis_i(point_gid,neighbor_point_lid,dim) - grad_basis_j(point_gid,neighbor_point_lid,dim));

                    // from neighbor point
                    conserve_check_lcl += 0.5*(grad_basis_i(neighbor_point_gid,neighbor_lid,dim) - grad_basis_j(neighbor_point_gid,neighbor_lid,dim));
                }
            }

        }, conserve_check);
        printf("conservation error = %f \n\n", conserve_check);



        printf("Testing sin(u) of vector field u = (x, y, z) \n\n");

        FOR_ALL(i, 0, num_points, {
            u(i, 0) = sin(point_positions(i, 0));        
            u(i, 1) = sin(point_positions(i, 1));
            u(i, 2) = sin(point_positions(i, 2));
        });
        u.update_device();

        DCArrayKokkos <double> div_sinx(num_points);
        div_sinx.set_values(0.0);

        DCArrayKokkos <double> div_fd_sinx(num_points);
        div_fd_sinx.set_values(0.0);

        FOR_ALL(i_gid, 0, num_points, {

            for(size_t j_lid = 0; j_lid<points_num_neighbors(i_gid); j_lid++){     

                size_t j_gid = points_in_point(i_gid, j_lid);
                size_t i_lid = reverse_neighbor_lid(i_gid, j_lid);

                if(i_gid != points_in_point(j_gid, i_lid)){
                    printf("CHECK: point i = %d, reverse map point i = %zu for j = %zu \n", i_gid, points_in_point(j_gid, i_lid), j_gid);
                }
                //printf("map check: edge points (%d,%zu), rev from j = %zu using i_lid = %zu \n", i_gid, j_gid,  points_in_point(j_gid, i_lid), i_lid);

                double g_ij[3];
                double g_ji[3];
                double sum[3];
                for (int dim=0; dim<3; ++dim) {
                    g_ij[dim] = grad_basis_i(i_gid,j_lid,dim);
                    g_ji[dim] = grad_basis_j(i_gid,j_lid,dim);
                    sum[dim] = g_ij[dim] + g_ji[dim];
                }
                double norm = sqrt(sum[0]*sum[0] + sum[1]*sum[1] + sum[2]*sum[2]);
                //printf("errors: norm=%f, g_ij = (%f, %f, %f), g_ji = (%f, %f, %f) \n", norm, g_ij[0], g_ij[1], g_ij[2], g_ji[0], g_ji[1], g_ji[2]);

                // conservative mesh-free FE
                double contrib = 0.0;
                for (int dim=0; dim<3; ++dim) {
                    contrib += 0.5*(g_ij[dim] - g_ji[dim]) * (u(j_gid, dim) + u(i_gid, dim));
                }
                div_sinx(i_gid) += vol(i_gid) * vol(j_gid) * contrib;

                // finite difference
                contrib = 0.0;
                for (int dim=0; dim<3; ++dim) {
                    contrib += g_ij[dim]*u(j_gid, dim)*vol(j_gid);
                }
                div_fd_sinx(i_gid) += contrib;

            } // end loop over neighbors

            div_sinx(i_gid) /= vol(i_gid);
            // remember: finite difference doesn't have the V_i on the right side, so no division
        });
        div_sinx.update_host();
        div_fd_sinx.update_host();



       







        printf("Writing VTK Graphics File \n\n");

        std::ofstream out("cloud.vtk");

        out << "# vtk DataFile Version 3.0\n";
        out << "3D point cloud\n";
        out << "ASCII\n";
        out << "DATASET POLYDATA\n";
        out << "POINTS " << num_points << " float\n";
        for (size_t point_gid = 0; point_gid < num_points; ++point_gid) {
            out << point_positions.host(point_gid,0) << " " 
                << point_positions.host(point_gid,1) << " " 
                << point_positions.host(point_gid,2) << "\n";
        }

        out << "\nPOINT_DATA " << num_points << "\n";
        out << "SCALARS error_div(x) float 1\n";
        out << "LOOKUP_TABLE default\n";
        for (size_t point_gid = 0; point_gid < num_points; ++point_gid) {
            out << div.host(point_gid)-3 << "\n";
        }

        out << "SCALARS error_div(sin(x)) float 1\n";
        out << "LOOKUP_TABLE default\n";
        for (size_t point_gid = 0; point_gid < num_points; ++point_gid) {
            // vec[0] = sin(x), dvec[0]/dx = cos(x)
            // vec[1] = sin(y), dvec[1]/dy = cos(y)
            // vec[2] = sin(z), dvec[2]/dz = cos(z)
            double val0 = cos(point_positions.host(point_gid,0));
            double val1 = cos(point_positions.host(point_gid,1));
            double val2 = cos(point_positions.host(point_gid,2));
            double exact_div = val0 + val1 + val2;
            out << (div_sinx.host(point_gid) - exact_div) << "\n";
        }

        printf("Finished \n\n");



    } // end of kokkos scope


    Kokkos::finalize();



    return 0;
    
} // end main
