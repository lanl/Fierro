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
#ifndef LIN_ALG_H
#define LIN_ALG_H

#include <iostream>  // std::cout etc.

#include "matar.h"

////////////////////////////
// Function Declarations  //
////////////////////////////

/* "small number" to avoid overflow in some cases */
#define TINY 1.e-30

/* ------------------------- */
/* LU decomposition function */
/* ------------------------- */

/*
LU-decomposition according to Crout's algorithm with pivoting.
Description:
int LU_decompos(
    Matrix &source_mat,  // matrix to invert
    Vector &indx,        // permutations
    int &parity,         // parity (+1 0r -1)
    const int n);        // matrix size

Parameters:
source_mat - source matrix (n x n) on input, destination on output;
n - the matrix size;
indx - integer array (size n) to remember permutations;
d - on output, contains +1 or -1 for even or odd permutations number.
vv - temporary array (size n).
Returns:
0 - the source matrix is singular (invalid for decomposition),
1 - if OK.
*/
template<typename T1, typename T2>
int LU_decompos(
    T1&  source_mat, // matrix to invert
    T2&  indx,       // permutations
    int& parity,     // parity (+1 0r -1)
    const int n);    // matrix size

/* ----------------------------- */
/* LU back substitution function */
/* ----------------------------- */
/*

Back substitution, using LU decomposed matrix.
Description:
void LU_backsub(
    Matrix &source_mat,  // input matrix
    Vector &indx,        // permutations
    Vector &coeffs,      // Least squares coefficents
    const int n);

Parameters:
source_mat - the matrix decomposed by Crout;
n - the matrix size;
indx - permutation order obtained by decomposition algorithm;
coeffs - the vector (size n) to be substituted on input, the result
of the substitution on output.
Note: a and indx are not modified by this routine and could be
used in multiple calls.
*/
template<typename T1, typename T2, typename T3>
void LU_backsub(
    T1& source_mat,  // input matrix
    T2& indx,        // permutations
    T3& coeffs,      // Least squares coefficents
    const int n);

/* ------------------- */
/* LU  invert function */
/* ------------------- */
/*
Invertation of matrix, using LU decomposed matrix.
Description:
void LU_invert(
    Matrix &source_mat,  // input matrix
    Vector &indx,        // permutations
    Matrix &inv_mat,     // inverse matrix
    Vector&col,          // tmp array
    int n)

Parameters:
source_mat - the matrix decomposed by Crout;
n - the matrix size;
indx - permutation order obtained by decomposition algorithm;
inv_mat - the destination matrix;
col - temporary array (size n).
Note: test for singularity has been already obtained on the
matrix decomposition, a and indx are not modified by this routine,
the routine uses multiple backsubstitutions (previous algorithm).
*/
template<typename T1, typename T2, typename T3, typename T4>
void LU_invert(
    T1& source_mat,  // input matrix
    T2& indx,        // permutations
    T3& inv_mat,     // inverse matrix
    T4& col,         // tmp array
    int n);          // size of matrix

/* ------------------------ */
/* LU  determinant function */
/* ------------------------ */
/*
Obtaining the matrix determinant, using LU-decomposed matrix
Description:
double LU_determ(
    Matrix &source_mat,  // input matrix
    Vector &indx,        // permutations
    const int parity,    // parity (+1 0r -1)
    const int n)         // size of source matrix

Parameters:
source_mat - the matrix decomposed by Crout;
n - the matrix size;
indx - permutation order obtained by decomposition algorithm;
parity - the parity sign (+1 or -1) obtained at decomposition.
Returns:
the determinant value. Note: non-zero (the matrix cannot be
singular, if decomposed properly); a, indx and parity are not modified
by this routine.
*/
template<typename T1, typename T2>
double LU_determ(
    T1& source_mat,     // input matrix
    T2& indx,           // permutations
    const int parity,   // parity (+1 0r -1)
    const int n);       // size of source matrix

////////////////////////////
// Function Definitions   //
////////////////////////////

/* ------------------------- */
/* LU decomposition function */
/* ------------------------- */
template<typename T1, typename T2>
int LU_decompos(
    T1&  source_mat, // matrix to invert
    T2&  indx,       // permutations
    int& parity,     // parity (+1 0r -1)
    const int n)     // matrix size
{
    int i, imax, k; // Indexing

    double big, sum, temp; // useful storage

    // double vv[n];       // temp arrary for solver

    DCArrayKokkos<double> vv = DCArrayKokkos<double>(n);

    parity = 1;
    /* search for the largest element in each row; save the scaling in the
    temporary array vv and return zero if the matrix is singular */
    FOR_ALL(i, 0, n, {
        big = 0.;
        for (int j = 0; j < n; j++) {
            if ((temp = fabs(source_mat(i, j))) > big) {
                big = temp;
            }
        }

        if (big == 0.) {
            return(0);
        }

        vv(i) = big;
    });

    /* the main loop for the Crout's algorithm */
    for (int j = 0; j < n; j++) {
        /* this is the part a) of the algorithm except for i==j */
        for (i = 0; i < j; i++) {
            sum = source_mat(i, j);

            for (k = 0; k < i; k++) {
                sum -= source_mat(i, k) * source_mat(k, j);
            }

            source_mat(i, j) = sum;
        }

        /* initialize for the search for the largest pivot element */
        big  = 0.;
        imax = j;

        /* this is the part a) for i==j and part b) for i>j + pivot search */
        for (i = j; i < n; i++) {
            sum = source_mat(i, j);

            for (k = 0; k < j; k++) {
                sum -= source_mat(i, k) * source_mat(k, j);
            }

            source_mat(i, j) = sum;

            /* is the figure of merit for the pivot better than the best so far? */
            if ((temp = vv(i) * fabs(sum)) >= big) {
                big = temp; imax = i;
            }
        }

        /* interchange rows, if needed, change parity and the scale factor */
        if (imax != j) {
            for (k = 0; k < n; k++) {
                temp = source_mat(imax, k);
                source_mat(imax, k) = source_mat(j, k);
                source_mat(j, k)    = temp;
            }

            parity   = -(parity);
            vv(imax) = vv(j);
        }

        /* store the index */
        indx(j) = imax;
        /* if the pivot element is zero, the matrix is singular but for some
        applications a tiny number is desirable instead */

        if (source_mat(j, j) == 0.) {
            source_mat(j, j) = TINY;
        }
        /* finally, divide by the pivot element */

        if (j < n - 1) {
            temp = 1. / source_mat(j, j);
            for (i = j + 1; i < n; i++) {
                source_mat(i, j) *= temp;
            }
        }
    }

    return(1);
}

/* ----------------------------- */
/* LU back substitution function */
/* ----------------------------- */

template<typename T1, typename T2, typename T3>
void LU_backsub(
    T1& source_mat,  // input matrix
    T2& indx,        // permutations
    T3& coeffs,      // Least squares coefficents
    const int n)     // Matrix size
{
    int    i, j, ip, ii = -1;
    double sum;

    /* First step of backsubstitution; the only wrinkle is to unscramble
    the permutation order. Note: the algorithm is optimized for a
    possibility of large amount of zeroes in b */

    for (i = 0; i < n; i++) {
        ip = indx(i);

        sum = coeffs(ip);
        coeffs(ip) = coeffs(i);

        if (ii >= 0) {
            for (j = ii; j < i; j++) {
                sum -= source_mat(i, j) * coeffs(j);
            }
        }
        else if (sum) {
            ii = i;             /* a nonzero element encounted */
        }
        coeffs(i) = sum;
    }

    /* the second step */
    for (i = n - 1; i >= 0; i--) {
        sum = coeffs(i);
        for (j = i + 1; j < n; j++) {
            sum -= source_mat(i, j) * coeffs(j);
        }

        coeffs(i) = sum / source_mat(i, i);
    }
}

/* ------------------ */
/* LU invert function */
/* ------------------ */
template<typename T1, typename T2, typename T3, typename T4>
void LU_invert(
    T1& source_mat,  // input matrix
    T2& indx,        // permutations
    T3& inv_mat,     // inverse matrix
    T4& col,         // tmp array
    int n)           // size of matrix
{
    for (int j = 0; j < n; j++) {

        // Initialize col to zero
        FOR_ALL(i, 0, n, {
            col(i) = 0.0;
        });

        col(j) = 1.0;
        LU_backsub(source_mat, indx, col, n);

        FOR_ALL(i, 0, n, {
            inv_mat(i, j) = col(i);
        });
    }
}

/* ----------------------- */
/* LU determinant function */
/* ----------------------- */
template<typename T1, typename T2>
double LU_determ(
    T1& source_mat,     // input matrix
    T2& indx,           // permutations
    const int parity,   // parity (+1 0r -1)
    const int n)        // size of source matrix
{
    int    j;
    double res = (double)(parity);

    for (j = 0; j < n; j++) {
        res *= source_mat(j, j);
    }

    return(res);
}

#endif // LIN_ALG_H
