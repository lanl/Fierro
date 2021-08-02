#ifndef SLAM_H
#define SLAM_H 

#include <iostream>  // std::cout etc.
#include <math.h> 

#include "utilities.h"
#include "matar.h"


using namespace utils;

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
    ViewCArray <real_t> &source_mat,  // matrix to invert
    ViewCArray <int> &indx,           // permutations
    int &parity,                        // parity (+1 0r -1)
    const int n);                       // matrix size                   


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
int LU_decompos(
    ViewCArray <real_t> &source_mat,  // matrix to invert
    ViewCArray <int> &indx,           // permutations
    int &parity,                        // parity (+1 0r -1)
    const int n);                       // matrix size  




/* ----------------------------- */
/* LU back substitution function */
/* ----------------------------- */
/* 

Back substitution, using LU decomposed matrix.
Description:
void LU_backsub(
    ViewCArray <real_t> &source_mat,  // input matrix
    ViewCArray <int> &indx,           // permutations
    ViewCArray <real_t> &coeffs,      // Least squares coefficents
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
void LU_backsub(
    ViewCArray <real_t> &source_mat,  // input matrix
    const ViewCArray <int> &indx,           // permutations
    ViewCArray <real_t> &b,           // Least squares coefficents
    const int n);




/* ------------------- */
/* LU  invert function */
/* ------------------- */
/* 
Invertation of matrix, using LU decomposed matrix.
Description:
void LU_invert(
    ViewCArray <real_t> &source_mat,  // input matrix
    ViewCArray <int> &indx,           // permutations
    ViewCArray <real_t> &inv_mat,     // inverse matrix
    ViewCArray <real_t> &col,         // tmp array
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
void LU_invert(
    ViewCArray <real_t> &source_mat,  // input matrix
    ViewCArray <int> &indx,           // permutations
    ViewCArray <real_t> &inv_mat, 	  // inverse matrix
    ViewCArray <real_t> &col, 		  // tmp array
    int n);



/* ------------------------ */
/* LU  determinant function */
/* ------------------------ */
/* 
Obtaining the matrix determinant, using LU-decomposed matrix
Description:
double LU_determ(
    ViewCArray <real_t> &source_mat,  // input matrix
    ViewCArray <int> &indx,           // permutations
    const int parity,                   // parity (+1 0r -1)
    const int n)                        // size of source matrix
        
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
double LU_determ(
    ViewCArray <real_t> &source_mat, 	// input matrix
    ViewCArray <int> &indx,           // permutations
    const int parity,					// parity (+1 0r -1)
    const int n);








#endif // SLAM_H 
