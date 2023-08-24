#pragma once
#include "common.h"

/* 
 * BLAS/LAPACK naming convention help
 *
 * GEMV  - general matrix-vector multiplication
 * GEMM  - general matrix-matrix multiplication
 * GETRF - general triangular (LU) factorization
 * GETRS - general triangular (LU) solution
 * STEV  - symmetric tridiagonal eigenvalues/vectors
 */

// Change definitions of routines
#define blas_gemv dgemv_
#define blas_gemm dgemm_       
#define lapack_getrf dgetrf_
#define lapack_getrs dgetrs_
#define lapack_stev dstev_

extern "C" {
  void dgemv_( 	
    const char *,
		int *, int *,
		Real *,
		Real *, int *,
		Real *, int *,
		Real *,
		Real *, int * 	
	);

  void dgemm_(
    const char *, const char *, 
    int	*, int	*, int	*, 
    Real *, Real *, int *, 
    Real *, int *, 
    Real *, Real *, int *
  );

  void dgetrf_(int*, int *, Real *, int	*, int *, int	*);

  void dgetrs_(
    const char *, 
    int *, int *, 
    Real *, int *, 
    int *, 
    Real *, int *, 
    int *);

  void dstev_(
    const char *,
		int *,
		Real *,
		Real *,
		Real *, int *,
		Real *,
		int *
	);
}
