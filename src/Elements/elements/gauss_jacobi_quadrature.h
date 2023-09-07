#pragma once

#include "jacobi_polynomials.h"
#include "matar_blas_lapack_interface.h"
#include "error.h"
#include "common.h"

void compute_gauss_jacobi_quadrature_rule(SizeType n, Real alpha, 
    Real beta, CArray<Real> &points, CArray<Real> &weights);
