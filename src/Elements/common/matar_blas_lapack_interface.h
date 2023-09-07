#pragma once

#include "c_blas_lapack_interface.h"

#include "error.h"
#include "common.h"

// Interface to selected BLAS routines for MATAR CArrays
namespace matar2blas {
  void transpose(const CArray<Real> &, CArray<Real> &);

  void matvec(const CArray<Real> &, 
      const CArray<Real> &, CArray<Real> &);

  void matmul(const CArray<Real> &, 
      const CArray<Real> &, CArray<Real> &);
}

// Interface to selected LAPACK routines for MATAR CArrays
namespace matar2lapack {
  void invert(const CArray<Real> &, CArray<Real> &);

  void eig_sym_tri(const CArray<Real> &diag, 
    const CArray<Real> &subdiag, CArray<Real> &eigvals, 
    CArray<Real> &eigvecs);
}
