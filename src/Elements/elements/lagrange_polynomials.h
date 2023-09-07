#pragma once

#include "common.h"

namespace lagrange {
  // Identifying singularities
  template <typename NumType>
  SizeType find_coincident_vertex(const SizeType &, const NumType *, 
      const NumType &);

  // Barycentric weights
  template <typename NumType>
  void compute_barycentric_weights(const SizeType &, const NumType *, 
      NumType *);

  // Lagrange polynomials (basis functions)
  template <typename NumType>
  NumType eval(const SizeType Nv, const SizeType i, const SizeType ic, 
      const NumType *Z, const NumType *w, const NumType X);

  template <typename NumType>
  NumType eval_der(const SizeType Nv, const SizeType n, 
      const SizeType i, const SizeType ic, const NumType *Z, 
      const NumType *w, const NumType X, NumType *c);

  // Lagrange interpolants (sums of products of bases and coefficients)
  template <typename NumType>
  NumType eval_interp(const SizeType Nv, const SizeType i, 
      const NumType *Z, const NumType *w, const NumType X, const NumType *c);

  template <typename NumType>
  NumType eval_der_interp(const SizeType Nv, const SizeType n, 
      const SizeType i, const NumType *Z, const NumType *w, const NumType X, 
      NumType *c);
}
