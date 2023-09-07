#pragma once

#include "common.h"

namespace legendre {
  // Legendre polynomials (basis functions)
  template <typename NumType>
  NumType eval(const int n, const NumType X);

  template <typename NumType>
  NumType eval_der(const int n, const int k, const NumType X);

  // Legendre approximations (sums of products of bases and coefficients)
  template <typename NumType>
  NumType eval_approx(const SizeType N, const NumType *c, const NumType X);

  template <typename NumType>
  NumType eval_der_approx(const SizeType N, const SizeType k, const NumType *c, 
      const NumType X);
}
