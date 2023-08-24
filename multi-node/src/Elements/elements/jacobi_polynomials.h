#pragma once

#include "common.h"

namespace jacobi {
  /* Recurrence relation parameters */
  inline Real a(Real alpha, Real beta, int n) {
    if (n == 1) return 0.5*(alpha + beta) + 1.0;
    return (2.0*Real(n) + alpha + beta - 1.0)*(2.0*Real(n) + alpha + beta)
        /(2.0*Real(n)*(Real(n) + alpha + beta));
  };

  inline Real b(Real alpha, Real beta, int n) {
    if (n == 1) return 0.5*(alpha - beta); 
    return (alpha*alpha - beta*beta)*(2*Real(n) + alpha + beta - 1.0)
        /(2.0*Real(n)*(Real(n) + alpha + beta)
            *(2.0*Real(n) + alpha + beta - 2.0));
  };

  inline Real c(Real alpha, Real beta, int n) {
    if (n == 1) return 0.0;
    return (Real(n) + alpha - 1.0)*(Real(n) + beta - 1.0)
       *(2.0*Real(n) + alpha + beta)
       /(Real(n)*(Real(n) + alpha + beta)
           *(2.0*Real(n) + alpha + beta - 2.0));
  };

  /* Polynomial and polynomial derivative evaluation */
  template <typename NumType>
  NumType eval(int n, Real alpha, Real beta, NumType X);

  template <typename NumType>
  NumType eval_der(const int n, const int k, const Real alpha, 
      const Real beta, const NumType X);
}
