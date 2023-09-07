#include "jacobi_polynomials.h"

namespace jacobi {
  template <typename NumType>
  NumType eval(const int n, const Real alpha, const Real beta, 
      const NumType X) {
    if (n == -1) return 0.0;
    if (n == 0) return 1.0;
    return (a(alpha, beta, n)*X + b(alpha, beta, n))
        *eval(n - 1, alpha, beta, X) 
        - c(alpha, beta, n)*eval(n - 2, alpha, beta, X);
  }

  template <typename NumType>
  NumType eval_der(const int n, const int k, const Real alpha, 
      const Real beta, const NumType X) {
    if (n == 0 || k > n) return 0.0;
    Real theta = Real(n) + Real(alpha + alpha) + 1.0;
    return std::pow(0.5, k)*std::tgamma(theta + Real(k))
        *eval(n - k, alpha + Real(k), beta + Real(k), X)/std::tgamma(theta);
  }

  // Explicit instantiations of template functions
  template Real eval(const int n, const Real alpha, const Real beta, 
      const Real X);
  template Complex eval(const int n, const Real alpha, const Real beta, 
      const Complex X);

  template Real eval_der(const int n, const int k, const Real alpha, 
      const Real beta, const Real X);
  template Complex eval_der(const int n, const int k, const Real alpha, 
      const Real beta, const Complex X);
}
