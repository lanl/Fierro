#include "legendre_polynomials.h"
#include "jacobi_polynomials.h"

namespace legendre {
  /*
   * Return evaluation of Legendre polynomial of specified order at given
   * coordinate
   *
   * Parameters
   * ----------
   * n : order of Legendre polynomial
   * X : coordinate in reference space [-1, 1]
   */
  template <typename NumType>
  NumType eval(const int n, const NumType X) {
    if (n == -1) return 0.0;
    if (n == 0) return 1.0;
    return (1.0/double(n))*((2.0*double(n) - 1.0)*X*eval(n - 1, X)
        - (double(n) - 1.0)*eval(n - 2, X));
  };

  /*
   * Return evaluation of Legendre polynomial of specified order at given
   * coordinate
   *
   * Parameters
   * ----------
   * n : order of Legendre polynomial
   * k : order of derivative
   * X : coordinate in reference space [-1, 1]
   */
  template <typename NumType>
  NumType eval_der(const int n, const int k, const NumType X) {
    return jacobi::eval_der(n, k, 0.0, 0.0, X);
  };

  /*
   * Return evaluation of Legendre polynomial approximation, which is the sum
   * of the products of Legendre polynomials and provided coefficients, at a
   * specified coordinate
   * 
   * Parameters
   * ----------
   * N : maximum polynomial order
   * c : coefficients
   * X : coordinate in reference space [-1, 1]
   */
  template <typename NumType>
  NumType eval_approx(const SizeType N, const NumType *c, const NumType X) {
    NumType sum = 0.0;
    for (SizeType j = 0; j < N; j++) {
      sum += c[j]*eval(j, X);
    }

    return sum;
  }

  /*
   * Return evaluation of derivative of Legendre polynomial approximation,
   * which is the sum of the products of the derivatives of the Legendre
   * polynomials and provided coefficients, at a specified coordinate
   * 
   * Parameters
   * ----------
   * N : maximum polynomial order
   * c : coefficients
   * X : coordinate in reference space [-1, 1]
   */
  template <typename NumType>
  NumType eval_der_approx(const SizeType N, const SizeType k, const NumType *c, 
      const NumType X) {
    NumType sum = 0.0;
    for (SizeType j = 0; j < N; j++) {
      sum += c[j]*eval_der(j, k, X);
    }

    return sum;
  }

  // Explicit instatiations of template functions
  template Real eval(const int n, const Real X);
  template Complex eval(const int n, const Complex X);

  template Real eval_der(const int n, const int k, const Real X);
  template Complex eval_der(const int n, const int k, const Complex X);

  template Real eval_approx(const SizeType N, const Real *c, const Real X);
  template Complex eval_approx(const SizeType N, const Complex *c, 
      const Complex X);

  template Real eval_der_approx(const SizeType N, const SizeType k, 
      const Real *c, const Real X);
  template Complex eval_der_approx(const SizeType N, const SizeType k, 
      const Complex *c, const Complex X);
}
