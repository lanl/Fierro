#include "gauss_jacobi_quadrature.h"

/* 
 * For a quadrature rule based on the interpolation of the integrand by a
 * n-degree Jacobi polynomial, compute the quadrature points and weights using
 * the algorithm presented in "Calculation of Guass Quadrature Rules" (Golub
 * and Welsh, 1969)
 */
void compute_gauss_jacobi_quadrature_rule(
    size_t n, Real alpha, Real beta, 
    CArray<Real> &points, CArray<Real> &weights) {
  assert(n > 0 and "Error: quadrature rule must have nonzero number of points");

  if (n > 1) {
    // Populate arrays corresponding to the diagonal and subdiagonal (or
    // superdiagonal) of the symmetric tridiagonal matrix
    CArray<Real> diag(n);
    CArray<Real> subdiag(n-1);

    for (int k = 0; k < n; k++) {
      Real 
      a_k = common::real(jacobi::a(alpha, beta, k+1)),
      b_k = common::real(jacobi::b(alpha, beta, k+1));

      diag(k) = -b_k/a_k;
    }

    for (int k = 0; k < n-1; k++) {
      Real 
      a_k   = common::real(jacobi::a(alpha, beta, k+1)),
      a_kp1 = common::real(jacobi::a(alpha, beta, k+2)),
      c_kp1 = common::real(jacobi::c(alpha, beta, k+2));

      subdiag(k) = std::sqrt(c_kp1/(a_k*a_kp1));
    }

    // Compute eigensolution of symmetric tridiagonal matrix
    CArray<Real> eigvals(n);
    CArray<Real> eigvecs(n,n);

    matar2lapack::eig_sym_tri(diag, subdiag, eigvals, eigvecs);

    // Compute zeroth moment (integral of weight function times 1); this integral
    // can be calculated analytically using the beta function (discovered this in
    // the SciPy source code); I implemented the beta function here in terms of
    // the gamma function
    Real mu0 = pow(2.0, common::real(alpha + beta + 1.0))
        *std::tgamma(common::real(alpha + 1.0))
        *std::tgamma(common::real(beta + 1.0))
        /std::tgamma(common::real(alpha + beta + 2.0));

    // Extract quadrature points from the eigenvalues and weights from the
    // magnitude of the eigenvector
    for (int j = 0; j < n; j++) {
      points(j) = eigvals(j); 

      Real q_0j = eigvecs(j,0);  // Note: eigvecs are along rows
      weights(j) = q_0j*q_0j*mu0;
    }
  } else if (n == 1) {
    points(0) = 0;
    weights(0) = 2;
  } 
}
