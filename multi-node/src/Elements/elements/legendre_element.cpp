#include "legendre_element.h"
#include "legendre_polynomials.h"

template <typename NumType>
LegendreElement<NumType>::LegendreElement(const SizeType order) : Np(order) {
  N = Np + 1;
  Ne = std::pow(N, Nd);

  rad[0] = N;
  rad[1] = N;
  rad[2] = N;

  // Allocate memory for intermediate coefficients
  C = new NumType[2*N];
}

template <typename NumType>
LegendreElement<NumType>::~LegendreElement() { }

/*
 * Return evaluation of basis function, which is the tensor product of Legendre
 * polynomials, at specified coordinates.
 *
 * Parameters
 * ----------
 * I : index of basis function
 * X : coordinates (in reference space)
 */
template <typename NumType>
NumType LegendreElement<NumType>::eval_basis(const SizeType I, const NumType *X) {
  // Decompose index of 3D tensor product basis function into indices of
  // Legendre polynomials
  common::base_10_to_mixed_radix(Nd, rad, I, ijk);

  // Evaluate Legendre polynomials
  NumType Pi = legendre::eval(int(ijk[0]), X[0]);
  NumType Pj = legendre::eval(int(ijk[1]), X[1]);
  NumType Pk = legendre::eval(int(ijk[2]), X[2]);

  return Pi*Pj*Pk;
}

/*
 * Evaluate gradient of basis function
 *
 * Parameters
 * ----------
 * I : index of basis function
 * X : coordinates (in reference space)
 *
 * Returns
 * -------
 * grad_phi : gradient of the basis function
 */
template <typename NumType>
void LegendreElement<NumType>::eval_grad_basis(const SizeType I, 
    const  NumType *X, NumType *grad_phi) {
  // Decompose index of 3D tensor product basis function into indices of
  // Legendre polynomials
  common::base_10_to_mixed_radix(Nd, rad, I, ijk);

  // Evaluate Legendre polynomials
  NumType Pi = legendre::eval(int(ijk[0]), X[0]);
  NumType Pj = legendre::eval(int(ijk[1]), X[1]);
  NumType Pk = legendre::eval(int(ijk[2]), X[2]);
  
  // Evaluate Legendre polynomial derivatives
  SizeType q = 1;  // order of derivative
  NumType dPi = legendre::eval_der(int(ijk[0]), q, X[0]);
  NumType dPj = legendre::eval_der(int(ijk[1]), q, X[1]);
  NumType dPk = legendre::eval_der(int(ijk[2]), q, X[2]);

  // Store partial derivatives in entries of gradient
  grad_phi[0] = dPi*Pj*Pk;
  grad_phi[1] = Pi*dPj*Pk;
  grad_phi[2] = Pi*Pj*dPk;
}

/*
 * Return evaluation of local function approximation, which is formed by the
 * sum of the products of tensor-product Legendre basis functions and the
 * provided coefficients, at specified coordinates
 *
 * Parameters
 * ----------
 * c : coefficients
 * X : coordinates (in reference space)
 *
 * Note
 * ----
 * The intermediate coefficients or workspace array is divided into 2 lines of
 * length equal to the number of vertices. In the dimension-by-dimension
 * approach used for the approximation evaluation, each line is filled by
 * evaluations at the previous dimensional level. The initial coefficients are
 * cycled through in the X-dimension loop to produce coefficients, stored in
 * the workspace array, to be used in the Y-dimension loop. The output from the
 * Y-dimension loop is coefficients, stored again in the workspace array, to be
 * used in the final Z-dimension evaluation. 
 */
template <typename NumType>
NumType LegendreElement<NumType>::eval_approx(const NumType *c, 
    const NumType *X) {
  for (int k = 0; k < N; k++) {
    for (int j = 0; j < N; j++) {
      // Collapse first dimension into coefficients for second dimension
      C[j] = legendre::eval_approx(N, &c[j*N+k*N*N], X[0]);
    }

    // Collapse second dimension into coefficients for third dimension
    C[N+k] = legendre::eval_approx(N, C, X[1]);
  }

  // Collapse third dimension into approximation evaluation
  return legendre::eval_approx(N, &C[N], X[2]);
}

/*
 * Evaluate gradient of the local function approximation, which is formed by
 * the sum of the products of tensor-product Legendre basis functions and the
 * provided coefficients, at specified coordinates
 *
 * Parameters
 * ----------
 * c : coefficients
 * X : coordinates (in reference space)
 *
 * Returns
 * -------
 * grad_f : gradient of the approximation
 */
template <typename NumType>
void LegendreElement<NumType>::eval_grad_approx(const NumType *c, 
    const NumType *X, NumType *grad_f) {
  SizeType q = 1;  // order of derivative
  for (int l = 0; l < Nd; l++) {
    for (int k = 0; k < N; k++) {
      for (int j = 0; j < N; j++) {
        // Collapse first dimension into coefficients for second dimension
        if (l == 0) {
          C[j] = legendre::eval_der_approx(N, 1, &c[j*N+k*N*N], X[0]);
        } else {
          C[j] = legendre::eval_approx(N, &c[j*N+k*N*N], X[0]);
        }
      }

      // Collapse second dimension into coefficients for third dimension
      if (l == 1) {
        C[N+k] = legendre::eval_der_approx(N, q, C, X[1]);
      } else {
        C[N+k] = legendre::eval_approx(N, C, X[1]);
      }
    }

    // Collapse third dimension into approximation evaluation
    if (l == 2) {
      grad_f[l] = legendre::eval_der_approx(N, q, &C[N], X[2]);
    } else {
      grad_f[l] = legendre::eval_approx(N, &C[N], X[2]);
    }
  }
}

// Explicit instantiation of template class
template class LegendreElement<Real>;
template class LegendreElement<Complex>;
