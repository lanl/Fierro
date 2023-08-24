#include "lagrange_element.h"
#include "lagrange_polynomials.h"

template <typename NumType>
LagrangeElement<NumType>::LagrangeElement(const SizeType order, 
    const NumType *vert_coords) : Np(order), Z(vert_coords) {
  N = Np + 1;
  Ne = std::pow(N, Nd);

  rad[0] = N;
  rad[1] = N;
  rad[2] = N;

  // Allocate memory for weights and compute
  w = new NumType[N];
  lagrange::compute_barycentric_weights(N, Z, w);

  // Allocate memory for intermediate coefficients
  C = new NumType[3*N];
}

template <typename NumType>
LagrangeElement<NumType>::~LagrangeElement() {
  delete [] w, C;
}

/*
 * Return evaluation of basis function, which is the tensor product of Lagrange
 * polynomials, at specified coordinates.
 *
 * Parameters
 * ----------
 * I : index of basis function
 * X : coordinates (in reference space)
 */
template <typename NumType>
NumType LagrangeElement<NumType>::eval_basis(const SizeType I, const NumType *X) {
  // Decompose index of 3D tensor product basis function into indices of
  // Lagrange polynomials
  common::base_10_to_mixed_radix(Nd, rad, I, ijk);

  // Check coincidence of coordinates with vertex coordinates
  SizeType ix = lagrange::find_coincident_vertex(N, Z, X[0]);
  SizeType iy = lagrange::find_coincident_vertex(N, Z, X[1]);
  SizeType iz = lagrange::find_coincident_vertex(N, Z, X[2]);

  // Evaluate Lagrange polynomials
  NumType li = lagrange::eval(N, ijk[0], ix, Z, w, X[0]);
  NumType lj = lagrange::eval(N, ijk[1], iy, Z, w, X[1]);
  NumType lk = lagrange::eval(N, ijk[2], iz, Z, w, X[2]);

  return li*lj*lk;
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
void LagrangeElement<NumType>::eval_grad_basis(const SizeType I, 
    const  NumType *X, NumType *grad_phi) {
  // Decompose index of 3D tensor product basis function into indices of
  // Lagrange polynomials
  common::base_10_to_mixed_radix(Nd, rad, I, ijk);

  // Check coincidence of coordinates with vertex coordinates
  SizeType ix = lagrange::find_coincident_vertex(N, Z, X[0]);
  SizeType iy = lagrange::find_coincident_vertex(N, Z, X[1]);
  SizeType iz = lagrange::find_coincident_vertex(N, Z, X[2]);

  // Evaluate Lagrange polynomials
  NumType li = lagrange::eval(N, ijk[0], ix, Z, w, X[0]);
  NumType lj = lagrange::eval(N, ijk[1], iy, Z, w, X[1]);
  NumType lk = lagrange::eval(N, ijk[2], iz, Z, w, X[2]);
  
  // Evaluate Lagrange polynomial derivatives
  SizeType n = 1;  // order of derivative
  NumType dli = lagrange::eval_der(N, n, ijk[0], ix, Z, w, X[0], C);
  NumType dlj = lagrange::eval_der(N, n, ijk[1], iy, Z, w, X[1], C);
  NumType dlk = lagrange::eval_der(N, n, ijk[2], iz, Z, w, X[2], C);

  // Store partial derivatives in entries of gradient
  grad_phi[0] = dli*lj*lk;
  grad_phi[1] = li*dlj*lk;
  grad_phi[2] = li*lj*dlk;
}

/*
 * Return evaluation of local function approximation, which is formed by
 * the sum of the products of tensor-product Lagrange basis functions and the
 * provided coefficients, at specified coordinates
 *
 * Parameters
 * ----------
 * c : coefficients
 * X : coordinates (in reference space)
 *
 * Note
 * ----
 * The intermediate coefficients or workspace array is divided into 3 lines of
 * length equal to the number of vertices. In the dimension-by-dimension
 * approach used for the approximation evaluation, each line is filled by
 * evaluations at the previous dimensional level. The initial coefficients are
 * cycled through in the X-dimension loop to produce coefficients, stored in
 * the workspace array, to be used in the Y-dimension loop. The output from the
 * Y-dimension loop is coefficients, stored again in the workspace array, to be
 * used in the final Z-dimension evaluation. This process, as described, only
 * requires 2 lines of workspace. The extra line of workspace is needed for
 * derivative evaluations, and so it is unused here.
 */
template <typename NumType>
NumType LagrangeElement<NumType>::eval_approx(const NumType *c, 
    const NumType *X) {
  // Check the coincidence of the coordinates with vertex coordinates
  SizeType ix = lagrange::find_coincident_vertex(N, Z, X[0]);
  SizeType iy = lagrange::find_coincident_vertex(N, Z, X[1]);
  SizeType iz = lagrange::find_coincident_vertex(N, Z, X[2]);

  for (int k = 0; k < N; k++) {
    for (int j = 0; j < N; j++) {
      // Collapse first dimension into coefficients for second dimension
      C[j] = lagrange::eval_interp(N, ix, Z, w, X[0], &c[j*N+k*N*N]);
    }

    // Collapse second dimension into coefficients for third dimension
    C[N+k] = lagrange::eval_interp(N, iy, Z, w, X[1], C);
  }

  // Collapse third dimension into interpolant evaluation
  return lagrange::eval_interp(N, iz, Z, w, X[2], &C[N]);
}

/*
 * Evaluate gradient of the local function approximation, which is formed by
 * the sum of the products of tensor-product Lagrange basis functions and the
 * provided coefficients, at specified coordinates
 *
 * Parameters
 * ----------
 * c : coefficients
 * X : coordinates (in reference space)
 *
 * Returns
 * -------
 * grad_f : gradient of the interpolant
 *
 * Note
 * ----
 * This routine uses a dimension-by-dimension approach as in the approximation
 * routine. See the documentation of that routine for details. This routine
 * differs from that one in that the first line of the workspace array is
 * reserved for the derivative evaluations. These derivative evaluations use
 * the first line to store divided differences.
 */
template <typename NumType>
void LagrangeElement<NumType>::eval_grad_approx(const NumType *c, 
    const NumType *X, NumType *grad_f) {
  const SizeType n = 1; // order of partial derivative

  // Check the coincidence of the coordinates with vertex coordinates
  SizeType ix = lagrange::find_coincident_vertex(N, Z, X[0]);
  SizeType iy = lagrange::find_coincident_vertex(N, Z, X[1]);
  SizeType iz = lagrange::find_coincident_vertex(N, Z, X[2]);

  for (int l = 0; l < Nd; l++) {
    for (int k = 0; k < N; k++) {
      for (int j = 0; j < N; j++) {
        // Collapse first dimension into coefficients for second dimension
        std::copy(c+j*N+k*N*N, c+j*N+k*N*N+N, C);  // load first line
        if (l == 0) {
          C[N+j] = lagrange::eval_der_interp(N, n, ix, Z, w, X[0], C);
        } else {
          C[N+j] = lagrange::eval_interp(N, ix, Z, w, X[0], C);
        }
      }

      // Collapse second dimension into coefficients for third dimension
      std::copy(C+N, C+2*N, C); // second line back into first line
      if (l == 1) {
        C[2*N+k] = lagrange::eval_der_interp(N, n, iy, Z, w, X[1], C);
      } else {
        C[2*N+k] = lagrange::eval_interp(N, iy, Z, w, X[1], C);
      }
    }

    // Collapse third dimension into interpolant evaluation
    std::copy(C+2*N, C+3*N, C); // third line back into first line
    if (l == 2) {
      grad_f[l] = lagrange::eval_der_interp(N, n, iz, Z, w, X[2], C);
    } else {
      grad_f[l] = lagrange::eval_interp(N, iz, Z, w, X[2], C);
    }
  }
}

/*
 * Evaluate the Jacobian of the spatial mapping
 *
 * Parameters
 * ----------
 * x : x coordinates in physical space
 * y : y coordinates in physical space
 * z : z coordinates in physical space
 * X : coordinates in reference space at which to evaluate
 * 
 * Returns
 * -------
 * J : Jacobian (stored in column-major order)
 */
template <typename NumType>
void LagrangeElement<NumType>::eval_jac(const NumType *x, const NumType *y, 
    const NumType *z, const NumType *X, NumType *J) {
  // Evaluate the gradient of x = x(X, Y, Z)
  this->eval_grad_approx(x, X, J);

  // Evaluate the gradient of y = y(X, Y, Z)
  this->eval_grad_approx(y, X, J+3);

  // Evaluate the gradient of z = z(X, Y, Z)
  this->eval_grad_approx(z, X, J+6);
}

/*
 * Return the determinant of the Jacobian of the spatial mapping
 *
 * Parameters
 * ----------
 * x : x coordinates in physical space
 * y : y coordinates in physical space
 * z : z coordinates in physical space
 * X : coordinates in reference space at which to evaluate
 */
template <typename NumType>
NumType LagrangeElement<NumType>::eval_det_jac(const NumType *x, 
    const NumType *y, const NumType *z, const NumType *X) {
  NumType J[9];
  this->eval_jac(x, y, z, X, J);

  return J[0]*(J[4]*J[8] - J[5]*J[7]) - J[1]*(J[3]*J[8] - J[5]*J[6]) 
      + J[2]*(J[3]*J[7] - J[4]*J[6]);
}

/* Evaluate the inverse of the Jacobian of the spatial mapping
 *
 * Parameters
 * ----------
 * x : x coordinates in physical space
 * y : y coordinates in physical space
 * z : z coordinates in physical space
 * X : coordinates in reference space at which to evaluate
 * 
 * Returns
 * -------
 * Jinv : Jacobian inverse (stored in column-major order)
 */
template <typename NumType>
void LagrangeElement<NumType>::eval_inv_jac(const NumType *x, const NumType *y, 
    const NumType *z, const NumType *X, NumType *Jinv) {
  NumType J[9];
  this->eval_jac(x, y, z, X, J);

  NumType d = 1.0/(this->eval_det_jac(x, y, z, X));

  Jinv[0] = d*(J[4]*J[8] - J[5]*J[7]);
  Jinv[1] = -d*(J[1]*J[8] - J[2]*J[7]);
  Jinv[2] = d*(J[1]*J[5] - J[2]*J[4]);

  Jinv[3] = -d*(J[3]*J[8] - J[5]*J[6]);
  Jinv[4] = d*(J[0]*J[8] - J[2]*J[6]);
  Jinv[5] = -d*(J[0]*J[5] - J[2]*J[3]);

  Jinv[6] = d*(J[3]*J[7] - J[4]*J[6]);
  Jinv[7] = -d*(J[0]*J[7] - J[1]*J[6]);
  Jinv[8] = d*(J[0]*J[4] - J[1]*J[3]);
}

// Explicit instantiation of template class
template class LagrangeElement<Real>;
template class LagrangeElement<Complex>;
