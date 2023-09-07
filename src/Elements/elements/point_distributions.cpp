#include "point_distributions.h"

/*
 * Fill an array with points equally spaced between end points (inclusive)
 *
 * Parameters
 * ----------
 * N  : number of points
 * Zl : left end point 
 * Zr : right end point 
 *
 * Returns
 * -------
 * Z : array of equispaced points
 */
template <typename NumType>
void equispaced_points(SizeType N, NumType &Zl, NumType &Zr, NumType *Z) {
  for (SizeType i = 0; i < N; i++)
    Z[i] = Zl + double(i)/double(N - 1)*(Zr - Zl);
}

/*
 * Fill an array with Chebyshev points of the second kind in the interval
 * defined by the specified end points. The symmetry-preserving technique from
 * Chebfun (chebtech1/chebpts.m) is used.
 *
 * Parameters
 * ----------
 * N  : number of points
 * Zl : left end point 
 * Zr : right end point 
 *
 * Returns
 * -------
 * Z : array of Chebyshev points
 */
template <typename NumType>
void chebyshev_points(SizeType N, NumType &Zl, NumType &Zr, NumType *Z) {
  // Evaluate the points using sine function to preserve symmetry
  NumType f = 0.5*M_PI/double(N);
  for (SizeType i = 0; i < N; i++) {
    int j = (-1*int(N) + 1) + 2*int(i);
    Z[i] = sin(f*double(j));
  }

  // Scale the points to fit the domain
  for (int i = 0; i < N; i++)
    Z[i] = 0.5*(1.0 - Z[i])*Zl + 0.5*(1.0 + Z[i])*Zr;
}

// Explicit instantations of template functions
template void equispaced_points(SizeType N, Real &Zl, Real &Zr, Real *Z);
template void equispaced_points(SizeType N, Complex &Zl, Complex &Zr, Complex *Z);

template void chebyshev_points(SizeType N, Real &Zl, Real &Zr, Real *Z);
template void chebyshev_points(SizeType N, Complex &Zl, Complex &Zr, Complex *Z);
