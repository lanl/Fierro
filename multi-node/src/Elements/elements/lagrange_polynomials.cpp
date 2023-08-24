#include "lagrange_polynomials.h"

namespace lagrange {
  /*
   * If the input coordinate is coincident with any of the input vertices,
   * return the index of the vertex with which it coincides
   *
   * It is typically bad practice to check if two floating point numbers, say x
   * and y, are equal by evaluating x - y == 0. But that is precisely what the
   * authors suggest to do in "Barycentric Lagrange Interpolation" (Berrut and
   * Trefethen, 2004). It is also what the authors of Nektar++ do in their
   * implementation of barycentric Lagrange interpolation. I'm not sure what to
   * do
   *
   * Parameters
   * ----------
   * N : number of vertices
   * z : vertex coordinates
   * x : coordinate to test
   */
  template <typename NumType>
  SizeType find_coincident_vertex(const SizeType &N, const NumType *z, 
      const NumType &x) {
    // The largest number an unsigned integer can be
    SizeType k = -1;

    // Adding to -1 in unsigned rolls over zero
    for (SizeType j = 0; j < N; j++) {
      //k += (j + 1)*SizeType(z[j] - x == 0.0);
      k += (j + 1)*SizeType(common::almost_equal(z[j], x));
    }

    return k;
  }

  /*
   * Calculation of barycentric weights of Lagrange interpolant vertices
   *
   * Parameters
   * ----------
   * N : number of vertices 
   * z : vertex coordinates
   *
   * Returns
   * -------
   * w : barycentric vertex weights
   *
   * Notes
   * -----
   * If one evaluates the interpolant using the barycentric formula of the
   * second kind, the barycentric weights may be scaled by an arbitrary scalar
   * (the maximum weight, for example). Since I am using the barycentric
   * formula of the first kind (also called the modified Lagrange formula), I
   * cannot and do not scale the weights 
   */
  template <typename NumType>
  void compute_barycentric_weights(const SizeType &N, const NumType *z, 
      NumType *w) {
    // Compute weights
    for (SizeType j = 0; j < N; j++) {
      w[j] = 1.0;
      for (SizeType k = 0; k < N; k++) {
        // TODO optimization: make branchless? Trefethen uses log o exp trick
        if (j != k) w[j] *= 1.0/(z[j] - z[k]);
      }
    }
  };

  /*
   * Return evaluation of Lagrange polynomial associated with a given vertex at
   * a given coordinate
   *
   * Parameters
   * ----------
   * Nv : number of vertices
   * i  : index of vertex
   * ic : index of vertex coincident with X (-1 if not coincident)
   * Z  : vertex coordinates
   * w  : barycentric weights
   * X  : coordinate at which to evaluate
   */
  template <typename NumType>
  NumType eval(const SizeType Nv, const SizeType i, const SizeType ic, 
      const NumType *Z, const NumType *w, const NumType X) {
    if (ic < Nv) return i == ic ? 1.0 : 0.0;

    // Evaluate nodal polynomial
    NumType L = 1.0;
    for (SizeType j = 0; j < Nv; j++) L *= (X - Z[j]);

    return L*w[i]/(X - Z[i]);
  }

  /* 
   * Return evaluation of derivative of Lagrange polynomial associated with a
   * given vertex at a given coordinate
   *
   * Currently implemented as a special case of interpolant derivative
   * evaluation. There is probably a better way to do this, but I haven't
   * gotten around to finding it and implementing it
   *
   * Parameters
   * ----------
   * Nv         : number of vertices
   * n          : order of derivative
   * i          : index of vertex
   * ic : index of vertex coincident with X (-1 if not coincident)
   * Z          : vertex coordinates
   * w          : barycentric weights
   * X          : coordinate at which to evaluate
   * c          : work array for intermediate coefficients
   */
  template <typename NumType>
  NumType eval_der(const SizeType Nv, const SizeType n, const SizeType i, 
      const SizeType ic, const NumType *Z, const NumType *w, 
      const NumType X, NumType *c) {
    // Initialize coefficients of interpolant to Kronecker deltas to create
    // interpolant equivalent to Lagrange polynomial at vertex
    for (SizeType j = 0; j < Nv; j++) {
      if (j != i) { c[j] = 0.0; } else { c[j] = 1.0; }
    }
    
    NumType dnl = eval_der_interp<NumType>(Nv, n, ic, Z, w, X, c);

    return dnl;
  }

  /* 
   * Return evaluation of Lagrange interpolant, which is the sum of the
   * products of Lagrange polynomials and provided coefficients, at specified
   * coordinates
   *
   * This implementation uses the modified Lagrange formula (barycentric
   * formula of the first kind) from "Barycentric Lagrange Interpolation"
   * (Berrut and Trefethen, 2004). 
   *
   * Parameters
   * ----------
   * Nv : number of vertices 
   * i  : index of vertex coincident with coordinate (-1 if not coincident)
   * Z  : vertex coordinates
   * w  : barycentric vertex weights
   * X  : coordinate at which to evaluate interpolant
   * c  : input coefficients (function values at vertices)
   */
  template <typename NumType>
  NumType eval_interp(const SizeType Nv, const SizeType i, const NumType *Z, 
      const NumType *w, const NumType X, const NumType *c) {
    NumType p;
    NumType L = 1.0;  // nodal polynomial evaluation

    // Evaluate the interpolant
    if (i < Nv) {  // coincident
      return c[i];
    } else {  // non-coincident
      // Loop over vertices
      p = 0.0;
      for (SizeType j = 0; j < Nv; j++) {
        // Contribution to interpolant evaluation
        p += w[j]*c[j]/(X - Z[j]);

        // Contribution to scalar (nodal polynomial evaluation)
        L *= (X - Z[j]);
      }

      // Apply scaling
      p *= L;
    }

    return p; 
  }

  /* 
   * Return evaluation of derivative of Lagrange interpolant, which is the sum
   * of the products of the derivatives of Lagrange polynomials and provided
   * coefficients, at specified coordinates
   *
   * This implementation is based on the formula for the derivative of a
   * rational interpolant given in "Some New Aspects of Rational Interpolation"
   * (Schneider and Werner, 1986).
   *
   * Parameters
   * ----------
   * Nv : number of vertices 
   * n  : order of derivative
   * i  : index of vertex coincident with coordinate (-1 if not coincident)
   * Z  : vertex coordinates
   * w  : barycentric vertex weights
   * X  : coordinate at which to eval interpolant
   * c  : input coefficients (will be modified)
   */
  template <typename NumType>
  NumType eval_der_interp(const SizeType Nv, const SizeType n, const SizeType i,
      const NumType *Z, const NumType *w, const NumType X, NumType *c) {
    // Interpolant and nodal polynomial evaluation
    NumType p = 0.0;  // interpolant evaluation
    NumType L = 1.0;  // nodal polynomial
    if (i < Nv) {  // coincident
      p = c[i];
    } else {  // non-coincident
      // Loop over vertices
      p = 0.0;
      for (SizeType j = 0; j < Nv; j++) {
        // Contribution to interpolant evaluation
        p += w[j]*c[j]/(X - Z[j]);

        // Contribution to scalar (nodal polynomial evaluation)
        L *= (X - Z[j]);
      }

      // Apply scaling
      p *= L;
    }

    // Evaluate derivatives, building up to specified order
    NumType dkp = p, dnp = 0.0;
    NumType M = 1.0;  // factorial(k)
    if (i < Nv) {  // coincident
      for (SizeType k = 1; k <= n; k++) {
        // Zero the output
        dnp = 0.0;

        for (SizeType j = 0; j < Nv; j++) {
          // Calculate divided difference and store in intermediate coefficients
          NumType sx = 1.0/(Z[i] - Z[j]);
          c[j] = sx*(dkp - c[j]);

          // Update output coefficient with contribution from vertex
          // TODO optimization: make branchless? Trefethen uses log o exp trick
          if (j != i) dnp += w[j]*c[j];
        }

        // Scale the output and copy for use in calculating next order
        M *= k;
        dnp *= -M/w[i];
        dkp = dnp;
      }
    } else {  // non-coincident
      for (SizeType k = 1; k <= n; k++) {
        // Zero the output
        dnp = 0.0;

        for (SizeType j = 0; j < Nv; j++) {
          // Calculate divided difference and store in intermediate coefficients
          NumType sx = 1.0/(X - Z[j]);
          c[j] = sx*(dkp - c[j]);

          // Update output coefficient with contribution from vertex
          dnp += w[j]*c[j]/(X - Z[j]);
        }

        // Scale the output and copy for use in calculating next order
        M *= k;
        dnp *= L*M;
        dkp = dnp;
      }
    }

    return dnp;
  }

  // Explicit instantiations of template functions
  template SizeType find_coincident_vertex(const SizeType &N, const Real *z, 
      const Real &x);
  template SizeType find_coincident_vertex(const SizeType &N, const Complex *z, 
      const Complex &x);

  template void compute_barycentric_weights(const SizeType &N, const Real *z, 
      Real *w);
  template void compute_barycentric_weights(const SizeType &N, const Complex *z, 
      Complex *w);

  template Real eval(const SizeType Nv, const SizeType i, const SizeType ic, 
      const Real *Z, const Real *w, const Real X);
  template Complex eval(const SizeType Nv, const SizeType i, const SizeType ic, 
      const Complex *Z, const Complex *w, const Complex X);

  template Real eval_der(const SizeType Nv, const SizeType n, const SizeType i, 
      const SizeType ic, const Real *Z, const Real *w, const Real X, Real *c);
  template Complex eval_der(const SizeType Nv, const SizeType n, 
      const SizeType i, const SizeType ic, const Complex *Z, const Complex *w, 
      const Complex X, Complex *c);

  template Real eval_interp(const SizeType Nv, const SizeType i, const Real *Z, 
      const Real *w, const Real X, const Real *c);
  template Complex eval_interp(const SizeType Nv, const SizeType i, 
      const Complex *Z, const Complex *w, const Complex X, const Complex *c);

  template Real eval_der_interp(const SizeType Nv, const SizeType n, 
      const SizeType i, const Real *Z, const Real *w, const Real X, 
      Real *c);
  template Complex eval_der_interp(const SizeType Nv, const SizeType n, 
      const SizeType i, const Complex *Z, const Complex *w, const Complex X, 
      Complex *c);
}
