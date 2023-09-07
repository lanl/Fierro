#include "lagrange_polynomials.h"

#include <iostream>
#include <iomanip>

// Polynomial function in three variables
NumType x(NumType X, NumType Y, NumType Z) {
  return 5.0 + 2.0*X*Y*Y + 3.0*X*Z + std::pow(X, 9);
}

// Polynomial function in three variables
NumType y(NumType X, NumType Y, NumType Z) {
  return 2.0 + 7.0*X*X*X + 5.0*Y + 9.0*Y*Z*Z*Z;
}

// Polynomial function in three variables
NumType z(NumType X, NumType Y, NumType Z) {
  return 2.0*X*Y + 2.0*X*Z;
}

// Determinant of the Jacobian based on three polynomial functions above
NumType J(NumType X, NumType Y, NumType Z) {
  NumType F[9];

  F[0] = 2.0*Y*Y + 3.0*Z + 9.0*std::pow(X, 8);
  F[1] = 21.0*X*X;
  F[2] = 2.0*Y + 2.0*Z;

  F[3] = 4.0*X*Y;
  F[4] = 5.0 + 9.0*Z*Z*Z;
  F[5] = 2.0*X;

  F[6] = 3.0*X;
  F[7] = 27.0*Y*Z*Z;
  F[8] = 2.0*X;

  return F[0]*(F[4]*F[8] - F[5]*F[7]) 
      - F[3]*(F[1]*F[8] - F[2]*F[7]) 
      + F[6]*(F[1]*F[5] - F[2]*F[4]);
}

/* Testing evaluation of Jacobian determinant using a convergence study
 *
 * Three polynomial functions (x, y, z) in three variables (X, Y, Z) are
 * considered together as a vector-valued function and their exact Jacobian
 * determinant is evaluated.
 *
 * Then successive orders of Lagrange tensor-product interpolants are created
 * and the Jacobian determinant approximated from these. Upon reaching the
 * highest order present in functions x, y, z, the approximation should match
 * the exact value.
 */
int main() {
  // Choose an evaluation point in the cube
  NumType X[3]; 
  X[0] = 0.5; 
  X[1] = 0.5; 
  X[2] = 0.5;

  // Evaluate the exact determinant of the Jacobian
  NumType Je = J(X[0], X[1], X[2]);

  // For a range of orders of interpolation, approximate the determinant
  std::cout << std::setw(20) << "p"
            << std::setw(20) << "J, exact"
            << std::setw(20) << "J, approx"
            << std::setw(20) << "error"
            << std::endl;
  std::cout << std::setw(20) << "---------"
            << std::setw(20) << "---------"
            << std::setw(20) << "---------"
            << std::setw(20) << "---------"
            << std::endl;
  for (int p = 1; p < 11; p++) {
    // Create a set of nodes in 1D
    SizeType N = p + 1;
    NumType Zl = -1.0, Zr = 1.0, Z[N];
    lagrange::equispaced_points(N, Zl, Zr, Z);
    
    // Compute barycentric weights
    NumType w[N];
    lagrange::compute_barycentric_weights(N, Z, w);

    // Determine the coefficients of the 3D interpolant
    NumType cx[N*N*N];
    NumType cy[N*N*N];
    NumType cz[N*N*N];
    for (int k = 0; k < N; k++) {
      for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
          cx[i + j*N + k*N*N] = x(Z[i], Z[j], Z[k]);
          cy[i + j*N + k*N*N] = y(Z[i], Z[j], Z[k]);
          cz[i + j*N + k*N*N] = z(Z[i], Z[j], Z[k]);
        }
      }
    }

    // Approximate the determinant of the Jacobian
    NumType co[3*N+1], Jp;
    lagrange::compute_jacobian_determinant(N, Z, w, cx, cy, cz, co, X, Jp);

    // Compute error
    NumType E = common::abs(Jp - Je);

    // Print results
    std::cout << std::setw(20) << p
              << std::setw(20) << Je 
              << std::setw(20) << Jp
              << std::setw(20) << E
              << std::endl;
  }

  return 0;
}
