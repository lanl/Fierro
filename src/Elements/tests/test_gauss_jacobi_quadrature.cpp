#include "common.h"
#include "gauss_jacobi_quadrature.h"

#include <iostream>
#include <iomanip>

using namespace std;

/** Test Gauss-Jacobi quadrature implementation */
int main() {
  cout.precision(15);

  const size_t n = 5;
  Real alpha = -0.5, beta = -0.5;
  CArray<Real> points(n); 
  CArray<Real> weights(n);
  compute_gauss_jacobi_quadrature_rule(n, alpha, beta, points, weights);

  cout << setw(30) << "points"
       << setw(30) << "weights"
       << endl;
  cout << setw(30) << "------"
       << setw(30) << "-------"
       << endl;
  for (int i = 0; i < n; i++) {
    cout << setw(30) << common::real(points(i))
         << setw(30) << common::real(weights(i))
         << endl;
  }

  return 0;
}
