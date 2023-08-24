#include "lagrange_polynomials.h"
#include "point_distributions.h"

#include <cstdlib>   // rand, RAND_MAX
#include <ctime>     // time
#include <iostream>

// There are 3 tests that cover the functionality of the Lagrange polynomial
// routines. These are described below

// Note: in the spirit of intellectual honesty, it should be disclaimed that
// the tolerances chosen in these tests were chosen so that the tests would
// pass when the results are "pretty damn close"; that is, the discrepancies
// between the results and the expectation or truth are presumed to come from
// finite precision errors. That may or may not be good enough. If in the
// course of using these routines elsewhere their numerical stabililty is
// called into question, these discrepancies and their sources may need to be
// investigated.

/* Test parameters */
template <typename NumType>
struct TestParams {
  SizeType Np;
  SizeType Nv;
  SizeType I;

  NumType Zl;
  NumType Zr;
  NumType X;

  NumType *c;
  NumType *Z;
  NumType *w;
  NumType *C;

  TestParams(SizeType order) {
    // Test parameters
    Np = order;
    Nv = Np + 1;  // number of vertices

    // Generate random array of coefficients between 0 and 1
    c = new NumType[Nv];
    for (SizeType i = 0; i < Nv; i++) {
      c[i] = Real(rand())/RAND_MAX;
    }

    // Generate a set of points between -1 and 1
    Zl = -1.0;
    Zr = 1.0;
    Z = new NumType[Nv];
    equispaced_points(Nv, Zl, Zr, Z);

    // Compute barycentric weights
    w = new NumType[Nv];
    lagrange::compute_barycentric_weights(Nv, Z, w);

    // Select a random coordinate between -1 and 1
    X = Zl + (Zr - Zl)*Real(rand())/Real(RAND_MAX);

    // Select a random vertex ID, between 0 and (Nv - 1)
    I = std::round((Nv - 1)*Real(rand())/Real(RAND_MAX));

    // Allocate workspace
    C = new NumType[Nv];
  };

  ~TestParams() {
    delete [] c, Z, w, C;
  }
};

/*
 * Test 1 
 * ------
 * This test checks the consistency of Lagrange polynomials (basis functions)
 * with the Lagrange interpolant. There are two ways to evaluate the sum of the
 * products of Lagrange polynomials and coefficients. One is to evaluate the
 * Lagrange polynomials individually, scale them by the corresponding
 * coefficients, and sum. The other is to use the interpolation routine. Either
 * way you do it, you should get the same answer; that is, they should be
 * consistent.
 *
 * Additionally, there are two cases of interest when using the Lagrange
 * polynomials in barycentric form, as we are here. The first case is when the
 * input coordinate is not coincident with any of the nodes. The second case is
 * when the input coordinate is coincident, which causes singularity to arise
 * in the barycentric form and must be handled differently.
 */
bool test1(TestParams<Real> &p) {
  Real tol = 1e-14;

  // Case 1: random coordinate, not coincident with any vertices
  bool case1 = false;

  Real sum1 = 0.0;

  for (SizeType i = 0; i < p.Nv; i++) {
    Real li = lagrange::eval(p.Nv, i, -1, p.Z, p.w, p.X);
    sum1 += p.c[i]*li;
  }

  Real sum2 = lagrange::eval_interp(p.Nv, -1, p.Z, p.w, p.X, p.c);

  Real rel_error = common::abs((sum1 - sum2)/sum1);

  if (rel_error < tol) {
    case1 = true;
  } else {
    std::cout << "Test 1, Case 1, error: " << rel_error << std::endl; 
  }

  // Case 2: random vertex, coincident case
  bool case2 = false;

  sum1 = 0.0;

  for (SizeType i = 0; i < p.Nv; i++) {
    Real li = lagrange::eval(p.Nv, i, p.I, p.Z, p.w, p.Z[p.I]);
    sum1 += p.c[i]*li;
  }

  sum2 = lagrange::eval_interp(p.Nv, p.I, p.Z, p.w, p.Z[p.I], p.c);

  rel_error = common::abs((sum1 - sum2)/sum1);

  if (rel_error < tol) {
    case2 = true;
  } else {
    std::cout << "Test 1, Case 2, error: " << rel_error << std::endl; 
  }

  return case1 && case2;
}


/*
 * Test 2
 * ------
 * This test checks the consistency of the (first) derivatives of the Lagrange
 * polynomials with the derivative of the Lagrange interpolant. As in Test 1,
 * there are two ways to obtain the derivative of the sum of the products of
 * the Lagrange polynomials and coefficients. The results from both should be
 * consistent.
 *
 * The test addresses the same two cases as in Test 1.
 */
bool test2(TestParams<Real> &p) {
  Real tol = 1e-10;

  // Case 1: random coordinate, not coincident with any vertices
  bool case1 = false;
  Real sum1 = 0.0;

  for (SizeType i = 0; i < p.Nv; i++) {
    Real dli = lagrange::eval_der(p.Nv, 1, i, -1, p.Z, p.w, p.X, p.C);
    sum1 += p.c[i]*dli;
  }

  std::copy(p.c, p.c+p.Nv, p.C);
  Real sum2 = lagrange::eval_der_interp(p.Nv, 1, -1, p.Z, p.w, p.X, p.C);

  Real rel_error = common::abs((sum1 - sum2)/sum1);

  if (rel_error < tol) {
    case1 = true;
  } else {
    std::cout << "Test 2, Case 1, error: " << rel_error << std::endl; 
  }

  // Case 2: random vertex, coincident case
  bool case2 = false;
  sum1 = 0.0;

  for (SizeType i = 0; i < p.Nv; i++) {
    Real dli = lagrange::eval_der(p.Nv, 1, i, p.I, p.Z, p.w, p.X, p.C);
    sum1 += p.c[i]*dli;
  }

  std::copy(p.c, p.c+p.Nv, p.C);
  sum2 = lagrange::eval_der_interp(p.Nv, 1, p.I, p.Z, p.w, p.X, p.C);

  rel_error = common::abs((sum1 - sum2)/sum1);

  if (rel_error < tol) {
    case2 = true;
  } else {
    std::cout << "Test 2, Case 2, error: " << rel_error << std::endl; 
  }

  return case1 && case2;
}

/*
 * Test 3
 * ------
 * This test checks the correctness of the interpolant derivative routine
 * against an approximation of the derivative based on the complex step method.
 *
 * The test addresses the same cases as in Test 1 and 2. In Case 1, we use the
 * complex step method. In Case 2, this is not possible since the interpolation
 * function, as implemented, is not algorithmically differentiable for
 * coincident coordinates. Therefore, we must verify the derivative by other
 * means. We use the formula for the exact first derivative from Berrut and
 * Trefethen, 2004
 */
bool test3(TestParams<Complex> &p) {
  Real tol = 1e-10;
  Real h = 1e-30;

  // Case 1: random coordinate, not coincident with any vertices
  bool case1 = false;

  std::copy(p.c, p.c+p.Nv, p.C);
  Real df = common::real(
      lagrange::eval_der_interp(p.Nv, 1, -1, p.Z, p.w, p.X, p.C));

  std::copy(p.c, p.c+p.Nv, p.C);
  Complex X = p.X + Complex(0, h);  // complex perturbation
  Complex f_cs = lagrange::eval_interp(p.Nv, -1, p.Z, p.w, X, p.C);
  Real df_cs = common::imag(f_cs)/h;

  Real rel_error = common::abs((df - df_cs)/df_cs);

  if (rel_error < tol) {
    case1 = true;
  } else {
    std::cout << "Test 3, Case 1, error: " << rel_error << std::endl; 
  }

  // Case 2: random vertex, coincident case
  bool case2 = false;

  std::copy(p.c, p.c+p.Nv, p.C);
  df = common::real(
      lagrange::eval_der_interp(p.Nv, 1, p.I, p.Z, p.w, p.X, p.C));

  Real df_exact = 0.0;
  for (SizeType i = 0; i < p.Nv; i++) {
    if (i == p.I) {
      Real dl = 0.0;
      for (SizeType j = 0; j < p.Nv; j++) {
        if (j != i) dl += common::real(-p.w[j]/p.w[i]/(p.Z[i] - p.Z[j]));
      }
      df_exact += common::real(p.c[i])*dl;
    } else {
      Real dl = common::real(p.w[i]/p.w[p.I]/(p.Z[p.I] - p.Z[i]));
      df_exact += common::real(p.c[i])*dl;
    }
  }

  rel_error = common::abs((df - df_exact)/df_exact);

  if (rel_error < tol) {
    case2 = true;
  } else {
    std::cout << "Test 3, Case 2, error: " << rel_error << std::endl; 
  }

  return case1 && case2;
}

int main() {
  std::cout.precision(15);

  std::cout << "TEST LAGRANGE POLYNOMIALS\n" 
            << "-------------------------" 
            << std::endl;

  // Initialize random seed
  srand(time(NULL));

  // Generate real-valued test parameters
  TestParams<Real> rp(8);

  // Run test 1
  bool pass1 = test1(rp);
  std::cout << "TEST 1 " << (pass1 ? "PASSED" : "FAILED") << std::endl;

  // Run test 2
  bool pass2 = test2(rp);
  std::cout << "TEST 2 " << (pass2 ? "PASSED" : "FAILED") << std::endl;

  // Generate complex-valued test parameters
  TestParams<Complex> ip(8);

  // Run test 3
  bool pass3 = test3(ip);
  std::cout << "TEST 3 " << (pass3 ? "PASSED" : "FAILED") << std::endl;

  std::cout << std::endl;
  std::cout << "PASSED " 
            << int(pass1) + int(pass2) + int(pass3) 
            << "/3" << std::endl;
  
  return 0;
}
