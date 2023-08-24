#include "legendre_polynomials.h"
#include "point_distributions.h"

#include <cstdlib>   // rand, RAND_MAX
#include <ctime>     // time
#include <iostream>

// There are 3 tests that cover the functionality of the Legendre polynomial
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
  SizeType N;

  NumType X;

  NumType *c;

  TestParams(SizeType order) {
    // Test parameters
    Np = order;
    N = Np + 1;  // number of polynomial terms

    // Generate random array of coefficients between 0 and 1
    c = new NumType[N];
    for (SizeType i = 0; i < N; i++) {
      c[i] = Real(rand())/RAND_MAX;
    }


    // Select a random coordinate between -1 and 1
    NumType Zl = -1.0;
    NumType Zr = 1.0;
    X = Zl + (Zr - Zl)*Real(rand())/Real(RAND_MAX);
  };

  ~TestParams() { delete[] c; }
};

/*
 * Test 1 
 * ------
 * This test checks the consistency of Legendre polynomials (basis functions)
 * with the Legendre approximation. There are two ways to evaluate the sum of
 * the products of Legendre polynomials and coefficients. One is to evaluate
 * the Legendre polynomials individually, scale them by the corresponding
 * coefficients, and sum. The other is to use the approximation routine. Either
 * way you do it, you should get the same answer; that is, they should be
 * consistent.
 */
bool test1(TestParams<Real> &p) {
  Real tol = 1e-14;

  bool pass = false;

  Real sum1 = 0.0;
  for (SizeType i = 0; i < p.N; i++) {
    Real Pi = legendre::eval(i, p.X);
    sum1 += p.c[i]*Pi;
  }

  Real sum2 = legendre::eval_approx(p.N, p.c, p.X);

  Real rel_error = common::abs((sum1 - sum2)/sum1);

  if (rel_error < tol) {
    pass = true;
  } else {
    std::cout << "Test 1, error: " << rel_error << std::endl; 
  }

  return pass;
}


/*
 * Test 2
 * ------
 * This test checks the consistency of the (first) derivatives of the Legendre
 * polynomials with the derivative of the Legendre approximation. As in Test 1,
 * there are two ways to obtain the derivative of the sum of the products of
 * the Legendre polynomials and coefficients. The results from both should be
 * consistent.
 */
bool test2(TestParams<Real> &p) {
  Real tol = 1e-10;

  bool pass = false;
  Real sum1 = 0.0;

  for (SizeType i = 0; i < p.N; i++) {
    Real dPi = legendre::eval_der(i, 1, p.X);
    sum1 += p.c[i]*dPi;
  }

  Real sum2 = legendre::eval_der_approx(p.N, 1, p.c, p.X);

  Real rel_error = common::abs((sum1 - sum2)/sum1);

  if (rel_error < tol) {
    pass = true;
  } else {
    std::cout << "Test 2, error: " << rel_error << std::endl; 
  }

  return pass;
}

/*
 * Test 3
 * ------
 * This test checks the correctness of the approximation derivative routine
 * against an approximation of the derivative based on the complex step method.
 */
bool test3(TestParams<Complex> &p) {
  Real tol = 1e-10;
  Real h = 1e-30;

  bool pass = false;

  Real df = common::real(
      legendre::eval_der_approx(p.N, 1, p.c, p.X));

  Complex X = p.X + Complex(0, h);  // complex perturbation
  Complex f_cs = legendre::eval_approx(p.N, p.c, X);
  Real df_cs = common::imag(f_cs)/h;

  Real rel_error = common::abs((df - df_cs)/df_cs);

  if (rel_error < tol) {
    pass = true;
  } else {
    std::cout << "Test 3, error: " << rel_error << std::endl; 
  }

  return pass;
}

int main() {
  std::cout.precision(15);

  std::cout << "TEST LEGENDRE POLYNOMIALS\n" 
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
