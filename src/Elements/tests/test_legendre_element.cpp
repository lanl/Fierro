#include "legendre_element.h"
#include "point_distributions.h"

#include <cstdlib>   // rand, RAND_MAX
#include <ctime>     // time
#include <iostream>

// There are 3 tests that cover the functionality of the Legendre element
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
  LegendreElement<NumType> *elem;

  NumType *c;

  NumType X[3];

  TestParams(SizeType order) {
    // Test parameters
    SizeType Np = order;
    SizeType N = order + 1;

    // Create Legendre element
    elem = new LegendreElement<NumType>(Np);

    // Generate random array of coefficients between 0 and 1
    c = new NumType[elem->Ne];
    for (SizeType i = 0; i < elem->Ne; i++) {
      //c[i] = Real(rand())/RAND_MAX;
      c[i] = i == 182 ? 1.0 : 0.0;
    }

    // Select coordinates between -1 and 1
    NumType Zl = -1.0;
    NumType Zr = 1.0;
    X[0] = Zl + (Zr - Zl)*(Real(rand())/RAND_MAX);
    X[1] = Zl + (Zr - Zl)*(Real(rand())/RAND_MAX);
    X[2] = Zl + (Zr - Zl)*(Real(rand())/RAND_MAX);
  };

  ~TestParams() { delete elem; }
};

/*
 * Test 1 
 * ------
 * This test checks the consistency of Legendre tensor-product basis functions
 * with the Legendre tensor-product approximation. There are two ways to evaluate
 * the sum of the products of Legendre tensor-product basis functions and
 * coefficients. One is to evaluate the basis functions individually, scale
 * them by the corresponding coefficients, and sum. The other is to use the
 * interpolation routine. Either way you do it, you should get the same answer;
 * that is, they should be consistent.
 */
bool test1(TestParams<Real> &p) {
  Real tol = 1e-14;

  bool pass = false;
  Real sum1 = 0.0;

  for (SizeType i = 0; i < p.elem->Ne; i++) {
    Real phi_i = p.elem->eval_basis(i, p.X);
    sum1 += p.c[i]*phi_i;
  }

  Real sum2 = p.elem->eval_approx(p.c, p.X);

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

  Real d0f = 0.0;
  Real d1f = 0.0;
  Real d2f = 0.0;

  for (SizeType i = 0; i < p.elem->Ne; i++) {
    Real grad_phi[3];
    p.elem->eval_grad_basis(i, p.X, grad_phi);
    d0f += p.c[i]*grad_phi[0];
    d1f += p.c[i]*grad_phi[1];
    d2f += p.c[i]*grad_phi[2];
  }

  Real grad_f[3];
  p.elem->eval_grad_approx(p.c, p.X, grad_f);

  Real rel_error0 = common::abs((d0f - grad_f[0])/d0f);
  Real rel_error1 = common::abs((d1f - grad_f[1])/d1f);
  Real rel_error2 = common::abs((d2f - grad_f[2])/d2f);

  if (rel_error0 < tol && rel_error1 < tol && rel_error2 < tol) {
    pass = true;
  } else {
    std::cout << "Test 2, errors: " 
              << rel_error0 << " " 
              << rel_error1 << " " 
              << rel_error2 << " " 
              << std::endl; 
  }

  return pass;
}

/*
 * Test 3
 * ------
 * This test checks the correctness of the basis gradient routine against an
 * approximation of the derivative based on the complex step method.
 */
bool test3(TestParams<Complex> &p) {
  // If the tensor-product basis is constant or linear, there is a risk that
  // the partials will be zero, which will cause this particular test to fail.
  // To ensure that this doesn't happen, let's choose a higher-order  basis
  // function to be tested 
  SizeType rad[3] = {p.elem->N, p.elem->N, p.elem->N};
  SizeType ijk[3] = {2, 2, 2};
  SizeType I = common::mixed_radix_to_base_10(3, rad, ijk);

  Real tol = 1e-10;
  Real h = 1e-30;

  bool pass = false;

  Complex Xc[3];
  Complex phi;

  std::copy(p.X, p.X+3, Xc);
  Xc[0] += Complex(0, h);
  phi = p.elem->eval_basis(I, Xc);
  Real d0phi = common::imag(phi)/h;

  std::copy(p.X, p.X+3, Xc);
  Xc[1] += Complex(0, h);
  phi = p.elem->eval_basis(I, Xc);
  Real d1phi = common::imag(phi)/h;

  std::copy(p.X, p.X+3, Xc);
  Xc[2] += Complex(0, h);
  phi = p.elem->eval_basis(I, Xc);
  Real d2phi = common::imag(phi)/h;

  Complex grad_phi[3];
  p.elem->eval_grad_basis(I, p.X, grad_phi);

  Real rel_error0 = common::abs((grad_phi[0] - d0phi)/d0phi);
  Real rel_error1 = common::abs((grad_phi[1] - d1phi)/d1phi);
  Real rel_error2 = common::abs((grad_phi[2] - d2phi)/d2phi);

  if (rel_error0 < tol && rel_error1 < tol && rel_error2 < tol) {
    pass = true;
  } else {
    std::cout << "Test 3, errors: " 
              << rel_error0 << " " 
              << rel_error1 << " " 
              << rel_error2 << " " 
              << std::endl; 
  }

  return pass;
}

int main() {
  std::cout.precision(15);

  std::cout << "TEST LEGENDRE ELEMENT\n" 
            << "---------------------" 
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
