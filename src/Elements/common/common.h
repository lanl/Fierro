#pragma once

#include "matar.h"

#include <cmath>
#include <complex>
#include <limits>

typedef double Real;
typedef std::complex<double> Complex;
typedef size_t SizeType;

#define NUM_EPS std::numeric_limits<Real>::epsilon()
#define NUM_MIN std::numeric_limits<Real>::min()
#define NUM_MAX std::numeric_limits<Real>::max()

namespace common {
  // Definitions used to ensure compatibility when switching from real number
  // type to complex number type to test derivative implementations via the
  // complex step method
  inline Real real(Real number) { return number; }
  inline Real real(Complex number) { return number.real(); }

  inline Real imag(Real number) { return 0.0; }
  inline Real imag(Complex number) { return number.imag(); }

  inline Real abs(Real number) { return std::abs(number); }
  inline Real abs(Complex number) { 
      return std::abs(number.real()); }

  // Check equality of two double precision floating point numbers
  template <typename NumType>
  inline bool almost_equal(NumType a, NumType b) {
    return common::abs(a - b) < 2.0*NUM_EPS;
  };

  // Converting between representations of array indices
  void base_10_to_mixed_radix(const SizeType Nb, const SizeType *b, 
      SizeType x, SizeType *y);
  SizeType mixed_radix_to_base_10(const SizeType Nb, const SizeType *b, 
      SizeType *x);

  // Encoding and decoding a requested partial derivative to and from an
  // unsigned integer
  SizeType encode_partial_derivative(const SizeType &nx, const SizeType &ny, 
      const SizeType &nz);
  void decode_partial_derivative(SizeType e, SizeType &nx, SizeType &ny, 
      SizeType &nz);
}
