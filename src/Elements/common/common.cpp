#include "common.h"

/*
 * Convert a number in its base 10 representation to its representation in
 * mixed-radix form for a provided set of radices
 *
 * Parameters
 * ----------
 * Nb : number of radices or bases
 * b  : an array of radices or bases
 * x  : an non-negative integer in base 10
 *
 * Returns
 * -------
 * y  : little-endian array of Nb entries, the mixed-radix form
 */
void common::base_10_to_mixed_radix(const SizeType Nb, const SizeType *b, 
    SizeType x, SizeType *y) {
  for (SizeType i = 0; i < Nb; i++) {
    y[i] = x % b[i];
    x /= b[i];
  }
};

/*
 * Convert a number in mixed-radix form to its base 10 representation for a
 * provided set of radices
 *
 * Parameters
 * ----------
 * Nb : number of radices or bases
 * b  : an array of radices or bases
 * x  : little-endian array of Nb entries, the mixed-radix form
 *
 * Returns
 * -------
 * y  : an non-negative integer in base 10
 */
SizeType common::mixed_radix_to_base_10(const SizeType Nb, const SizeType *b, 
    SizeType *x) {
  SizeType y = 0;
  SizeType z = 1;
  for (SizeType i = 0; i < Nb; i++) {
    y += x[i]*z; 
    z *= b[i];
  }

  return y;
};

/*
 * Encoding the information specifying a partial derivative of a function in
 * three variables. The information is the order of the partial derivative in
 * each of the three variables (x, y, z).
 *
 * The encoding is simple. A value of 1 in one of the first three bits
 * corresponds to a first-order derivative with respect to x, y, or z,
 * depending on whether it is the first (x), second (y), or third (z) bit.
 * The scheme generalizes to higher orders by shifting the bits up three
 * places. 
 *
 * Parameters
 * ----------
 * nx : order of the partial derivative with respect to x
 * ny : order of the partial derivative with respect to y
 * nz : order of the partial derivative with respect to z
 *
 * Returns
 * -------
 * e : encoding
 */
SizeType common::encode_partial_derivative(const SizeType &nx, 
    const SizeType &ny, const SizeType &nz) {
  SizeType e = (nx > 0) << (0 + 3*(nx - 1)) | (ny > 0) << (1 + 3*(ny - 1)) 
      | (nz > 0) << (2 + 3*(nz - 1));

  return e;
};

/*
 * Decoding the partial derivative information encoding above
 *
 * Parameters
 * ----------
 * e : encoding
 *
 * Returns
 * ----------
 * nx : order of the partial derivative with respect to x
 * ny : order of the partial derivative with respect to y
 * nz : order of the partial derivative with respect to z
 */
void common::decode_partial_derivative(SizeType e, SizeType &nx, SizeType &ny, 
    SizeType &nz) {
  nx = 0, ny = 0, nz = 0; 
  SizeType p = 1;  // order of derivative

  while (e > 0) {
    nx += ((e & 1) > 0)*p;
    ny += ((e & 2) > 0)*p;
    nz += ((e & 4) > 0)*p;
    e >>= 3;
    p += 1;
  }
};
