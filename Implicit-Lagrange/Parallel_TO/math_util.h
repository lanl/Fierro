#ifndef MATH_UTILITY_H
#define MATH_UTILITY_H  

#include "utilities.h"

using namespace utils;

namespace MathUtility{

  int LUPDecompose(real_t **A, int N, real_t Tol, int *P);
  void LUPSolve(real_t **A, int *P, real_t *b, int N, real_t *x);
};

#endif // end HEADER_H
