#pragma once

#include "definitions.h"

using namespace utils;

template<typename T>
KOKKOS_FUNCTION
T PowIntExpo(T base, int exponent);

KOKKOS_FUNCTION
int invert_matrix(double *matrix, int n);

KOKKOS_FUNCTION
int solve_linear_system(double* A, double* b, int n);

// Optimized pow for integer exponents only
template<typename T>
KOKKOS_FUNCTION
T PowIntExpo(T base, int exponent) {
    if (exponent == 0)
        return 1;
    else if (exponent == 1)
        return base;

    T result = 1;
    bool negativeExponent = (exponent < 0);
    exponent = abs(exponent);

    while (exponent > 0) {
        if (exponent & 1)
            result *= base;
        base *= base;
        exponent >>= 1;
    }

    return negativeExponent ? 1.0 / result : result;
}


