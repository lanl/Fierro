#pragma once

#include "definitions.h"

using namespace utils;

template<typename T>
KOKKOS_FUNCTION
T PowIntExpo(T base, int exponent);

KOKKOS_FUNCTION
void swap_rows(double *matrix, int row1, int row2, int n);
 
KOKKOS_FUNCTION
void scale_row(double *matrix, int row, double factor, int n);

KOKKOS_FUNCTION
void add_rows(double *matrix, int src_row, int dest_row, double factor, int n);

template<int n>
KOKKOS_FUNCTION
int invert_matrix(double *matrix);

KOKKOS_FUNCTION
int solve_linear_system(double* A, double* b, int n);



template<int n>
KOKKOS_FUNCTION
int invert_matrix(double *matrix)
{
    int error_flag=0;
    double identity[n*n];
    if (identity == NULL) {
        printf("Error: Failed to allocate memory for identity matrix.\n");
        return 1;
    }

    // Initialize identity matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            *(identity + i * n + j) = (i == j) ? 1.0 : 0.0;
        }
    }

    // Perform Gauss-Jordan elimination
    for (int i = 0; i < n; i++) {
        if (*(matrix + i * n + i) == 0) {
            //printf("Error: Matrix is not invertible.\n");
            return error_flag;
        }

        double pivot = *(matrix + i * n + i);
        scale_row(matrix, i, 1.0 / pivot, n);
        scale_row(identity, i, 1.0 / pivot, n);

        for (int j = 0; j < n; j++) {
            if (j != i) {
                double factor = -*(matrix + j * n + i);
                add_rows(matrix, i, j, factor, n);
                add_rows(identity, i, j, factor, n);
            }
        }
    }

    // Copy the inverted matrix to the input matrix pointer
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            *(matrix + i * n + j) = *(identity + i * n + j);
        }
    }

    error_flag = 0;
    return error_flag;
}
// End New Gauss-Jordan elimination matrix inverse functions


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