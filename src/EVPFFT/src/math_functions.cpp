#include <stdio.h>
#include "math_functions.h"


// New Gauss-Jordan elimination matrix inverse functions 
KOKKOS_FUNCTION
void swap_rows(double *matrix, int row1, int row2, int n) {
    for (int i = 0; i < n; i++) {
        double temp = *(matrix + row1 * n + i);
        *(matrix + row1 * n + i) = *(matrix + row2 * n + i);
        *(matrix + row2 * n + i) = temp;
    }
}

KOKKOS_FUNCTION
void scale_row(double *matrix, int row, double factor, int n) {
    for (int i = 0; i < n; i++) {
        *(matrix + row * n + i) *= factor;
    }
}

KOKKOS_FUNCTION
void add_rows(double *matrix, int src_row, int dest_row, double factor, int n) {
    for (int i = 0; i < n; i++) {
        *(matrix + dest_row * n + i) += *(matrix + src_row * n + i) * factor;
    }
}

KOKKOS_FUNCTION
int invert_matrix(double *matrix, int n) {
    int error_flag=0;
    double *identity = new double[n*n];
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
            delete [] identity;
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

    delete [] identity;
    error_flag = 0;
    return error_flag;
}
// End New Gauss-Jordan elimination matrix inverse functions

KOKKOS_FUNCTION
int solve_linear_system(double* A, double* b, int n) {
    // Function to solve linear system using Gaussian elimination
    
    int error_flag=1;
    
    for (int i = 0; i < n; i++) {
        // Find the pivot row
        int pivotRow = i;
        for (int j = i + 1; j < n; j++) {
            if (fabs(A[j * n + i]) > fabs(A[pivotRow * n + i])) {
                pivotRow = j;
            }   
        }   
        
        // Swap the current row with the pivot row using the provided function
        if (pivotRow != i) {
            swap_rows(A, i, pivotRow, n);
            
            // C-style swap for the b vector elements
            double temp = b[i];
            b[i] = b[pivotRow];
            b[pivotRow] = temp;
        }   
        
        // Check for singular or nearly singular matrix
        if (A[i * n + i] == 0.0) {
            return error_flag; // Matrix is singular or nearly singular
        }   
        
        // Make the diagonal element 1
        double pivot = A[i * n + i];
        for (int j = i; j < n; j++) {
            A[i * n + j] /= pivot;
        }   
        b[i] /= pivot;
        
        // Eliminate other rows
        for (int j = 0; j < n; j++) {
            if (j != i) { 
                double factor = A[j * n + i];
                for (int k = i; k < n; k++) {
                    A[j * n + k] -= factor * A[i * n + k];
                }   
                b[j] -= factor * b[i];
            }   
        }   
    }   
    
    error_flag = 0;
    return error_flag; // Matrix is non-singular
}
