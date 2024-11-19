#ifndef LIN_ALG
#define LIN_ALG

#include "state.h"
#include "mesh.h"
#include "rdh.h"
#define TINY 1.e-16
#include <Kokkos_Core.hpp>

#include <Kokkos_Macros.hpp>
#include <cmath> 

#ifndef PI
#define PI 3.14159265358979323846
#endif

// LU Decomposition function
KOKKOS_INLINE_FUNCTION
int lu_decomp(double source_mat[50][50], int indx[50], int &parity, const int n) {
    double vv[50];
    parity = 1;

    // Search for the largest element in each row
    for (int i = 0; i < n; i++) {
        double big = 0.0;
        for (int j = 0; j < n; j++) {
            double temp = fabs(source_mat[i][j]);
            if (temp > big) {
                big = temp;
            }
        }
        if (big == 0.0) {
            return -1; // Singular matrix, decomposition fails
        }
        vv[i] = 1.0 / big;
    }

    // Perform Gaussian elimination with partial pivoting
    for (int k = 0; k < n; k++) {
        double big = 0.0;
        int imax = k;
        for (int i = k; i < n; i++) {
            double temp = vv[i] * fabs(source_mat[i][k]);
            if (temp > big) {
                big = temp;
                imax = i;
            }
        }

        if (k != imax) {
            for (int j = 0; j < n; j++) {
                double temp = source_mat[imax][j];
                source_mat[imax][j] = source_mat[k][j];
                source_mat[k][j] = temp;
            }
            parity = -parity;
            vv[imax] = vv[k];
        }

        indx[k] = imax;
        if (fabs(source_mat[k][k]) < TINY) {
            source_mat[k][k] = (source_mat[k][k] < 0) ? -TINY : TINY;
        }

        double diag_inv = 1.0 / source_mat[k][k];
        for (int i = k + 1; i < n; i++) {
            double factor = source_mat[i][k] * diag_inv;
            source_mat[i][k] = factor;
            for (int j = k + 1; j < n; j++) {
                source_mat[i][j] -= factor * source_mat[k][j];
            }
        }
    }

    return 0; // Assuming successful decomposition
}


// LU Inversion function
KOKKOS_INLINE_FUNCTION
void lu_invert(double lu_mtx[50][50], Kokkos::View<double**> mtx_inv, double col[50], int indx[50], int n, int col_index) {
    double x[50];

    // Forward substitution (Ly = b)
    for (int i = 0; i < n; i++) {
        int ip = indx[i];
        double sum = col[ip];
        col[ip] = col[i];
        for (int j = 0; j < i; j++) {
            sum -= lu_mtx[i][j] * x[j];
        }
        x[i] = sum;
    }

    // Backward substitution (Ux = y)
    for (int i = n - 1; i >= 0; i--) {
        double sum = x[i];
        for (int j = i + 1; j < n; j++) {
            sum -= lu_mtx[i][j] * x[j];
        }
        x[i] = sum / lu_mtx[i][i];
    }

    // Copy solution into the appropriate column of the inverse matrix
    for (int i = 0; i < n; i++) {
        mtx_inv(i, col_index) = x[i];
    }
}


KOKKOS_INLINE_FUNCTION
void invert_matrix(Kokkos::View<double**> mtx, Kokkos::View<double**> mtx_inv, int size) {
    double col[50]; // Use standard arrays instead of Kokkos::View
    int index[50];
    double lu_mtx[50][50];

    int parity;

    // Copy mtx to lu_mtx for LU decomposition
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            lu_mtx[i][j] = mtx(i, j);
        }
    }

    // LU decomposition
    if (lu_decomp(lu_mtx, index, parity, size) != 0) {
        // If decomposition fails, return early
        return;
    }

    // Invert matrix using LU decomposition
    for (int j = 0; j < size; j++) {
        // Set up the column of the identity matrix
        for (int i = 0; i < size; i++) {
            col[i] = (i == j) ? 1.0 : 0.0;
        }
        
        // Solve for each column of the inverse matrix
        lu_invert(lu_mtx, mtx_inv, col, index, size, j);
    }
}


KOKKOS_INLINE_FUNCTION
void get_spectrum(const double matrix[3][3], double eigenvalues[3], double eigenvectors[3][3]) {
    // Compute the eigenvalues and eigenvectors of a symmetric 3x3 matrix

    // Copy matrix elements to local variables for easier reference
    double d11 = matrix[0][0];
    double d12 = matrix[0][1]; // Since the matrix is symmetric, matrix[1][0] == matrix[0][1]
    double d13 = matrix[0][2];
    double d22 = matrix[1][1];
    double d23 = matrix[1][2];
    double d33 = matrix[2][2];

    // Scaling to improve numerical stability
    double d_max = fabs(d11);
    d_max = fmax(d_max, fabs(d22));
    d_max = fmax(d_max, fabs(d33));
    d_max = fmax(d_max, fabs(d12));
    d_max = fmax(d_max, fabs(d13));
    d_max = fmax(d_max, fabs(d23));

    double mult = 1.0;
    if (d_max > 0.0) {
        mult = 1.0 / d_max;
    }

    // Scale the matrix entries
    d11 *= mult;  d22 *= mult;  d33 *= mult;
    d12 *= mult;  d13 *= mult;  d23 *= mult;

    // Center the matrix by subtracting the mean of the diagonal elements
    double aa = (d11 + d22 + d33) / 3.0;  // Mean of the diagonal elements
    double c1 = d11 - aa;
    double c2 = d22 - aa;
    double c3 = d33 - aa;

    // Compute Q and R for the characteristic equation
    double Q = (2.0 * (d12*d12 + d13*d13 + d23*d23) + c1*c1 + c2*c2 + c3*c3) / 6.0;
    double R = (c1 * (d23*d23 - c2*c3) + d12 * (d12 * c3 - 2.0 * d13 * d23) + d13 * d13 * c2) / 2.0;

    double eigenvalue1, eigenvalue2, eigenvalue3;

    if (Q <= 1e-12) {
        // All eigenvalues are equal
        eigenvalue1 = eigenvalue2 = eigenvalue3 = aa / mult;
        // Eigenvectors are the standard basis vectors
        eigenvectors[0][0] = 1.0; eigenvectors[0][1] = 0.0; eigenvectors[0][2] = 0.0;
        eigenvectors[1][0] = 0.0; eigenvectors[1][1] = 1.0; eigenvectors[1][2] = 0.0;
        eigenvectors[2][0] = 0.0; eigenvectors[2][1] = 0.0; eigenvectors[2][2] = 1.0;

        eigenvalues[0] = eigenvalue1;
        eigenvalues[1] = eigenvalue2;
        eigenvalues[2] = eigenvalue3;
    } else {
        double sqrtQ = sqrt(Q);
        double R_div_sqrtQ3 = R / (sqrtQ * Q);

        // Ensure R_div_sqrtQ3 is within the valid range of acos
        R_div_sqrtQ3 = fmax(-1.0, fmin(1.0, R_div_sqrtQ3));

        // Compute the three eigenvalues using trigonometric identities
        double theta = acos(R_div_sqrtQ3) / 3.0;
        double sqrtQ2 = 2.0 * sqrtQ;
        eigenvalue1 = aa + sqrtQ2 * cos(theta);
        eigenvalue2 = aa + sqrtQ2 * cos(theta + (2.0 * PI / 3.0));
        eigenvalue3 = aa + sqrtQ2 * cos(theta + (4.0 * PI / 3.0));

        // Scale the eigenvalues back
        eigenvalue1 /= mult;
        eigenvalue2 /= mult;
        eigenvalue3 /= mult;

        // Assign the eigenvalues to the array
        eigenvalues[0] = eigenvalue1;
        eigenvalues[1] = eigenvalue2;
        eigenvalues[2] = eigenvalue3;

        // Compute eigenvectors for each eigenvalue
        for (int idx = 0; idx < 3; ++idx) {
            double lambda = eigenvalues[idx]; // Eigenvalues are already unscaled

            // Form the matrix (A - lambda * I)
            double A[3][3];
            A[0][0] = matrix[0][0] - lambda;
            A[0][1] = matrix[0][1];
            A[0][2] = matrix[0][2];
            A[1][0] = matrix[1][0];
            A[1][1] = matrix[1][1] - lambda;
            A[1][2] = matrix[1][2];
            A[2][0] = matrix[2][0];
            A[2][1] = matrix[2][1];
            A[2][2] = matrix[2][2] - lambda;

            // Compute the eigenvector corresponding to the eigenvalue
            double v[3];

            // Use the cross product of two rows of A to get a vector in the null space
            double r0[3] = {A[0][0], A[0][1], A[0][2]};
            double r1[3] = {A[1][0], A[1][1], A[1][2]};
            double r2[3] = {A[2][0], A[2][1], A[2][2]};

            // Compute the cross product of any two rows
            double cp1[3] = {
                r0[1]*r1[2] - r0[2]*r1[1],
                r0[2]*r1[0] - r0[0]*r1[2],
                r0[0]*r1[1] - r0[1]*r1[0]
            };

            double norm_cp1 = sqrt(cp1[0]*cp1[0] + cp1[1]*cp1[1] + cp1[2]*cp1[2]);

            if (norm_cp1 > 1e-6) {
                // Use cp1 as the eigenvector
                v[0] = cp1[0];
                v[1] = cp1[1];
                v[2] = cp1[2];
            } else {
                // If cp1 is too small, use the cross product of another pair
                double cp2[3] = {
                    r0[1]*r2[2] - r0[2]*r2[1],
                    r0[2]*r2[0] - r0[0]*r2[2],
                    r0[0]*r2[1] - r0[1]*r2[0]
                };
                double norm_cp2 = sqrt(cp2[0]*cp2[0] + cp2[1]*cp2[1] + cp2[2]*cp2[2]);

                if (norm_cp2 > 1e-6) {
                    v[0] = cp2[0];
                    v[1] = cp2[1];
                    v[2] = cp2[2];
                } else {
                    // Last resort, use the cross product of the last pair
                    double cp3[3] = {
                        r1[1]*r2[2] - r1[2]*r2[1],
                        r1[2]*r2[0] - r1[0]*r2[2],
                        r1[0]*r2[1] - r1[1]*r2[0]
                    };
                    v[0] = cp3[0];
                    v[1] = cp3[1];
                    v[2] = cp3[2];
                }
            }

            // Normalize the eigenvector
            double norm_v = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
            if (norm_v > 1e-12) {
                v[0] /= norm_v;
                v[1] /= norm_v;
                v[2] /= norm_v;
            } else {
                // If the norm is zero (should not happen), set eigenvector to zero vector
                v[0] = v[1] = v[2] = 0.0;
            }

            // Assign the eigenvector
            eigenvectors[idx][0] = v[0];
            eigenvectors[idx][1] = v[1];
            eigenvectors[idx][2] = v[2];
        }

        // Sort the eigenvalues and corresponding eigenvectors in ascending order
        for (int i = 0; i < 2; ++i) {
            int min_idx = i;
            for (int j = i + 1; j < 3; ++j) {
                if (eigenvalues[j] < eigenvalues[min_idx]) {
                    min_idx = j;
                }
            }
            if (min_idx != i) {
                // Swap eigenvalues
                double temp_val = eigenvalues[i];
                eigenvalues[i] = eigenvalues[min_idx];
                eigenvalues[min_idx] = temp_val;
                // Swap eigenvectors
                for (int k = 0; k < 3; ++k) {
                    double temp_vec = eigenvectors[i][k];
                    eigenvectors[i][k] = eigenvectors[min_idx][k];
                    eigenvectors[min_idx][k] = temp_vec;
                }
            }
        }
    }
}


KOKKOS_INLINE_FUNCTION
void compute_Fro_norm(const double matrix[3][3],
                     double &norm){
        
		norm = 0.0;
		double temp = 0.0;

		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++){
				temp += abs( matrix[i][j]*matrix[i][j] );
			}// j
		}// i
		
        norm = sqrt(temp);

}// end compute Fro norm

KOKKOS_INLINE_FUNCTION
void compute_l2_norm(const double vector[3],
                     double &norm){
    	norm = 0.0;

        norm = sqrt( vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2] );

}// end compute l2 norm

KOKKOS_INLINE_FUNCTION
void symmetrize_matrix(const double matrix[3][3],
						double sym_matrix[3][3]){
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			sym_matrix[i][j] = 0.0;
		}// j
	}// i

	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			sym_matrix[i][j] = 0.5*(matrix[i][j] + matrix[j][i]);
		}// j
	}// i

}// end symmetrix_matrix

#endif
