#ifndef LIN_ALG
#define LIN_ALG


#include <cmath>
#include <cstdio>
#include <algorithm> // For std::swap
#include "state.h"
#include "mesh.h"
#include "rdh.h"
#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>

#define TINY 1.e-16

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

// Swap function
KOKKOS_INLINE_FUNCTION
void Swap(double &a, double &b) {
    double tmp = a;
    a = b;
    b = tmp;
}


// Auxiliary normalization function to normalize a 3D vector by avoiding overflow/underflow
KOKKOS_INLINE_FUNCTION
void Vec_normalize3_aux(const double &x1, const double &x2, const double &x3,
                        double &n1, double &n2, double &n3) {
    double t, r;

    const double m = fabs(x1);
    if (m > 0.0) {
        r = x2 / m;
        t = 1.0 + r * r;
        r = x3 / m;
        t = sqrt(1.0 / (t + r * r));
        n1 = copysign(t, x1);
        t /= m;
        n2 = x2 * t;
        n3 = x3 * t;
    } else {
        n1 = n2 = n3 = 0.0;
    }
}

// Main normalization function for a 3D vector
KOKKOS_INLINE_FUNCTION
void Vec_normalize3(const double &x1, const double &x2, const double &x3,
                    double &n1, double &n2, double &n3) {
    // Identify the largest component to prevent overflow
    if (fabs(x1) >= fabs(x2) && fabs(x1) >= fabs(x3)) {
        if (x1 != 0.0) {
            Vec_normalize3_aux(x1, x2, x3, n1, n2, n3);
        } else {
            n1 = n2 = n3 = 0.0;
        }
    } else if (fabs(x2) >= fabs(x3)) {
        Vec_normalize3_aux(x2, x1, x3, n2, n1, n3);
    } else {
        Vec_normalize3_aux(x3, x1, x2, n3, n1, n2);
    }
}

// Compute the eigensystem of a 2x2 symmetric matrix
KOKKOS_INLINE_FUNCTION
void Eigensystem2S(const double &d12, double &d1, double &d2,
                   double &c, double &s) {
    const double sqrt_1_eps = sqrt(1.0 / std::numeric_limits<double>::epsilon());
    if (d12 == 0.0) {
        c = 1.0;
        s = 0.0;
    } else {
        double t;
        const double zeta = (d2 - d1) / (2.0 * d12);
        const double azeta = fabs(zeta);
        if (azeta < sqrt_1_eps) {
            t = copysign(1.0 / (azeta + sqrt(1.0 + zeta * zeta)), zeta);
        } else {
            t = copysign(0.5 / azeta, zeta);
        }
        c = 1.0 / sqrt(1.0 + t * t);
        s = c * t;
        t *= d12;
        d1 -= t;
        d2 += t;
    }
}


// Compute a vector in the "near"-kernel of a 2x2 matrix
KOKKOS_INLINE_FUNCTION
bool KernelVector2G(const int &mode,
                    double &d1, double &d12, double &d21, double &d2) {
    // Find a vector (z1, z2) in the "near"-kernel of the matrix
    double n1 = fabs(d1) + fabs(d21);
    double n2 = fabs(d2) + fabs(d12);

    bool swap_columns = (n2 > n1);
    double mu;

    if (!swap_columns) {
        if (n1 == 0.0) {
            return true;
        }
        if (mode == 0) {
            if (fabs(d1) > fabs(d21)) {
                Swap(d1, d21);
                Swap(d12, d2);
            }
        } else {
            if (fabs(d1) < fabs(d21)) {
                Swap(d1, d21);
                Swap(d12, d2);
            }
        }
    } else {
        if (mode == 0) {
            if (fabs(d12) > fabs(d2)) {
                Swap(d1, d2);
                Swap(d12, d21);
            } else {
                Swap(d1, d12);
                Swap(d21, d2);
            }
        } else {
            if (fabs(d12) < fabs(d2)) {
                Swap(d1, d2);
                Swap(d12, d21);
            } else {
                Swap(d1, d12);
                Swap(d21, d2);
            }
        }
    }

    n1 = hypot(d1, d21);
    if (d21 != 0.0) {
        mu = copysign(n1, d1);
        double temp_n1 = -d21 * (d21 / (d1 + mu));
        d1 = mu;
        if (fabs(temp_n1) <= fabs(d21)) {
            double n1_over_d21 = temp_n1 / d21;
            mu = (2.0 / (1.0 + n1_over_d21 * n1_over_d21)) * (n1_over_d21 * d12 + d2);
            d2 -= mu;
            d12 -= mu * n1_over_d21;
        } else {
            double d21_over_n1 = d21 / temp_n1;
            mu = (2.0 / (1.0 + d21_over_n1 * d21_over_n1)) * (d12 + d21_over_n1 * d2);
            d2 -= mu * d21_over_n1;
            d12 -= mu;
        }
    }

    mu = -d12 / d1;
    n2 = 1.0 / (1.0 + fabs(mu));

    if (fabs(d1) <= n2 * fabs(d2)) {
        d2 = 0.0;
        d1 = 1.0;
    } else {
        d2 = n2;
        d1 = mu * n2;
    }

    if (swap_columns) {
        Swap(d1, d2);
    }

    return false;
}

// Auxiliary function to compute the kernel vector for a 3x3 matrix
KOKKOS_INLINE_FUNCTION
int KernelVector3G_aux(const int &mode,
                       double &d1, double &d2, double &d3,
                       double &c12, double &c13, double &c23,
                       double &c21, double &c31, double &c32) {
    int kdim;
    double mu, n1, n2, n3, s1, s2, s3;

    s1 = hypot(c21, c31);
    n1 = hypot(d1, s1);

    if (s1 != 0.0) {
        mu = copysign(n1, d1);
        n1 = -s1 * (s1 / (d1 + mu));
        d1 = mu;

        double max_n = fmax(fabs(n1), fmax(fabs(c21), fabs(c31)));
        double s_norm = (max_n != 0.0) ? max_n : 1.0;
        s2 = c21 / s_norm;
        s3 = c31 / s_norm;
        mu = 2.0 / (1.0 + s2 * s2 + s3 * s3);

        n2 = mu * (c12 + s2 * d2 + s3 * c32);
        n3 = mu * (c13 + s2 * c23 + s3 * d3);

        c12 -= n2;
        d2 -= s2 * n2;
        c32 -= s3 * n2;
        c13 -= n3;
        c23 -= s2 * n3;
        d3 -= s3 * n3;
    }

    if (KernelVector2G(mode, d2, c23, c32, d3)) {
        d2 = c12 / d1;
        d3 = c13 / d1;
        d1 = 1.0;
        kdim = 2;
    } else {
        d1 = -(c12 * d2 + c13 * d3) / d1;
        kdim = 1;
    }

    Vec_normalize3(d1, d2, d3, d1, d2, d3);

    return kdim;
}

// Compute a vector in the "near"-kernel of a symmetric 3x3 matrix
KOKKOS_INLINE_FUNCTION
int KernelVector3S(const int &mode, const double &d12,
                   const double &d13, const double &d23,
                   double &d1, double &d2, double &d3) {
    double c12 = d12, c13 = d13, c23 = d23;
    double c21, c31, c32;
    int col, row;

    double n1 = fabs(d1) + fabs(c12) + fabs(c13);
    double n2 = fabs(d2) + fabs(c12) + fabs(c23);
    double n3 = fabs(d3) + fabs(c13) + fabs(c23);

    if (n1 >= n3) {
        col = (n1 >= n2) ? 1 : 2;
    } else {
        col = (n2 >= n3) ? 2 : 3;
    }

    switch (col) {
        case 1:
            if (n1 == 0.0) {
                return 3;
            }
            break;
        case 2:
            if (n2 == 0.0) {
                return 3;
            }
            Swap(c13, c23);
            Swap(d1, d2);
            break;
        case 3:
            if (n3 == 0.0) {
                return 3;
            }
            Swap(c12, c23);
            Swap(d1, d3);
            break;
    }

    if (mode == 0) {
        if (fabs(d1) <= fabs(c13)) {
            row = (fabs(d1) <= fabs(c12)) ? 1 : 2;
        } else {
            row = (fabs(c12) <= fabs(c13)) ? 2 : 3;
        }
    } else {
        if (fabs(d1) >= fabs(c13)) {
            row = (fabs(d1) >= fabs(c12)) ? 1 : 2;
        } else {
            row = (fabs(c12) >= fabs(c13)) ? 2 : 3;
        }
    }

    switch (row) {
        case 1:
            c21 = c12;
            c31 = c13;
            c32 = c23;
            break;
        case 2:
            c21 = d1;
            c31 = c13;
            c32 = c23;
            d1 = c12;
            c12 = d2;
            d2 = d1;
            c13 = c23;
            c23 = c31;
            break;
        case 3:
            c21 = c12;
            c31 = d1;
            c32 = c12;
            d1 = c13;
            c12 = c23;
            c13 = d3;
            d3 = d1;
            break;
    }

    int kdim = KernelVector3G_aux(mode, d1, d2, d3, c12, c13, c23, c21, c31, c32);

    switch (col) {
        case 2:
            Swap(d1, d2);
            break;
        case 3:
            Swap(d1, d3);
            break;
    }

    return kdim;
}

// Implement Reduce3S
KOKKOS_INLINE_FUNCTION
int Reduce3S(const int &mode,
             double &d1, double &d2, double &d3,
             double &d12, double &d13, double &d23,
             double &z1, double &z2, double &z3,
             double &v1, double &v2, double &v3,
             double &g) {
    int k;
    double s, w1, w2, w3;

    if (mode == 0) {
        if (fabs(z1) <= fabs(z3)) {
            k = (fabs(z1) <= fabs(z2)) ? 1 : 2;
        } else {
            k = (fabs(z2) <= fabs(z3)) ? 2 : 3;
        }
    } else {
        if (fabs(z1) >= fabs(z3)) {
            k = (fabs(z1) >= fabs(z2)) ? 1 : 2;
        } else {
            k = (fabs(z2) >= fabs(z3)) ? 2 : 3;
        }
    }

    switch (k) {
        case 2:
            Swap(d13, d23);
            Swap(d1, d2);
            Swap(z1, z2);
            break;
        case 3:
            Swap(d12, d23);
            Swap(d1, d3);
            Swap(z1, z3);
            break;
    }

    s = hypot(z2, z3);

    if (s == 0.0) {
        v1 = v2 = v3 = 0.0;
        g = 1.0;
    } else {
        g = copysign(1.0, z1);
        v1 = -s * (s / (z1 + g)); // v1 = z1 - g

        // Normalize (v1, z2, z3) by max-norm to avoid overflow/underflow
        double max_v = fmax(fabs(v1), fmax(fabs(z2), fabs(z3)));
        v1 /= max_v;
        v2 = z2 / max_v;
        v3 = z3 / max_v;
        g = 2.0 / (v1 * v1 + v2 * v2 + v3 * v3);

        // Compute w = g * A * v
        w1 = g * (d1 * v1 + d12 * v2 + d13 * v3);
        w2 = g * (d12 * v1 + d2 * v2 + d23 * v3);
        w3 = g * (d13 * v1 + d23 * v2 + d3 * v3);

        // w := w - (g/2)(v^t w) v
        s = (g / 2.0) * (v1 * w1 + v2 * w2 + v3 * w3);
        w1 -= s * v1;
        w2 -= s * v2;
        w3 -= s * v3;

        // Update matrix entries
        d1 -= 2.0 * v1 * w1;
        d2 -= 2.0 * v2 * w2;
        d23 -= v2 * w3 + v3 * w2;
        d3 -= 2.0 * v3 * w3;
    }

    switch (k) {
        case 2:
            Swap(z1, z2);
            break;
        case 3:
            Swap(z1, z3);
            break;
    }

    return k;
}
// Utility function used to get the scaling factor for normalizing matrix entries
KOKKOS_INLINE_FUNCTION
void GetScalingFactor(const double &d_max, double &mult) {
    int d_exp;

    if (d_max > 0.0) {
        mult = frexp(d_max, &d_exp);
        
        // Handle edge case where d_max has the maximum possible exponent
        if (d_exp == std::numeric_limits<double>::max_exponent) {
            mult *= std::numeric_limits<double>::radix;
        }
        mult = d_max / mult;
    } else {
        mult = 1.0;
    }

    // `mult = 2^d_exp` ensures that `d_max / mult` is in the interval [0.5, 1).
}


// Gram-Schmidt Orthogonalization to maintain orthogonality of eigenvectors (columns)
KOKKOS_INLINE_FUNCTION
void GramSchmidtOrthogonalize(double vecs[3][3]) {
    for(int i = 1; i < 3; ++i) {
        for(int j = 0; j < i; ++j) {
            double dot = vecs[0][j] * vecs[0][i] +
                         vecs[1][j] * vecs[1][i] +
                         vecs[2][j] * vecs[2][i];
            vecs[0][i] -= dot * vecs[0][j];
            vecs[1][i] -= dot * vecs[1][j];
            vecs[2][i] -= dot * vecs[2][j];
        }
        // Normalize the column after orthogonalizing
        double norm = sqrt(vecs[0][i] * vecs[0][i] +
                           vecs[1][i] * vecs[1][i] +
                           vecs[2][i] * vecs[2][i]);
        if(norm > 1e-12){
            vecs[0][i] /= norm;
            vecs[1][i] /= norm;
            vecs[2][i] /= norm;
        }
    }
}

// Jacobi Rotation for eigenvalue computation
KOKKOS_INLINE_FUNCTION
void ApplyJacobiRotation(double A[3][3], double V[3][3], int p, int q) {
    if (A[p][q] != 0.0) {
        double tau = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
        double t = (tau > 0.0 ? 1.0 : -1.0) / (fabs(tau) + sqrt(1.0 + tau * tau));
        double c = 1.0 / sqrt(1.0 + t * t);
        double s = t * c;

        double app = A[p][p];
        double aqq = A[q][q];
        double apq = A[p][q];

        // Update the matrix A with rotation
        A[p][p] = c * c * app + s * s * aqq - 2.0 * s * c * apq;
        A[q][q] = s * s * app + c * c * aqq + 2.0 * s * c * apq;
        A[p][q] = A[q][p] = 0.0;

        // Update the other elements
        for (int r = 0; r < 3; ++r) {
            if (r != p && r != q) {
                double arp = A[r][p];
                double arq = A[r][q];
                A[r][p] = A[p][r] = c * arp - s * arq;
                A[r][q] = A[q][r] = c * arq + s * arp;
            }
        }

        // Update the eigenvector matrix V with rotation
        for (int r = 0; r < 3; ++r) {
            double vrp = V[r][p];
            double vrq = V[r][q];
            V[r][p] = c * vrp - s * vrq;
            V[r][q] = c * vrq + s * vrp;
        }
    }
}

// Jacobi Method for Eigenvalue and Eigenvector Computation
KOKKOS_INLINE_FUNCTION
void get_spectrum(const double data[3][3], double eigenvalues[3], double eigenvectors[3][3]) {
    // Initialize the matrix A as a copy of data
    double A[3][3] = {
        {data[0][0], data[0][1], data[0][2]},
        {data[1][0], data[1][1], data[1][2]},
        {data[2][0], data[2][1], data[2][2]}
    };

    // Initialize eigenvectors as the identity matrix
    eigenvectors[0][0] = 1.0; eigenvectors[0][1] = 0.0; eigenvectors[0][2] = 0.0;
    eigenvectors[1][0] = 0.0; eigenvectors[1][1] = 1.0; eigenvectors[1][2] = 0.0;
    eigenvectors[2][0] = 0.0; eigenvectors[2][1] = 0.0; eigenvectors[2][2] = 1.0;

    // Jacobi method iteration parameters
    const int MAX_ITER = 1000;
    const double EPS = 1e-15;

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        // Find the largest off-diagonal element in A
        int p = 0, q = 1;
        double max_offdiag = fabs(A[0][1]);

        for (int i = 0; i < 3; ++i) {
            for (int j = i + 1; j < 3; ++j) {
                if (fabs(A[i][j]) > max_offdiag) {
                    max_offdiag = fabs(A[i][j]);
                    p = i;
                    q = j;
                }
            }
        }

        // If the off-diagonal elements are small enough, we can stop
        if (max_offdiag < EPS) {
            break;
        }

        // Apply Jacobi rotation to eliminate A[p][q]
        ApplyJacobiRotation(A, eigenvectors, p, q);
    }

    // Extract the eigenvalues from the diagonal of A
    for (int i = 0; i < 3; ++i) {
        eigenvalues[i] = A[i][i];
    }

    // Sort eigenvalues and corresponding eigenvectors in ascending order
    for (int i = 0; i < 2; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            if (eigenvalues[i] > eigenvalues[j]) {
                // Swap eigenvalues
                Swap(eigenvalues[i], eigenvalues[j]);
                // Swap eigenvectors
                for (int k = 0; k < 3; ++k) {
                    Swap(eigenvectors[k][i], eigenvectors[k][j]);
                }
            }
        }
    }

    // Orthogonalize the eigenvectors using Gram-Schmidt to maintain accuracy
    GramSchmidtOrthogonalize(eigenvectors);

    // Normalize all eigenvectors to ensure unit length
    for (int i = 0; i < 3; ++i) {
        Vec_normalize3(eigenvectors[0][i], eigenvectors[1][i], eigenvectors[2][i],
                      eigenvectors[0][i], eigenvectors[1][i], eigenvectors[2][i]);
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

KOKKOS_INLINE_FUNCTION
void antisymmetrize_matrix(const double matrix[3][3],
						double sym_matrix[3][3]){
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			sym_matrix[i][j] = 0.0;
		}// j
	}// i

	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			sym_matrix[i][j] = 0.5*(matrix[i][j] - matrix[j][i]);
		}// j
	}// i

}// end antisymmetrix_matrix

#endif
