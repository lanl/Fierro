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

// Normalize a 3D vector
KOKKOS_INLINE_FUNCTION
void normalize(double &x, double &y, double &z) {
    double norm = sqrt(x*x + y*y + z*z);
    if (norm > 1e-12) {
        x /= norm;
        y /= norm;
        z /= norm;
    }
}

// Normalize a 3D vector (array version)
KOKKOS_INLINE_FUNCTION
void normalize(double v[3]) {
    normalize(v[0], v[1], v[2]);
}

// Normalize a 3D vector and return the result
KOKKOS_INLINE_FUNCTION
void Vec_normalize3(double x, double y, double z, double &nx, double &ny, double &nz) {
    double norm = sqrt(x*x + y*y + z*z);
    if (norm > 1e-12) {
        nx = x / norm;
        ny = y / norm;
        nz = z / norm;
    } else {
        nx = ny = nz = 0.0;
    }
}

// Eigensystem of a 2x2 symmetric matrix
KOKKOS_INLINE_FUNCTION
void Eigensystem2S(const double &d12, double &d1, double &d2,
                   double &c, double &s) {
    const double sqrt_1_eps = sqrt(1.0 / std::numeric_limits<double>::epsilon());
    if (d12 == 0.0) {
        c = 1.0;
        s = 0.0;
    } else {
        // "The Symmetric Eigenvalue Problem", B. N. Parlett, pp.189-190
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

// Utility function used in KernelVector3G_aux
KOKKOS_INLINE_FUNCTION
bool KernelVector2G(const int &mode,
                    double &d1, double &d12, double &d21, double &d2) {
    // Find a vector (z1,z2) in the "near"-kernel of the matrix
    // |  d1  d12 |
    // | d21   d2 |
    // using QR factorization.
    // The vector (z1,z2) is returned in (d1,d2). Return 'true' if the matrix
    // is zero without setting (d1,d2).
    // Note: in the current implementation |z1| + |z2| = 1.

    // l1-norms of the columns
    double n1 = fabs(d1) + fabs(d21);
    double n2_col = fabs(d2) + fabs(d12);

    bool swap_columns = (n2_col > n1);
    double mu;

    if (!swap_columns) {
        if (n1 == 0.0) {
            return true;
        }

        if (mode == 0) { // eliminate the larger entry in the column
            if (fabs(d1) > fabs(d21)) {
                Swap(d1, d21);
                Swap(d12, d2);
            }
        } else { // eliminate the smaller entry in the column
            if (fabs(d1) < fabs(d21)) {
                Swap(d1, d21);
                Swap(d12, d2);
            }
        }
    } else {
        // n2_col > n1, swap columns 1 and 2
        if (mode == 0) { // eliminate the larger entry in the column
            if (fabs(d12) > fabs(d2)) {
                Swap(d1, d2);
                Swap(d12, d21);
            } else {
                Swap(d1, d12);
                Swap(d21, d2);
            }
        } else { // eliminate the smaller entry in the column
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
        // v = (n1, n2)^t,  |v| = 1
        // Q = I - 2 v v^t,  Q (d1, d21)^t = (mu, 0)^t
        mu = copysign(n1, d1);
        double temp_n1 = -d21 * (d21 / (d1 + mu)); // temp_n1 = d1 - mu
        d1 = mu;
        // Normalize (temp_n1, d21) to avoid overflow/underflow
        if (fabs(temp_n1) <= fabs(d21)) {
            // (temp_n1, d21) <-- (temp_n1 / d21, 1)
            double n1_over_d21 = temp_n1 / d21;
            mu = (2.0 / (1.0 + n1_over_d21 * n1_over_d21)) * (n1_over_d21 * d12 + d2);
            d2 = d2 - mu;
            d12 = d12 - mu * n1_over_d21;
        } else {
            // (temp_n1, d21) <-- (1, d21 / temp_n1)
            double d21_over_n1 = d21 / temp_n1;
            mu = (2.0 / (1.0 + d21_over_n1 * d21_over_n1)) * (d12 + d21_over_n1 * d2);
            d2 = d2 - mu * d21_over_n1;
            d12 = d12 - mu;
        }
    }

    // Solve:
    // | d1 d12 | | z1 | = | 0 |
    // |  0  d2 | | z2 |   | 0 |

    // Evaluate z2 at t = t1
    mu = -d12 / d1;
    double n2 = 1.0 / (1.0 + fabs(mu));

    // Check if |d1| <= |d2| * z2
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

    return false; // Indicates the matrix is not zero
}

// Implementation of KernelVector3G_aux
KOKKOS_INLINE_FUNCTION
int KernelVector3G_aux(const int &mode,
                       double &d1, double &d2, double &d3,
                       double &c12, double &c13, double &c23,
                       double &c21, double &c31, double &c32) {
    int kdim;
    double mu, n1, n2, n3, s1, s2, s3;

    s1 = hypot(c21, c31);
    n1 = hypot(d1, s1);

    if (s1 != 0.) {
        // v = (n1, c21, c31)
        mu = copysign(n1, d1);
        n1 = -s1*(s1/(d1 + mu)); // = d1 - mu
        d1 = mu;

        // Normalize (n1, c21, c31) to avoid overflow/underflow
        double max_n = std::max({fabs(n1), fabs(c21), fabs(c31)});
        s1 = n1 / max_n;
        s2 = c21 / max_n;
        s3 = c31 / max_n;
        mu = 2.0 / (s1*s1 + s2*s2 + s3*s3);

        // Compute n2 and n3
        n2  = mu * (c12 + s2*d2  + s3*c32);
        n3  = mu * (c13 + s2*c23 + s3*d3);

        // Update matrix entries
        c12 = c12 -    n2;
        d2  = d2  - s2*n2;
        c32 = c32 - s3*n2;
        c13 = c13 -    n3;
        c23 = c23 - s2*n3;
        d3  = d3  - s3*n3;
    }

    // Solve:
    // |  d2 c23 | | z2 | = | 0 |
    // | c32  d3 | | z3 |   | 0 |
    if (KernelVector2G(mode, d2, c23, c32, d3) == 2) {
        // Have two solutions:
        // two vectors in the kernel are P (-c12/d1, 1, 0)^t and
        // P (-c13/d1, 0, 1)^t where P is the permutation matrix swapping
        // entries 1 and col.

        // A vector orthogonal to both these vectors is P (1, c12/d1, c13/d1)^t
        d2 = c12/d1;
        d3 = c13/d1;
        d1 = 1.;
        kdim = 2;
    } else {
        // solve for z1:
        d1 = -(c12*d2 + c13*d3)/d1;
        kdim = 1;
    }

    Vec_normalize3(d1, d2, d3, d1, d2, d3);

    return kdim;
}

// Implement KernelVector3S
KOKKOS_INLINE_FUNCTION
int KernelVector3S(const int &mode, const double &d12,
                   const double &d13, const double &d23,
                   double &d1, double &d2, double &d3) {
    // Find a unit vector (z1,z2,z3) in the "near"-kernel of the matrix
    // |  d1  d12  d13 |
    // | d12   d2  d23 |
    // | d13  d23   d3 |
    // using QR factorization.
    // The vector (z1,z2,z3) is returned in (d1,d2,d3).
    // Returns the dimension of the kernel, kdim, but never zero.

    double c12 = d12, c13 = d13, c23 = d23;
    double c21, c31, c32;
    int col, row;

    // l1-norms of the columns:
    double n1 = fabs(d1) + fabs(c12) + fabs(c13);
    double n2 = fabs(d2) + fabs(c12) + fabs(c23);
    double n3 = fabs(d3) + fabs(c13) + fabs(c23);

    // Column pivoting: choose the column with the largest norm
    if (n1 >= n3) {
        col = (n1 >= n2) ? 1 : 2;
    } else {
        col = (n2 >= n3) ? 2 : 3;
    }

    switch (col) {
        case 1:
            if (n1 == 0.0) {
                return 3; // Zero matrix
            }
            break;
        case 2:
            if (n2 == 0.0) {
                return 3; // Zero matrix
            }
            Swap(c13, c23);
            Swap(d1, d2);
            break;
        case 3:
            if (n3 == 0.0) {
                return 3; // Zero matrix
            }
            Swap(c12, c23);
            Swap(d1, d3);
            break;
    }

    // Row pivoting depending on 'mode'
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
        double max_v = std::max({fabs(v1), fabs(z2), fabs(z3)});
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

        // Off-diagonal entries should be zero (for debugging)
        // d12 = d12 - v1 * w2 - v2 * w1; // Should be zero
        // d13 = d13 - v1 * w3 - v3 * w1; // Should be zero
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

// Main function to compute eigenvalues and eigenvectors
KOKKOS_INLINE_FUNCTION
void get_spectrum(const double data[3][3], double eigenvalues[3], double eigenvectors[3][3]) {
    // Initialize eigenvalues and eigenvectors to zero
    for (int i = 0; i < 3; ++i) {
        eigenvalues[i] = 0.0;
        for (int j = 0; j < 3; ++j) {
            eigenvectors[i][j] = 0.0;
        }
    }

    // Extract matrix elements
    double d11 = data[0][0];
    double d12 = data[0][1];
    double d13 = data[0][2];
    double d22 = data[1][1];
    double d23 = data[1][2];
    double d33 = data[2][2];

    double mult;
    {
        double d_max = fabs(d11);
        if (d_max < fabs(d22)) { d_max = fabs(d22); }
        if (d_max < fabs(d33)) { d_max = fabs(d33); }
        if (d_max < fabs(d12)) { d_max = fabs(d12); }
        if (d_max < fabs(d13)) { d_max = fabs(d13); }
        if (d_max < fabs(d23)) { d_max = fabs(d23); }

        mult = (d_max > 0.0) ? d_max : 1.0;
    }

    // Scale the matrix entries
    d11 /= mult;  d22 /= mult;  d33 /= mult;
    d12 /= mult;  d13 /= mult;  d23 /= mult;

    double aa = (d11 + d22 + d33)/3;  // aa = tr(A)/3
    double c1 = d11 - aa;
    double c2 = d22 - aa;
    double c3 = d33 - aa;

    double Q = (2*(d12*d12 + d13*d13 + d23*d23) + c1*c1 + c2*c2 + c3*c3)/6;
    double R = (c1*(d23*d23 - c2*c3)+ d12*(d12*c3 - 2*d13*d23) + d13*d13*c2)/2;

    const double epsilon = 1e-12;

    if (Q <= epsilon) {
        eigenvalues[0] = eigenvalues[1] = eigenvalues[2] = aa * mult;
        eigenvectors[0][0] = 1.0; eigenvectors[0][1] = 0.0; eigenvectors[0][2] = 0.0;
        eigenvectors[1][0] = 0.0; eigenvectors[1][1] = 1.0; eigenvectors[1][2] = 0.0;
        eigenvectors[2][0] = 0.0; eigenvectors[2][1] = 0.0; eigenvectors[2][2] = 1.0;
    } else {
        double sqrtQ = sqrt(Q);
        double sqrtQ3 = Q*sqrtQ;
        double r;
        if (fabs(R) >= sqrtQ3) {
            r = (R < 0.0) ? 2*sqrtQ : -2*sqrtQ;
        } else {
            R = R / sqrtQ3;
            double theta = acos(R);
            r = -2*sqrtQ*cos(theta/3);
        }

        aa += r;
        c1 = d11 - aa;
        c2 = d22 - aa;
        c3 = d33 - aa;

        const int mode = 0;

        // Find a unit vector z = (z1,z2,z3) in the "near"-kernel of
        //  |  c1  d12  d13 |
        //  | d12   c2  d23 | = A - aa*I
        //  | d13  d23   c3 |
        // This vector is also an eigenvector for A corresponding to aa.
        // The vector z overwrites (c1,c2,c3).
        int kdim = KernelVector3S(mode, d12, d13, d23, c1, c2, c3);

        if (kdim == 3) {
            // 'aa' is a triple eigenvalue
            eigenvalues[0] = eigenvalues[1] = eigenvalues[2] = aa * mult;
            eigenvectors[0][0] = 1.0; eigenvectors[0][1] = 0.0; eigenvectors[0][2] = 0.0;
            eigenvectors[1][0] = 0.0; eigenvectors[1][1] = 1.0; eigenvectors[1][2] = 0.0;
            eigenvectors[2][0] = 0.0; eigenvectors[2][1] = 0.0; eigenvectors[2][2] = 1.0;
            return;
        }

        // Using the eigenvector c=(c1,c2,c3) transform A into
        //                   | d11   0   0 |
        // A <-- Q P A P Q = |  0  d22 d23 |
        //                   |  0  d23 d33 |
        double v1, v2, v3, g;
        int k = Reduce3S(mode, d11, d22, d33, d12, d13, d23, c1, c2, c3, v1, v2, v3, g);

        // find the eigenvalues and eigenvectors for
        // | d22 d23 |
        // | d23 d33 |
        double c, s;
        Eigensystem2S(d23, d22, d33, c, s);

        // Assign the eigenvalues
        eigenvalues[0] = aa * mult;
        eigenvalues[1] = d22 * mult;
        eigenvalues[2] = d33 * mult;

        // Assign eigenvectors
        // First eigenvector is c = (c1, c2, c3)
        normalize(c1, c2, c3);
        eigenvectors[0][0] = c1;
        eigenvectors[0][1] = c2;
        eigenvectors[0][2] = c3;

        // Second and third eigenvectors
        double v[3];
        v[0] = v1;
        v[1] = v2;
        v[2] = v3;

        double evec2[3], evec3[3];
        evec2[0] = -v[0]*(g*(v[1]*c - v[2]*s));
        evec2[1] = c - v[1]*(g*(v[1]*c - v[2]*s));
        evec2[2] = -s - v[2]*(g*(v[1]*c - v[2]*s));

        evec3[0] = -v[0]*(g*(v[1]*s + v[2]*c));
        evec3[1] = s - v[1]*(g*(v[1]*s + v[2]*c));
        evec3[2] = c - v[2]*(g*(v[1]*s + v[2]*c));

        switch (k) {
            case 2:
                Swap(evec2[0], evec2[1]);
                Swap(evec3[0], evec3[1]);
                break;

            case 3:
                Swap(evec2[0], evec2[2]);
                Swap(evec3[0], evec3[2]);
        }

        // Normalize eigenvectors
        normalize(evec2);
        normalize(evec3);

        // Assign eigenvectors
        eigenvectors[1][0] = evec2[0];
        eigenvectors[1][1] = evec2[1];
        eigenvectors[1][2] = evec2[2];

        eigenvectors[2][0] = evec3[0];
        eigenvectors[2][1] = evec3[1];
        eigenvectors[2][2] = evec3[2];

        // Now, sort the eigenvalues and eigenvectors in ascending order
        for (int i = 0; i < 2; ++i) {
            int min_idx = i;
            for (int j = i + 1; j < 3; ++j) {
                if (eigenvalues[j] < eigenvalues[min_idx]) {
                    min_idx = j;
                }
            }
            if (min_idx != i) {
                // Swap eigenvalues
                Swap(eigenvalues[i], eigenvalues[min_idx]);
                // Swap eigenvectors
                for (int k = 0; k < 3; ++k) {
                    Swap(eigenvectors[i][k], eigenvectors[min_idx][k]);
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
