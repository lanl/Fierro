// // #include <iostream>
// // #include <stdio.h>
// // #include <iomanip> // Include this header for std::fixed and std::setprecision
// // #include <stdlib.h>
// // #include <Kokkos_Core.hpp>
// // #include <sys/stat.h>

// #include "ref_elem.h"
// #include "mesh.h"
// #include "state.h"

// #include "matar.h"

// using namespace mtr;

// #define TINY 1.e-16

// // void print_mat(DCArrayKokkos<double> &mat){


// //     mat.update_host();

// //     std::cout << std::fixed; // Use fixed-point notation
// //     std::cout << std::setprecision(4); // Set the precision to 2 decimal places

// //     std::cout<<"**************************************************"<<std::endl;

// //     // ************ WARNING: CUDA MAY NOT LIKE THIS  ************ 
// //     RUN({

// //         for(int i = 0; i < mat.dims(0); i++){
// //             std::cout<<std::endl;
// //             for(int j = 0; j < mat.dims(1); j++){
// //                 std::cout<<mat.host(i,j)<<"  ";
// //             }
// //             // std::cout<<std::endl;
// //         }
// //         std::cout<<std::endl;
// //     });  // end RUN
// //     std::cout<<std::endl;
// //     std::cout<<"**************************************************"<<std::endl;
// //     std::cout<<std::endl;
// // }

// KOKKOS_FUNCTION
// int lu_decomp(
//     ViewCArrayKokkos<double> &source_mat,
//     CArrayKokkos<int> &indx,
//     int &parity, 
//     const int n)
// {


//     CArrayKokkos<double> vv(n);// = DCArrayKokkos<double>(n);

//     parity = 1;
//     /* search for the largest element in each row; save the scaling in the
//     temporary array vv and return zero if the matrix is singular */

//     // int i, imax, k; // Indexing
//     // double big, sum, temp; // useful storage

//     CArrayKokkos <double> extra_double_var(2);  // big1 = extra_double_var(0), temp1 = extra_double_var(1)
//     CArrayKokkos <int> extra_int_var(2);

    
//     extra_int_var(0) = 1; // singular
//     extra_int_var(1) = 1; // parity
    
//     for(int i = 0; i < n; i++){
//         double *big1 = &extra_double_var(0);
//         double *temp1 = &extra_double_var(1);

//         *big1 = 0.;
//         for(int j = 0; j < n; j++) if((*temp1=fabs(source_mat(i,j))) > *big1) *big1=*temp1;
        
//         if (*big1 == 0.) {
//             extra_int_var(0) = 0;
//         }
//         vv(i) = *big1;
//     }

//     // Return if singular

//     // if (extra_int_var.host(0) == 0) return(0);
//     //extra_int_var.update_host();
//     if(extra_int_var(0) == 0) return(0);//.host(0) == 0) return(0);

//     /* the main loop for the Crout's algorithm */
//     CArrayKokkos <int> int_vars(0);
//     CArrayKokkos <double> double_vars(2);
//     int *imax = &int_vars(0);

//     double *big = &double_vars(0);
//     double *temp = &double_vars(1);

//     for (int j = 0; j < n; j++) {
//         /* this is the part a) of the algorithm except for i==j */
//         for (int i = 0; i < j; i++) {
//             double sum = source_mat(i, j);

//             for (int k = 0; k < i; k++) {
//                 sum -= source_mat(i, k) * source_mat(k, j);
//             }
//             source_mat(i, j) = sum;
//         }

//         /* initialize for the search for the largest pivot element */
//         *big  = 0.;
//         *imax = j;

//         /* this is the part a) for i==j and part b) for i>j + pivot search */
//         for (int i = j; i < n; i++) {
//             double sum = source_mat(i, j);

//             for (int k = 0; k < j; k++) {
//                 sum -= source_mat(i, k) * source_mat(k, j);
//             }

//             source_mat(i, j) = sum;

//             /* is the figure of merit for the pivot better than the best so far? */
//             if ((*temp = vv(i) * fabs(sum)) >= *big) {
//                 *big = *temp; *imax = i;
//             }
//         }

//         /* interchange rows, if needed, change parity and the scale factor */
//         if (*imax != j) {
//             for (int k = 0; k < n; k++) {
//                 *temp = source_mat(*imax, k);
//                 source_mat(*imax, k) = source_mat(j, k);
//                 source_mat(j, k)    = *temp;
//             }

//             extra_int_var(1)   = -(extra_int_var(1));
//             vv(*imax) = vv(j);
//         }

//         /* store the index */
//         indx(j) = *imax;
//         /* if the pivot element is zero, the matrix is singular but for some
//         applications a tiny number is desirable instead */

//         if (source_mat(j, j) == 0.) {
//             source_mat(j, j) = TINY;
//         }
//         /* finally, divide by the pivot element */

//         if (j < n - 1) {
//             *temp = 1. / source_mat(j, j);
//             for (int i = j + 1; i < n; i++) {
//                 source_mat(i, j) *= *temp;
//             }
//         }
//     }

//     //extra_int_var.update_host();
//     parity = extra_int_var(1);

//     return(extra_int_var(0));
// }


// KOKKOS_FUNCTION
// void lu_backsub(
//     ViewCArrayKokkos<double> &mat, 
//     CArrayKokkos<int> &indx,
//     CArrayKokkos<double> &col_vec,
//     const int not_used)
// {
//     // int j, ip, ii = -1;
//     // double sum;
    
//     CArrayKokkos <int> int_vars(3);
//     CArrayKokkos <double> double_vars(1);

//     int_vars(0) = -1;
//     int_vars(1) = -1;
//     int_vars(2) = -1;
//     double_vars(0) = 0.0;

//     /* First step of backsubstitution; the only wrinkle is to unscramble
//     the permutation order. Note: the algorithm is optimized for a
//     possibility of large amount of zeroes in b */

//     int n = mat.dims(0);
//     int *j = &int_vars(0);
//     int *ip = &int_vars(1);
//     int *ii = &int_vars(2);

//     double *sum = &double_vars(0);

//     for (int i = 0; i < n; i++) {
//         *ip = indx(i);

//         *sum = col_vec(*ip);
//         col_vec(*ip) = col_vec(i);

//         if (*ii >= 0) {
//             for (*j = *ii; *j < i; (*j)++) {
//                 *sum -= mat(i, *j) * col_vec(*j);
//             }
//         }
//         else if (*sum) {
//             *ii = i;             /* a nonzero element encounted */
//         }
//         col_vec(i) = *sum;
//     }

//     /* the second step */
//     for (int i = n - 1; i >= 0; i--) {
//         *sum = col_vec(i);

//         // std::cout<<"i = "<<i<<std::endl;
//         // std::cout<<"j = "<<*j<<std::endl;

//         for (*j = i + 1; *j < n; (*j)++) {

//             *sum -= mat(i, *j) * col_vec(*j);
//         }

//         col_vec(i) = *sum / mat(i, i);
//     }

// }

// KOKKOS_FUNCTION
// void lu_invert(
//     ViewCArrayKokkos<double> &mat, 
//     ViewCArrayKokkos<double> &mat_inv, 
//     CArrayKokkos<double> &col_vec,
//     CArrayKokkos<int> &indx,
//     int n)
// {
//     //int n = mat.dims(0);

//     for (int j = 0; j < n; j++) {
//         // Initialize col to zero
//         for(int i = 0; i < n; i++) {
//             col_vec(i) = 0.0;
//         }

//         col_vec(j) = 1.0;
//         lu_backsub(mat, indx, col_vec, n);

//         for(int i = 0; i < n; i++) {
//             mat_inv(i, j) = col_vec(i);
//         }
//     }
// }


// int main(int argc, char* argv[])
// {

//     Kokkos::initialize();
//     {

//         std::cout << "**** Testing Linear Algebra **** " << std::endl;


//         int matrix_size = 8;


//         // Create and initialize matrix
//         DCArrayKokkos<double> MAT = DCArrayKokkos<double>(matrix_size, matrix_size, "input_matrix");
//         DCArrayKokkos<double> MAT_LU = DCArrayKokkos<double>(matrix_size, matrix_size, "input_matrix_LU_decomp");
//         DCArrayKokkos<double> MAT_INV = DCArrayKokkos<double>(matrix_size, matrix_size, "input_matrix_inverse");

//         FOR_ALL (i, 0, MAT.dims(0),
//                  j, 0, MAT.dims(1),{
            
//             if(i > j){
//                 MAT(i,j) = (double)i - (double)j;
//                 MAT_LU(i,j) = (double)i - (double)j;
//                 MAT_INV(i,j) = (double)i - (double)j;
//             }
//             if(i < j){
//                 MAT(i,j) = (double)i + (double)j + 1.0;
//                 MAT_LU(i,j) = (double)i + (double)j + 1.0;
//                 MAT_INV(i,j) = (double)i + (double)j + 1.0;
//             }

//             if(i == j){
//                 MAT(i,j) = (double)i + (double)j + 1.0;
//                 MAT_LU(i,j) = (double)i + (double)j + 1.0;
//                 MAT_INV(i,j) = (double)i + (double)j + 1.0;
//             }
//         });



//         // Print matrix
//         std::cout<<"Printing MAT"<<std::endl;

//         // ************ WARNING: CUDA MAY NOT LIKE THIS  ************ 
//         print_mat(MAT);

//         // Needed for LU decomp
//         DCArrayKokkos<double> col_vec = DCArrayKokkos<double>(matrix_size, "least_squares_coeffs");
//         DCArrayKokkos<int> indx(matrix_size, "permutations");
//         int parity = 0;
        
//         FOR_ALL (i, 0, matrix_size, {
//             indx(i) = 0;
//             col_vec(i) = 0.0;
//         });

//         int singular = 1;
//         singular = lu_decomp_test(MAT_LU, indx, parity, matrix_size); 

//         if(singular == 0) std::cout<<"WARNING: SINGULAR MATRIX"<<std::endl;


//         // std::cout<<"Printing MAT_LU"<<std::endl;
//         // print_mat(MAT_LU);

//         lu_invert_test(MAT_LU, MAT_INV, col_vec, indx);

//         std::cout<<"Printing MAT_INV"<<std::endl;
//         print_mat(MAT_INV);


//         DCArrayKokkos<double> MAT_TEST = DCArrayKokkos<double>(matrix_size, matrix_size, "test_matrix");
//         RUN({
//             int n = MAT_TEST.dims(0);
//             for(int i=0; i<n; i++){
//                 for(int j=0; j<n; j++){
//                     MAT_TEST(i,j) = 0.0;
//                 }
//             }

        
//             for(int i=0; i<n; i++){
//                 for(int j=0; j<n; j++){
//                     for(int k=0; k<n; k++){
//                         MAT_TEST(i,k) += MAT(i,j)*MAT_INV(j,k);
//                     }
//                 }
//             }
//         });

//         std::cout<<"Printing MAT_TEST: MAT*MAT_INV"<<std::endl;
//         // ************ WARNING: CUDA MAY NOT LIKE THIS  ************ 
//         print_mat(MAT_TEST);

//     }
//     Kokkos::finalize();

    
//     return 0;
// }

#include "state.h"
#include "mesh.h"
#define TINY 1.e-16
#include <Kokkos_Core.hpp>

// KOKKOS_INLINE_FUNCTION
// void invert_matrix(Kokkos::View<double**> mtx, Kokkos::View<double**> mtx_inv, int size);
// KOKKOS_INLINE_FUNCTION
// int lu_decomp(Kokkos::View<double**> source_mat, Kokkos::View<int*> indx, int &parity, const int n);
// KOKKOS_INLINE_FUNCTION
// void lu_invert(Kokkos::View<double**> lu_mtx, Kokkos::View<double**> mtx_inv, Kokkos::View<double*> col, Kokkos::View<int*> indx, int n);

// // / Matrix inversion function using LU decomposition
// KOKKOS_INLINE_FUNCTION
// void invert_matrix(Kokkos::View<double**> mtx, Kokkos::View<double**> mtx_inv, int size) {
//     Kokkos::View<double*> col("col", size);
//     Kokkos::View<int*> index("index", size);
//     Kokkos::View<double**> lu_mtx("lu_mtx", size, size);

//     int parity;

//     // Copy mtx to lu_mtx for LU decomposition
//     for (int i = 0; i < size; i++) {
//         for (int j = 0; j < size; j++) {
//             lu_mtx(i, j) = mtx(i, j);
//         }
//     }

//     // LU decomposition
//     lu_decomp(lu_mtx, index, parity, size);

//     // Invert matrix using LU decomposition
//     for (int j = 0; j < size; j++) {
//         // Set up the column of the identity matrix
//         for (int i = 0; i < size; i++) {
//             col(i) = (i == j) ? 1.0 : 0.0;
//         }
        
//         // Solve for each column of the inverse matrix
//         lu_invert(lu_mtx, mtx_inv, col, index, size);
//     }
// }

// // LU Decomposition function
// KOKKOS_INLINE_FUNCTION
// int lu_decomp(Kokkos::View<double**> source_mat, Kokkos::View<int*> indx, int &parity, const int n) {
//     Kokkos::View<double*> vv("vv", n);
//     parity = 1;

//     // Search for the largest element in each row
//     for (int i = 0; i < n; i++) {
//         double big = 0.0;
//         for (int j = 0; j < n; j++) {
//             double temp = fabs(source_mat(i, j));
//             if (temp > big) {
//                 big = temp;
//             }
//         }
//         vv(i) = 1.0 / big;
//     }

//     // Perform Gaussian elimination with partial pivoting
//     for (int k = 0; k < n; k++) {
//         double big = 0.0;
//         int imax = k;
//         for (int i = k; i < n; i++) {
//             double temp = vv(i) * fabs(source_mat(i, k));
//             if (temp > big) {
//                 big = temp;
//                 imax = i;
//             }
//         }

//         if (k != imax) {
//             for (int j = 0; j < n; j++) {
//                 double temp = source_mat(imax, j);
//                 source_mat(imax, j) = source_mat(k, j);
//                 source_mat(k, j) = temp;
//             }
//             parity = -parity;
//             vv(imax) = vv(k);
//         }

//         indx(k) = imax;
//         if (source_mat(k, k) == 0.0) {
//             source_mat(k, k) = TINY;
//         }

//         double diag_inv = 1.0 / source_mat(k, k);
//         for (int i = k + 1; i < n; i++) {
//             double factor = source_mat(i, k) * diag_inv;
//             source_mat(i, k) = factor;
//             for (int j = k + 1; j < n; j++) {
//                 source_mat(i, j) -= factor * source_mat(k, j);
//             }
//         }
//     }

//     return 0; // Assuming successful decomposition
// }

// // LU Inversion function
// KOKKOS_INLINE_FUNCTION
// void lu_invert(Kokkos::View<double**> lu_mtx, Kokkos::View<double**> mtx_inv, Kokkos::View<double*> col, Kokkos::View<int*> indx, int n) {
//     Kokkos::View<double*> x("x", n);

//     // Forward substitution (Ly = b)
//     for (int i = 0; i < n; i++) {
//         int ip = indx(i);
//         double sum = col(ip);
//         col(ip) = col(i);
//         for (int j = 0; j < i; j++) {
//             sum -= lu_mtx(i, j) * x(j);
//         }
//         x(i) = sum;
//     }

//     // Backward substitution (Ux = y)
//     for (int i = n - 1; i >= 0; i--) {
//         double sum = x(i);
//         for (int j = i + 1; j < n; j++) {
//             sum -= lu_mtx(i, j) * x(j);
//         }
//         x(i) = sum / lu_mtx(i, i);
//     }

//     // Copy solution into the appropriate column of the inverse matrix
//     for (int i = 0; i < n; i++) {
//         mtx_inv(i, col.extent(0) - 1) = x(i);
//     }
// }
// KOKKOS_INLINE_FUNCTION
// void invert_matrix(Kokkos::View<double**> mtx, Kokkos::View<double**> mtx_inv, int size);
// KOKKOS_INLINE_FUNCTION
// int lu_decomp(double source_mat[50][50], int indx[50], int &parity, const int n);
// KOKKOS_INLINE_FUNCTION
// void lu_invert(double lu_mtx[50][50], Kokkos::View<double**> mtx_inv, double col[50], int indx[50], int n, int col_index);


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

// The following spectral routines only work for 3D
KOKKOS_INLINE_FUNCTION
void min_eigenvalue_eigenvector(ViewCArrayKokkos<double> &matrix, double &min_eigenvalue, double eigenvector[3]) {
    // Inverse power iteration method to find the smallest eigenvalue and eigenvector
    const int max_iter = 100;
    const double tolerance = 1e-6;
    double b_k[3] = {1.0, 1.0, 1.0}; // Initial guess for eigenvector
    double b_k1[3];

    // Normalize initial vector
    double norm = sqrt(b_k[0] * b_k[0] + b_k[1] * b_k[1] + b_k[2] * b_k[2]);
    for (int i = 0; i < 3; i++) {
        b_k[i] /= norm;
    }

    for (int iter = 0; iter < max_iter; iter++) {
        // Solve system: b_k1 = inverse(matrix) * b_k
        double temp_vec[3] = {0.0, 0.0, 0.0};
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                temp_vec[i] += matrix(i, j) * b_k[j];
            }
        }
        for (int i = 0; i < 3; i++) {
            b_k1[i] = temp_vec[i];
        }

        // Calculate the norm of b_k1
        norm = sqrt(b_k1[0] * b_k1[0] + b_k1[1] * b_k1[1] + b_k1[2] * b_k1[2]);

        // Normalize b_k1 to get the next iteration of b_k
        for (int i = 0; i < 3; i++) {
            b_k1[i] /= norm;
        }

        // Check for convergence
        double diff = 0.0;
        for (int i = 0; i < 3; i++) {
            diff += fabs(b_k1[i] - b_k[i]);
        }
        if (diff < tolerance) {
            break;
        }

        // Update b_k for the next iteration
        for (int i = 0; i < 3; i++) {
            b_k[i] = b_k1[i];
        }
    }

    // Rayleigh quotient to find the corresponding eigenvalue
    min_eigenvalue = 0.0;
    for (int i = 0; i < 3; i++) {
        double temp = 0.0;
        for (int j = 0; j < 3; j++) {
            temp += matrix(i, j) * b_k[j];
        }
        min_eigenvalue += b_k[i] * temp;
    }

    // Copy eigenvector
    for (int i = 0; i < 3; i++) {
        eigenvector[i] = b_k[i];
    }
}

KOKKOS_INLINE_FUNCTION
void max_eigenvalue_eigenvector(CArrayKokkos<double> &matrix, double &max_eigenvalue, double eigenvector[3]) {
    // Power iteration method to find the dominant eigenvalue and eigenvector
    const int max_iter = 100;
    const double tolerance = 1e-6;
    double b_k[3] = {1.0, 1.0, 1.0}; // Initial guess for eigenvector
    double b_k1[3];

    // Normalize initial vector
    double norm = sqrt(b_k[0] * b_k[0] + b_k[1] * b_k[1] + b_k[2] * b_k[2]);
    for (int i = 0; i < 3; i++) {
        b_k[i] /= norm;
    }

    for (int iter = 0; iter < max_iter; iter++) {
        // Multiply matrix by current vector: b_k1 = matrix * b_k
        for (int i = 0; i < 3; i++) {
            b_k1[i] = 0.0;
            for (int j = 0; j < 3; j++) {
                b_k1[i] += matrix(i, j) * b_k[j];
            }
        }

        // Calculate the norm of b_k1
        norm = sqrt(b_k1[0] * b_k1[0] + b_k1[1] * b_k1[1] + b_k1[2] * b_k1[2]);

        // Normalize b_k1 to get the next iteration of b_k
        for (int i = 0; i < 3; i++) {
            b_k1[i] /= norm;
        }

        // Check for convergence
        double diff = 0.0;
        for (int i = 0; i < 3; i++) {
            diff += fabs(b_k1[i] - b_k[i]);
        }
        if (diff < tolerance) {
            break;
        }

        // Update b_k for the next iteration
        for (int i = 0; i < 3; i++) {
            b_k[i] = b_k1[i];
        }
    }

    // Rayleigh quotient to find the corresponding eigenvalue
    max_eigenvalue = 0.0;
    for (int i = 0; i < 3; i++) {
        double temp = 0.0;
        for (int j = 0; j < 3; j++) {
            temp += matrix(i, j) * b_k[j];
        }
        max_eigenvalue += b_k[i] * temp;
    }

    // Copy eigenvector
    for (int i = 0; i < 3; i++) {
        eigenvector[i] = b_k[i];
    }
}

// #include <iostream>
// #include <cmath>

// void test_min_eigenvalue_eigenvector() {
//     // Test 1: Diagonal Matrix
//     CArrayKokkos<double> matrix_1(3, 3);
//     matrix_1(0, 0) = 1.0;
//     matrix_1(0, 1) = 0.0;
//     matrix_1(0, 2) = 0.0;
//     matrix_1(1, 0) = 0.0;
//     matrix_1(1, 1) = 2.0;
//     matrix_1(1, 2) = 0.0;
//     matrix_1(2, 0) = 0.0;
//     matrix_1(2, 1) = 0.0;
//     matrix_1(2, 2) = 3.0;

//     double min_eigenvalue_1;
//     double eigenvector_1[3];

//     min_eigenvalue_eigenvector(matrix_1, min_eigenvalue_1, eigenvector_1);

//     std::cout << "Test 1 (Diagonal Matrix) - Minimum Eigenvalue: " << min_eigenvalue_1 << std::endl;
//     std::cout << "Eigenvector: [" << eigenvector_1[0] << ", " << eigenvector_1[1] << ", " << eigenvector_1[2] << "]" << std::endl;

//     // Test 2: Symmetric Matrix
//     CArrayKokkos<double> matrix_2(3, 3);
//     matrix_2(0, 0) = 4.0;
//     matrix_2(0, 1) = 1.0;
//     matrix_2(0, 2) = 1.0;
//     matrix_2(1, 0) = 1.0;
//     matrix_2(1, 1) = 3.0;
//     matrix_2(1, 2) = 0.0;
//     matrix_2(2, 0) = 1.0;
//     matrix_2(2, 1) = 0.0;
//     matrix_2(2, 2) = 2.0;

//     double min_eigenvalue_2;
//     double eigenvector_2[3];

//     min_eigenvalue_eigenvector(matrix_2, min_eigenvalue_2, eigenvector_2);

//     std::cout << "Test 2 (Symmetric Matrix) - Minimum Eigenvalue: " << min_eigenvalue_2 << std::endl;
//     std::cout << "Eigenvector: [" << eigenvector_2[0] << ", " << eigenvector_2[1] << ", " << eigenvector_2[2] << "]" << std::endl;

//     // Test 3: Non-Symmetric Matrix
//     CArrayKokkos<double> matrix_3(3, 3);
//     matrix_3(0, 0) = 2.0;
//     matrix_3(0, 1) = -1.0;
//     matrix_3(0, 2) = 0.0;
//     matrix_3(1, 0) = -1.0;
//     matrix_3(1, 1) = 2.0;
//     matrix_3(1, 2) = -1.0;
//     matrix_3(2, 0) = 0.0;
//     matrix_3(2, 1) = -1.0;
//     matrix_3(2, 2) = 2.0;

//     double min_eigenvalue_3;
//     double eigenvector_3[3];

//     min_eigenvalue_eigenvector(matrix_3, min_eigenvalue_3, eigenvector_3);

//     std::cout << "Test 3 (Non-Symmetric Matrix) - Minimum Eigenvalue: " << min_eigenvalue_3 << std::endl;
//     std::cout << "Eigenvector: [" << eigenvector_3[0] << ", " << eigenvector_3[1] << ", " << eigenvector_3[2] << "]" << std::endl;

//     // Test 4: Identity Matrix
//     CArrayKokkos<double> matrix_4(3, 3);
//     matrix_4(0, 0) = 1.0;
//     matrix_4(0, 1) = 0.0;
//     matrix_4(0, 2) = 0.0;
//     matrix_4(1, 0) = 0.0;
//     matrix_4(1, 1) = 1.0;
//     matrix_4(1, 2) = 0.0;
//     matrix_4(2, 0) = 0.0;
//     matrix_4(2, 1) = 0.0;
//     matrix_4(2, 2) = 1.0;

//     double min_eigenvalue_4;
//     double eigenvector_4[3];

//     min_eigenvalue_eigenvector(matrix_4, min_eigenvalue_4, eigenvector_4);

//     std::cout << "Test 4 (Identity Matrix) - Minimum Eigenvalue: " << min_eigenvalue_4 << std::endl;
//     std::cout << "Eigenvector: [" << eigenvector_4[0] << ", " << eigenvector_4[1] << ", " << eigenvector_4[2] << "]" << std::endl;
// }