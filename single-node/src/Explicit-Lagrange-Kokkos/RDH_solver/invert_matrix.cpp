// #include "linear_algebra.h"
// #include "state.h"
// #include "mesh.h"
// #define TINY 1.e-16
// KOKKOS_FUNCTION
// void invert_matrix( ViewCArrayKokkos <double> &mtx_inv,
//                     ViewCArrayKokkos <double> &matrix,
//                    const mesh_t &mesh,
//                    const int size){
    
//     int singular = 0;
//     int parity = 0;
//     CArrayKokkos <double> col(size);
//     CArrayKokkos <int> index(size);

//     for (int i = 0; i < size; i++){
//         col(i) = 0.0;
//         index(i) = 0.0;
//     }

//     int permutation;

//     for(int i = 0; i <  size; i++){
//         for (int j =0; j < size; j++){
//             mtx_inv(i,j) = 0.0;
//         }
//     }
   
//     ViewCArrayKokkos <double> lu_mtx(&matrix(0,0), size, size);
//     ViewCArrayKokkos <double> view_mtx_inv(&mtx_inv(0,0), size, size);

//     lu_decomp(lu_mtx, index, parity, size);

//     lu_invert(lu_mtx, view_mtx_inv, col, index, size);

// }// end invert_matrix


// #include <Kokkos_Core.hpp>

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


