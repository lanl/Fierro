#include "linear_algebra.h"
#include "state.h"
#include "mesh.h"

KOKKOS_FUNCTION
void invert_matrix(CArrayKokkos <double> &mtx_inv,
                   CArrayKokkos <double> &matrix,
                   const mesh_t &mesh,
                   const int size){
    
    int singular = 0;
    int parity = 0;
    CArrayKokkos <double> col(size);
    CArrayKokkos <int> index(size);

    for (int i = 0; i < size; i++){
        col(i) = 0.0;
        index(i) = 0.0;
    }

    int permutation;

    for(int i = 0; i <  size; i++){
        for (int j =0; j < size; j++){
            mtx_inv(i,j) = 0.0;
        }
    }
   
    ViewCArrayKokkos <double> lu_mtx(&matrix(0,0), size, size);
    ViewCArrayKokkos <double> view_mtx_inv(&mtx_inv(0,0), size, size);

    lu_decomp(lu_mtx, index, parity, size);

    lu_invert(lu_mtx, view_mtx_inv, col, index, size);

}// end invert_matrix


