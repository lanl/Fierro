#include <matar.h>
#include <stdio.h>
#include <math.h>

using namespace mtr; // matar namespace

/*
  Functions or subroutines that will be called from fortran should include "_" at the end of name. 
  Example: Given a C++ definded subroutine as "subroutineName_", 
           it should be called in fortran as "call subroutineName(...)"
  Also the functions or subroutines should be decleared with extern "C"
*/

extern "C" void kokkos_initialize_();
extern "C" void kokkos_finalize_();

extern "C" void square_array_elements_(double *array, int *nx, int *ny);
extern "C" void sum_array_elements_(double *array, int *nx, int *ny, double *sum_of_elements);

void kokkos_initialize_()
{
    Kokkos::initialize();
}

void kokkos_finalize_()
{
    Kokkos::finalize();
}

void square_array_elements_(double *array, int *nx, int *ny)
{
    // define private copys of nx and ny
    // this enables kokkos to copy stack variables
    // if used in kokkos kernal
    int nx_ = *nx;
    int ny_ = *ny;

    // create DViewFMatrixKokkos since array is fortran allocated
    auto array_2D_dual_view = DViewFMatrixKokkos <double> (array, nx_, ny_);

    DO_ALL(j, 1, ny_,
           i, 1, nx_, {
        array_2D_dual_view(i,j) = pow(array_2D_dual_view(i,j), 2);
    });

    array_2D_dual_view.update_host();
}

void sum_array_elements_(double *array, int *nx, int *ny, double *sum_of_elements)
{
    // define private copys of nx and ny
    int nx_ = *nx;
    int ny_ = *ny;

    // create DViewFMatrixKokkos since array is fortran allocated
    auto array_2D_dual_view = DViewFMatrixKokkos <double> (array, nx_, ny_);

    double global_sum;
    double local_sum;
    DO_REDUCE_SUM(j, 1, ny_,
                  i, 1, nx_, 
                  local_sum, {
        local_sum += array_2D_dual_view(i,j);
    }, global_sum);

    // update sum_of_elements memory location
    *sum_of_elements = global_sum;
}
