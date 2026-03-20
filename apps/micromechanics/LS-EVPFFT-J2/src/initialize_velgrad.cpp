#include "evpfft.h"
#include "utilities.h"
#include "Profiler.h"
#include "reduction_data_structures.h"

void EVPFFT::initialize_velgrad()
{
  Profiler profiler(__FUNCTION__);

  const size_t n = 9; // for dvelgradavg (3,3)
  ArrayType <real_t, n> all_reduce;

  real_t dvelgradavg_[3*3];

  // create views of thread private arrays
  ViewMatrixTypeReal dvelgradavg(dvelgradavg_,3,3);
 
  // average of correction in current
  Kokkos::parallel_reduce(
    Kokkos::MDRangePolicy<Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER>>({1,1,1}, {npts3+1,npts2+1,npts1+1}),
    KOKKOS_CLASS_LAMBDA(const int k, const int j, const int i, ArrayType <real_t,n> & loc_reduce) {

    // averages packed in array
    int ic;
    real_t dum;

    ic = -1;

    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        dum = 0.0;
        for (int kk = 1; kk <= 3; kk++) {
          dum += work(ii,kk,i,j,k)*defgradinv(kk,jj,i,j,k);
        }
        ic = ic + 1;
        loc_reduce.array[ic] += dum * wgtc(i,j,k);
      }
    }

  }, all_reduce);
  Kokkos::fence(); // needed to prevent race condition

  MPI_Allreduce(MPI_IN_PLACE, all_reduce.array, all_reduce.size, MPI_REAL_T, MPI_SUM, mpi_comm);

  // averages unpacked from array
  int ic;
  ic = -1;
  
  for (int ii = 1; ii <= 3; ii++) {
    for (int jj = 1; jj <= 3; jj++) {
      ic = ic + 1;
      dvelgradavg(ii,jj) = all_reduce.array[ic];
    }
  }

  MatrixTypeRealDual dudot_dX_avg(3,3);
  
  for (int ii = 1; ii <= 3; ii++) {
    for (int jj = 1; jj <= 3; jj++) {
      dudot_dX_avg.host(ii,jj) = 0.0;
      for (int kk = 1; kk <= 3; kk++) {
        dudot_dX_avg.host(ii,jj) += (dvelgradmacro(ii,kk) - dvelgradavg(ii,kk))*defgradinvavgc_inv(kk,jj);
      }
    }
  }

  // update device
  dudot_dX_avg.update_device();

  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    real_t dum;

    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        dum = 0.0;
        for (int kk = 1; kk <= 3; kk++) {
          dum += (work(ii,kk,i,j,k) + dudot_dX_avg(ii,kk))*defgradinv(kk,jj,i,j,k);
        }
        velgrad(ii,jj,i,j,k) += dum;
      }
    }
  }); // end FOR_ALL_CLASS

}
