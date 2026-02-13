#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "utilities.h"
#include "Profiler.h"
#include "math_functions.h"
#include "reduction_data_structures.h"

void EVPFFT::Cauchy_to_PK1()
{
  Profiler profiler(__FUNCTION__);

  const size_t n = 9; // for sgavg (3,3)
  ArrayType <real_t, n> all_reduce;

  MatrixTypeRealDual sgavg(3,3);
 
  Kokkos::parallel_reduce(
    Kokkos::MDRangePolicy<Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER>>({1,1,1}, {npts3+1,npts2+1,npts1+1}),
    KOKKOS_CLASS_LAMBDA(const int k, const int j, const int i, ArrayType <real_t,n> & loc_reduce) {

    // averages packed in array
    int ic;
    ic = -1;

    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        ic = ic + 1;
        loc_reduce.array[ic] += sg(ii,jj,i,j,k) * wgtc(i,j,k);
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
      sgavg.host(ii,jj) = all_reduce.array[ic];
    }
  }
  sgavg.update_device();

  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {


    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        sgPK1(ii,jj,i,j,k) = 0.0;
        for (int kk = 1; kk <= 3; kk++) {
          sgPK1(ii,jj,i,j,k) += (sg(ii,kk,i,j,k) - sgavg(ii,kk))*defgradinv(jj,kk,i,j,k)*detF(i,j,k);
        }
      }
    }
  
  }); // end FOR_ALL_CLASS

}
