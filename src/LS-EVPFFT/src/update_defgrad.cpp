#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "utilities.h"
#include "Profiler.h"
#include "matrix_exp.h"
#include "determinant33.h"
#include "math_functions.h"
#include "reduction_data_structures.h"

void EVPFFT::update_defgrad()
{
  Profiler profiler(__FUNCTION__);

  const size_t n = 1 + 9 + 9 + 9; // for detFavg (1), defgradinvavgc (9), velgradavg (9), and defgradavg (9)
  ArrayType <real_t, n> all_reduce;

  real_t detFavg;
  real_t velgradavg_[3*3];
  real_t defgradinvavgc_[3*3];

  // create views of thread private arrays
  ViewMatrixTypeReal velgradavg(velgradavg_,3,3);
  ViewMatrixTypeReal defgradinvavgc(defgradinvavgc_,3,3);

  Kokkos::parallel_reduce(
    Kokkos::MDRangePolicy<Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER>>({1,1,1}, {npts3+1,npts2+1,npts1+1}),
    KOKKOS_CLASS_LAMBDA(const int k, const int j, const int i, ArrayType <real_t,n> & loc_reduce) {

    real_t detFtmp;

    // thread private arrays
    real_t defgradold_[3*3];
    real_t defgradnew_[3*3];
    real_t defgradinvold_[3*3];
    real_t disgradinc_[3*3];
    real_t expdisgradinc_[3*3];
    real_t defgradinvtmp_[3*3];
    real_t velgradold_[3*3];
    real_t defgradincinv_[3*3];
    
    // create views of thread private arrays
    ViewMatrixTypeReal defgradold(defgradold_,3,3);
    ViewMatrixTypeReal defgradnew(defgradnew_,3,3);
    ViewMatrixTypeReal defgradinvold(defgradinvold_,3,3);
    ViewMatrixTypeReal disgradinc(disgradinc_,3,3);
    ViewMatrixTypeReal expdisgradinc(expdisgradinc_,3,3);
    ViewMatrixTypeReal defgradinvtmp(defgradinvtmp_,3,3);
    ViewMatrixTypeReal velgradold(velgradold_,3,3);
    ViewMatrixTypeReal defgradincinv(defgradincinv_,3,3);

    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        defgradold(ii,jj) = defgrad(ii,jj,i,j,k);
        defgradinvold(ii,jj) = defgradinv(ii,jj,i,j,k);
        disgradinc(ii,jj) = velgrad(ii,jj,i,j,k)*tdot;
      }
    }

    // matrix exponential of increment
    matrix_exp(disgradinc.pointer(), expdisgradinc.pointer());

    // update 
    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        defgradnew(ii,jj) = 0.0;
        for (int kk = 1; kk <= 3; kk++) {
          defgradnew(ii,jj) += expdisgradinc(ii,kk)*defgrad(kk,jj,i,j,k);
        }
      }
    }

    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        defgrad(ii,jj,i,j,k) = defgradnew(ii,jj);
      }
    }

    determinant33(defgradnew.pointer(), detFtmp);

    if (detFtmp < 0.0) {
      printf(" -> WARNING: detF = %E24.14E in voxel %d %d %d\n", detFtmp, i, j, k);
      printf("elem = %d\n",elem_id);
    }

    detF(i,j,k) = detFtmp;

    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        defgradinvtmp(ii,jj) = defgradnew(ii,jj);
      }
    }

    invert_matrix <3> (defgradinvtmp.pointer());

    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        defgradinv(ii,jj,i,j,k) = defgradinvtmp(ii,jj);
      }
    }

    wgtc(i,j,k) = wgt*detFtmp;

    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        disgrad(ii,jj,i,j,k) = - defgradinvtmp(ii,jj);
      }
      disgrad(jj,jj,i,j,k) = disgrad(jj,jj,i,j,k) + 1.0;
    }

    // update velgrad configuration
    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        velgradold(ii,jj) = velgrad(ii,jj,i,j,k);
        defgradincinv(ii,jj) = expdisgradinc(ii,jj);
      }
    }

    invert_matrix <3> (defgradincinv.pointer());

    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        velgrad(ii,jj,i,j,k) = 0.0;
        for (int kk = 1; kk <= 3; kk++) {
          velgrad(ii,jj,i,j,k) += velgradold(ii,kk)*defgradincinv(kk,jj);
        }
      }
    }

    // averages packed in array
    int ic;
    ic = -1;
    
    ic = ic + 1;
    loc_reduce.array[ic] += detFtmp*wgt;

    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        ic = ic + 1;
        loc_reduce.array[ic] += defgradinvtmp(ii,jj) * wgt * detFtmp;
      }
    }

    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        ic = ic + 1;
        loc_reduce.array[ic] += velgrad(ii,jj,i,j,k) * wgt * detFtmp;
      }
    }

    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        ic = ic + 1;
        loc_reduce.array[ic] += defgrad(ii,jj,i,j,k) * wgt;
      }
    }
  }, all_reduce);
  Kokkos::fence(); // needed to prevent race condition

  MPI_Allreduce(MPI_IN_PLACE, all_reduce.array, all_reduce.size, MPI_REAL_T, MPI_SUM, mpi_comm);

  // averages unpacked from array
  int ic;
  ic = -1;
  
  ic = ic + 1;
  detFavg = all_reduce.array[ic];

  for (int ii = 1; ii <= 3; ii++) {
    for (int jj = 1; jj <= 3; jj++) {
      ic = ic + 1;
      defgradinvavgc(ii,jj) = all_reduce.array[ic]/detFavg;
    }
  }

  for (int ii = 1; ii <= 3; ii++) {
    for (int jj = 1; jj <= 3; jj++) {
      ic = ic + 1;
      velgradavg(ii,jj) = all_reduce.array[ic]/detFavg;
    }
  }

  for (int ii = 1; ii <= 3; ii++) {
    for (int jj = 1; jj <= 3; jj++) {
      ic = ic + 1;
      defgradavg(ii,jj) = all_reduce.array[ic];
    }
  }

  invert_matrix <3> (defgradinvavgc.pointer());
  for (int ii = 1; ii <= 3; ii++) {
    for (int jj = 1; jj <= 3; jj++) {
      defgradinvavgc_inv(ii,jj) = defgradinvavgc(ii,jj);
    }
  }

  // weight normalization and correction to velgrad
  MatrixTypeRealDual dudot_dXini_avg(3,3);

  for (int jj = 1; jj <= 3; jj++) {
    for (int ii = 1; ii <= 3; ii++) {
      dudot_dXini_avg.host(ii,jj) = 0.0;
      for (int kk = 1; kk <= 3; kk++) {
        dudot_dXini_avg.host(ii,jj) += (velgradmacroactual(ii,kk) - velgradavg(ii,kk))*defgradinvavgc_inv(kk,jj);
      }
    }
  }
  dudot_dXini_avg.update_device();

  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    wgtc(i,j,k) = wgtc(i,j,k)/detFavg;


    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        for (int kk = 1; kk <= 3; kk++) {
          velgrad(ii,jj,i,j,k) += dudot_dXini_avg(ii,kk)*defgradinv(kk,jj,i,j,k);
        }
      }
    }
  
  }); // end FOR_ALL_CLASS

}
