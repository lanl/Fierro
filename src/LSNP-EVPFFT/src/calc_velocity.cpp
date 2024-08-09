#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "utilities.h"
#include "Profiler.h"
#include "math_functions.h"
#include "reduction_data_structures.h"
#include "mod_frequency.h"

void EVPFFT::calc_velocity()
{
  Profiler profiler(__FUNCTION__);

  const size_t n = 9; // for velgradrefavg (9)
  ArrayType <real_t, n> all_reduce;

  MatrixTypeRealDual velgradrefavg(3,3);
  
  const real_t twopi = 8.*ATAN(1.0);

  Kokkos::parallel_reduce(
    Kokkos::MDRangePolicy<Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER>>({1,1,1}, {npts3+1,npts2+1,npts1+1}),
    KOKKOS_CLASS_LAMBDA(const int k, const int j, const int i, ArrayType <real_t,n> & loc_reduce) {

    if (ibc == 0) {
      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          velgradref(ii,jj,i,j,k) = 0.0;
          for (int kk = 1; kk <= 3; kk++) {
            velgradref(ii,jj,i,j,k) += (velgrad(ii,kk,i,j,k))*defgrad(kk,jj,i,j,k);
          }
        }
      }
    }

    // averages packed in array
    int ic;
    ic = -1;

    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        ic = ic + 1;
        loc_reduce.array[ic] += velgradref(ii,jj,i,j,k) * wgt;
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
      velgradrefavg.host(ii,jj) = all_reduce.array[ic];
    }
  }
  velgradrefavg.update_device();

  for (int ii = 1; ii <= 3; ii++) {
    for (int jj = 1; jj <= 3; jj++) {

      // prep for forward FFT
      FOR_ALL_CLASS(k, 1, npts3+1,
                    j, 1, npts2+1,
                    i, 1, npts1+1, {
        data(i,j,k) = velgradref(ii,jj,i,j,k);
      }); // end FOR_ALL_CLASS
      Kokkos::fence();

#if defined USE_FFTW || USE_MKL
      data.update_host();

      // perform forward FFT
      fft->forward(data.host_pointer(), (std::complex<double>*) data_cmplx.host_pointer());
      data_cmplx.update_device();
#else
      // perform forward FFT
      fft->forward(data.device_pointer(), (std::complex<double>*) data_cmplx.device_pointer());
#endif

      // write result to ouput
      FOR_ALL_CLASS(k, 1, npts3_cmplx+1,
                    j, 1, npts2_cmplx+1,
                    i, 1, npts1_cmplx+1, {
        work(ii,jj,i,j,k)   = data_cmplx(1,i,j,k);
        workim(ii,jj,i,j,k) = data_cmplx(2,i,j,k);
      }); // end FOR_ALL_CLASS
      Kokkos::fence();

    } // end for jj
  } // end for ii

  FOR_ALL_CLASS(k, 1, npts3_cmplx+1,
                j, 1, npts2_cmplx+1,
                i, 1, npts1_cmplx+1, {

    real_t xkxk;
    real_t krot_ckrot;

    // thread private arrays
    real_t xk_[3];
    real_t velhatr_[3];
    real_t velhatim_[3];
    real_t krot_re_[3];
    real_t krot_im_[3];

    // create views of thread private arrays
    ViewMatrixTypeReal xk(xk_,3);
    ViewMatrixTypeReal velhatr(velhatr_,3);
    ViewMatrixTypeReal velhatim(velhatim_,3);
    ViewMatrixTypeReal krot_re(krot_re_,3);
    ViewMatrixTypeReal krot_im(krot_im_,3);

    xk(1) = twopi*xk_gb(i);
    xk(2) = twopi*yk_gb(j);
    xk(3) = twopi*zk_gb(k);

    if (igamma == 0) {

      xkxk = xk(1)*xk(1) + xk(2)*xk(2) + xk(3)*xk(3);
  
      //if (xkxk > 0.0 && 
      //  !(i + local_start1_cmplx == npts1_g/2+1 || 
      //    j + local_start2_cmplx == npts2_g/2+1 || 
      //    (npts3_g > 1 && k + local_start3_cmplx == npts3_g/2+1))) {
      if (xkxk > 0.0) {

        for (int ii = 1; ii <= 3; ii++) {
          velhatr(ii) = 0.0;
          velhatim(ii) = 0.0;
          for (int jj = 1; jj <= 3; jj++) {
            velhatr(ii) += workim(ii,jj,i,j,k)*xk(jj)/xkxk;
            velhatim(ii) += - work(ii,jj,i,j,k)*xk(jj)/xkxk;
          }
        }
  
      } else {
  
        for (int ii = 1; ii <= 3; ii++) {
          velhatr(ii) = 0.0;
          velhatim(ii) = 0.0;
        }
      }

    } else if (igamma == 1) {

      mod_frequency(xk.pointer(), krot_re.pointer(), krot_im.pointer());

      krot_ckrot = POW2(krot_re(1)) + POW2(krot_im(1)) + POW2(krot_re(2)) +
          POW2(krot_im(2)) + POW2(krot_re(3)) + POW2(krot_im(3));
      // if (abs(krot_ckrot) >= 1.0e-15 && 
      //   !(i + local_start1_cmplx == npts1_g/2+1 || 
      //     j + local_start2_cmplx == npts2_g/2+1 || 
      //     (npts3_g > 1 && k + local_start3_cmplx == npts3_g/2+1))) {
      if (abs(krot_ckrot) >= 1.0e-15) {

        for (int ii = 1; ii <= 3; ii++) {
          velhatr(ii) = 0.0;
          velhatim(ii) = 0.0;
          for (int jj = 1; jj <= 3; jj++) {
            velhatr(ii) += (work(ii,jj,i,j,k)*krot_re(jj) + workim(ii,jj,i,j,k)*krot_im(jj))/krot_ckrot;
            velhatim(ii) += (workim(ii,jj,i,j,k)*krot_re(jj) - work(ii,jj,i,j,k)*krot_im(jj))/krot_ckrot;
          }
        }

      } else {

        for (int ii = 1; ii <= 3; ii++) {
          velhatr(ii) = 0.0;
          velhatim(ii) = 0.0;
        }
      }

    }

    for (int ii = 1; ii <= 3; ii++) {
      work(ii,1,i,j,k) = velhatr(ii);
      workim(ii,1,i,j,k) = velhatim(ii);
    }

  }); // end FOR_ALL_CLASS

  //for (int i = 1; i <= npts1_cmplx; i++) {
  //for (int j = 1; j <= npts2_cmplx; j++) {
  //for (int k = 1; k <= npts3_cmplx; k++) {
  //  printf("%d %d %d %24.14E\n", i, j, k,  work(1,1,i,j,k));
  //  printf("%d %d %d %24.14E\n", i, j, k,  workim(1,1,i,j,k));
  //}
  //}
  //}
  //exit(1);

  for (int ii = 1; ii <= 3; ii++) {

    // prep for backward FFT
    FOR_ALL_CLASS(k, 1, npts3_cmplx+1,
                  j, 1, npts2_cmplx+1,
                  i, 1, npts1_cmplx+1, {
      data_cmplx(1,i,j,k) = work(ii,1,i,j,k);
      data_cmplx(2,i,j,k) = workim(ii,1,i,j,k);
    }); // end FOR_ALL_CLASS
    Kokkos::fence();

#if defined USE_FFTW || USE_MKL
    data_cmplx.update_host();

    // perform backward FFT
    fft->backward((std::complex<double>*) data_cmplx.host_pointer(), data.host_pointer());
    data.update_device();
#else
    // perform backward FFT
    fft->backward((std::complex<double>*) data_cmplx.device_pointer(), data.device_pointer());
#endif

    // write result to ouput
    FOR_ALL_CLASS(k, 1, npts3+1,
                  j, 1, npts2+1,
                  i, 1, npts1+1, {
      velocity(ii,i,j,k) = data(i,j,k);
    }); // end FOR_ALL_CLASS
    Kokkos::fence();

  } // end for ii

  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    // thread private arrays
    real_t X_ref_[3];

    // create views of thread private arrays
    ViewMatrixTypeReal X_ref(X_ref_,3);

    if (igamma == 0) {

      X_ref(1) = float(i + local_start1);
      X_ref(2) = float(j + local_start2);
      X_ref(3) = float(k + local_start3);

      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          velocity(ii,i,j,k) += velgradrefavg(ii,jj)*X_ref(jj);
        }
      }

    } else if (igamma == 1) {

      X_ref(1) = float(i + local_start1) - 0.5;
      X_ref(2) = float(j + local_start2) - 0.5;
      X_ref(3) = float(k + local_start3) - 0.5;
      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          velocity(ii,i,j,k) += velgradrefavg(ii,jj)*X_ref(jj);
        }
      }

    }

  }); // end FOR_ALL_CLASS
  Kokkos::fence();

}
