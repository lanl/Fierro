#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "utilities.h"
#include "Profiler.h"
#include "math_functions.h"
#include "reduction_data_structures.h"
#include "mod_frequency.h"

void EVPFFT::update_grid_velgrad()
{
  Profiler profiler(__FUNCTION__);

  const size_t n = 9; // for velgradrefavg (9)
  ArrayType <real_t, n> all_reduce;

  MatrixTypeRealDual velgradrefavg(3,3);
  
  const real_t twopi = 8.*ATAN(1.0);

  Kokkos::parallel_reduce(
    Kokkos::MDRangePolicy<Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER>>({1,1,1}, {npts3+1,npts2+1,npts1+1}),
    KOKKOS_CLASS_LAMBDA(const int k, const int j, const int i, ArrayType <real_t,n> & loc_reduce) {

    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        velgradref(ii,jj,i,j,k) = 0.0;
        for (int kk = 1; kk <= 3; kk++) {
          velgradref(ii,jj,i,j,k) += (velgrad(ii,kk,i,j,k))*defgrad(kk,jj,i,j,k);
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
  
      if (xkxk > 0.0 && 
        !(i + local_start1_cmplx == npts1_g/2+1 || 
          j + local_start2_cmplx == npts2_g/2+1 || 
          (npts3_g > 1 && k + local_start3_cmplx == npts3_g/2+1))) {
  
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
      if (abs(krot_ckrot) >= 1.0e-15 && 
        !(i + local_start1_cmplx == npts1_g/2+1 || 
          j + local_start2_cmplx == npts2_g/2+1 || 
          (npts3_g > 1 && k + local_start3_cmplx == npts3_g/2+1))) {

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
      work(ii,1,i,j,k) = data(i,j,k);
    }); // end FOR_ALL_CLASS
    Kokkos::fence();

  } // end for ii

  FOR_ALL_CLASS(k, 1, npts3+2,
                j, 1, npts2+2,
                i, 1, npts1+2, {

    int iv1;
    int iv2;
    int jv1;
    int jv2;
    int kv1;
    int kv2;
    int in;
    int jn;
    int kn;

    // thread private arrays
    real_t dvnode_[3];
    real_t Xnode_ref_[3];

    // create views of thread private arrays
    ViewMatrixTypeReal dvnode(dvnode_,3);
    ViewMatrixTypeReal Xnode_ref(Xnode_ref_,3);

    Xnode_ref(1) = 1.0*(i + local_start1) - 0.5;
    Xnode_ref(2) = 1.0*(j + local_start2) - 0.5;
    Xnode_ref(3) = 1.0*(k + local_start3) - 0.5;

    if (igamma == 0) {

      if (i == 1) {
        iv1 = npts1;
        iv2 = 1;
      } else if (i == npts1 + 1) {
        iv1 = npts1;
        iv2 = 1;
      } else {
        iv1 = i - 1;
        iv2 = i;
      }
  
      if (j == 1) {
        jv1 = npts2;
        jv2 = 1;
      } else if (j == npts2 + 1) {
        jv1 = npts2;
        jv2 = 1;
      } else {
        jv1 = j - 1;
        jv2 = j;
      }
  
      if (k == 1) {
        kv1 = npts3;
        kv2 = 1;
      } else if (k == npts3 + 1) {
        kv1 = npts3;
        kv2 = 1;
      } else {
        kv1 = k - 1;
        kv2 = k;
      }

      for (int ii = 1; ii <= 3; ii++) {
        dvnode(ii) = 0.125*(work(ii,1,iv1,jv1,kv1) + work(ii,1,iv2,jv1,kv1) + 
         work(ii,1,iv1,jv2,kv1) + work(ii,1,iv2,jv2,kv1) +
         work(ii,1,iv1,jv1,kv2) + work(ii,1,iv2,jv1,kv2) + 
         work(ii,1,iv1,jv2,kv2) + work(ii,1,iv2,jv2,kv2));
      }

    } else if (igamma == 1) {

      if (i == npts1 + 1) {
        in = 1;
      } else {
        in = i;
      }

      if (j == npts2 + 1) {
        jn = 1;
      } else {
        jn = j;
      }

      if (k == npts3 + 1) {
        kn = 1;
      } else {
        kn = k;
      }

      for (int ii = 1; ii <= 3; ii++) {
        dvnode(ii) = work(ii,1,in,jn,kn);
      }

    }

    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        dvnode(ii) += velgradrefavg(ii,jj)*Xnode_ref(jj);
      }
    }

    for (int ii = 1; ii <= 3; ii++) {
      xnode(ii,i,j,k) += dvnode(ii)*tdot;
    }

  }); // end FOR_ALL_CLASS
  Kokkos::fence();
  
}
