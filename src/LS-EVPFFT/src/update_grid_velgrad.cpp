#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "utilities.h"
#include "Profiler.h"
#include "math_functions.h"
#include "reduction_data_structures.h"

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
        loc_reduce.array[ic] += velgradref(ii,jj,i,j,k) * wgtc(i,j,k);
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

    // thread private arrays
    real_t xk_[3];
    real_t velhatr_[3];
    real_t velhatim_[3];

    // create views of thread private arrays
    ViewMatrixTypeReal xk(xk_,3);
    ViewMatrixTypeReal velhatr(velhatr_,3);
    ViewMatrixTypeReal velhatim(velhatim_,3);

    xk(1) = twopi*xk_gb(i);
    xk(2) = twopi*yk_gb(j);
    xk(3) = twopi*zk_gb(k);

    xkxk = xk(1)*xk(1) + xk(2)*xk(2) + xk(3)*xk(3);

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

    for (int ii = 1; ii <= 3; ii++) {
      work(ii,1,i,j,k) = velhatr(ii);
      workim(ii,1,i,j,k) = velhatim(ii);
    }

  }); // end FOR_ALL_CLASS

  MatrixTypeRealDual workG(3,3,npts1_g,npts2_g,npts3_g);

  int jj;
  jj = 1;
  for (int ii = 1; ii <= 3; ii++) {

    // prep for backward FFT
    FOR_ALL_CLASS(k, 1, npts3_cmplx+1,
                  j, 1, npts2_cmplx+1,
                  i, 1, npts1_cmplx+1, {
      data_cmplx(1,i,j,k) = work(ii,jj,i,j,k);
      data_cmplx(2,i,j,k) = workim(ii,jj,i,j,k);
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
      workG(ii,jj,i+local_start1,j+local_start2,k+local_start3) = data(i,j,k);
    }); // end FOR_ALL_CLASS
    Kokkos::fence();

  } // end for ii
  workG.update_host();

  MPI_Allreduce(MPI_IN_PLACE, workG.host_pointer(), 3*3*npts1_g*npts2_g*npts3_g, MPI_REAL_T, MPI_SUM, mpi_comm);

  FOR_ALL_CLASS(k, 1, npts3+2,
                j, 1, npts2+2,
                i, 1, npts1+2, {

    int iv1;
    int iv2;
    int jv1;
    int jv2;
    int kv1;
    int kv2;

    // thread private arrays
    real_t dvnode_[3];
    real_t Xnode_ref_[3];

    // create views of thread private arrays
    ViewMatrixTypeReal dvnode(dvnode_,3);
    ViewMatrixTypeReal Xnode_ref(Xnode_ref_,3);

    if ((i+local_start1) == 1) {
      iv1 = npts1_g;
      iv2 = 1;
    } else if ((i+local_start1) == npts1_g + 1) {
      iv1 = npts1_g;
      iv2 = 1;
    } else {
      iv1 = i+local_start1 - 1;
      iv2 = i+local_start1;
    }

    if ((j+local_start2) == 1) {
      jv1 = npts2_g;
      jv2 = 1;
    } else if ((j+local_start2) == npts2_g + 1) {
      jv1 = npts2_g;
      jv2 = 1;
    } else {
      jv1 = j+local_start2 - 1;
      jv2 = j+local_start2;
    }

    if ((k+local_start3) == 1) {
      kv1 = npts3_g;
      kv2 = 1;
    } else if ((k+local_start3) == npts3_g + 1) {
      kv1 = npts3_g;
      kv2 = 1;
    } else {
      kv1 = k+local_start3 - 1;
      kv2 = k+local_start3;
    }

    Xnode_ref(1) = float(i) - 0.5;
    Xnode_ref(2) = float(j) - 0.5;
    Xnode_ref(3) = float(k) - 0.5;

    for (int ii = 1; ii <= 3; ii++) {
      dvnode(ii) = 0.125*(workG(ii,jj,iv1,jv1,kv1) + workG(ii,jj,iv2,jv1,kv1) + 
       workG(ii,jj,iv1,jv2,kv1) + workG(ii,jj,iv2,jv2,kv1) +
       workG(ii,jj,iv1,jv1,kv2) + workG(ii,jj,iv2,jv1,kv2) + 
       workG(ii,jj,iv1,jv2,kv2) + workG(ii,jj,iv2,jv2,kv2));
    }

    for (int ii = 1; ii <= 3; ii++) {
      xnode(ii,i,j,k) += dvnode(ii)*tdot;
      for (int jj = 1; jj <= 3; jj++) {
        xnode(ii,i,j,k) += velgradrefavg(ii,jj)*Xnode_ref(jj)*tdot;
      }
    }

  }); // end FOR_ALL_CLASS
  Kokkos::fence();
  
}
