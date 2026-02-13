#include "evpfft.h"
#include "utilities.h"
#include "Profiler.h"

void EVPFFT::forward_fft()
{
  Profiler profiler(__FUNCTION__);

  for (int ii = 1; ii <= 3; ii++) {
    for (int jj = 1; jj <= 3; jj++) {

      // prep for forward FFT
      FOR_ALL_CLASS(k, 1, npts3+1,
                    j, 1, npts2+1,
                    i, 1, npts1+1, {
        data(i,j,k) = sgPK1(ii,jj,i,j,k);
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

}


void EVPFFT::backward_fft()
{
  Profiler profiler(__FUNCTION__);

  for (int ii = 1; ii <= 3; ii++) {
    for (int jj = 1; jj <= 3; jj++) {

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
        work(ii,jj,i,j,k) = data(i,j,k);
      }); // end FOR_ALL_CLASS
      Kokkos::fence();

    } // end for jj
  } // end for ii

}
