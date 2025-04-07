#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "utilities.h"


void EVPFFT::calc_wfhat()
{
  FOR_ALL_CLASS(kz, 1, npts3+1,
                ky, 1, npts2+1,
                kx, 1, npts1+1, {

    int kx_g;
    int ky_g;
    int kz_g;
    kx_g = kx + local_start1;
    ky_g = ky + local_start2;
    kz_g = kz + local_start3;

    real_t xn2;
    
    real_t x_[3];

    ViewMatrixTypeReal x(x_,3);

    if (kx_g <= npts1_g/2){
      x(1) = 1.0*(kx_g - 1);
    } else if (kx_g > npts1_g/2){
      x(1) = 1.0*(kx_g - npts1_g - 1);
    } else {
      x(1) = 0.0;
    }

    if (ky_g <= npts2_g/2){
      x(2) = 1.0*(ky_g - 1);
    } else if (ky_g > npts2_g/2){
      x(2) = 1.0*(ky_g - npts2_g - 1);
    } else {
      x(2) = 0.0;
    }

    if (kz_g <= npts3_g/2){
      x(3) = 1.0*(kz_g - 1);
    } else if (kz_g > npts3_g/2){
      x(3) = 1.0*(kz_g - npts3_g - 1);
    } else {
      x(3) = 0.0;
    }

    xn2 = (POW2(x(1)) + POW2(x(2)) + POW2(x(3)));
    if (xn2 > 0.0) {
      data(kx,ky,kz) = 1.0/POW2(xn2);
    } else {
      data(kx,ky,kz) = 0.0;
    }

  });

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
    wfhat_re(i,j,k) = data_cmplx(1,i,j,k);
    wfhat_im(i,j,k) = data_cmplx(2,i,j,k);
  }); // end FOR_ALL_CLASS
  Kokkos::fence();

}
