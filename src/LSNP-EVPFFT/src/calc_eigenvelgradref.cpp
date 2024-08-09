#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "utilities.h"
#include "mod_frequency.h"

void EVPFFT::calc_eigenvelgradref()
{

  const real_t twopi = 8.*ATAN(1.0);
  real_t indf_[npts1*npts2*npts3];
  real_t wfn_[npts1*npts2*npts3];

  // create views of thread private arrays
  ViewMatrixTypeReal indf(indf_,npts1,npts2,npts3);
  ViewMatrixTypeReal wfn(wfn_,npts1,npts2,npts3);

  // indicator function for velocity (1.0 if applied)
  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    int i_g;
    int j_g;
    int k_g;
    int i_gs;
    int j_gs;
    int k_gs;

    i_g = i + local_start1;
    j_g = j + local_start2;
    k_g = k + local_start3;
    i_gs = i_g - dnpts1_g;
    j_gs = j_g - dnpts2_g;
    k_gs = k_g - dnpts3_g;

    if ( (i_gs == 1 || i_gs == npts1_g - 2*dnpts1_g + 1) && (j_gs >= 1 && j_gs <= npts1_g-2*dnpts1_g+1 && k_gs >= 1 && k_gs <= npts3_g-2*dnpts3_g+1) ||
         (j_gs == 1 || j_gs == npts2_g - 2*dnpts2_g + 1) && (i_gs >= 1 && i_gs <= npts2_g-2*dnpts2_g+1 && k_gs >= 1 && k_gs <= npts3_g-2*dnpts3_g+1) ||
         (k_gs == 1 || k_gs == npts2_g - 2*dnpts2_g + 1) && (i_gs >= 1 && i_gs <= npts2_g-2*dnpts2_g+1 && j_gs >= 1 && j_gs <= npts2_g-2*dnpts2_g+1)){

      indf(i,j,k) = 1.0;
 
    } else if (i_g == 1 || j_g == 1 || k_g == 1) {

      indf(i,j,k) = 1.0;

    } else {

      indf(i,j,k) = 0.0;

    }
    data(i,j,k) = indf(i,j,k);

  });

  // normalizing function for interpolation 
#if defined USE_FFTW || USE_MKL
  data.update_host();

  // perform forward FFT
  fft->forward(data.host_pointer(), (std::complex<double>*) data_cmplx.host_pointer());
  data_cmplx.update_device();
#else
  // perform forward FFT
  fft->forward(data.device_pointer(), (std::complex<double>*) data_cmplx.device_pointer());
#endif

  // convolution
  FOR_ALL_CLASS(k, 1, npts3_cmplx+1,
                j, 1, npts2_cmplx+1,
                i, 1, npts1_cmplx+1, {
    data_cmplx(1,i,j,k) = data_cmplx(1,i,j,k)*wfhat_re(i,j,k);
    data_cmplx(2,i,j,k) = data_cmplx(2,i,j,k)*wfhat_re(i,j,k);
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
    wfn(i,j,k) = 1.0/data(i,j,k);
  }); // end FOR_ALL_CLASS
  Kokkos::fence();

  // applied velocity interpolation
  for (int ii = 1; ii <= 3; ii++) {

    // prep for forward FFT
    FOR_ALL_CLASS(k, 1, npts3+1,
                  j, 1, npts2+1,
                  i, 1, npts1+1, {
      data(i,j,k) = velapp_node(ii,i,j,k)*indf(i,j,k);
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

    // convolution
    FOR_ALL_CLASS(k, 1, npts3_cmplx+1,
                  j, 1, npts2_cmplx+1,
                  i, 1, npts1_cmplx+1, {
      data_cmplx(1,i,j,k) = data_cmplx(1,i,j,k)*wfhat_re(i,j,k);
      data_cmplx(2,i,j,k) = data_cmplx(2,i,j,k)*wfhat_re(i,j,k);
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
      if (indf(i,j,k) < 1.0) {
        velapp_node(ii,i,j,k) = data(i,j,k)*wfn(i,j,k);
      }
    }); // end FOR_ALL_CLASS
    Kokkos::fence();

  } // end for ii

  // FFT of applied velocity and gradient in Fourier space
  for (int ii = 1; ii <= 3; ii++) {

    // prep for forward FFT
    FOR_ALL_CLASS(k, 1, npts3+1,
                  j, 1, npts2+1,
                  i, 1, npts1+1, {
      data(i,j,k) = velapp_node(ii,i,j,k);
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

    // gradient
    FOR_ALL_CLASS(k, 1, npts3_cmplx+1,
                  j, 1, npts2_cmplx+1,
                  i, 1, npts1_cmplx+1, {

      // thread private arrays
      real_t xk_[3];
      real_t krot_re_[3];
      real_t krot_im_[3];

      // create views of thread private arrays
      ViewMatrixTypeReal xk(xk_,3);
      ViewMatrixTypeReal krot_re(krot_re_,3);
      ViewMatrixTypeReal krot_im(krot_im_,3);

      xk(1) = twopi*xk_gb(i);
      xk(2) = twopi*yk_gb(j);
      xk(3) = twopi*zk_gb(k);
      mod_frequency(xk.pointer(), krot_re.pointer(), krot_im.pointer());

      
      // if (!(i + local_start1_cmplx == npts1_g/2+1 || 
      //       j + local_start2_cmplx == npts2_g/2+1 || 
      //       (npts3_g > 1 && k + local_start3_cmplx == npts3_g/2+1))) {

        for (int jj = 1; jj <= 3; jj++) {
          work(ii,jj,i,j,k)   = data_cmplx(1,i,j,k)*krot_re(jj) - data_cmplx(2,i,j,k)*krot_im(jj);
          workim(ii,jj,i,j,k) = data_cmplx(1,i,j,k)*krot_im(jj) + data_cmplx(2,i,j,k)*krot_re(jj);
        }
        
//      } else {
//
//        for (int jj = 1; jj <= 3; jj++) {
//          work(ii,jj,i,j,k)   = 0.0;
//          workim(ii,jj,i,j,k) = 0.0;
//        }
//
//      }
        
    }); // end FOR_ALL_CLASS
    Kokkos::fence();

  } // end for ii

  // gradient to real space
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
        eigenvelgradref(ii,jj,i,j,k) = data(i,j,k);
      }); // end FOR_ALL_CLASS
      Kokkos::fence();

    } // end for jj
  } // end for ii


//  // print for each cpu
//  for (int k = 1; k <= npts3; k++) {
//  for (int j = 1; j <= npts2; j++) {
//  for (int i = 1; i <= npts1; i++) {
//
//    int i_g;
//    int j_g;
//    int k_g;
//  
//    i_g = i + local_start1;
//    j_g = j + local_start2;
//    k_g = k + local_start3;
//    printf("%d %d %d %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E\n",i_g,j_g,k_g,
//     eigenvelgradref(1,1,i,j,k),eigenvelgradref(2,1,i,j,k),eigenvelgradref(3,1,i,j,k),
//     eigenvelgradref(1,2,i,j,k),eigenvelgradref(2,2,i,j,k),eigenvelgradref(3,2,i,j,k),
//     eigenvelgradref(1,3,i,j,k),eigenvelgradref(2,3,i,j,k),eigenvelgradref(3,3,i,j,k));
//    printf("%d %d %d %16.8E %16.8E %16.8E\n",i_g,j_g,k_g,velapp_node(1,i,j,k),velapp_node(2,i,j,k),velapp_node(3,i,j,k));
//
//  }
//  }
//  }
//
//  exit(1);

}
