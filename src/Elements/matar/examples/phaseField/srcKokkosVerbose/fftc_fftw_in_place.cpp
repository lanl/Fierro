#ifdef IN_PLACE_FFT
#ifdef HAVE_OPENMP

#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <complex.h>

#include <fftw3.h>
//#ifdef FFTW_OMP
#include <omp.h>
//#endif

// ----------------------------------------------------------------------
// FFTW

static void
fft_forward_fftw(double *data, int nn[3])
{
  static fftw_plan plan;
  if (!plan) {   
    plan = fftw_plan_dft_3d(nn[2], nn[1], nn[0],
			    (fftw_complex *) data, (fftw_complex *) data,
			    FFTW_FORWARD, FFTW_ESTIMATE);

  }
  
  fftw_execute_dft(plan, (fftw_complex *) data, (fftw_complex *) data);
}

static void
fft_backward_fftw(double *data, int nn[3])
{
  static fftw_plan plan;
  if (!plan) {
    plan = fftw_plan_dft_3d(nn[2], nn[1], nn[0],
			    (fftw_complex *) data, (fftw_complex *) data,
			    FFTW_BACKWARD, FFTW_ESTIMATE);
  }

  fftw_execute_dft(plan, (fftw_complex *) data, (fftw_complex *) data);
}

void fftc_fftw_in_place_(double data[], int nn[], int *ndim, int *isign)
{
  if (*isign == -1) {
    fft_forward_fftw(data, nn);
  } else {
    fft_backward_fftw(data, nn);
  }
}

void fftc_fftw_init_in_place_(void)
{
//#ifdef FFTW_OMP
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
//#endif
}

#endif
#endif
