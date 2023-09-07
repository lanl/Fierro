#ifdef OUT_OF_PLACE_FFT
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
fft_forward_fftw(double *input, double *output, int nn[3])
{
  static fftw_plan plan;
  if (!plan) {   
    plan = fftw_plan_dft_r2c_3d(nn[0], nn[1], nn[2],
			    (double *) input, (fftw_complex *) output, 
                            FFTW_ESTIMATE);

  }
  
  fftw_execute_dft_r2c(plan, (double *) input, (fftw_complex *) output);
}

static void
fft_backward_fftw(double *input, double *output, int nn[3])
{
  static fftw_plan plan;
  if (!plan) {
    plan = fftw_plan_dft_c2r_3d(nn[0], nn[1], nn[2],
			    (fftw_complex *) input, (double *) output,
			    FFTW_ESTIMATE);
  }

  fftw_execute_dft_c2r(plan, (fftw_complex *) input, (double *) output);
}

void fftc_fftw_out_of_place_(double input[], double output[], int nn[], int *ndim, int *isign)
{
  if (*isign == -1) {
    fft_forward_fftw(input, output, nn);
  } else {
    fft_backward_fftw(input, output, nn);
  }
}

void fftc_fftw_init_out_of_place_(void)
{
//#ifdef FFTW_OMP
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
//#endif
}

#endif
#endif
