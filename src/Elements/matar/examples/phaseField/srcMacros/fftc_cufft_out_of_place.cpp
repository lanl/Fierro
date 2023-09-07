#ifdef OUT_OF_PLACE_FFT
#ifdef HAVE_CUDA

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <complex.h>

#include <cufft.h> 

// ----------------------------------------------------------------------
// CUFFT

static void
fft_cufft_forward(double *input, double *output, int nn[3])
{
  int stride = 2 * nn[0] * nn[1] * nn[2];
  int rc, i;

  static cufftHandle planD2Z;
//  typedef cuComplex cufftComplex;
  typedef cuDoubleComplex cufftDoubleComplexi;
  typedef double cufftDoubleReal;
  if (!planD2Z) {
    cufftPlan3d(&planD2Z, nn[0], nn[1], nn[2], CUFFT_D2Z);
  }

//#pragma acc data copy(data[0:batch*stride])
  {
//      printf("data1 %p\n", data);
//#pragma acc host_data use_device(data)
    {
//      printf("data2 %p\n", data);
      rc = cufftExecD2Z(planD2Z, (cufftDoubleReal *) input,
			(cufftDoubleComplex *) output);
      assert(rc == CUFFT_SUCCESS);
    }
  }
}

static void
fft_cufft_backward(double *input, double *output, int nn[3])
{
  int stride = 2 * nn[0] * nn[1] * nn[2];
  int rc, i;
  
  static cufftHandle planZ2D;
//  typedef cuComplex cufftComplex;
  typedef cuDoubleComplex cufftDoubleComplex;
  typedef double cufftDoubleReal;
  if (!planZ2D) {
    cufftPlan3d(&planZ2D, nn[0], nn[1], nn[2], CUFFT_Z2D);
  }

//#pragma acc data copy(data[0:batch*stride])
  {
//#pragma acc host_data use_device(data)
    rc = cufftExecZ2D(planZ2D, (cufftDoubleComplex *) input,
		      (cufftDoubleReal *) output);
    assert(rc == CUFFT_SUCCESS);
  }
}

// ----------------------------------------------------------------------

void fftc_cufft_out_of_place_(double input[], double output[], int nn[], int *ndim, int *isign)
{
  //assert(*ndim == 3);
  if (*isign == -1) {
    fft_cufft_forward(input, output, nn);
  } else {
    fft_cufft_backward(input, output, nn);
  }
}

void fftc_cufft_init_out_of_place_(void)
{
}


#endif
#endif
