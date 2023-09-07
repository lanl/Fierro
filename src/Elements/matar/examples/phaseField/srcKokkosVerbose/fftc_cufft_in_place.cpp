#ifdef IN_PLACE_FFT
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
fft_cufft_forward(double *data, int nn[3])
{
  int stride = 2 * nn[0] * nn[1] * nn[2];
  int rc, i;

  static cufftHandle planZ2Z;
//  typedef cuComplex cufftComplex;
  typedef cuDoubleComplex cufftDoubleComplex;
  if (!planZ2Z) {
    cufftPlan3d(&planZ2Z, nn[0], nn[1], nn[2], CUFFT_Z2Z);
  }

//#pragma acc data copy(data[0:batch*stride])
  {
//      printf("data1 %p\n", data);
//#pragma acc host_data use_device(data)
    {
//      printf("data2 %p\n", data);
      rc = cufftExecZ2Z(planZ2Z, (cufftDoubleComplex *) data,
			(cufftDoubleComplex *) data,
			CUFFT_FORWARD);
      assert(rc == CUFFT_SUCCESS);
    }
  }
}

static void
fft_cufft_backward(double *data, int nn[3])
{
  int stride = 2 * nn[0] * nn[1] * nn[2];
  int rc, i;
  
  static cufftHandle planZ2Z;
//  typedef cuComplex cufftComplex;
  typedef cuDoubleComplex cufftDoubleComplex;
  if (!planZ2Z) {
    cufftPlan3d(&planZ2Z, nn[0], nn[1], nn[2], CUFFT_Z2Z);
  }

//#pragma acc data copy(data[0:batch*stride])
  {
//#pragma acc host_data use_device(data)
    rc = cufftExecZ2Z(planZ2Z, (cufftDoubleComplex *) data,
		      (cufftDoubleComplex *) data,
		      CUFFT_INVERSE);
    assert(rc == CUFFT_SUCCESS);
  }
}

// ----------------------------------------------------------------------

void fftc_cufft_in_place_(double data[], int nn[], int *ndim, int *isign)
{
  //assert(*ndim == 3);
  if (*isign == -1) {
    fft_cufft_forward(data, nn);
  } else {
    fft_cufft_backward(data, nn);
  }
}

void fftc_cufft_init_in_place_(void)
{
}


#endif
#endif
