#pragma once

#include <assert.h>

#ifdef USE_FFTW
#include "fftw3.h"
#elif USE_CUFFT
#include <cufft.h>
#endif

template <typename T>
class FFT3D_R2C
{    
public:
    const size_t N1;
    const size_t N2;
    const size_t N3;
    int result;
#ifdef USE_FFTW
    fftwf_plan planf_forward;
    fftwf_plan planf_backward;
    fftw_plan plan_forward;
    fftw_plan plan_backward;
#elif USE_CUFFT
    cufftHandle plan_forward;
    cufftHandle plan_backward;
#endif
    FFT3D_R2C(size_t N1_, size_t N2_, size_t N3_);
    void make_plan(T* input, T*output);
    void forward(T* input, T* output);
    void backward(T* input, T* output);

private:
    void make_plan_(float* input, float* output);
    void make_plan_(double* input, double* output);
    void forward_(float* input, float* output);
    void forward_(double* input, double* output);
    void backward_(float* input, float* output);
    void backward_(double* input, double* output);
};

template <typename T>
FFT3D_R2C<T>::FFT3D_R2C(size_t N1_, size_t N2_, size_t N3_)
    : N1(N1_)
    , N2(N2_)
    , N3(N3_)
{
}

template <typename T>
void FFT3D_R2C<T>::make_plan(T* input, T* output)
{
    make_plan_(input, output);
}

template <typename T>
void FFT3D_R2C<T>::forward(T* input, T* output)
{
    forward_(input, output);
}

template <typename T>
void FFT3D_R2C<T>::backward(T* input, T* output)
{
    backward_(input, output);
}

template <typename T>
void FFT3D_R2C<T>::make_plan_(float* input, float* output)
{
#ifdef USE_FFTW
    planf_forward = fftwf_plan_dft_r2c_3d(N1, N2, N3, (float *) input, (fftwf_complex *) output, FFTW_ESTIMATE);
    planf_backward = fftwf_plan_dft_c2r_3d(N1, N2, N3, (fftwf_complex *) output, (float *) input, FFTW_ESTIMATE);
#elif USE_CUFFT
    cufftPlan3d(&plan_forward, N1, N2, N3, CUFFT_R2C);
    cufftPlan3d(&plan_backward, N1, N2, N3, CUFFT_C2R);
#endif
}

template <typename T>
void FFT3D_R2C<T>::make_plan_(double* input, double* output)
{
#ifdef USE_FFTW
    plan_forward = fftw_plan_dft_r2c_3d(N1, N2, N3, (double *) input, (fftw_complex *) output, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_c2r_3d(N1, N2, N3, (fftw_complex *) output, (double *) input, FFTW_ESTIMATE);
#elif USE_CUFFT
    cufftPlan3d(&plan_forward, N1, N2, N3, CUFFT_D2Z);
    cufftPlan3d(&plan_backward, N1, N2, N3, CUFFT_Z2D);
#endif
}

template <typename T>
void FFT3D_R2C<T>::forward_(float* input, float* output)
{
#ifdef USE_FFTW
    fftwf_execute_dft_r2c(planf_forward, (float *) input, (fftwf_complex *) output);
#elif USE_CUFFT
    result = cufftExecR2C(plan_forward, (cufftReal *) input, (cufftComplex *) output);
    assert(result == CUFFT_SUCCESS);
#endif
}

template <typename T>
void FFT3D_R2C<T>::forward_(double* input, double* output)
{
#ifdef USE_FFTW
    fftw_execute_dft_r2c(plan_forward, (double *) input, (fftw_complex *) output);
#elif USE_CUFFT
    result = cufftExecD2Z(plan_forward, (cufftDoubleReal *) input, (cufftDoubleComplex *) output);
    assert(result == CUFFT_SUCCESS);
#endif
}

template <typename T>
void FFT3D_R2C<T>::backward_(float* input, float* output)
{
#ifdef USE_FFTW
    fftwf_execute_dft_c2r(planf_backward, (fftwf_complex *) input, (float *) output);
#elif USE_CUFFT
    result = cufftExecC2R(plan_backward, (cufftComplex *) input, (cufftReal *) output);
    assert(result == CUFFT_SUCCESS);
#endif
}

template <typename T>
void FFT3D_R2C<T>::backward_(double* input, double* output)
{
if (!plan_backward) printf("Not Good!\n");
#ifdef USE_FFTW
    fftw_execute_dft_c2r(plan_backward, (fftw_complex *) input, (double *) output);
#elif USE_CUFFT
    result = cufftExecZ2D(plan_backward, (cufftDoubleComplex *) input, (cufftDoubleReal *) output);
    assert(result == CUFFT_SUCCESS);
#endif
}

