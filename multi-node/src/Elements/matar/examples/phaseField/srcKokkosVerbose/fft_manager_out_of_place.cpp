#ifdef OUT_OF_PLACE_FFT

#include "fft_manager_out_of_place.h"

FFTManagerOutOfPlace::FFTManagerOutOfPlace(int * nn)
{
    nn_ = nn;
    nx_  = nn_[0];
    ny_  = nn_[1];
    nz_  = nn_[2];
    nz21_ = nz_/2 + 1;

    // initialize fft
    #ifdef HAVE_CUDA 
        fftc_cufft_init_out_of_place_();
    #else
        fftc_fftw_init_out_of_place_();
    #endif
}




void FFTManagerOutOfPlace::perform_forward_fft(double *input, double *output)
{
    // this function performs forward fft on "input" array and 
    // writes the result to "output" array.
    // it calls the appropriate function to perform the forward out-of-place fft
    // either using OPENMP or CUDA.

    // perform foward fft
    isign_ = -1;
    #ifdef HAVE_CUDA
        fftc_cufft_out_of_place_(input, output, nn_, &ndim_, &isign_);
    #else
        fftc_fftw_out_of_place_(input, output, nn_, &ndim_, &isign_);
    #endif

}

void FFTManagerOutOfPlace::perform_backward_fft(double *input, double *output)
{
    // this function performs backward fft on "input" array and 
    // writes the result to "output" array.
    // it calls the appropriate function to perform the backward out-of-place fft
    // either using OPENMP or CUDA.

    // perform backward fft
    isign_ = 1;
    #ifdef HAVE_CUDA
        fftc_cufft_out_of_place_(input, output, nn_, &ndim_, &isign_);
    #else
        fftc_fftw_out_of_place_(input, output, nn_, &ndim_, &isign_);
    #endif
}



#endif
