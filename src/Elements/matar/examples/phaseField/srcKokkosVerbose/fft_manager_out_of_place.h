#ifdef OUT_OF_PLACE_FFT

#pragma once

#include "matar.h"

using namespace mtr; // matar namespace

class FFTManagerOutOfPlace
{
    private:
        int* nn_;
        int  nx_;
        int  ny_;
        int  nz_;
        int  nz21_;
        int  ndim_;
        int  isign_;

    public:
        FFTManagerOutOfPlace(int * nn);
        void perform_forward_fft(double *input, double *output);
        void perform_backward_fft(double *input, double *output);

};



#ifdef HAVE_CUDA
void fftc_cufft_init_out_of_place_();
void fftc_cufft_out_of_place_(double input[], double output[], int nn[], int *ndim, int *isign);
#else
void fftc_fftw_init_out_of_place_();
void fftc_fftw_out_of_place_(double input[], double output[], int nn[], int *ndim, int *isign);
#endif


#endif
