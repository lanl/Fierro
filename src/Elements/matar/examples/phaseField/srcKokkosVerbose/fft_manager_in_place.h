#ifdef IN_PLACE_FFT

#pragma once

#include "matar.h"

using namespace mtr; // matar namespace

class FFTManagerInPlace
{
    private:
        int* nn_;
        int  nx_;
        int  ny_;
        int  nz_;
        int  ndim_;
        int  isign_;
        CArrayKokkos<double> data_;

    public:
        FFTManagerInPlace(int * nn);
        void perform_forward_fft(double *input, double *output);
        void perform_backward_fft(double *input, double *output);

        void prep_for_forward_fft_(double *input);
        void get_forward_fft_result_(double *output);
        void prep_for_backward_fft_(double *input);
        void get_backward_fft_result_(double *output);

};


#ifdef HAVE_CUDA
void fftc_cufft_init_in_place_();
void fftc_cufft_in_place_(double data[], int nn[], int *ndim, int *isign);
#else
void fftc_fftw_init_in_place_();
void fftc_fftw_in_place_(double data[], int nn[], int *ndim, int *isign);
#endif


#endif
