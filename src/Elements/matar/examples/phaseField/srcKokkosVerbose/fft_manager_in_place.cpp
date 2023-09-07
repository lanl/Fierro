#ifdef IN_PLACE_FFT

# include "fft_manager_in_place.h"


FFTManagerInPlace::FFTManagerInPlace(int * nn)
{
    nn_ = nn;
    nx_  = nn_[0];
    ny_  = nn_[1];
    nz_  = nn_[2];
    data_ = CArrayKokkos<double>(nx_, ny_, nz_, 2);

    // calculate ndim
    ndim_ = 0;
    for (int i = 0; i < 3; i++) {
        if (nn_[i] > 1) ++ndim_;
    }

    // initialize fft
    #ifdef HAVE_CUDA 
        fftc_cufft_init_in_place_();
    #else
        fftc_fftw_init_in_place_();
    #endif
}


void FFTManagerInPlace::prep_for_forward_fft_(double *input)
{
    // this function writes the data in "input" array to "data_" array
    // in order to ready "data_" for in-place forward fft.

    // create view of input
    auto input_view = ViewCArrayKokkos <double> (input, nx_, ny_, nz_);

    // write input to data for in-place forward fft
    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0}, {nx_, ny_, nz_}),
        KOKKOS_CLASS_LAMBDA(const int i, const int j, const int k){
            data_(i,j,k,0) = input_view(i,j,k);
            data_(i,j,k,1) = 0.0; 
    });

}


void FFTManagerInPlace::get_forward_fft_result_(double *output)
{
    // this function writes the result of in-place forward fft
    // in "data_" array into "output" array.

    // create view of output
    auto output_view = ViewCArrayKokkos <double> (output, nx_, ny_, nz_, 2);

    // write data to output after in-place fft
    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0}, {nx_, ny_, nz_}),
        KOKKOS_CLASS_LAMBDA(const int i, const int j, const int k){
            output_view(i,j,k,0) = data_(i,j,k,0);
            output_view(i,j,k,1) = data_(i,j,k,1);

    });

}


void FFTManagerInPlace::prep_for_backward_fft_(double *input)
{
    // this function writes the data in "input" array to "data_" array
    // in order to ready "data_" for in-place backward fft.

    // create view of input
    auto input_view = ViewCArrayKokkos <double> (input, nx_, ny_, nz_, 2);

    // write input to data for in-place fft
    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0}, {nx_, ny_, nz_}),
        KOKKOS_CLASS_LAMBDA(const int i, const int j, const int k){
            data_(i,j,k,0) = input_view(i,j,k,0);
            data_(i,j,k,1) = input_view(i,j,k,1);
    });

}


void FFTManagerInPlace::get_backward_fft_result_(double *output)
{
    // this function writes the result of in-place backward fft
    // in "data_" array into "output" array.

    // create view of output
    auto output_view = ViewCArrayKokkos <double> (output, nx_, ny_, nz_);

    // write data to output after in-place fft
    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0}, {nx_, ny_, nz_}),
        KOKKOS_CLASS_LAMBDA(const int i, const int j, const int k){
            output_view(i,j,k) = data_(i,j,k,0);
    });

}


void FFTManagerInPlace::perform_forward_fft(double *input, double *output)
{
    // this function performs forward fft on "input" array and 
    // writes the result to "output" array.
    // it calls the appropriate function to perform the forward in-place fft
    // either using OPENMP or CUDA.

    // prep for forward fft
    prep_for_forward_fft_(input);

    // perform foward fft
    isign_ = -1;
    #ifdef HAVE_CUDA
        fftc_cufft_in_place_(data_.pointer(), nn_, &ndim_, &isign_);
    #else
        fftc_fftw_in_place_(data_.pointer(), nn_, &ndim_, &isign_);
    #endif

    // get result after performing foward fft
    get_forward_fft_result_(output);
}


void FFTManagerInPlace::perform_backward_fft(double *input, double *output)
{
    // this function performs backward fft on "input" array and 
    // writes the result to "output" array.
    // it calls the appropriate function to perform the backward in-place fft
    // either using OPENMP or CUDA.

    // prep for backward fft
    prep_for_backward_fft_(input);

    // perform backward fft
    isign_ = 1;
    #ifdef HAVE_CUDA
        fftc_cufft_in_place_(data_.pointer(), nn_, &ndim_, &isign_);
    #else
        fftc_fftw_in_place_(data_.pointer(), nn_, &ndim_, &isign_);
    #endif

    // get result after performing backward fft
    get_backward_fft_result_(output);
}

#endif
