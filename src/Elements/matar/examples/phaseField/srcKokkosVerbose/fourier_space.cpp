#include <iostream>

#include "fourier_space.h"

FourierSpace::FourierSpace(int* nn, double* delta)
{
    // initialize class data
    nn_    = nn;
    nx_    = nn[0];
    ny_    = nn[1];
    nz_    = nn[2];
    delta_ = delta;
    dx_    = delta[0];
    dy_    = delta[1];
    dz_    = delta[2];

    // kx_, ky_, and kz_ initialization
    kx_ = CArrayKokkos<double>(nx_);
    ky_ = CArrayKokkos<double>(ny_);
#ifdef IN_PLACE_FFT
    kz_ = CArrayKokkos<double>(nz_);
#elif OUT_OF_PLACE_FFT
    nz21_ = nz_/2 + 1;
    kz_ = CArrayKokkos<double>(nz21_);
#endif
    

    // set values of kx_, ky_, and kz_
    set_kx_ky_kz_();

}


void FourierSpace::set_kx_ky_kz_()
{
    // calculate kx_
    Kokkos::parallel_for(
        Kokkos::RangePolicy<>(0, nx_),
        KOKKOS_CLASS_LAMBDA(const int i){
            int ti;
            ti = i;
            if (ti > nx_/2) ti = ti - nx_;
            kx_(i) = (float(ti) * twopi_) / (nx_ * dx_);
    });


    // calculate ky_
    Kokkos::parallel_for(
        Kokkos::RangePolicy<>(0, ny_),
        KOKKOS_CLASS_LAMBDA(const int j){
            int tj;
            tj = j;
            if (tj > ny_/2) tj = tj - ny_;
            ky_(j) = (float(tj) * twopi_) / (ny_ * dy_);
    });


    // calculate kz_ for in-place-fft
#ifdef IN_PLACE_FFT
    Kokkos::parallel_for(
        Kokkos::RangePolicy<>(0, nz_),
        KOKKOS_CLASS_LAMBDA(const int k){
            int tk;
            tk = k;
            if (tk > nz_/2) tk = tk - nz_;
            kz_(k) = (float(tk) * twopi_) / (nz_ * dz_);
    });
#elif OUT_OF_PLACE_FFT
    Kokkos::parallel_for(
        Kokkos::RangePolicy<>(0, nz21_),
        KOKKOS_CLASS_LAMBDA(const int k){
            int tk;
            tk = k;
            kz_(k) = (float(tk) * twopi_) / (nz_ * dz_);
    }); 
#endif
}


CArrayKokkos<double>& FourierSpace::get_kx()
{
    return kx_;
}


CArrayKokkos<double>& FourierSpace::get_ky()
{
    return ky_;
}


CArrayKokkos<double>& FourierSpace::get_kz()
{
    return kz_;
}


