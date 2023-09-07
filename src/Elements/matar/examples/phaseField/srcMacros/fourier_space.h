#pragma once

#include "matar.h"

using namespace mtr; // matar namespace

class FourierSpace
{
    private:
        int* nn_;
        int nx_;
        int ny_;
        int nz_;
#ifdef OUT_OF_PLACE_FFT
        int nz21_;
#endif
        double* delta_;
        double dx_;
        double dy_;
        double dz_;
        const double pi_ = 3.141592653589793238463;
        const double twopi_ = 2.0*pi_;
        CArrayKokkos<double> kx_; 
        CArrayKokkos<double> ky_;
        CArrayKokkos<double> kz_;  
  
    public:
        FourierSpace(int* nn, double* delta);
        CArrayKokkos<double>& get_kx();
        CArrayKokkos<double>& get_ky();
        CArrayKokkos<double>& get_kz();

        void set_kx_ky_kz_();
};
