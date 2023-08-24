#pragma once

#include "sim_parameters.h"
#include "matar.h"

using namespace mtr; // matar namespace

class CHFourierSpectralSolver
{
    private:
        // simulation parameters
        int*    nn_;
        int     nn_img_[3];
        int     nx_;
        int     ny_;
        int     nz_;
#ifdef OUT_OF_PLACE_FFT
        int nz21_;
#endif
        int     ndim_;
        double* delta_;
        double  dx_;
        double  dy_;
        double  dz_;
        double  dt_;
        double  M_;
        double  kappa_;
        
        // arrays needed by solver 
        CArrayKokkos<double> comp_img_;
        CArrayKokkos<double> dfdc_img_;
        CArrayKokkos<double> kpow2_;
        CArrayKokkos<double> denominator_;
   
    public:
        CHFourierSpectralSolver(SimParameters &sp);
        void set_kpow2_();
        void set_denominator_();
        void time_march(DCArrayKokkos<double> &comp, CArrayKokkos<double> &dfdc);
};
