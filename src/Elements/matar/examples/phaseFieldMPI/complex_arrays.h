#pragma once

#include "sim_parameters.h"
#include "matar.h"
#include "fourier_space.h"
#include "heffte_fft.h"

using namespace mtr; // matar namespace

class ComplexArrays
{
public:
// arrays needed by solver 
DCArrayKokkos<double> comp_img;
DCArrayKokkos<double> dfdc_img;
DCArrayKokkos<double> kpow2;
CArrayKokkos<double> denominator;
FourierSpace fs;

ComplexArrays(const SimParameters & sp, const std::array<int,3> & nn_img, const std::array<int,3> & start_index);
void set_kpow2();
void set_denominator(const SimParameters & sp);
};
