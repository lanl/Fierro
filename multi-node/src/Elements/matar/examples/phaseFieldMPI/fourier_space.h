#pragma once

#include "matar.h"
#include <string>

using namespace mtr; // matar namespace

class FourierSpace
{
public:
const double pi = 3.141592653589793238463;
const double twopi = 2.0*pi;
CArrayKokkos<double> kx; 
CArrayKokkos<double> ky;
CArrayKokkos<double> kz; 
  
FourierSpace(const std::array<int,3> & glob_nn_real, 
             const std::array<int,3> & loc_nn_cmplx, 
             const std::array<int,3> & loc_start_index, 
             const std::array<double,3> & delta);
void set_kx_ky_kz(const std::array<int,3> & glob_nn_real, 
                  const std::array<int,3> & loc_nn_cmplx, 
                  const std::array<int,3> & loc_start_index, 
                  const std::array<double,3> & delta);
};
