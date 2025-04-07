#include "mod_frequency.h"

KOKKOS_FUNCTION
void mod_frequency(real_t *xk_, real_t *kmod_re_, real_t *kmod_im_)
{

  ViewMatrixTypeReal xk(xk_, 3);
  ViewMatrixTypeReal kmod_re(kmod_re_, 3);
  ViewMatrixTypeReal kmod_im(kmod_im_, 3);

  real_t prefactor_re;
  real_t prefactor_im;

  prefactor_re = -0.25*sin(xk(1)) - 0.25*cos(xk(2))*sin(xk(1)) - 0.25*cos(xk(3))*sin(xk(1)) -
    0.25*cos(xk(2))*cos(xk(3))*sin(xk(1)) - 0.25*sin(xk(2)) - 0.25*cos(xk(1))*sin(xk(2)) -
    0.25*cos(xk(3))*sin(xk(2)) - 0.25*cos(xk(1))*cos(xk(3))*sin(xk(2)) - 0.25*sin(xk(3)) -
    0.25*cos(xk(1))*sin(xk(3)) - 0.25*cos(xk(2))*sin(xk(3)) -
    0.25*cos(xk(1))*cos(xk(2))*sin(xk(3)) + 0.25*sin(xk(1))*sin(xk(2))*sin(xk(3));
  prefactor_im = 0.25 + 0.25*cos(xk(1)) + 0.25*cos(xk(2)) + 0.25*cos(xk(1))*cos(xk(2)) + 
    0.25*cos(xk(3)) + 0.25*cos(xk(1))*cos(xk(3)) + 0.25*cos(xk(2))*cos(xk(3)) + 
    0.25*cos(xk(1))*cos(xk(2))*cos(xk(3)) - 0.25*sin(xk(1))*sin(xk(2)) - 
    0.25*cos(xk(3))*sin(xk(1))*sin(xk(2)) - 0.25*sin(xk(1))*sin(xk(3)) - 
    0.25*cos(xk(2))*sin(xk(1))*sin(xk(3)) - 0.25*sin(xk(2))*sin(xk(3)) - 
    0.25*cos(xk(1))*sin(xk(2))*sin(xk(3));
  kmod_re(1) = tan(xk(1)/2.0)*prefactor_re;
  kmod_re(2) = tan(xk(2)/2.0)*prefactor_re;
  kmod_re(3) = tan(xk(3)/2.0)*prefactor_re;
  kmod_im(1) = tan(xk(1)/2.0)*prefactor_im;
  kmod_im(2) = tan(xk(2)/2.0)*prefactor_im;
  kmod_im(3) = tan(xk(3)/2.0)*prefactor_im;

  
}
