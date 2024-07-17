#include "euler.h"

KOKKOS_FUNCTION
void euler(int iopt, real_t &ph, real_t &th, real_t &tm, real_t *a_)
{
  ViewMatrixTypeReal a(a_,3,3);
  real_t sth,sph,cph,cth,stm,ctm; 
  const real_t pi = 4.*ATAN(1.0);

  //  CALCULATE THE EULER ANGLES ASSOCIATED WITH THE TRANSFORMATION
  //  MATRIX A(I,J) IF IOPT=1 AND VICEVERSA IF IOPT=2
  //  A(i,j) TRANSFORMS FROM SYSTEM sa TO SYSTEM ca.
  //  ph,th,om ARE THE EULER ANGLES OF ca REFERRED TO sa.

  if (iopt == 1) {
    th = ACOS(a(3,3));
    if(ABS(a(3,3)) >= 0.9999) {
      tm = 0.0;
      ph = ATAN2(a(1,2),a(1,1));
    } else {
      sth = SIN(th);
      tm  = ATAN2(a(1,3)/sth,a(2,3)/sth);
      ph  = ATAN2(a(3,1)/sth,-a(3,2)/sth);
    }
    th = th*180.0/pi;
    ph = ph*180.0/pi;
    tm = tm*180.0/pi;
  } else if (iopt == 2) {
    sph = SIN(ph);
    cph = COS(ph);
    sth = SIN(th);
    cth = COS(th);
    stm = SIN(tm);
    ctm = COS(tm);
    a(1,1) =  ctm*cph-sph*stm*cth;
    a(2,1) = -stm*cph-sph*ctm*cth;
    a(3,1) =  sph*sth;
    a(1,2) =  ctm*sph+cph*stm*cth;
    a(2,2) = -sph*stm+cph*ctm*cth;
    a(3,2) = -sth*cph;
    a(1,3) =  sth*stm;
    a(2,3) =  ctm*sth;
    a(3,3) =  cth;
  } // end if (iopt == 1)

}
