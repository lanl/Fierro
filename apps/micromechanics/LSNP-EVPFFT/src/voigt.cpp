#include "voigt.h"

void voigt(real_t *c2_voigt_, real_t *c4_voigt_, int iopt)
{

  ViewMatrixTypeReal c2_voigt(c2_voigt_, 6,6);
  ViewMatrixTypeReal c4_voigt(c4_voigt_, 3,3,3,3);

  int ijv_[6*2] = {1,2,3,2,1,1,1,2,3,3,3,2};
  real_t f_[6*6];
  ViewMatrixTypeInt  ijv(ijv_,6,2);
  ViewMatrixTypeReal f(f_, 6,6);
  int i1;
  int i2;
  int j1;
  int j2;
  
  //
  if(iopt == 1) {
    for (int i = 1; i <= 6; i++) {
      i1 = ijv(i,1);
      i2 = ijv(i,2);
      for (int j = 1; j <= 6; j++) {
        j1 = ijv(j,1);
        j2 = ijv(j,2);
        c4_voigt(i1,i2,j1,j2) = c2_voigt(i,j);
        c4_voigt(i2,i1,j1,j2) = c2_voigt(i,j);
        c4_voigt(i1,i2,j2,j1) = c2_voigt(i,j);
        c4_voigt(i2,i1,j2,j1) = c2_voigt(i,j);
      }
    }
  }  // end if (iopt == 1)
  
  //
  if (iopt == 2) {
    for (int i = 1; i <= 6; i++) {
      i1 = ijv(i,1);
      i2 = ijv(i,2);
      for (int j = 1; j <= 6; j++) {
        j1 = ijv(j,1);
        j2 = ijv(j,2);
        c2_voigt(i,j) = c4_voigt(i1,i2,j1,j2);
      }
    }
  } // end if (iopt == 2)
  
  //
  if (iopt == 3) {
    for (int i = 1; i <= 6; i++) {
      for (int j = 1; j <= 6; j++) {
        f(i,j) = 1.0;
        if (i > 3) f(i,j) = 0.5;
        if (j > 3) f(i,j) = 0.5*f(i,j);
      }
    }
    
    for (int i = 1; i <= 6; i++) {
      i1 = ijv(i,1);
      i2 = ijv(i,2);
      for (int j = 1; j <= 6; j++) {
        j1 = ijv(j,1);
        j2 = ijv(j,2);
        c4_voigt(i1,i2,j1,j2) = f(i,j)*c2_voigt(i,j);
        c4_voigt(i2,i1,j1,j2) = f(i,j)*c2_voigt(i,j);
        c4_voigt(i1,i2,j2,j1) = f(i,j)*c2_voigt(i,j);
        c4_voigt(i2,i1,j2,j1) = f(i,j)*c2_voigt(i,j);
      }
    }
  } // end if (iopt == 3)
  
  //
  if (iopt == 4) {
    for (int i = 1; i <= 6; i++) {
      i1 = ijv(i,1);
      i2 = ijv(i,2);
      for (int j = 1; j <= 6; j++) {
        j1 = ijv(j,1);
        j2 = ijv(j,2);
        if (i < 3) {
          c4_voigt(i1,i2,j1,j2) = c2_voigt(i,j);
          c4_voigt(i2,i1,j1,j2) = c2_voigt(i,j);
          c4_voigt(i1,i2,j2,j1) = c2_voigt(i,j);
          c4_voigt(i2,i1,j2,j1) = c2_voigt(i,j);
        } else {
          c4_voigt(i1,i2,j1,j2) =  c2_voigt(i,j);
          c4_voigt(i2,i1,j1,j2) = -c2_voigt(i,j);
          c4_voigt(i1,i2,j2,j1) =  c2_voigt(i,j);
          c4_voigt(i2,i1,j2,j1) = -c2_voigt(i,j);
        } // end if (i < 3)
      } // end for j
    } // end for i
  } // end if (iopt == 4)
  
}
