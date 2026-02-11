#include "orient.h"

KOKKOS_FUNCTION
void orient(real_t *a_, real_t *c_)
{
  ViewMatrixTypeReal a(a_,3,3);
  ViewMatrixTypeReal c(c_,3,3);

  real_t snorm;
  real_t snorm1;

  // thread private arrays
  real_t th2_[3*3];
  real_t v_[3];
  real_t vbar_[3];
  real_t th_[3*3];
  real_t rot_[3*3];
  real_t anew_[3*3];

  // create views of thread private arrays
  ViewMatrixTypeReal th2(th2_,3,3);
  ViewMatrixTypeReal v(v_,3);
  ViewMatrixTypeReal vbar(vbar_,3);
  ViewMatrixTypeReal th(th_,3,3);
  ViewMatrixTypeReal rot(rot_,3,3);
  ViewMatrixTypeReal anew(anew_,3,3);
  

//     BUILD ROTATION TENSOR BASED ON RODRIGUES FORMULA

  v(1) = c(3,2);
  v(2) = c(1,3);
  v(3) = c(2,1);
  snorm = SQRT(v(1)*v(1)+v(2)*v(2)+v(3)*v(3));
  snorm1 = TAN(snorm/2.0);
  if (snorm > 1.0E-06) goto label_97;
  snorm = 1.0;

label_97 : {
  for (int i = 1; i <= 3; i++) {
    vbar(i) = snorm1*v(i)/snorm;
  }
} // end label_97

  snorm = vbar(1)*vbar(1)+vbar(2)*vbar(2)+vbar(3)*vbar(3);
  th(3,2) =  vbar(1);
  th(1,3) =  vbar(2);
  th(2,1) =  vbar(3);
  th(2,3) = -vbar(1);
  th(3,1) = -vbar(2);
  th(1,2) = -vbar(3);

  for (int i = 1; i <= 3; i++) {
    th(i,i) = 0.0;
  }

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      th2(i,j) = 0.0;
      for (int k = 1; k <= 3; k++) {
        th2(i,j) += th(i,k)*th(k,j);
      }
    }
  }

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      rot(i,j) = (i/j)*(j/i)+2.0*(th(i,j)+th2(i,j))/(1.0+snorm);
    }
  }

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      anew(i,j) = 0.0;
      for (int k = 1; k <= 3; k++) {
        anew(i,j) += rot(i,k)*a(k,j);
      }
    }
  }

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      a(i,j) = anew(i,j);
    }
  }

}
