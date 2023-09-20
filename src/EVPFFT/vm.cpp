#include "vm.h"


real_t vm(real_t *dtensor_)
{
  /* Returens VonMises strain */

  ViewMatrixTypeReal dtensor(dtensor_,3,3);
  MatrixTypeRealHost dt(3,3);

  real_t trace;
  real_t result;

  trace = dtensor(1,1) + dtensor(2,2) + dtensor(3,3);

  for (int j = 1; j <= 3; j++) {
    for (int i = 1; i <= 3; i++) {
      dt(i,j) = (dtensor(i,j)+dtensor(j,i))/2.0 - 
                (i/j)*(j/i)*trace/3.0;
    }
  }

  result = 0.0;
  for (int j = 1; j <= 3; j++) {
    for (int i = 1; i <= 3; i++) {
      result += POW2(dt(i,j));
    }
  }
  result = sqrt(2.0/3.0*result);

  return result;
}

real_t vm_stress(real_t *stress) {
    /* Returns VonMises stress */

    real_t von_mises = sqrtf((((stress[0] - stress[4])*(stress[0] - stress[4])) 
                       +        ((stress[4] - stress[8])*(stress[4] - stress[8]))
                       +        ((stress[8] - stress[0])*(stress[8] - stress[0]))
                       +   6*(pow(stress[1], 2.0)+pow(stress[5],2.0)+pow(stress[2],2.0)))/2.0);
    return von_mises;
}

