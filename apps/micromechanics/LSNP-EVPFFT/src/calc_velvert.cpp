#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "math_functions.h"

void EVPFFT::calc_velvert()
{

  for (int inode = 1; inode <= 8; inode++) {
    for (int ii = 1; ii <= 3; ii++) {
      velvert(ii,inode) = 0.0;
      for (int jj = 1; jj <= 3; jj++) {
        velvert(ii,inode) += udot(ii,jj)*xvert_current(jj,inode);
      }
    }
  }

}
