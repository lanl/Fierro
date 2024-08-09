#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "utilities.h"

void EVPFFT::update_grid()
{

  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    x_grid(1,i,j,k) += velocity(1,i,j,k)*tdot;
    x_grid(2,i,j,k) += velocity(2,i,j,k)*tdot;
    x_grid(3,i,j,k) += velocity(3,i,j,k)*tdot;

  }); // end FOR_ALL_CLASS

  // TODO: pass to xnode array for igamma = 0 and 1

}
