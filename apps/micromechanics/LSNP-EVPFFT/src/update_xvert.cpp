#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "math_functions.h"

void EVPFFT::update_xvert()
{

  for (int inode = 1; inode <= 8; inode++) {
    for (int ii = 1; ii <= 3; ii++) {
      xvert_current(ii,inode) += velvert(ii,inode)*tdot;
    }
  }

}
