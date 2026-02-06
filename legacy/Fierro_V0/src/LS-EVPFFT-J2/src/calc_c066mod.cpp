#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "utilities.h"
#include "Profiler.h"
#include "math_functions.h"
#include "determinant33.h"

void EVPFFT::calc_c066mod()
{
  Profiler profiler(__FUNCTION__);

  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    real_t dum;

    // thread private arrays
    real_t c0mod_[3*3*3*3];
    real_t c0modsym_[3*3*3*3];
    real_t c066modtmp_[6*6];

    // create views of thread private arrays
    ViewMatrixTypeReal c0mod(c0mod_,3,3,3,3);
    ViewMatrixTypeReal c0modsym(c0modsym_,3,3,3,3);
    ViewMatrixTypeReal c066modtmp(c066modtmp_,6,6);


    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        for (int kk = 1; kk <= 3; kk++) {
          for (int ll = 1; ll <= 3; ll++) {
            dum = 0.0;
            for (int m = 1; m <= 3; m++) {
              for (int n = 1; n <= 3; n++) {
                dum += c0(ii,m,kk,n)*defgrad(jj,m,i,j,k)*defgrad(ll,n,i,j,k);
              }
            }
            c0mod(ii,jj,kk,ll) = dum/detF(i,j,k);
          }
        }
      }
    }

    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        for (int kk = 1; kk <= 3; kk++) {
          for (int ll = 1; ll <= 3; ll++) {
            c0modsym(ii,jj,kk,ll) = 0.25*(c0mod(ii,jj,kk,ll) + c0mod(jj,ii,kk,ll) +
                                  c0mod(ii,jj,ll,kk) + c0mod(jj,ii,ll,kk));
          }
        }
      }
    }

    cb.chg_basis_4(c066modtmp.pointer(), c0modsym.pointer(), 4, 6, cb.B_basis_device_pointer());

    for (int ii = 1; ii <= 6; ii++) {
      for (int jj = 1; jj <= 6; jj++) {
        c066mod(ii,jj,i,j,k) = c066modtmp(ii,jj);
      }
    }

  }); // end FOR_ALL_CLASS

}
