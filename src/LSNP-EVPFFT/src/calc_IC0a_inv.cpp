#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "utilities.h"
#include "math_functions.h"


void EVPFFT::calc_IC0a_inv()
{

  real_t s66_[6*6];
  ViewMatrixTypeReal s66(s66_,6,6);

  for (int ii = 1; ii <= 6; ii++) {
    for (int jj = 1; jj <= 6; jj++) {
      s66(ii,jj) = dF6_dP6.host(ii,jj)/frame_c;
    }
  }

  real_t IC0a66_inv_[6*6];
  ViewMatrixTypeReal IC0a66_inv(IC0a66_inv_,6,6);

  for (int ii = 1; ii <= 6; ii++) {
    for (int jj = 1; jj <= 6; jj++) {
      IC0a66_inv(ii,jj) = 0.0;
      for (int kk = 1; kk <= 6; kk++) {
        IC0a66_inv(ii,jj) += c066(ii,kk)*s66(kk,jj);
      }
    }
  }

  for (int ii = 1; ii <= 6; ii++) {
    IC0a66_inv(ii,ii) += 1.0;
  }

  invert_matrix <6> (IC0a66_inv.pointer());
  cb.chg_basis_3(IC0a66_inv.pointer(), IC0a_inv.host_pointer(), 3, 6, cb.B_basis_host_pointer());

  real_t s_frame_[3*3*3*3];
  ViewMatrixTypeReal s_frame(s_frame_,3,3,3,3);
  cb.chg_basis_3(s66.pointer(), s_frame.pointer(), 3, 6, cb.B_basis_host_pointer());

  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    if (iframe(i,j,k) == 1) {
      real_t dum;
      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          dum = 0.0;
          for (int kk = 1; kk <= 3; kk++) {
            for (int ll = 1; ll <= 3; ll++) {
              dum += s_frame(ii,jj,kk,ll)*sgPK1(kk,ll,i,j,k);
            }
          }
          Cinv_sgPK1old(ii,jj,i,j,k) = dum;
        }
      }
    }
  });

  IC0a_inv.update_device();
  Cinv_sgPK1old.update_device();

}
