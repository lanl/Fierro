#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "utilities.h"
#include "Profiler.h"
#include "math_functions.h"
#include "determinant33.h"
#include "defgrad_dcmp.h"

void EVPFFT::update_defgrade()
{
  Profiler profiler(__FUNCTION__);

  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    real_t detdefgradeinc;

    // thread private arrays
    real_t defgradeinvold_[3*3];
    real_t FpFiniinv_[3*3];
    real_t defgradeinc_[3*3];
    real_t Finc_[3*3];
    real_t Vinc_[3*3];
    real_t Rinc_[3*3];

    // create views of thread private arrays
    ViewMatrixTypeReal defgradeinvold(defgradeinvold_,3,3);
    ViewMatrixTypeReal FpFiniinv(FpFiniinv_,3,3);
    ViewMatrixTypeReal defgradeinc(defgradeinc_,3,3);
    ViewMatrixTypeReal Finc(Finc_,3,3);
    ViewMatrixTypeReal Vinc(Vinc_,3,3);
    ViewMatrixTypeReal Rinc(Rinc_,3,3);

    int iph;

    iph = jphase(i,j,k);

    if (igas(iph) == 0) {

      for (int jj = 1; jj <= 3; jj++) {
        for (int ii = 1; ii <= 3; ii++) {
          defgradeinvold(ii,jj) = defgrade(ii,jj,i,j,k);
        }
      }

      invert_matrix <3> (defgradeinvold.pointer());

      for (int jj = 1; jj <= 3; jj++) {
        for (int ii = 1; ii <= 3; ii++) {
          FpFiniinv(ii,jj) = 0.0;
          for (int kk = 1; kk <= 3; kk++) {
            FpFiniinv(ii,jj) += defgradp(ii,kk,i,j,k)*defgradini(kk,jj,i,j,k);
          }
        }
      }

      invert_matrix <3> (FpFiniinv.pointer());

      // update 
      for (int jj = 1; jj <= 3; jj++) {
        for (int ii = 1; ii <= 3; ii++) {
          defgrade(ii,jj,i,j,k) = 0.0;
          for (int kk = 1; kk <= 3; kk++) {
            defgrade(ii,jj,i,j,k)  += defgrad(ii,kk,i,j,k)*FpFiniinv(kk,jj);
          }
        }
      }

      for (int jj = 1; jj <= 3; jj++) {
        for (int ii = 1; ii <= 3; ii++) {
          defgradeinc(ii,jj) = 0.0;
          for (int kk = 1; kk <= 3; kk++) {
            defgradeinc(ii,jj)  += defgrade(ii,kk,i,j,k)*defgradeinvold(kk,jj);
          }
          Finc(ii,jj) = defgradeinc(ii,jj);
        }
      }

      defgrad_dcmp(Finc.pointer(), Vinc.pointer(), Rinc.pointer());
      
      // determinant33(defgradeinc.pointer(), detdefgradeinc);

      for (int jj = 1; jj <= 3; jj++) {
        for (int ii = 1; ii <= 3; ii++) {
          sgt(ii,jj,i,j,k) = 0.0;
          for (int kk = 1; kk <= 3; kk++) {
            for (int ll = 1; ll <= 3; ll++) {
              // sgt(ii,jj,i,j,k) += defgradeinc(ii,kk)*defgradeinc(jj,ll)*sg(kk,ll,i,j,k)/detdefgradeinc;
              sgt(ii,jj,i,j,k) += Rinc(ii,kk)*Rinc(jj,ll)*sg(kk,ll,i,j,k);
            }
          }
        }
      }

    } // end if (igas(iph) == 0)

  }); // end FOR_ALL_CLASS

}
