#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "utilities.h"
#include "Profiler.h"
#include "matrix_exp.h"
#include "math_functions.h"

void EVPFFT::update_defgradp()
{
  Profiler profiler(__FUNCTION__);

  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    // thread private arrays
    real_t aa_[3*3];
    real_t velgradp_[3*3];
    real_t disgradincp_[3*3];
    real_t expdisgradincp_[3*3];
    real_t defgradpnew_[3*3];
    real_t dnsa_[3];
    real_t dbsa_[3];
    real_t defgradeinv_[3*3];

    // create views of thread private arrays
    ViewMatrixTypeReal aa(aa_,3,3);
    ViewMatrixTypeReal velgradp(velgradp_,3,3);
    ViewMatrixTypeReal disgradincp(disgradincp_,3,3);
    ViewMatrixTypeReal expdisgradincp(expdisgradincp_,3,3);
    ViewMatrixTypeReal defgradpnew(defgradpnew_,3,3);
    ViewMatrixTypeReal dnsa(dnsa_,3);
    ViewMatrixTypeReal dbsa(dbsa_,3);
    ViewMatrixTypeReal defgradeinv(defgradeinv_,3,3);

    int iph;

    iph = jphase(i,j,k);

    if (igas(iph) == 0) {

      if (iJ2(iph) == 0) { 

        for (int jj = 1; jj <= 3; jj++) {
          for (int ii = 1; ii <= 3; ii++) {
            aa(ii,jj) = ag(ii,jj,i,j,k);
            velgradp(ii,jj) = 0.0;
          }
        }
  
        // plastic velocity gradient in intermediate configuration
        for (int is = 1; is <= nsyst(iph); is++) {
          for (int ii = 1; ii <= 3; ii++) {
            dnsa(ii) = 0.0;
            dbsa(ii) = 0.0;
            for (int jj = 1; jj <= 3; jj++) {
              dnsa(ii) += aa(ii,jj) * dnca(jj,is,iph);
              dbsa(ii) += aa(ii,jj) * dbca(jj,is,iph);
            }
          }
  
          for (int ii = 1; ii <= 3; ii++) {
            for (int jj = 1; jj <= 3; jj++) {
              velgradp(ii,jj) += dbsa(ii) * dnsa(jj) * gamdot(is,i,j,k);
            }
          }
        } // end for is

      } else { // else for if (iJ2(jph) == 0)

        for (int ii = 1; ii <= 3; ii++) {
          for (int jj = 1; jj <= 3; jj++) {
            defgradeinv(ii,jj) = defgrade(ii,jj,i,j,k);
          }
        }
        invert_matrix <3> (defgradeinv.pointer());

        for (int ii = 1; ii <= 3; ii++) {
          for (int jj = 1; jj <= 3; jj++) {
            velgradp(ii,jj) = 0.0;
            for (int kk = 1; kk <= 3; kk++) {
              for (int ll = 1; ll <= 3; ll++) {
                velgradp(ii,jj) += defgradeinv(ii,kk)*edotp(kk,ll,i,j,k)*defgrade(ll,jj,i,j,k);
              }
            }
          }
        }

      } // end if (iJ2(jph) == 0)

      // increment
      for (int jj = 1; jj <= 3; jj++) {
        for (int ii = 1; ii <= 3; ii++) {
          disgradincp(ii,jj) = velgradp(ii,jj)*tdot;
        }
      }

      // matrix exponential of increment
      matrix_exp(disgradincp.pointer(), expdisgradincp.pointer());

      // update 
      for (int jj = 1; jj <= 3; jj++) {
        for (int ii = 1; ii <= 3; ii++) {
          defgradpnew(ii,jj) = 0.0;
          for (int kk = 1; kk <= 3; kk++) {
            defgradpnew(ii,jj) += expdisgradincp(ii,kk)*defgradp(kk,jj,i,j,k);
          }
        }
      }

      for (int jj = 1; jj <= 3; jj++) {
        for (int ii = 1; ii <= 3; ii++) {
          defgradp(ii,jj,i,j,k) = defgradpnew(ii,jj);
        }
      }


    } // end if (igas(iph) == 0) 

  }); // end FOR_ALL_CLASS

}
