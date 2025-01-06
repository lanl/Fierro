#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "utilities.h"
#include "Profiler.h"
#include "math_functions.h"
#include "determinant33.h"
#include "defgrad_dcmp.h"

void EVPFFT::update_el_stiff()
{
  Profiler profiler(__FUNCTION__);

  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    real_t detFe;
    real_t dum;

    // thread private arrays
    real_t aa_[3*3];
    real_t aat_[3*3];
    real_t defgradetmp_[3*3];
    real_t caux3333_[3*3*3*3];
    real_t caux66_[6*6];
    real_t F_[3*3];
    real_t V_[3*3];
    real_t R_[3*3];

    // create views of thread private arrays
    ViewMatrixTypeReal aa(aa_,3,3);
    ViewMatrixTypeReal aat(aat_,3,3);
    ViewMatrixTypeReal defgradetmp(defgradetmp_,3,3);
    ViewMatrixTypeReal caux3333(caux3333_,3,3,3,3);
    ViewMatrixTypeReal caux66(caux66_,6,6);
    ViewMatrixTypeReal F(F_,3,3);
    ViewMatrixTypeReal V(V_,3,3);
    ViewMatrixTypeReal R(R_,3,3);

    int iph;

    iph = jphase(i,j,k);

    if (igas(iph) == 0) {

      for (int jj = 1; jj <= 3; jj++) {
        for (int ii = 1; ii <= 3; ii++) {
          F(ii,jj) = defgrade(ii,jj,i,j,k);
        }
      }
      defgrad_dcmp(F.pointer(), V.pointer(), R.pointer());
      
      for (int jj = 1; jj <= 3; jj++) {
        for (int ii = 1; ii <= 3; ii++) {
          aa(ii,jj) = 0.0;
          for (int kk = 1; kk <= 3; kk++) {
            // aa(ii,jj) += defgrade(ii,kk,i,j,k)*ag(kk,jj,i,j,k);
            aa(ii,jj) += R(ii,kk)*ag(kk,jj,i,j,k);
          }
        }
      }

      for (int jj = 1; jj <= 3; jj++) {
        for (int ii = 1; ii <= 3; ii++) {
          aat(ii,jj) = aa(jj,ii);
          defgradetmp(ii,jj) = defgrade(ii,jj,i,j,k);
        }
      }

      // determinant33(defgradetmp.pointer(), detFe);

      for (int i1 = 1; i1 <= 3; i1++) {
        for (int j1 = 1; j1 <= 3; j1++) {
          for (int k1 = 1; k1 <= 3; k1++) {
            for (int l1 = 1; l1 <= 3; l1++) {
            
              dum = 0.0;
              
              for (int i2 = 1; i2 <= 3; i2++) {
                for (int j2 = 1; j2 <= 3; j2++) {
                  for (int k2 = 1; k2 <= 3; k2++) {
                    for (int l2 = 1; l2 <= 3; l2++) {
                    
                      dum += aa(i1,i2) *
                             aa(j1,j2) *
                             cc(i2,j2,k2,l2,iph) * 
                             aat(k2,k1) * 
                             aat(l2,l1);
                             
                    }
                  }
                }
              }
                            
              // caux3333(i1,j1,k1,l1) = dum/detFe;
              caux3333(i1,j1,k1,l1) = dum;
             
            }
          }
        }
      }
            
      cb.chg_basis_4(caux66.pointer(), caux3333.pointer(), 4, 6, cb.B_basis_device_pointer());

      for (int ii = 1; ii <= 6; ii++) {
        for (int jj = 1; jj <= 6; jj++) {
          cg66(ii,jj,i,j,k) = caux66(ii,jj);
        }
      }

    } // end if (igas(iph) == 0) 

  }); // end FOR_ALL_CLASS

}
