#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "utilities.h"

void EVPFFT::update_schmid()
{

  FOR_ALL_CLASS(kk, 1, npts3+1,
                jj, 1, npts2+1,
                ii, 1, npts1+1, {

    for (int is = 1; is <= NSYSMX; is++) {
      for (int j = 1; j <= 5; j++) {
        sch(j,is,ii,jj,kk)    = 0.0;
#ifdef NON_SCHMID_EFFECTS
        schnon(j,is,ii,jj,kk) = 0.0;
#endif
      } // end for j
    } // end for is
  });


  FOR_ALL_CLASS(kk, 1, npts3+1,
                jj, 1, npts2+1,
                ii, 1, npts1+1, {

    // thread private arrays
    real_t aux5_[5];
    real_t aux33_[3*3]; 
    real_t aux33r_[3*3];

    // create views of thread private arrays
    ViewMatrixTypeReal aux5    (aux5_,5);
    ViewMatrixTypeReal aux33   (aux33_,3,3);
    ViewMatrixTypeReal aux33r  (aux33r_,3,3);
    

    int jph;
    jph = jphase(ii,jj,kk);

    if (igas(jph) == 0) {
      for (int is = 1; is <= nsyst(jph); is++) {

        for (int j = 1; j <= 5; j++) {
          aux5(j) = schca(j,is,jph);
        } // end for j

        cb.chg_basis_1(aux5.pointer(), aux33.pointer(), 1, 5, cb.B_basis_device_pointer());

        for (int i = 1; i <= 3; i++) {
          for (int j = 1; j <= 3; j++) {
            aux33r(i,j) = 0.0;
            for (int i1 = 1; i1 <= 3; i1++) {
              for (int j1 = 1; j1 <= 3; j1++) {
                aux33r(i,j) += ag(i,i1,ii,jj,kk)*
                               ag(j,j1,ii,jj,kk)*aux33(i1,j1);
              } // end for j1
            } // end for i1
          } // end for j
        } // end for i

        cb.chg_basis_2(aux5.pointer(), aux33r.pointer(), 2, 5, cb.B_basis_device_pointer());

        for (int j = 1; j <= 5; j++) {
          sch(j,is,ii,jj,kk) = aux5(j);
        } // end for j

#ifdef NON_SCHMID_EFFECTS
          for (int j = 1; j <= 5; j++) {
            aux5(j) = schcnon(j,is,jph);
          } // end for j

          cb.chg_basis_1(aux5.pointer(), aux33.pointer(), 1, 5, cb.B_basis_device_pointer());

          for (int i = 1; i <= 3; i++) {
            for (int j = 1; j <= 3; j++) {
              aux33r(i,j) = 0.0;
              for (int i1 = 1; i1 <= 3; i1++) {
                for (int j1 = 1; j1 <= 3; j1++) {
                  aux33r(i,j) += ag(i,i1,ii,jj,kk)*
                                 ag(j,j1,ii,jj,kk)*aux33(i1,j1);
                } // end for j1
              } // end for i1
            } // end for j
          } // end for i

          cb.chg_basis_2(aux5.pointer(), aux33r.pointer(), 2, 5, cb.B_basis_device_pointer());
          
          for (int j = 1; j <= 5; j++) {
            schnon(j,is,ii,jj,kk) = aux5(j);
          } // end for j
#endif

      } // end for is 
    }  // end if (igas(jph) == 0)

  }); // end FOR_ALL_CLASS

}
