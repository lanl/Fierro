#include "evpfft.h"
#include "utilities.h"
#include "Profiler.h"

void EVPFFT::kinhard_param()
{
  Profiler profiler(__FUNCTION__);

  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    int jph;

    // thread private arrays
    real_t sc5_[5];
    real_t sg5_[5];
    real_t sgx_[3*3];

    // create views of thread private arrays
    ViewMatrixTypeReal sc5(sc5_,5);
    ViewMatrixTypeReal sg5(sg5_,5);
    ViewMatrixTypeReal sgx(sgx_,3,3);

    jph = jphase(i,j,k);

    if (igas(jph) == 0) {

      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          sgx(ii,jj) = sg(ii,jj,i,j,k);
        }
      }

      cb.chg_basis_2(sg5.pointer(), sgx.pointer(), 2, 5, cb.B_basis_device_pointer());

      for (int is = 1; is <= nsyst(jph); is++) {
        for (int jj = 1; jj <= 5; jj++) {
           sc5(jj) = sch(jj,is,i,j,k);
        }

#ifdef NON_SCHMID_EFFECTS       
        for (int jj = 1; jj <= 5; jj++) {
          sc5(jj) = schnon(jj,is,i,j,k);
        }
#endif

        xkin(is,i,j,k) = sc5(1)*sg5(1) +
                            sc5(2)*sg5(2) +
                            sc5(3)*sg5(3) + 
                            sc5(4)*sg5(4) + 
                            sc5(5)*sg5(5);
      } // end for is

    } //end if (igas(jph) == 0)

  }); // end FOR_ALL_CLASS

}
