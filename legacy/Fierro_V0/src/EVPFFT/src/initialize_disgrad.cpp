#include "evpfft.h"
#include "utilities.h"
#include "Profiler.h"

void EVPFFT::initialize_disgrad()
{
  Profiler profiler(__FUNCTION__);

  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        disgrad(ii,jj,i,j,k) += ddisgradmacro(ii,jj) + 
                                    work(ii,jj,i,j,k);
      }
    }
  }); // end FOR_ALL_CLASS

}
