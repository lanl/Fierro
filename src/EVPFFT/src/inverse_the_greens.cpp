#include "evpfft.h"
#include "math_functions.h"
#include "utilities.h"
#include "Profiler.h"

void EVPFFT::inverse_the_greens()
{
  Profiler profiler(__FUNCTION__);

  FOR_ALL_CLASS(kzz, 1, npts3_cmplx+1,
                kyy, 1, npts2_cmplx+1,
                kxx, 1, npts1_cmplx+1, {

    real_t xknorm;

    // thread private arrays
    real_t xk_[3];
    real_t ddisgrad_[3*3];
    real_t ddisgradim_[3*3];
    real_t a_[3*3];
    real_t g1_[3*3*3*3];

    // create views of thread private arrays
    ViewMatrixTypeReal xk(xk_,3);
    ViewMatrixTypeReal ddisgrad(ddisgrad_,3,3);
    ViewMatrixTypeReal ddisgradim(ddisgradim_,3,3);
    ViewMatrixTypeReal a(a_,3,3);
    ViewMatrixTypeReal g1(g1_,3,3,3,3);

    xk(1) = xk_gb(kxx);
    xk(2) = yk_gb(kyy);
    xk(3) = zk_gb(kzz);

    xknorm = sqrt( xk(1)*xk(1) + 
                   xk(2)*xk(2) + 
                   xk(3)*xk(3) );

    if (xknorm != 0.0) {
      for (int i = 1; i <= 3; i++) {
        xk(i) = xk(i) / xknorm;
      } // end for i
    } // end if (xknorm != 0.0)

    if ( kxx + local_start1_cmplx == npts1_g/2+1 || 
         kyy + local_start2_cmplx == npts2_g/2+1 || 
         (npts3_g > 1 && kzz + local_start3_cmplx == npts3_g/2+1) ) {
      for (int l = 1; l <= 3; l++) {
        for (int k = 1; k <= 3; k++) {
          for (int j = 1; j <= 3; j++) {
            for (int i = 1; i <= 3; i++) {
              g1(i,j,k,l) = -s0(i,j,k,l);
            }
          }
        }
      }
    } else {
      // TODO: optimize indexing of this loop
      for (int i = 1; i <= 3; i++) {
        for (int k = 1; k <= 3; k++) {
          a(i,k) = 0.0;
          for (int j = 1; j <= 3; j++) {
            for (int l = 1; l <= 3; l++) {
              a(i,k) += c0(i,j,k,l)*xk(j)*xk(l);
            }
          }
        }
      }

      invert_matrix <3> (a.pointer());

      // TODO: optimize indexing of this loop 
      for (int p = 1; p <= 3; p++) {
        for (int qq = 1; qq <= 3; qq++) {
          for (int i = 1; i <= 3; i++) {
            for (int j = 1; j <= 3; j++) {
              g1(p,qq,i,j) = -a(p,i)*xk(qq)*xk(j);
            }
          }
        }
      }

    } // end if ( kxx == npts1/2 || kyy == npts2/2 || (npts3 > 1 && kzz == npts3/2) )

    if ( kxx + local_start1_cmplx == 1 &&
         kyy + local_start2_cmplx == 1 && 
         kzz + local_start3_cmplx == 1 ) {
      for (int j = 1; j <= 3; j++) {
        for (int i = 1; i <= 3; i++) {
          ddisgrad(i,j)   = 0.0;
          ddisgradim(i,j) = 0.0;
        }
      }
    } else {
      // TODO: optimize indexing of this loop
      for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
          ddisgrad(i,j)   = 0.0;
          ddisgradim(i,j) = 0.0;
          for (int k = 1; k <= 3; k++) {
            for (int l = 1; l <= 3; l++) {
              ddisgrad(i,j)   += g1(i,j,k,l) * work(k,l,kxx,kyy,kzz);
              ddisgradim(i,j) += g1(i,j,k,l) * workim(k,l,kxx,kyy,kzz);
            }
          }
        }
      }
    } // end if ( kxx == 1 && kyy == 1 && kzz == 1 )

    for (int j = 1; j <= 3; j++) {
      for (int i = 1; i <= 3; i++) { 
        work(i,j,kxx,kyy,kzz)   = ddisgrad(i,j);
        workim(i,j,kxx,kyy,kzz) = ddisgradim(i,j);
      }
    }

  }); // end FOR_ALL_CLASS

}
