#include "evpfft.h"
#include "math_functions.h"
#include "utilities.h"
#include "Profiler.h"
#include "mod_frequency.h"

void EVPFFT::inverse_the_greens()
{
  Profiler profiler(__FUNCTION__);

  FOR_ALL_CLASS(kzz, 1, npts3_cmplx+1,
                kyy, 1, npts2_cmplx+1,
                kxx, 1, npts1_cmplx+1, {

    real_t xknorm;
    real_t prefactor_re;
    real_t prefactor_im;

    // thread private arrays
    real_t xk_[3];
    real_t krot_re_[3];
    real_t krot_im_[3];
    real_t ddisgrad_[3*3];
    real_t ddisgradim_[3*3];
    real_t a_[3*3];
    real_t D_dft_[3*3];
    real_t g1_[3*3*3*3];

    // create views of thread private arrays
    ViewMatrixTypeReal xk(xk_,3);
    ViewMatrixTypeReal krot_re(krot_re_,3);
    ViewMatrixTypeReal krot_im(krot_im_,3);
    ViewMatrixTypeReal ddisgrad(ddisgrad_,3,3);
    ViewMatrixTypeReal ddisgradim(ddisgradim_,3,3);
    ViewMatrixTypeReal a(a_,3,3);
    ViewMatrixTypeReal D_dft(D_dft_,3,3);
    ViewMatrixTypeReal g1(g1_,3,3,3,3);

    xk(1) = xk_gb(kxx);
    xk(2) = yk_gb(kyy);
    xk(3) = zk_gb(kzz);

    if (igamma == 0) {

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

    } else if (igamma == 1) { 

      // TODO: this can be precomputed and stored in D_dft array (not sure if memory is issue)
      xk(1) =  xk(1)*2.0*pi;
      xk(2) =  xk(2)*2.0*pi;
      xk(3) =  xk(3)*2.0*pi;

      mod_frequency(xk.pointer(), krot_re.pointer(), krot_im.pointer());

      if (abs(POW2(krot_re(1)) + POW2(krot_im(1)) + POW2(krot_re(2)) + POW2(krot_im(2)) + 
          POW2(krot_re(3)) + POW2(krot_im(3))) < 1.0e-15) {

        for (int l = 1; l <= 3; l++) {
          for (int k = 1; k <= 3; k++) {
            for (int j = 1; j <= 3; j++) {
              for (int i = 1; i <= 3; i++) {
                g1(i,j,k,l) = 0.0;
              }
            }
          }
        }
      } else {

        for (int i = 1; i <= 3; i++) {
          for (int j = 1; j <= 3; j++) {
            D_dft(i,j) = krot_re(i)*krot_re(j) + krot_im(i)*krot_im(j);
          }
        }

        // TODO: optimize indexing of this loop
        for (int i = 1; i <= 3; i++) {
          for (int k = 1; k <= 3; k++) {
            a(i,k) = 0.0;
            for (int j = 1; j <= 3; j++) {
              for (int l = 1; l <= 3; l++) {
                a(i,k) += c0(i,j,k,l)*D_dft(j,l);
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
                g1(p,qq,i,j) = -a(p,i)*D_dft(qq,j);
              }
            }
          }
        }

      }

    }// end if (igamma)


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
