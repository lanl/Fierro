#include "evpfft.h"
#include "math_functions.h"
#include "utilities.h"
#include "Profiler.h"
#include "reduction_data_structures.h"
#include "mod_frequency.h"

void EVPFFT::calc_Goperr0()
{
  Profiler profiler(__FUNCTION__);

  const size_t n = 9;
  const real_t twopi = 8.*ATAN(1.0);
  ArrayType <real_t, n> all_reduce;

  FOR_ALL_CLASS(kzz, 1, npts3_cmplx+1,
                kyy, 1, npts2_cmplx+1,
                kxx, 1, npts1_cmplx+1, {

    real_t xknorm;

    // thread private arrays
    real_t xk_[3];
    real_t a_[3*3];
    real_t krot_re_[3];
    real_t krot_im_[3];
    real_t D_dft_[3*3];

    // create views of thread private arrays
    ViewMatrixTypeReal xk(xk_,3);
    ViewMatrixTypeReal a(a_,3,3);
    ViewMatrixTypeReal krot_re(krot_re_,3);
    ViewMatrixTypeReal krot_im(krot_im_,3);
    ViewMatrixTypeReal D_dft(D_dft_,3,3);

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
  
      if ( (kxx + local_start1_cmplx == npts1_g/2+1 && npts1_g % 2 == 0) || 
           (kyy + local_start2_cmplx == npts2_g/2+1 && npts2_g % 2 == 0) || 
           (npts3_g > 1 && (kzz + local_start3_cmplx == npts3_g/2+1 && npts3_g % 2 == 0)) ) {
         for (int j = 1; j <= 3; j++) {
           for (int i = 1; i <= 3; i++) {
             a(i,j) = 0.0;
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
  
      } // end if ( kxx == npts1/2 || kyy == npts2/2 || (npts3 > 1 && kzz == npts3/2) )

    } else if (igamma == 1) { 

      // TODO: this can be precomputed and stored in D_dft array (not sure if memory is issue)
      xk(1) =  xk(1)*twopi;
      xk(2) =  xk(2)*twopi;
      xk(3) =  xk(3)*twopi;

      mod_frequency(xk.pointer(), krot_re.pointer(), krot_im.pointer());

      if (abs(POW2(krot_re(1)) + POW2(krot_im(1)) + POW2(krot_re(2)) + POW2(krot_im(2)) + 
          POW2(krot_re(3)) + POW2(krot_im(3))) < 1.0e-15) {
        for (int j = 1; j <= 3; j++) {
          for (int i = 1; i <= 3; i++) {
            a(i,j) = 0.0;
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

      }

    }

    // for (int i = 1; i <= 3; i++) {
    //   for (int j = 1; j <= 3; j++) {
    //     printf(" a %d %d  = %24.14E \n", i,j,a(i,j));
    //   }
    // }

    // TODO: optimize indexing of this loop 
    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
        Ghat(i,j,kxx,kyy,kzz)   = a(i,j);
      }
    }

  }); // end FOR_ALL_CLASS


  for (int ii = 1; ii <= 3; ii++) {
    for (int jj = 1; jj <= 3; jj++) {
      FOR_ALL_CLASS(kzz, 1, npts3_cmplx+1,
                    kyy, 1, npts2_cmplx+1,
                    kxx, 1, npts1_cmplx+1, {

        real_t xknorm;
    
        // thread private arrays
        real_t xk_[3];
        real_t krot_re_[3];
        real_t krot_im_[3];
        real_t D_dft_[3*3];
    
        // create views of thread private arrays
        ViewMatrixTypeReal xk(xk_,3);
        ViewMatrixTypeReal krot_re(krot_re_,3);
        ViewMatrixTypeReal krot_im(krot_im_,3);
        ViewMatrixTypeReal D_dft(D_dft_,3,3);
    
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
  
          if ( (kxx + local_start1_cmplx == npts1_g/2+1 && npts1_g % 2 == 0) || 
               (kyy + local_start2_cmplx == npts2_g/2+1 && npts2_g % 2 == 0) || 
               (npts3_g > 1 && (kzz + local_start3_cmplx == npts3_g/2+1 && npts3_g % 2 == 0)) ) {
             for (int i = 1; i <= 3; i++) {
               for (int j = 1; j <= 3; j++) {
                 work(i,j,kxx,kyy,kzz) = -s0(i,ii,j,jj);
                 workim(i,j,kxx,kyy,kzz) = 0.0;
               }
             }
          } else {
             for (int i = 1; i <= 3; i++) {
               for (int j = 1; j <= 3; j++) {
                 work(i,j,kxx,kyy,kzz) = -Ghat(i,j,kxx,kyy,kzz)*xk(ii)*xk(jj);
                 workim(i,j,kxx,kyy,kzz) = 0.0;
               }
             }
          }
          // printf(" g1 %d %d %d  = %24.14E \n", kxx,kyy,kzz,work(1,1,kxx,kyy,kzz));

        } else if (igamma == 1) { 

          // TODO: this can be precomputed and stored in D_dft array (not sure if memory is issue)
          xk(1) =  xk(1)*twopi;
          xk(2) =  xk(2)*twopi;
          xk(3) =  xk(3)*twopi;
    
          mod_frequency(xk.pointer(), krot_re.pointer(), krot_im.pointer());

          if (abs(POW2(krot_re(1)) + POW2(krot_im(1)) + POW2(krot_re(2)) + POW2(krot_im(2)) + 
              POW2(krot_re(3)) + POW2(krot_im(3))) < 1.0e-15) {
            for (int j = 1; j <= 3; j++) {
              for (int i = 1; i <= 3; i++) {
                work(i,j,kxx,kyy,kzz) = 0.0;
                workim(i,j,kxx,kyy,kzz) = 0.0;
              }
            }
             
          } else {
    
            for (int i = 1; i <= 3; i++) {
              for (int j = 1; j <= 3; j++) {
                D_dft(i,j) = krot_re(i)*krot_re(j) + krot_im(i)*krot_im(j);
              }
            }
    
            for (int i = 1; i <= 3; i++) {
              for (int j = 1; j <= 3; j++) {
                work(i,j,kxx,kyy,kzz) = -Ghat(i,j,kxx,kyy,kzz)*D_dft(ii,jj);
                workim(i,j,kxx,kyy,kzz) = 0.0;
              }
            }
    
          }

        }

      }); // end FOR_ALL_CLASS

      backward_fft();

      Kokkos::parallel_reduce(
        Kokkos::MDRangePolicy<Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER>>({1,1,1}, {npts3+1,npts2+1,npts1+1}),
        KOKKOS_CLASS_LAMBDA(const int k, const int j, const int i, ArrayType <real_t,n> & loc_reduce) {
 
        // printf(" local_start1, local_start2, local_start3 %d %d %d \n", local_start1, local_start2, local_start3);
 
        if ( i + local_start1 == 1 &&
             j + local_start2 == 1 && 
             k + local_start3 == 1 ) {
          // printf(" i + local_start1 %d \n", i + local_start1);
          // printf(" j + local_start2 %d \n", j + local_start2);
          // printf(" k + local_start3 %d \n", k + local_start3);
          int ic;
          ic = -1;
          for (int kk = 1; kk <= 3; kk++) {
            for (int ll = 1; ll <= 3; ll++) {
              ic = ic + 1;
              loc_reduce.array[ic]= work(kk,ll,i,j,k);
            }
          }
        }
         
      }, all_reduce);
      Kokkos::fence(); // needed to prevent race condition

      MPI_Allreduce(MPI_IN_PLACE, all_reduce.array, all_reduce.size, MPI_REAL_T, MPI_SUM, mpi_comm);
      int ic;
      ic = -1;     
      for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
          ic = ic + 1;
          Goperr0(i,ii,j,jj) = all_reduce.array[ic];
          // printf(" g0 %d %d %d %d  = %24.14E \n", i,ii,j,jj,Goperr0(i,ii,j,jj));
        }
      }
    }
  }
  // exit(1);
}
