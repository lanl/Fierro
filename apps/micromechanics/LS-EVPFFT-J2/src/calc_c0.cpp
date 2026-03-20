#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "utilities.h"
#include "Profiler.h"
#include "matrix_exp.h"
#include "determinant33.h"
#include "math_functions.h"
#include "reduction_data_structures.h"

void EVPFFT::calc_c0()
{
  Profiler profiler(__FUNCTION__);

  const size_t n = 6*6; // for dP6_dF6(6,6)
  ArrayType <real_t, n> all_reduce;

  real_t dP6_dF6_[6*6];

  // create views of thread private arrays
  ViewMatrixTypeReal dP6_dF6(dP6_dF6_,6,6);

  Kokkos::parallel_reduce(
    Kokkos::MDRangePolicy<Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER>>({1,1,1}, {npts3+1,npts2+1,npts1+1}),
    KOKKOS_CLASS_LAMBDA(const int k, const int j, const int i, ArrayType <real_t,n> & loc_reduce) {

    int jph;
    int isign;

    real_t sgvm;
    real_t mag;

    // thread private arrays
    real_t sg6_[6];
    real_t sgaux_[3*3];
    real_t dedotp66_[6*6];
    real_t sg66_[6*6];
    real_t dsg6_de6_[6*6];
    real_t dPg6_dF6_[6*6];
    real_t dsg_de_[3*3*3*3];
    real_t dPg_dF_[3*3*3*3];

    // create views of thread private arrays
    ViewMatrixTypeReal sg6(sg6_,6);
    ViewMatrixTypeReal sgaux(sgaux_,3,3);
    ViewMatrixTypeReal dedotp66(dedotp66_,6,6);
    ViewMatrixTypeReal sg66(sg66_,6,6);
    ViewMatrixTypeReal dsg6_de6(dsg6_de6_,6,6);
    ViewMatrixTypeReal dPg6_dF6(dPg6_dF6_,6,6);
    ViewMatrixTypeReal dsg_de(dsg_de_,3,3,3,3);
    ViewMatrixTypeReal dPg_dF(dPg_dF_,3,3,3,3);

    jph = jphase(i,j,k);

    if (igas(jph) == 0) {

      // TODO: optimize indexing of this loop
      for (int ii = 1; ii <= 6; ii++) {
        for (int jj = 1; jj <= 6; jj++) {
          sg66(ii,jj) = cg66(ii,jj,i,j,k);
        }
      }

      invert_matrix <6> (sg66.pointer());
      
      // TODO: optimize indexing of this loop
      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          sgaux(ii,jj)     = sg(ii,jj,i,j,k);
        }
      }

      cb.chg_basis_2(sg6.pointer(), sgaux.pointer(), 2, 6, cb.B_basis_device_pointer());

      if (iJ2(jph) == 0) { 

        //     GET RESOLVED SHEAR STRESSES 'rss' AND SHEAR RATES 'gamdot'.
        //     SIGN(GAMDOT)=SIGN(RSS).
        //     NRS CAN BE EVEN OR ODD.
        //     RSS1 IS ALWAYS > 0 AND IS USED TO CALCULATE VISCOUS COMPLIANCE.
    
        for (int is = 1; is <= nsyst(jph); is++) {
#ifdef NON_SCHMID_EFFECTS
          rss(is,i,j,k) = schnon(1,is,,i,j,k)*sg6(1) + 
                    schnon(2,is,i,j,k)*sg6(2) + 
                    schnon(3,is,i,j,k)*sg6(3) + 
                    schnon(4,is,i,j,k)*sg6(4) + 
                    schnon(5,is,i,j,k)*sg6(5);
#else
          rss(is,i,j,k) = sch(1,is,i,j,k)*sg6(1) +
                    sch(2,is,i,j,k)*sg6(2) + 
                    sch(3,is,i,j,k)*sg6(3) + 
                    sch(4,is,i,j,k)*sg6(4) + 
                    sch(5,is,i,j,k)*sg6(5);
#endif      isign = 1;
          if ( (rss(is,i,j,k)-xkin(is,i,j,k)) < 0.0 ) {
            isign = 2;
          }
       
          rss(is,i,j,k) = (rss(is,i,j,k)-xkin(is,i,j,k))/crss(is,isign,i,j,k);
    
#ifdef TWO_SIGN_SLIP_SYSTEMS
          if ( rss(is,i,j,k) < 0.0 ) {
            rss(is,i,j,k) = 0.0;
          }
#endif
  
#ifdef TWO_SIGN_SLIP_SYSTEMS
          rss1(is,i,j,k) = gamd0(is,jph) * nrs(is,jph) * ABS(PowIntExpo(rss(is,i,j,k),(nrs(is,jph)-1))) / crss(is,isign,i,j,k);
          rss2(is,i,j,k) = gamd0(is,jph) * ABS(PowIntExpo(rss(is,i,j,k),nrs(is,jph))) * COPYSIGN(1.0,rss(is,i,j,k));
#else
          rss1(is,i,j,k) = gamd0(is,jph) * nrs(is,jph) * ABS(PowIntExpo(rss(is,i,j,k),(nrs(is,jph)-1))) / crss(is,isign,i,j,k);
          rss2(is,i,j,k) = gamd0(is,jph) * ABS(PowIntExpo(rss(is,i,j,k),nrs(is,jph))) * COPYSIGN(1.0,rss(is,i,j,k));
#endif
  
          gamdot(is,i,j,k) = rss2(is,i,j,k);
        } // end for is
    
        for (int jj = 1; jj <= 5; jj++) {
          for (int ii = 1; ii <= 5; ii++) {
    
            dedotp66(ii,jj) = 0.0;
    
            for (int k1 = 1; k1 <= nsyst(jph); k1++) {         
#ifdef NON_SCHMID_EFFECTS   
              dedotp66(ii,jj) += sch(ii,k1,i,j,k)*schnon(jj,k1,i,j,k)*rss1(k1,i,j,k);
#else
              dedotp66(ii,jj) += sch(ii,k1,i,j,k)*sch(jj,k1,i,j,k)*rss1(k1,i,j,k);
#endif
            } // end for k1
  
          } // end for ii
        } // end for jj

        for (int ii = 1; ii <= 6; ii++) {
          dedotp66(ii,6) = 0.0;
          dedotp66(6,ii) = 0.0;
        }

      } else { // else for if (iJ2(jph) == 0)

        sgvm = 0.0;
        for (int ii = 1; ii <= 5; ii++) {     
          sgvm += POW2(sg6(ii));
        } // end for ii
        sgvm = sqrt(sgvm);

        mag = 1.5*edotp0_j2(jph)/sigma0gr(i,j,k)*PowIntExpo(sgvm/sigma0(i,j,k),(nrs_j2(jph)-1));
        for (int jj = 1; jj <= 5; jj++) {
          for (int ii = 1; ii <= 5; ii++) {
    
            dedotp66(ii,jj) = mag*((nrs_j2(jph)-1)*sg6(ii)/sgvm*sg6(jj)/sgvm + (ii/jj)*(jj/ii));

          } // end for ii
        } // end for jj
        
        for (int ii = 1; ii <= 6; ii++) {
          dedotp66(ii,6) = 0.0;
          dedotp66(6,ii) = 0.0;
        }

      } // end if (iJ2(jph) == 0)

      for (int ii = 1; ii <= 6; ii++) {
        for (int jj = 1; jj <= 6; jj++) {
          dsg6_de6(ii,jj) = dedotp66(ii,jj) + sg66(ii,jj)/tdot;
        }
      }
      invert_matrix <6> (dsg6_de6.pointer());

      cb.chg_basis_3(dsg6_de6.pointer(), dsg_de.pointer(), 3, 6, cb.B_basis_device_pointer());

      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          for (int kk = 1; kk <= 3; kk++) {
            for (int ll = 1; ll <=3; ll++) {
              dPg_dF(ii,jj,kk,ll) = 0.0;
              for (int p = 1; p <= 3; p++) {
                for (int r = 1; r <= 3; r++) {
                  dPg_dF(ii,jj,kk,ll) += dsg_de(ii,p,kk,r)*defgradinv(jj,p,i,j,k)*defgradinv(ll,r,i,j,k)*detF(i,j,k);
                }
              }
            }
          }
        }
      }

      cb.chg_basis_4(dPg6_dF6.pointer(), dPg_dF.pointer(), 4, 6, cb.B_basis_device_pointer());

      // averages packed in array
      int ic;
      ic = -1;
      
      for (int ii = 1; ii <= 6; ii++) {
        for (int jj = 1; jj <= 6; jj++) {
          ic = ic + 1;
          loc_reduce.array[ic] += dPg6_dF6(ii,jj) * wgt;
        }
      }

    }

  }, all_reduce);
  Kokkos::fence(); // needed to prevent race condition

  MPI_Allreduce(MPI_IN_PLACE, all_reduce.array, all_reduce.size, MPI_REAL_T, MPI_SUM, mpi_comm);

  // averages unpacked from array
  int ic;
  ic = -1;

  for (int ii = 1; ii <= 6; ii++) {
    for (int jj = 1; jj <= 6; jj++) {
      ic = ic + 1;
      dP6_dF6(ii,jj) = all_reduce.array[ic];
    }
  }

  real_t s066_[6*6];
  ViewMatrixTypeReal s066(s066_,6,6);
  for (int ii = 1; ii <= 6; ii++) {
    for (int jj = 1; jj <= 6; jj++) {
      c066.host(ii,jj) = dP6_dF6(ii,jj);
      s066(ii,jj) = dP6_dF6(ii,jj);
    }
  }
  invert_matrix <6> (s066.pointer());

  cb.chg_basis_3(c066.host_pointer(), c0.host_pointer(), 3, 6, cb.B_basis_host_pointer());
  cb.chg_basis_3(s066.pointer(), s0.host_pointer(), 3, 6, cb.B_basis_host_pointer());

  // update_device
  c066.update_device();
  c0.update_device();
  s0.update_device();

}
