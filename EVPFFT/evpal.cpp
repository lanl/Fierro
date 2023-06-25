#include "evpfft.h"
#include "inverse.h"
#include "utilities.h"
#include "reduction_data_structures.h"

#if 0
// custom data structure for reduction in function evpal
struct EvpalReduce
{
  real_t errs = 0.0;
  real_t erre = 0.0;
  real_t nr_iter_ave_voxels = 0.0;

  KOKKOS_INLINE_FUNCTION    // Default constructor
  EvpalReduce(){}

  KOKKOS_INLINE_FUNCTION    // Copy constructor
  EvpalReduce(const EvpalReduce & rhs) {
    this->errs = rhs.errs;
    this->erre = rhs.erre;
    this->nr_iter_ave_voxels = rhs.nr_iter_ave_voxels;
  }

  KOKKOS_INLINE_FUNCTION    // add operator
  EvpalReduce& operator+=(const EvpalReduce & rhs) {
    this->errs += rhs.errs;
    this->erre += rhs.erre;
    this->nr_iter_ave_voxels += rhs.nr_iter_ave_voxels;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION    // volatile add operator
  void operator+=(const volatile EvpalReduce & rhs) volatile {
    this->errs += rhs.errs;
    this->erre += rhs.erre;
    this->nr_iter_ave_voxels += rhs.nr_iter_ave_voxels; 
  }

};
#endif

void EVPFFT::evpal(int imicro)
{

  // Reduction is performed on errs, erre
  // Note: errs = all_reduce[0]
  //       erre = all_reduce[1]
  //
#if BUILD_EVPFFT_FIERRO
  // create space to perform reduction for dsde.
  const size_t n = 2+36; // 2 for errs, erre. 36 for M66
  ArrayType <real_t, n> all_reduce;
#else
  const size_t n = 2; // 2 for errs, erre
  ArrayType <real_t, n> all_reduce;
#endif

  Kokkos::parallel_reduce(
    Kokkos::MDRangePolicy<Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER>>({1,1,1}, {npts3+1,npts2+1,npts1+1}),
    KOKKOS_CLASS_LAMBDA(const int k, const int j, const int i, ArrayType <real_t,n> & loc_reduce) {

    int jph;
    int itmaxal;
    int iter1;
    int isign;

    real_t sgnorm;
    real_t dgnorm;
    real_t erroral;
    real_t erral;
    real_t dsgnorm1;
    real_t dsgnorm2;
    real_t ddgnorm;
    real_t errald;

    // thread private arrays
    real_t xlambda_[3*3];
    real_t xlambda6_[6];
    real_t sgaux_[3*3];
    real_t sg6_[6];
    real_t sg6old_[6];
    real_t eptaux_[3*3];
    real_t ept6_[6];
    real_t edotpaux_[3*3];
    real_t edotp6_[6];
    real_t dedotp66_[6*6];
    real_t strainaux_[3*3];
    real_t strain6_[6];
    real_t sg66_[6*6];
    real_t strainceq_[3*3];
    real_t strainceq6_[6];
    real_t res_[6];
    real_t xjacobinv_[6*6];
    real_t rss_[NSYSMX];
    real_t rss1_[NSYSMX];
    real_t rss2_[NSYSMX];
    real_t sc_[5*NSYSMX];
    real_t taux_[NSYSMX*2];
    real_t nrsx_[NSYSMX];
#ifdef NON_SCHMID_EFFECTS
    real_t scnon_[5*NSYSMX];
#endif
    real_t xkinaux_[NSYSMX];

    // create views of thread private arrays
    ViewMatrixTypeReal xlambda(xlambda_,3,3);
    ViewMatrixTypeReal xlambda6(xlambda6_,6);
    ViewMatrixTypeReal sgaux(sgaux_,3,3);
    ViewMatrixTypeReal sg6(sg6_,6);
    ViewMatrixTypeReal sg6old(sg6old_,6);
    ViewMatrixTypeReal eptaux(eptaux_,3,3);
    ViewMatrixTypeReal ept6(ept6_,6);
    ViewMatrixTypeReal edotpaux(edotpaux_,3,3);
    ViewMatrixTypeReal edotp6(edotp6_,6);
    ViewMatrixTypeReal dedotp66(dedotp66_,6,6);
    ViewMatrixTypeReal strainaux(strainaux_,3,3);
    ViewMatrixTypeReal strain6(strain6_,6);
    ViewMatrixTypeReal sg66(sg66_,6,6);
    ViewMatrixTypeReal strainceq(strainceq_,3,3);
    ViewMatrixTypeReal strainceq6(strainceq6_,6);
    ViewMatrixTypeReal res(res_,6);
    ViewMatrixTypeReal xjacobinv(xjacobinv_,6,6);
    ViewMatrixTypeReal rss(rss_,NSYSMX);
    ViewMatrixTypeReal rss1(rss1_,NSYSMX);
    ViewMatrixTypeReal rss2(rss2_,NSYSMX);
    ViewMatrixTypeReal sc(sc_,5,NSYSMX);
    ViewMatrixTypeReal taux(taux_,NSYSMX,2);
    ViewMatrixTypeReal nrsx(nrsx_,NSYSMX);
#ifdef NON_SCHMID_EFFECTS
    ViewMatrixTypeReal scnon(scnon_,5,NSYSMX);
#endif
    ViewMatrixTypeReal xkinaux(xkinaux_,NSYSMX);


    jph = jphase(i,j,k);

    if (igas(jph) == 0) {

      // TODO: optimize indexing of this loop
      for (int ii = 1; ii <= 6; ii++) {
        for (int jj = 1; jj <= 6; jj++) {
          sg66(ii,jj) = cg66(ii,jj,i,j,k);
        }
      }

#ifdef LU_MATRIX_INVERSE
      lu_inverse(sg66.pointer(), 6);
#elif GJE_MATRIX_INVERSE
      inverse_gj(sg66.pointer(), 6);
#endif

      // TODO: optimize indexing of this loop
      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          xlambda(ii,jj)   = sg(ii,jj,i,j,k);
          sgaux(ii,jj)     = sg(ii,jj,i,j,k);
          eptaux(ii,jj)    = ept(ii,jj,i,j,k);
          strainaux(ii,jj) = (disgrad(ii,jj,i,j,k) + disgrad(jj,ii,i,j,k)) / 2.0;
        }
      }

      cb.chg_basis_2(xlambda6.pointer(), xlambda.pointer(), 2, 6, cb.B_basis_device_pointer());
      cb.chg_basis_2(sg6.pointer(), sgaux.pointer(), 2, 6, cb.B_basis_device_pointer());
      cb.chg_basis_2(ept6.pointer(), eptaux.pointer(), 2, 6, cb.B_basis_device_pointer());
      cb.chg_basis_2(strain6.pointer(), strainaux.pointer(), 2, 6, cb.B_basis_device_pointer());

      sgnorm = 0.0;
      dgnorm = 0.0;
      // TODO: optimize indexing of this loop
      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          sgnorm += xlambda(ii,jj)*xlambda(ii,jj);
          dgnorm += strainaux(ii,jj)*strainaux(ii,jj);
        }
      }
      sgnorm = SQRT(sgnorm);
      dgnorm = SQRT(dgnorm);

      erroral = 0.0000001;
      erral   = 2.0*erroral;
      itmaxal = MAX_ITER_NR; //100
      iter1 = 0;

      while (iter1 < itmaxal && erral > erroral) {
      //for (int counter = 1; counter <= MAX_ITER_NR; counter++) {
        //if (iter1 < itmaxal && erral > erroral) {
          iter1 += 1;
          for (int iii = 1; iii <= 6; iii++) {
            sg6old(iii) = sg6(iii);
          }

          if (ithermo != 1 || imicro != 1) {
         
            for (int is = 1; is <= nsyst(jph); is++) {
              nrsx(is) = nrs(is,jph);

              taux(is,1) = crss(is,1,i,j,k);
              taux(is,2) = crss(is,2,i,j,k);

              xkinaux(is) = xkin(is,i,j,k);

              for (int jj = 1; jj <= 5; jj++) {
                sc(jj,is) = sch(jj,is,i,j,k);
              } // end for jj

#ifdef NON_SCHMID_EFFECTS
              for (int jj = 1; jj <= 5; jj++) {
                scnon(jj,is) = schnon(jj,is,i,j,k);
              }
#endif

            } // end for is

            //     GET RESOLVED SHEAR STRESSES 'rss' AND SHEAR RATES 'gamdot'.
            //     SIGN(GAMDOT)=SIGN(RSS).
            //     NRS CAN BE EVEN OR ODD.
            //     RSS1 IS ALWAYS > 0 AND IS USED TO CALCULATE VISCOUS COMPLIANCE.

            for (int is = 1; is <= nsyst(jph); is++) {
#ifdef NON_SCHMID_EFFECTS
              rss(is) = scnon(1,is)*sg6(1) + 
                        scnon(2,is)*sg6(2) + 
                        scnon(3,is)*sg6(3) + 
                        scnon(4,is)*sg6(4) + 
                        scnon(5,is)*sg6(5);
#else
              rss(is) = sc(1,is)*sg6(1) +
                        sc(2,is)*sg6(2) + 
                        sc(3,is)*sg6(3) + 
                        sc(4,is)*sg6(4) + 
                        sc(5,is)*sg6(5);
#endif
              isign = 1;
              if ( (rss(is)-xkinaux(is)) < 0.0 ) {
                isign = 2;
              }
   
              rss(is) = (rss(is)-xkinaux(is))/taux(is,isign);

#ifdef TWO_SIGN_SLIP_SYSTEMS
              if ( rss(is) < 0.0 ) {
                rss(is) = 0.0;
              }
#endif

#ifdef TWO_SIGN_SLIP_SYSTEMS
              rss1(is) = gamd0(is,jph) * nrsx(is) * ABS(POW(rss(is),(nrsx(is)-1))) / taux(is,isign);
              rss2(is) = gamd0(is,jph) * ABS(POW(rss(is),nrsx(is))) * COPYSIGN(1.0,rss(is));
#else
              rss1(is) = gamd0(is,jph) * nrsx(is) * ABS(POW(rss(is),(nrsx(is)-1))) / taux(is,isign);
              rss2(is) = gamd0(is,jph) * ABS(POW(rss(is),nrsx(is))) * COPYSIGN(1.0,rss(is));
#endif

              gamdot(is,i,j,k) = rss2(is);
            } // end for is

            for (int ii = 1; ii <= 5; ii++) {
              edotp6(ii) = 0.0;
              for (int k1 = 1; k1 <= nsyst(jph); k1++) {
                edotp6(ii) += sc(ii,k1) * rss2(k1);
              }
            }
            edotp6(6) = 0.0;

            for (int jj = 1; jj <= 5; jj++) {
              for (int ii = 1; ii <= 5; ii++) {

                dedotp66(ii,jj) = 0.0;

                for (int k1 = 1; k1 <= nsyst(jph); k1++) {         
#ifdef NON_SCHMID_EFFECTS   
                  dedotp66(ii,jj) += sc(ii,k1)*scnon(jj,k1)*rss1(k1);
#else
                  dedotp66(ii,jj) += sc(ii,k1)*sc(jj,k1)*rss1(k1);
#endif
                } // end for k1

              } // end for ii
            } // end for jj

            for (int ii = 1; ii <= 6; ii++) {
              dedotp66(ii,6) = 0.0;
              dedotp66(6,ii) = 0.0;
            }

          } else { // else for if (ithermo != 1 || imicro != 1)
            for (int ii = 1; ii <= 6; ii++) {
              edotp6(ii) = 0.0;
            } // end for ii

          } // end for else (ithermo != 1 || imicro != 1)


          for (int ii = 1; ii <= 6; ii++) {
            if (ithermo == 1 && imicro == 1) {
              strainceq6(ii) = ept6(ii);
            } else {
              strainceq6(ii) = ept6(ii)+edotp6(ii)*tdot;
            }
            for (int jj = 1; jj <= 6; jj++) {
              strainceq6(ii) += sg66(ii,jj)*sg6(jj);
            } // end for jj
          } // end for ii

          cb.chg_basis_1(strainceq6.pointer(), strainceq.pointer(), 1, 6, cb.B_basis_device_pointer());

          for (int ii = 1; ii <= 6; ii++) {
            res(ii) = sg6(ii) - xlambda6(ii);
            for (int jj = 1; jj <= 6; jj++) {
              res(ii) += c066(ii,jj) * (strainceq6(jj) - strain6(jj));
            }
          }

          for (int ii = 1; ii <= 6; ii++) {
            for (int jj = 1; jj <= 6; jj++) {

              xjacobinv(ii,jj) = (ii/jj)*(jj/ii);

              for (int kk = 1; kk <= 6; kk++) {
                if (ithermo == 1 && imicro == 1) {
                  xjacobinv(ii,jj) += c066(ii,kk)*sg66(kk,jj);
                } else {
                  xjacobinv(ii,jj) += c066(ii,kk)*(sg66(kk,jj)+dedotp66(kk,jj)*tdot);
                }
              } // end for kk

            } // end for jj
          } // end for ii

#ifdef LU_MATRIX_INVERSE
          lu_inverse(xjacobinv.pointer(), 6);
#elif GJE_MATRIX_INVERSE
          inverse_gj(xjacobinv.pointer(), 6);
#endif

          // TODO: optimize indexing of this loop
          for (int ii = 1; ii <= 6; ii++) {
            for (int jj = 1; jj <= 6; jj++) {
              sg6(ii) += -xjacobinv(ii,jj)*res(jj);
            } // end for jj
          } // end for ii

          dsgnorm1 = 0.0;
          dsgnorm2 = 0.0;
          ddgnorm  = 0.0;

          for (int ii = 1; ii <= 6; ii++) {     
            dsgnorm1 += POW2(sg6(ii)-sg6old(ii));
            dsgnorm2 += POW2(sg6(ii)-xlambda6(ii));
            ddgnorm  += POW2(strainceq6(ii)-strain6(ii));
          } // end for ii

          erral  = dsgnorm1 / sgnorm;
          errald = ddgnorm / dgnorm;

        //} // end if (iter1 < itmaxal && erral > erroral)
      //} // end for counter
      } // end while (iter1 < itmaxal && erral > erroral)

      cb.chg_basis_1(sg6.pointer(), sgaux.pointer(), 1, 6, cb.B_basis_device_pointer());
      cb.chg_basis_1(edotp6.pointer(), edotpaux.pointer(), 1, 6, cb.B_basis_device_pointer());

      // TODO: optimize indexing of this loop
      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          sg(ii,jj,i,j,k)    = sgaux(ii,jj);
          edotp(ii,jj,i,j,k) = edotpaux(ii,jj);
        } // end for jj
      } // end for ii

      loc_reduce.array[0] += dsgnorm2*wgt;
      loc_reduce.array[1] += ddgnorm*wgt;

    } else { // else for if (igas(jph) == 0)

      // TODO: optimize indexing of this loop
      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          sg(ii,jj,i,j,k) = 0.0;  
          edotp(ii,jj,i,j,k) = 0.0;    //  ????
        } // end for jj
      } // end for ii

    } // end if (igas(jph) == 0)

#if BUILD_EVPFFT_FIERRO
    real_t dsde_[36];
    ViewMatrixTypeReal dsde(dsde_,6,6);
    if (igas(jph) == 0) {
      for (int ii = 1; ii <= 6; ii++) {
        for (int jj = 1; jj <= 6; jj++) {
          dsde(ii,jj) = sg66(ii,jj) + dedotp66(ii,jj)*tdot;
        }
      }
#ifdef LU_MATRIX_INVERSE
      lu_inverse(dsde.pointer(), 6);
#elif GJE_MATRIX_INVERSE
      inverse_gj(dsde.pointer(), 6);
#endif
      ViewMatrixTypeReal dsde_avg(&loc_reduce.array[2],6,6);
      for (int ii = 1; ii <= 6; ii++) {
        for (int jj = 1; jj <= 6; jj++) {
          dsde_avg(ii,jj) += dsde(ii,jj) * wgt;
        }
      }
    //} else { // else for if (igas(jph) == 0)
    } // end if (igas(jph) == 0)
#endif


  }, all_reduce);
  Kokkos::fence(); // needed to prevent race condition

  //update host
  gamdot.update_host();

  errs = all_reduce.array[0];
  erre = all_reduce.array[1];

#if BUILD_EVPFFT_FIERRO
  // copy dsde_avg to M66
  for (int i = 0; i < 36; i++) {
    M66.pointer()[i] = (&all_reduce.array[2])[i];
  }
#endif

}
