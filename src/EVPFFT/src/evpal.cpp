#include "evpfft.h"
#include "math_functions.h"
#include "utilities.h"
#include "reduction_data_structures.h"
#include "Profiler.h"

void EVPFFT::evpal(int imicro)
{
  Profiler profiler(__FUNCTION__);

#if BUILD_EVPFFT_FIERRO
  // create space to perform reduction for dsde and cg66.
  const size_t n = 3+36+36+9; // 2 for errs, erre, iter. 36 for dsde_avg. 36 for cg66_avg. 9 edotp_avg
  ArrayType <real_t, n> all_reduce;
#else
  const size_t n = 3; // 2 for errs, erre
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
      itmaxal = 100;
      iter1 = 0;

      while (iter1 < itmaxal && erral > erroral) {
          iter1 += 1;
          for (int iii = 1; iii <= 6; iii++) {
            sg6old(iii) = sg6(iii);
          }

          if (ithermo != 1 || imicro != 1) {
         
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
#endif
              isign = 1;
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

            for (int ii = 1; ii <= 5; ii++) {
              edotp6(ii) = 0.0;
              for (int k1 = 1; k1 <= nsyst(jph); k1++) {
                edotp6(ii) += sch(ii,k1,i,j,k) * rss2(k1,i,j,k);
              }
            }
            edotp6(6) = 0.0;

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

#if 0  
          int error_flag = invert_matrix <6> (xjacobinv.pointer());

          // TODO: optimize indexing of this loop
          for (int ii = 1; ii <= 6; ii++) {
            for (int jj = 1; jj <= 6; jj++) {
              sg6(ii) += -xjacobinv(ii,jj)*res(jj);
            } // end for jj
          } // end for ii
#endif

          // Calculate new stress by solving the system -[J][delt_sg] = [R]
          int error_flag = solve_linear_system(xjacobinv.pointer(), res.pointer(), 6);
          //if (error_flag != 0) {
          //  printf("Failed solve linear system\n");
          //  exit(1);
          //}
        
          //sg6 = sg6 - res
          for (int ii = 1; ii <= 6; ii++) {
            sg6(ii) -= res(ii);
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
      loc_reduce.array[2] += iter1*wgt;

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
    ViewMatrixTypeReal cg66_avg(&loc_reduce.array[3],6,6);
    ViewMatrixTypeReal dedotp66_avg(&loc_reduce.array[39],6,6);
    for (int ii = 1; ii <= 6; ii++) {
      for (int jj = 1; jj <= 6; jj++) {
        cg66_avg(ii,jj) += cg66(ii,jj,i,j,k) * wgt;
        dedotp66_avg(ii,jj) += dedotp66(ii,jj) * wgt;
      }
    }

    ViewMatrixTypeReal edotp_avg_loc(&loc_reduce.array[75],3,3);
    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        edotp_avg_loc(ii,jj) += edotp(ii,jj,i,j,k) * wgt;
      }
    }
#endif

  }, all_reduce);
  Kokkos::fence(); // needed to prevent race condition

  //update host
  gamdot.update_host();

  MPI_Allreduce(MPI_IN_PLACE, all_reduce.array, all_reduce.size, MPI_REAL_T, MPI_SUM, mpi_comm);

  errs = all_reduce.array[0];
  erre = all_reduce.array[1];
  avg_nr_iter = all_reduce.array[2];

#if BUILD_EVPFFT_FIERRO
  ViewMatrixTypeReal cg66_avg_view (&all_reduce.array[3],6,6);
  ViewMatrixTypeReal dedotp66_avg_view (&all_reduce.array[39],6,6);
  ViewMatrixTypeReal edotp_avg_view (&all_reduce.array[75],3,3); 

  for (int ii = 1; ii <= 6; ii++) {
    for (int jj = 1; jj <= 6; jj++) {
      cg66_avg(ii,jj) = cg66_avg_view(ii,jj);
      sg66_avg(ii,jj) = cg66_avg_view(ii,jj);
      dedotp66_avg(ii,jj) = dedotp66_avg_view(ii,jj);
    }
  }

  invert_matrix <6> (sg66_avg.pointer());

  // copy edotp_avg_view into edotp_avg
  for (int ii = 1; ii <= 3; ii++) {
    for (int jj = 1; jj <= 3; jj++) {
      edotp_avg(ii,jj) = edotp_avg_view(ii,jj);
    }
  }

// endif for if BUILD_EVPFFT_FIERRO
#endif

}
