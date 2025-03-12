#include "evpfft.h"
#include "orient.h"
#include "reduction_data_structures.h"

#include "utilities.h"



void EVPFFT::update_orient()
{

  // reduction is performed on rslbar and rlcbar
  // Note: rslbar = all_reduce[0]
  //       rlcbar = all_reduce[1]
  ArrayType <real_t, 2> all_reduce;


  Kokkos::parallel_reduce(
    Kokkos::MDRangePolicy<Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER>>({1,1,1}, {npts3+1,npts2+1,npts1+1}),
    KOKKOS_CLASS_LAMBDA(const int k, const int j, const int i, ArrayType <real_t, 2> & loc_reduce) {

    int iph;

    // thread private arrays
    real_t aa_[3*3];
    real_t distor_[3*3];
    real_t dnsa_[3];
    real_t dbsa_[3];
    real_t rotslip_[3*3];
    real_t rotloc_[3*3];
    real_t rot_[3*3];

    // create views of thread private arrays
    ViewMatrixTypeReal aa(aa_,3,3);
    ViewMatrixTypeReal distor(distor_,3,3);
    ViewMatrixTypeReal dnsa(dnsa_,3);
    ViewMatrixTypeReal dbsa(dbsa_,3);
    ViewMatrixTypeReal rotslip(rotslip_,3,3);
    ViewMatrixTypeReal rotloc(rotloc_,3,3);
    ViewMatrixTypeReal rot(rot_,3,3);


    iph = jphase(i,j,k);

    if (igas(iph) == 0) {
//
//      LOCAL ROTATION RATE: ANTISYM(VELGRAD)
//

      for (int jj = 1; jj <= 3; jj++) {
        for (int ii = 1; ii <= 3; ii++) {
          rotloc(ii,jj) = (velgrad(ii,jj,i,j,k) - velgrad(jj,ii,i,j,k)) / 2.0;
        }
      }
//
//      SLIP ROTATION RATE
//

      for (int jj = 1; jj <= 3; jj++) {
        for (int ii = 1; ii <= 3; ii++) {
          aa(ii,jj) = ag(ii,jj,i,j,k);
          distor(ii,jj) = 0.0;
        }
      }

      for (int is = 1; is <= nsyst(iph); is++) {
        for (int ii = 1; ii <= 3; ii++) {
            dnsa(ii) = 0.0;
            dbsa(ii) = 0.0;
          for (int jj = 1; jj <= 3; jj++) {
            dnsa(ii) += aa(ii,jj) * dnca(jj,is,iph);
            dbsa(ii) += aa(ii,jj) * dbca(jj,is,iph);
          }
        }

        for (int ii = 1; ii <= 3; ii++) {
          for (int jj = 1; jj <= 3; jj++) {
            distor(ii,jj) += dbsa(ii) * dnsa(jj) * gamdot(is,i,j,k);
          }
        }
      } // end for is

      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          rotslip(ii,jj) = (distor(ii,jj)-distor(jj,ii)) / 2.0;
        }
      }

//
//      AVERAGE ROTATION RATE
//
      loc_reduce.array[0] += SQRT(POW2(rotslip(3,2)) + POW2(rotslip(1,3)) + POW2(rotslip(2,1))) * wgtc(i,j,k);
      loc_reduce.array[1] += SQRT(POW2(rotloc(3,2)) + POW2(rotloc(1,3)) + POW2(rotloc(2,1))) * wgtc(i,j,k);

//
//      TOTAL ROTATION
//
      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          rot(ii,jj) = (tomtot(ii,jj)+rotloc(ii,jj)-rotslip(ii,jj)) * tdot;
        }
      }

//
//      REORIENTATION
//
      orient(aa.pointer(), rot.pointer());

//
//      UPDATE ORIENTATION MATRIX
//
      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          ag(ii,jj,i,j,k) = aa(ii,jj);
        }
      }

    } // end if (igas(iph) == 0)

  }, all_reduce);
  Kokkos::fence(); // needed to prevent race condition

  MPI_Allreduce(MPI_IN_PLACE, all_reduce.array, all_reduce.size, MPI_REAL_T, MPI_SUM, mpi_comm);

  real_t rslbar = all_reduce.array[0];
  real_t rlcbar = all_reduce.array[1];

#ifndef ABSOLUTE_NO_OUTPUT
  if (0 == my_rank) {
    printf("\n");
    printf(" AVERAGE PLASTIC ROTATION = %24.14E\n", rslbar);
    printf(" AVERAGE LOCAL ROTATION = %24.14E\n", rlcbar);
  }
#endif

}
