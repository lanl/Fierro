#include "evpfft.h"
#include "utilities.h"
#include "reduction_data_structures.h"
#include "Profiler.h"

void EVPFFT::get_smacro()
{
  Profiler profiler(__FUNCTION__);
  
  // Note: reduction is performed on scauav(3,3) and scauav1(3,3)
  //       all_reduce[0] = scauav
  //       all_reduce[1] = scauav1
  ArrayOfArrayType <2, real_t, 3*3> all_reduce;

  Kokkos::parallel_reduce(
    Kokkos::MDRangePolicy<Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER>>({1,1,1}, {npts3+1,npts2+1,npts1+1}),
    KOKKOS_CLASS_LAMBDA(const int k, const int j, const int i,
                  ArrayOfArrayType <2, real_t, 3*3> & loc_reduce) {

      ViewMatrixTypeReal scauav_loc  (&loc_reduce.array[0].array[0], 3, 3);
      ViewMatrixTypeReal scauav1_loc (&loc_reduce.array[1].array[0], 3, 3);

      for (int jj = 1; jj <= 3; jj++) {
        for (int ii = 1; ii <= 3; ii++) {

          scauav_loc(ii,jj) += sg(ii,jj,i,j,k) * wgtc(i,j,k);

          if (jphase(i,j,k) == 1) {
            scauav1_loc(ii,jj) += sg(ii,jj,i,j,k) * wgtc(i,j,k)/wph1;
          }

        } // end for ii
      } // end for jj

  }, all_reduce);
  Kokkos::fence(); // needed to prevent race condition on scauav and scauav1

  ViewMatrixTypeReal scauav_loc  (&all_reduce.array[0].array[0], 3, 3);
  ViewMatrixTypeReal scauav1_loc (&all_reduce.array[1].array[0], 3, 3);
  for (int jj = 1; jj <= 3; jj++) {
    for (int ii = 1; ii <= 3; ii++) {
      scauav(ii,jj) = scauav_loc(ii,jj);
      scauav1(ii,jj) = scauav1_loc(ii,jj);
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, scauav.pointer(), scauav.size(), MPI_REAL_T, MPI_SUM, mpi_comm);
  MPI_Allreduce(MPI_IN_PLACE, scauav1.pointer(), scauav1.size(), MPI_REAL_T, MPI_SUM, mpi_comm);

//
//    MIXED BC
//

  int ii;
  int jj;
  int kk;
  int ll;
  int ijv_[6*2] = {1,2,3,2,1,1,1,2,3,3,3,2};
  ViewMatrixTypeInt  ijv(ijv_,6,2);
  MatrixTypeRealHost sav6(6);
  MatrixTypeRealHost sav5(5);

  for (int i = 1; i <= 6; i++) {

    ii = ijv(i,1);
    jj = ijv(i,2);

    dvelgradmacro(ii,jj) = 0.0;

    if (idsim(i) == 0) {
      for (int k = 1; k <= 6; k++) {
        kk = ijv(k,1);
        ll = ijv(k,2);

        dvelgradmacro(ii,jj) += s0.host(ii,jj,kk,ll) * iscau(k) *
                                     (scauchy(kk,ll) - scauav(kk,ll));
      } // end for k
    } // end if (idsim(i) == 0)

  } // end for i

  for (int jj = 1; jj <= 3; jj++) {
    for (int ii = 1; ii <= 3; ii++) {
      dvelgradmacroacum(ii,jj) += dvelgradmacro(ii,jj);
    }
  }

  cb.chg_basis_2(sav6.pointer(), scauav.pointer(), 2, 6, cb.B_basis_host_pointer());

  for (int ii = 1; ii <= 5; ii++) {
    sav5(ii) = sav6(ii);
  }

  cb.chg_basis_1(sav5.pointer(), sdeviat.pointer(), 1, 5, cb.B_basis_host_pointer());

  svm = 0.0;
  for (int j = 1; j <= 3; j++) {
    for (int i = 1; i <= 3; i++) {
      svm += sdeviat(i,j) * sdeviat(i,j);
    }
  }
  svm = sqrt(3.0/2.0*svm);

  cb.chg_basis_2(sav6.pointer(), scauav1.pointer(), 2, 6, cb.B_basis_host_pointer());

  for (int ii = 1; ii <= 5; ii++) {
    sav5(ii) = sav6(ii);
  }

  cb.chg_basis_1(sav5.pointer(), sdeviat.pointer(), 1, 5, cb.B_basis_host_pointer());

  svm1 = 0.0;
  for (int j = 1; j <= 3; j++) {
    for (int i = 1; i <= 3; i++) {
      svm1 += sdeviat(i,j) * sdeviat(i,j);
    }
  }
  svm1 = sqrt(3.0/2.0*svm1);

}
