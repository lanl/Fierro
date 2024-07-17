#include "evpfft.h"
#include "vm.h"
#include "reduction_data_structures.h"
#include "utilities.h"
#include "Profiler.h"

void EVPFFT::step_update_velgrad()
{
  Profiler profiler(__FUNCTION__);

  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {
    
    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        velgrad(ii,jj,i,j,k) = velgradmacro(ii,jj);
      } // end for ii
    } // end for jj
 
  }); // end FOR_ALL_CLASS

}

void EVPFFT::step_set_dvelgradmacro_and_dvelgradmacroacum_to_zero()
{
  for (int j = 1; j <= 3; j++) {
    for (int i = 1; i <= 3; i++) {
      dvelgradmacro(i,j) = 0.0;
      dvelgradmacroacum(i,j)  = 0.0;
    }
  }
}


void EVPFFT::step_update()
{
  Profiler profiler(__FUNCTION__);

  //  velgrad, which contained disgrad at t, is updated
  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        sgt(ii,jj,i,j,k) = sg(ii,jj,i,j,k);
      } // end for ii
    } // end for jj

  }); // end FOR_ALL_CLASS

  for (int jj = 1; jj <= 3; jj++) {
    for (int ii = 1; ii <= 3; ii++) {
      velgradmacroactual(ii,jj) = velgradmacro.host(ii,jj) + dvelgradmacroacum(ii,jj);
      disgradmacro.host(ii,jj) = disgradmacrot(ii,jj) + velgradmacroactual(ii,jj) * tdot;
    } // end for ii
  } // end for jj


  // TOTAL (EL+PL) VM
  evm = vm(disgradmacro.host.pointer());
  dvm = vm(velgradmacroactual.pointer());
  for (int jj = 1; jj <= 3; jj++) {
    for (int ii = 1; ii <= 3; ii++) {
      disgradmacrot(ii,jj) = disgradmacro.host(ii,jj);
    } // end for ii
  } // end for jj

}


void EVPFFT::step_vm_calc()
{
  Profiler profiler(__FUNCTION__);

#if 0
  //     INITIAL GUESS OF DISGRADMACRO AT t+dt ALWAYS ELASTIC   
  for (int jj = 1; jj <= 3; jj++) {
    for (int ii = 1; ii <= 3; ii++) {
      disgradmacro.host(ii,jj) = disgradmacrot(ii,jj) + udot.host(ii,jj) * tdot;
    } // end for ii
  } // end for jj
  // update device
  disgradmacro.update_device();
#endif

  // Note: reduction is performed on epav(3,3) and edotpav(3,3) 
  //       epav = all_reduce[0]
  //       edotpav = all_reduce[1]
  ArrayOfArrayType <2, real_t, 3*3> all_reduce;

  Kokkos::parallel_reduce(
  Kokkos::MDRangePolicy<Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER>>({1,1,1}, {npts3+1,npts2+1,npts1+1}),
  KOKKOS_CLASS_LAMBDA(const int k, const int j, const int i, 
                ArrayOfArrayType <2, real_t, 3*3> & loc_reduce) {

    ViewMatrixTypeReal epav_loc (&loc_reduce.array[0].array[0], 3, 3);
    ViewMatrixTypeReal edotpav_loc (&loc_reduce.array[1].array[0], 3, 3);

    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        ept(ii,jj,i,j,k) += edotp(ii,jj,i,j,k) * tdot;
      } // end for ii
    } // end for jj

//
//     PLASTIC VM
//
    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        epav_loc(ii,jj) += ept(ii,jj,i,j,k) * wgtc(i,j,k);
        edotpav_loc(ii,jj) += edotp(ii,jj,i,j,k) * wgtc(i,j,k);
      }
    }

  }, all_reduce);
  Kokkos::fence(); // needed to prevent race condition

  ViewMatrixTypeReal epav_loc (&all_reduce.array[0].array[0], 3, 3);
  ViewMatrixTypeReal edotpav_loc (&all_reduce.array[1].array[0], 3, 3);
  for (int jj = 1; jj <= 3; jj++) {
    for (int ii = 1; ii <= 3; ii++) {
      epav(ii,jj) = epav_loc(ii,jj);
      edotpav(ii,jj) = edotpav_loc(ii,jj);  
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, epav.pointer(), epav.size(), MPI_REAL_T, MPI_SUM, mpi_comm);
  MPI_Allreduce(MPI_IN_PLACE, edotpav.pointer(), edotpav.size(), MPI_REAL_T, MPI_SUM, mpi_comm);

  evmp = 0.0;
  dvmp = 0.0;
  for (int jj = 1; jj <= 3; jj++) {
    for (int ii = 1; ii <= 3; ii++) {
      evmp += POW2(epav(ii,jj));
      dvmp += POW2(edotpav(ii,jj));
    }
  }
  evmp = sqrt(2.0/3.0*evmp);
  dvmp = sqrt(2.0/3.0*dvmp);
//printf("\t%26.16E, \t%26.16E", evmp, dvmp);
//PAUSE;

}


void EVPFFT::step_texture_rve_update()
{
  Profiler profiler(__FUNCTION__);

  //   VELMAX
  velmax(1) = dsim(1,1)*delt(1)*(npts1-1);
  velmax(2) = dsim(2,2)*delt(2)*(npts2-1);
  velmax(3) = dsim(3,3)*delt(3)*(npts3-1);

  //   UPDATE ORIENTATIONS
  update_orient();

  //   UPDATE DELT
  delt(1) = (delt(1)*(npts1-1)+velmax(1)*tdot)/(npts1-1);
  delt(2) = (delt(2)*(npts2-1)+velmax(2)*tdot)/(npts2-1);
  if (npts3 > 1) {
     delt(3)=(delt(3)*(npts3-1)+velmax(3)*tdot)/(npts3-1);
  }

}


