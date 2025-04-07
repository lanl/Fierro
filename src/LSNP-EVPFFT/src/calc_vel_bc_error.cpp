#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "utilities.h"
#include "math_functions.h"
#include "reduction_data_structures.h"

void EVPFFT::calc_vel_bc_error()
{

  const size_t n = 7; // for velocities (7)
  ArrayType <real_t, n> all_reduce;

  real_t count;
  real_t velbound_avg_[3];
  real_t velboundapp_avg_[3];

  // create views of thread private arrays
  ViewMatrixTypeReal velbound_avg(velbound_avg_,3);
  ViewMatrixTypeReal velboundapp_avg(velboundapp_avg_,3);

  Kokkos::parallel_reduce(
    Kokkos::MDRangePolicy<Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER>>({1,1,1}, {npts3+1,npts2+1,npts1+1}),
    KOKKOS_CLASS_LAMBDA(const int k, const int j, const int i, ArrayType <real_t,n> & loc_reduce) {

    int i_g;
    int j_g;
    int k_g;
    int i_gs;
    int j_gs;
    int k_gs;

    i_g = i + local_start1;
    j_g = j + local_start2;
    k_g = k + local_start3;
    i_gs = i_g - dnpts1_g;
    j_gs = j_g - dnpts2_g;
    k_gs = k_g - dnpts3_g;

    if ( (i_gs == 1 || i_gs == npts1_g - 2*dnpts1_g + 1) && (j_gs >= 1 && j_gs <= npts1_g-2*dnpts1_g+1 && k_gs >= 1 && k_gs <= npts3_g-2*dnpts3_g+1) ||
         (j_gs == 1 || j_gs == npts2_g - 2*dnpts2_g + 1) && (i_gs >= 1 && i_gs <= npts2_g-2*dnpts2_g+1 && k_gs >= 1 && k_gs <= npts3_g-2*dnpts3_g+1) ||
         (k_gs == 1 || k_gs == npts2_g - 2*dnpts2_g + 1) && (i_gs >= 1 && i_gs <= npts2_g-2*dnpts2_g+1 && j_gs >= 1 && j_gs <= npts2_g-2*dnpts2_g+1)){

      // averages packed in array
      int ic;
      ic = -1;
  
      ic = ic + 1;
      loc_reduce.array[ic] += 1.0;
      for (int ii = 1; ii <= 3; ii++) {
        ic = ic + 1;
        loc_reduce.array[ic] += velocity(ii,i,j,k);
        ic = ic + 1;
        loc_reduce.array[ic] += velapp_node(ii,i,j,k);
      }

    }

  }, all_reduce);
  Kokkos::fence(); // needed to prevent race condition

  MPI_Allreduce(MPI_IN_PLACE, all_reduce.array, all_reduce.size, MPI_REAL_T, MPI_SUM, mpi_comm);

  // averages unpacked from array
  int ic;
  ic = -1;

  ic = ic + 1;
  count = all_reduce.array[ic];
  for (int ii = 1; ii <= 3; ii++) {
    ic = ic + 1;
    velbound_avg(ii) = all_reduce.array[ic]/count;
    ic = ic + 1;
    velboundapp_avg(ii) = all_reduce.array[ic]/count;
  }

  // printf("velbound_avg %16.8E %16.8E %16.8E\n", velbound_avg(1), velbound_avg(2), velbound_avg(3));
  // printf("velboundapp_avg %16.8E %16.8E %16.8E\n", velboundapp_avg(1), velboundapp_avg(2), velboundapp_avg(3));

  real_t veln_avg;
  real_t dveln_avg;

  Kokkos::parallel_reduce(
    Kokkos::MDRangePolicy<Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER>>({1,1,1}, {npts3+1,npts2+1,npts1+1}),
    KOKKOS_CLASS_LAMBDA(const int k, const int j, const int i, ArrayType <real_t,n> & loc_reduce) {

    int i_g;
    int j_g;
    int k_g;
    int i_gs;
    int j_gs;
    int k_gs;

    i_g = i + local_start1;
    j_g = j + local_start2;
    k_g = k + local_start3;
    i_gs = i_g - dnpts1_g;
    j_gs = j_g - dnpts2_g;
    k_gs = k_g - dnpts3_g;

    if ( (i_gs == 1 || i_gs == npts1_g - 2*dnpts1_g + 1) && (j_gs >= 1 && j_gs <= npts1_g-2*dnpts1_g+1 && k_gs >= 1 && k_gs <= npts3_g-2*dnpts3_g+1) ||
         (j_gs == 1 || j_gs == npts2_g - 2*dnpts2_g + 1) && (i_gs >= 1 && i_gs <= npts2_g-2*dnpts2_g+1 && k_gs >= 1 && k_gs <= npts3_g-2*dnpts3_g+1) ||
         (k_gs == 1 || k_gs == npts2_g - 2*dnpts2_g + 1) && (i_gs >= 1 && i_gs <= npts2_g-2*dnpts2_g+1 && j_gs >= 1 && j_gs <= npts2_g-2*dnpts2_g+1)){

      real_t veln;
      real_t dveln;

      real_t dvel_[3];
      ViewMatrixTypeReal dvel(dvel_,3);

      veln = SQRT(POW2(velapp_node(1,i,j,k) - velboundapp_avg(1)) + POW2(velapp_node(2,i,j,k) - velboundapp_avg(2)) + 
               POW2(velapp_node(3,i,j,k) - velboundapp_avg(3)));

      dvel(1) = (velapp_node(1,i,j,k) - velboundapp_avg(1)) - (velocity(1,i,j,k) - velbound_avg(1));
      dvel(2) = (velapp_node(2,i,j,k) - velboundapp_avg(2)) - (velocity(2,i,j,k) - velbound_avg(2));
      dvel(3) = (velapp_node(3,i,j,k) - velboundapp_avg(3)) - (velocity(3,i,j,k) - velbound_avg(3));

      dveln = SQRT(POW2(dvel(1)) + POW2(dvel(2)) + POW2(dvel(3)));

      // averages packed in array
      loc_reduce.array[0] += veln;
      loc_reduce.array[1] += dveln;

    }

  }, all_reduce);
  Kokkos::fence(); // needed to prevent race condition

  MPI_Allreduce(MPI_IN_PLACE, all_reduce.array, all_reduce.size, MPI_REAL_T, MPI_SUM, mpi_comm);

  veln_avg = all_reduce.array[0]/count;
  dveln_avg = all_reduce.array[1]/count;
  err_disp_bc_vel = dveln_avg/veln_avg;
#ifndef ABSOLUTE_NO_OUTPUT
  if (0 == my_rank) {
    printf(" err_disp_bc_vel %16.8E\n", err_disp_bc_vel);
  }
#endif

}
