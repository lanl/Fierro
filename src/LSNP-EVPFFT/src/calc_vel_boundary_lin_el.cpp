#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "utilities.h"
#include "math_functions.h"
#include "reduction_data_structures.h"

void EVPFFT::calc_vel_boundary_lin_el()
{

  const size_t n = 4; // for velappavg (9)
  ArrayType <real_t, 4> all_reduce;

  // FOR_ALL_CLASS(i, 1, npts1+1,
  //               j, 1, npts2+1,
  //               k, 1, npts3+1, {

  Kokkos::parallel_reduce(
    Kokkos::MDRangePolicy<Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER>>({1,1,1}, {npts3+1,npts2+1,npts1+1}),
    KOKKOS_CLASS_LAMBDA(const int k, const int j, const int i, ArrayType <real_t,n> & loc_reduce) {

    int i_g;
    int j_g;
    int k_g;

    i_g = i + local_start1 - dnpts1_g;
    j_g = j + local_start2 - dnpts2_g;
    k_g = k + local_start3 - dnpts3_g;

    if ( (i_g == 1 || i_g == npts1_g - 2*dnpts1_g + 1) && (j_g >= 1 && j_g <= npts1_g-2*dnpts1_g+1 && k_g >= 1 && k_g <= npts3_g-2*dnpts3_g+1) ||
         (j_g == 1 || j_g == npts2_g - 2*dnpts2_g + 1) && (i_g >= 1 && i_g <= npts2_g-2*dnpts2_g+1 && k_g >= 1 && k_g <= npts3_g-2*dnpts3_g+1) ||
         (k_g == 1 || k_g == npts2_g - 2*dnpts2_g + 1) && (i_g >= 1 && i_g <= npts2_g-2*dnpts2_g+1 && j_g >= 1 && j_g <= npts2_g-2*dnpts2_g+1)){

    // if ((i_g == 1 || j_g == 1 || k_g == 1 || 
    //     i_g == npts1_g - 2*dnpts1_g + 1 || j_g == npts2_g - 2*dnpts2_g + 1 || 
    //     k_g == npts3_g - 2*dnpts3_g + 1) && iframe(i,j,k) == 0){

      real_t g;
      real_t h;
      real_t r;

      real_t phi_[8];

      ViewMatrixTypeReal phi(phi_,8);

      g = 2.0*(i_g - 1)/(npts1_g - 2*dnpts1_g) - 1.0;
      h = 2.0*(j_g - 1)/(npts2_g - 2*dnpts2_g) - 1.0;
      r = 2.0*(k_g - 1)/(npts3_g - 2*dnpts3_g) - 1.0;

      phi(1) = 0.125*(1.0 - g)*(1.0 - h)*(1.0 - r);
      phi(2) = 0.125*(1.0 + g)*(1.0 - h)*(1.0 - r);
      phi(3) = 0.125*(1.0 + g)*(1.0 + h)*(1.0 - r);
      phi(4) = 0.125*(1.0 - g)*(1.0 + h)*(1.0 - r);
      phi(5) = 0.125*(1.0 - g)*(1.0 - h)*(1.0 + r);
      phi(6) = 0.125*(1.0 + g)*(1.0 - h)*(1.0 + r);
      phi(7) = 0.125*(1.0 + g)*(1.0 + h)*(1.0 + r);
      phi(8) = 0.125*(1.0 - g)*(1.0 + h)*(1.0 + r);

      for (int ii = 1; ii <= 3; ii++) {
        velapp_node(ii,i,j,k) = 0.0;
        for (int jj = 1; jj <= 8; jj++) {
          velapp_node(ii,i,j,k) += velvert(ii,jj)*phi(jj);
        }
      }

      // averages packed in array
      int ic;
      ic = -1;
      for (int ii = 1; ii <= 3; ii++) {
        ic = ic + 1;
        loc_reduce.array[ic] += velapp_node(ii,i,j,k);
      }
      ic = ic + 1;
      loc_reduce.array[ic] += 1.0;

    }

//  });
  }, all_reduce);
  Kokkos::fence(); // needed to prevent race condition

  MPI_Allreduce(MPI_IN_PLACE, all_reduce.array, all_reduce.size, MPI_REAL_T, MPI_SUM, mpi_comm);

  real_t velappavg_[3];
  ViewMatrixTypeReal velappavg(velappavg_,3);
  int ic;
  ic = -1;
  for (int ii = 1; ii <=3; ii++) {
    ic = ic + 1;
    velappavg(ii) = all_reduce.array[ic];
  }
  ic = ic + 1;
  for (int ii = 1; ii <= 3; ii++) {
    velappavg(ii) = velappavg(ii)/all_reduce.array[ic];
  }

  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    int i_g;
    int j_g;
    int k_g;

    i_g = i + local_start1 - dnpts1_g;
    j_g = j + local_start2 - dnpts2_g;
    k_g = k + local_start3 - dnpts3_g;

    if ( (i_g == 1 || i_g == npts1_g - 2*dnpts1_g + 1) && (j_g >= 1 && j_g <= npts1_g-2*dnpts1_g+1 && k_g >= 1 && k_g <= npts3_g-2*dnpts3_g+1) ||
         (j_g == 1 || j_g == npts2_g - 2*dnpts2_g + 1) && (i_g >= 1 && i_g <= npts2_g-2*dnpts2_g+1 && k_g >= 1 && k_g <= npts3_g-2*dnpts3_g+1) ||
         (k_g == 1 || k_g == npts2_g - 2*dnpts2_g + 1) && (i_g >= 1 && i_g <= npts2_g-2*dnpts2_g+1 && j_g >= 1 && j_g <= npts2_g-2*dnpts2_g+1)){

    // if ((i_g == 1 || j_g == 1 || k_g == 1 || 
    //     i_g == npts1_g - 2*dnpts1_g + 1 || j_g == npts2_g - 2*dnpts2_g + 1 || 
    //     k_g == npts3_g - 2*dnpts3_g + 1) && iframe(i,j,k) == 0){

      for (int ii = 1; ii <= 3; ii++) {
        velapp_node(ii,i,j,k) -= velappavg(ii);
      }
    }
  });

}
