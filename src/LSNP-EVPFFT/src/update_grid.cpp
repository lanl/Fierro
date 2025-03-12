#include <iostream>
#include <stdio.h>
#include "evpfft.h"
#include "utilities.h"

void EVPFFT::update_grid()
{

  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    x_grid(1,i,j,k) += velocity(1,i,j,k)*tdot;
    x_grid(2,i,j,k) += velocity(2,i,j,k)*tdot;
    x_grid(3,i,j,k) += velocity(3,i,j,k)*tdot;

  }); // end FOR_ALL_CLASS

  // TODO: pass to xnode array for igamma = 0 and 1

  MatrixTypeRealDual workG(3,npts1_g,npts2_g,npts3_g);

  // write result to ouput
  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {
    for (int ii = 1; ii <= 3; ii++) {
      // workG(ii,i+local_start1,j+local_start2,k+local_start3) = x_grid(ii,i,j,k);
      workG(ii,i+local_start1,j+local_start2,k+local_start3) = velocity(ii,i,j,k);
    }
  }); // end FOR_ALL_CLASS
  Kokkos::fence();

  workG.update_host();
  MPI_Allreduce(MPI_IN_PLACE, workG.host_pointer(), 3*npts1_g*npts2_g*npts3_g, MPI_REAL_T, MPI_SUM, mpi_comm);

  FOR_ALL_CLASS(k, 1, npts3+2,
                j, 1, npts2+2,
                i, 1, npts1+2, {

    if (igamma == 0) {

      int iv1;
      int iv2;
      int jv1;
      int jv2;
      int kv1;
      int kv2;

      if ((i+local_start1) == 1) {
        iv1 = npts1_g;
        iv2 = 1;
      } else if ((i+local_start1) == npts1_g + 1) {
        iv1 = npts1_g;
        iv2 = 1;
      } else {
        iv1 = i+local_start1 - 1;
        iv2 = i+local_start1;
      }
  
      if ((j+local_start2) == 1) {
        jv1 = npts2_g;
        jv2 = 1;
      } else if ((j+local_start2) == npts2_g + 1) {
        jv1 = npts2_g;
        jv2 = 1;
      } else {
        jv1 = j+local_start2 - 1;
        jv2 = j+local_start2;
      }
  
      if ((k+local_start3) == 1) {
        kv1 = npts3_g;
        kv2 = 1;
      } else if ((k+local_start3) == npts3_g + 1) {
        kv1 = npts3_g;
        kv2 = 1;
      } else {
        kv1 = k+local_start3 - 1;
        kv2 = k+local_start3;
      }
  
      for (int ii = 1; ii <= 3; ii++) {
        xnode(ii,i,j,k) += 0.125*(workG(ii,iv1,jv1,kv1) + workG(ii,iv2,jv1,kv1) + 
         workG(ii,iv1,jv2,kv1) + workG(ii,iv2,jv2,kv1) +
         workG(ii,iv1,jv1,kv2) + workG(ii,iv2,jv1,kv2) + 
         workG(ii,iv1,jv2,kv2) + workG(ii,iv2,jv2,kv2))*tdot;
      }
  
    } else if (igamma == 1) {

      if (ibc == 0) {
        printf("needs to be fixed since velgradrefavg is no longer zero so we cannot jsut copy velocity.");
      }

      int ind;
      int jnd;
      int knd;

      if ((i+local_start1) == npts1_g + 1) {
        ind = 1;
      } else {
        ind = i;
      }
      if ((j+local_start2) == npts2_g + 1) {
        jnd = 1;
      } else {
        jnd = j;
      }
      if ((k+local_start3) == npts3_g + 1) {
        knd = 1;
      } else {
        knd = k;
      }

      for (int ii = 1; ii <= 3; ii++) {
        xnode(ii,i,j,k) += workG(ii,ind,jnd,knd)*tdot;
      }

    }

  }); // end FOR_ALL_CLASS
  Kokkos::fence();

}
