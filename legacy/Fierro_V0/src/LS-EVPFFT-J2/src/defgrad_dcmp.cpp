#include "defgrad_dcmp.h"
#include "math_functions.h"
#include "utilities.h"

KOKKOS_FUNCTION
real_t frobenius_norm(const real_t* A, int m, int n);

KOKKOS_FUNCTION
void defgrad_dcmp(real_t *F_, real_t *V_, real_t *R_)
{

  ViewMatrixTypeReal F (F_,3,3);
  ViewMatrixTypeReal V (V_,3,3);
  ViewMatrixTypeReal R (R_,3,3);

  real_t A_curr_  [3*3];
  real_t A_next_  [3*3];
  real_t A_inv_   [3*3];
  real_t A_trans_ [3*3];
  real_t R_inv_   [3*3];
  ViewMatrixTypeReal A_curr (A_curr_,3,3);
  ViewMatrixTypeReal A_next (A_next_,3,3);
  ViewMatrixTypeReal A_inv  (A_inv_,3,3);
  ViewMatrixTypeReal A_trans  (A_trans_,3,3);
  ViewMatrixTypeReal R_inv  (R_inv_,3,3);

  const int iter_max = 50;
  const real_t tol = 1.0e-12;
  real_t conv_err = 1;

  // Initialize A_curr as F before starting the iteration
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      A_curr(i,j) = F(i,j);
      A_next(i,j) = 0.0;
    }
  }

  // Do iterations
  int iter;
  for (iter = 1; iter <= iter_max; iter++)
  {
    // Store A_curr in A_inv before taking inverse
    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
        A_inv(i,j) = A_curr(i,j);
      }
    }

    invert_matrix <3> (A_inv.pointer());

    // Transpose A_inv
    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
        A_trans(i,j) = A_inv(j,i);
      }
    }

    // Set scaling factor
    real_t c;
    if (conv_err < 10e-2) {
      c = 1.0;
    } else {
      // c = 1.0;
      c = SQRT(frobenius_norm(A_inv.pointer(),3,3) / 
               frobenius_norm(A_curr.pointer(),3,3));
    }

    // Update A_next using the iterative formula
    real_t c_inv = 1.0 / c;
    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
        A_next(i,j) = 0.5 * ( c * A_curr(i,j) + c_inv * A_trans(i,j) ) ;
      }
    }

    // Check for convergence
    conv_err = 0.0;
    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
        real_t temp = A_next(i,j) - A_curr(i,j);
        conv_err += temp * temp;
      }
    }
    if (conv_err < tol) {
      break;
    }

    // Update A_curr for the next iteration
    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
        A_curr(i,j) = A_next(i,j);
      }
    }
  }
  // printf("\nTotal iter = %d\n", iter);

  if (iter > iter_max) {
    printf("Max number of iterations exceeded in %s\n", __FUNCTION__);
    exit(1);
  }

  // The resulting matrix A_next is the rotation matrix R
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      R(i,j) = A_next(i,j);
    }
  }

  // Calculate V using the relationship F = VR, V = F * R^{-1}
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      R_inv(i,j) = R(i,j);
    }
  }
  invert_matrix <3> (R_inv.pointer());

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      V(i,j) = 0.0;
      for (int k = 1; k <= 3; k++) {
        V(i,j) += F(i,k) * R_inv(k,j);
      }   
    }   
  }

  return;
}


KOKKOS_FUNCTION
real_t frobenius_norm(const real_t* A, int m, int n) {
  real_t norm = 0.0;

  for (int i = 0; i < m * n; i++) {
    norm += A[i] * A[i];
  }

  return sqrt(norm);
}