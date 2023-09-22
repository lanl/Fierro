#include "chg_basis.h"
#include "utilities.h"

ChgBasis::ChgBasis()
  : B_basis_ (3,3,6)
{
  set_b_basis();
}

void ChgBasis::set_b_basis()
{

// *** CALCULATES BASIS TENSORS B(N)
  for (int k = 1; k <= 6; k++) {
    for (int j = 1; j <= 3; j++) {
      for (int i = 1; i <= 3; i++) {
        B_basis_.host(i,j,k) = 0.0;
      }
    }
  }
  
  B_basis_.host(1,1,2) = -RSQ6_;
  B_basis_.host(2,2,2) = -RSQ6_;
  B_basis_.host(3,3,2) =  real_t(2.0)*RSQ6_;

  B_basis_.host(1,1,1) = -RSQ2_;
  B_basis_.host(2,2,1) =  RSQ2_;

  B_basis_.host(2,3,3) =  RSQ2_;
  B_basis_.host(3,2,3) =  RSQ2_;

  B_basis_.host(1,3,4) =  RSQ2_;
  B_basis_.host(3,1,4) =  RSQ2_;

  B_basis_.host(1,2,5) =  RSQ2_;
  B_basis_.host(2,1,5) =  RSQ2_;

  B_basis_.host(1,1,6) =  RSQ3_;
  B_basis_.host(2,2,6) =  RSQ3_;
  B_basis_.host(3,3,6) =  RSQ3_;

  // update device
  B_basis_.update_device();
  Kokkos::fence();

}


KOKKOS_FUNCTION
void ChgBasis::chg_basis_1(real_t *CE2_, real_t *C2_, int IOPT, int KDIM, real_t *B_basis_ptr) const
{
  assert(IOPT == 1 && "ERROR: IOPT should be 1 for call to chg_basis_1");

  ViewMatrixTypeReal CE2 (CE2_, KDIM);
  ViewMatrixTypeReal C2  (C2_, 3,3);
  ViewMatrixTypeReal B_basis (B_basis_ptr, 3,3,6);

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      C2(i,j) = 0.0;
      for (int n = 1; n <= KDIM; n++) {
        C2(i,j) += CE2(n) * B_basis(i,j,n);
      }
    }
  }

}


KOKKOS_FUNCTION
void ChgBasis::chg_basis_2(real_t *CE2_, real_t *C2_, int IOPT, int KDIM, real_t *B_basis_ptr) const
{
  assert(IOPT == 2 && "ERROR: IOPT should be 2 for call to chg_basis_2");

  ViewMatrixTypeReal CE2 (CE2_, KDIM);
  ViewMatrixTypeReal C2  (C2_, 3,3);
  ViewMatrixTypeReal B_basis (B_basis_ptr, 3,3,6);

  for (int n = 1; n <= KDIM; n++) {
    CE2(n) = 0.0;
    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
        CE2(n) += C2(i,j) * B_basis(i,j,n);
      }
    }
  }

}


KOKKOS_FUNCTION
void ChgBasis::chg_basis_3(real_t *CE4_, real_t *C4_, int IOPT, int KDIM, real_t *B_basis_ptr) const
{
  assert(IOPT == 3 && "ERROR: IOPT should be 3 for call to chg_basis_3");

  ViewMatrixTypeReal CE4 (CE4_, KDIM,KDIM);
  ViewMatrixTypeReal C4  (C4_, 3,3,3,3);
  ViewMatrixTypeReal B_basis (B_basis_ptr, 3,3,6);

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      for (int k = 1; k <= 3; k++) {
        for (int l = 1; l <= 3; l++) {
          C4(i,j,k,l) = 0.0;
          for (int n = 1; n <= KDIM; n++) {
            for (int m = 1; m <= KDIM; m++) {
              C4(i,j,k,l) += CE4(n,m) * B_basis(i,j,n) * B_basis(k,l,m);
            }
          }
        }
      }
    }
  }

}


KOKKOS_FUNCTION
void ChgBasis::chg_basis_4(real_t *CE4_, real_t *C4_, int IOPT, int KDIM, real_t *B_basis_ptr) const
{
  assert(IOPT == 4 && "ERROR: IOPT should be 4 for call to chg_basis_4");

  ViewMatrixTypeReal CE4 (CE4_, KDIM,KDIM);
  ViewMatrixTypeReal C4  (C4_, 3,3,3,3);
  ViewMatrixTypeReal B_basis (B_basis_ptr, 3,3,6);

  for (int n = 1; n <= KDIM; n++) {
    for (int m = 1; m <= KDIM; m++) {
      CE4(n,m) = 0.0;
      for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
          for (int k = 1; k <= 3; k++) {
            for (int l = 1; l <= 3; l++) {
              CE4(n,m) += C4(i,j,k,l) * B_basis(i,j,n) * B_basis(k,l,m);
            }
          }
        }
      }
    }
  }

}

KOKKOS_FUNCTION
real_t* ChgBasis::B_basis_device_pointer() const
{
  return B_basis_.device_pointer();
}

real_t* ChgBasis::B_basis_host_pointer() const
{
  return B_basis_.host_pointer();
}
