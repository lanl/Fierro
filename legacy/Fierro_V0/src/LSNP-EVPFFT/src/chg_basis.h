

// ************************************************************************
//     SUBROUTINE CHG_BASIS    --->   VERSION 19/JUL/01
//
//     (modif. 06/FEB/98 - same convention as SELFPOLY - C.N.T.)
//     (modif. 16/JUN/99 - same convention as Maudlin  - C.N.T.)
//     (modif. 10/MAY/01 - KDIM version - R.L.)
//
//     KDIM=5 or 6, FOR DEVIATORIC or DEV+HYDROST TENSORS, RESPECTIVELY.
//     IOPT=0: DEFINES A BASIS OF 6 SECOND ORDER TENSORS B(N).
//     IOPT=1: CALCULATES SECOND ORDER TENSOR 'C2' AS AN EXPANSION IN TERMS
//             OF VECTOR COMPONENTS CE2(KDIM) AND THE BASIS TENSORS B(KDIM).
//     IOPT=2: CALCULATES COMPONENTS OF C2 AS A VECTOR CE2(KDIM).
//     IOPT=3: CALCULATES FOURTH ORDER TENSOR 'C4' AS AN EXPANSION IN TERMS
//             OF MATRIX COMPONENTS CE4(K,K) AND THE BASIS TENSORS B(KDIM).
//     IOPT=4: CALCULATES MATRIX COMPONENTS CE4(K,K) OF TENSOR 'C4'.
// **************************************************************************

#pragma once

#include "definitions.h"

using namespace utils;

class ChgBasis
{
  private:
    const real_t SQR2_ = 1.41421356237309;
    const real_t RSQ2_ = 0.70710678118654744;
    const real_t RSQ3_ = 0.57735026918962584;
    const real_t RSQ6_ = 0.40824829046386304;

    MatrixTypeRealDual B_basis_;
        
  public:
    ChgBasis();

    void set_b_basis();

    KOKKOS_FUNCTION
    void chg_basis_1(real_t *CE2_, real_t *C2_, int IOPT, int KDIM, real_t *B_basis_ptr) const;

    KOKKOS_FUNCTION
    void chg_basis_2(real_t *CE2_, real_t *C2_, int IOPT, int KDIM, real_t *B_basis_ptr) const;

    KOKKOS_FUNCTION
    void chg_basis_3(real_t *CE4_, real_t *C4_, int IOPT, int KDIM, real_t *B_basis_ptr) const;

    KOKKOS_FUNCTION
    void chg_basis_4(real_t *CE4_, real_t *C4_, int IOPT, int KDIM, real_t *B_basis_ptr) const;

    KOKKOS_FUNCTION
    real_t* B_basis_device_pointer() const;

    real_t* B_basis_host_pointer() const;
};
