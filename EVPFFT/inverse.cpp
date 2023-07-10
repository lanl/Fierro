#include <stdio.h>
#include "inverse.h"

//
// Optimized pow for integer exponents only
//
KOKKOS_FUNCTION
double optimizedPow(double base, int exponent) {
    if (exponent == 0)
        return 1.0;
    else if (exponent == 1)
        return base;

    double result = 1.0;
    bool negativeExponent = (exponent < 0);
    exponent = abs(exponent);

    while (exponent > 0) {
        if (exponent & 1)
            result *= base;
        base *= base;
        exponent >>= 1;
    }

    return negativeExponent ? 1.0 / result : result;
}


//
// New Gauss-Jordan elimination matrix inverse functions 
//
KOKKOS_FUNCTION
void swap_rows(double *matrix, int row1, int row2, int n) {
    for (int i = 0; i < n; i++) {
        double temp = *(matrix + row1 * n + i);
        *(matrix + row1 * n + i) = *(matrix + row2 * n + i);
        *(matrix + row2 * n + i) = temp;
    }
}

KOKKOS_FUNCTION
void scale_row(double *matrix, int row, double factor, int n) {
    for (int i = 0; i < n; i++) {
        *(matrix + row * n + i) *= factor;
    }
}

KOKKOS_FUNCTION
void add_rows(double *matrix, int src_row, int dest_row, double factor, int n) {
    for (int i = 0; i < n; i++) {
        *(matrix + dest_row * n + i) += *(matrix + src_row * n + i) * factor;
    }
}

KOKKOS_FUNCTION
int invert_matrix(double *matrix, int n) {
    double *identity = new double[n*n];
    if (identity == NULL) {
        printf("Error: Failed to allocate memory for identity matrix.\n");
        return 0;
    }

    // Initialize identity matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            *(identity + i * n + j) = (i == j) ? 1.0 : 0.0;
        }
    }

    // Perform Gauss-Jordan elimination
    for (int i = 0; i < n; i++) {
        if (*(matrix + i * n + i) == 0) {
            //printf("Error: Matrix is not invertible.\n");
            delete [] identity;
            return 0;
        }

        double pivot = *(matrix + i * n + i);
        scale_row(matrix, i, 1.0 / pivot, n);
        scale_row(identity, i, 1.0 / pivot, n);

        for (int j = 0; j < n; j++) {
            if (j != i) {
                double factor = -*(matrix + j * n + i);
                add_rows(matrix, i, j, factor, n);
                add_rows(identity, i, j, factor, n);
            }
        }
    }

    // Copy the inverted matrix to the input matrix pointer
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            *(matrix + i * n + j) = *(identity + i * n + j);
        }
    }

    delete [] identity;
    return 1;
}
// End New Gauss-Jordan elimination matrix inverse functions



//
//  inverse_gj
//
KOKKOS_FUNCTION
void inverse_gj(real_t *a_, int n)
{
  invert_matrix(a_, n);
#if 0
  ViewMatrixTypeReal a(a_,n,n);
  real_t             tmp;
  real_t             fac;

  for (int i = 1; i <= n; i++) { // loop 1
    fac = 1.0 / a(i,i); // This would become inverse if done on blocks
    a(i,i) = 1.0;
    for (int j = 1; j <= n; j++) {
      a(j,i) = a(j,i)*fac; 
    }

    // NOTE:These loops don't shrink as i increases so load should
    // be balanced?
    for (int k = 1; k <= i-1; k++) {
      tmp = a(i,k);
      a(i,k) = 0.0;
      for (int j = 1; j <= n; j++) {
        a(j,k) += -a(j,i)*tmp;
      }
    }
    for (int k = i+1; k <= n; k++) {
      tmp = a(i,k);
      a(i,k) = 0.0;
      for (int j = 1; j <= n; j++) {
        a(j,k) += -a(j,i)*tmp;
      }
    } 
  } // end loop 1
#endif
}

//
//  lu_inverse
//
KOKKOS_FUNCTION
void lu_inverse(real_t *a_, int n)
{
  real_t              *y_ = new real_t[n*n];
  int                 *indx_ = new int[n];
  ViewMatrixTypeReal  y(y_,n,n);
  ViewMatrixTypeInt   indx(indx_,n);
  ViewMatrixTypeReal  a(a_,n,n);
  real_t              d = 0.0;
  int                 isingular = 0;

  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      y(i,j) = 0.0;
    }
    y(i,i) = 1.0;
  }

  ludcmp(a_,n,n,indx_,d,isingular);

  if (isingular == 0) {
    for (int j = 1; j <= n; j++) {
      lubksb(a_,n,n,indx_,&y(1,j));
    }

    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= n; j++) {
        a(i,j) = y(i,j);
      }
    }
  }

  delete [] y_;
  delete [] indx_;

  return;
}

//
//  ludcmp
//
KOKKOS_FUNCTION
void ludcmp(real_t *a_, int n, int np, int *indx_, real_t d, int &isingular)
{
  ViewMatrixTypeReal  a(a_,n,n);
  ViewMatrixTypeInt   indx(indx_,n);
  int                 imax;
  const int           nmax = 500;
  real_t              vv_[nmax];
  ViewMatrixTypeReal  vv(vv_,nmax);
  real_t              aamax;
  real_t              dum;
  real_t              sum;

  d = 1.0;
  for (int i = 1; i <= n; i++) { // loop 12
    aamax = 0.0;
    for (int j = 1; j <= n; j++) { // loop 11
      if (abs(a(i,j)) > aamax) aamax = abs(a(i,j));
    } // end loop 11
    if (aamax == 0.0) {
      //printf("singular matrix in ludcmp\n");
      //PAUSE;
    } // end if (aamax == 0.0)
    if (aamax == 0.0) {
      isingular = 1;
      return;
    }
    vv(i) = 1.0 / aamax;
  } // end loop 12

  for (int j = 1; j <= n; j++) { // loop 19

    for (int i = 1; i <= j-1; i++) { // loop 14
      sum = a(i,j);
      for (int k = 1; k <= i-1; k++) { // loop 13
        sum += - a(i,k)*a(k,j);
      } // end loop 13
      a(i,j) = sum;
    } // end loop 14
    aamax = 0.0;

    for (int i = j; i <= n; i++) { // loop 16
      sum = a(i,j);
      for (int k = 1; k <= j-1; k++) { // loop 15
        sum += - a(i,k)*a(k,j);
      } // end loop 15
      a(i,j) = sum;
      dum = vv(i)*abs(sum);
      if (dum >= aamax) {
        imax = i;
        aamax = dum;
      } // end if (dum >= aamax)
    } // end loop 16

    if (j != imax) {
      for (int k = 1; k <= n; k++) { // loop 17
        dum = a(imax,k);
        a(imax,k) = a(j,k);
        a(j,k) = dum;
      } // end loop 17
      d = -d;
      vv(imax) = vv(j);
    } // end if (j != imax)
    indx(j) = imax;

    //if (a(j,j) == 0.0) a(j,j) = TINY;
    if (aamax == 0.0) {
      isingular = 1;
      return;
    }

    if (j != n) {
      dum = 1.0 / a(j,j);
      for (int i = j+1; i <= n; i++) { // loop 18
        a(i,j) = a(i,j)*dum;
      } // end loop 18
    } // end if (j != n)

  } // end loop 19

  isingular = 0;
}

//
//  lubksb
//
KOKKOS_FUNCTION
void lubksb(real_t *a_, int n, int np, int *indx_, real_t *b_)
{
  ViewMatrixTypeReal  a(a_,n,n);
  ViewMatrixTypeInt   indx(indx_,n);
  ViewMatrixTypeReal  b(b_,n);
  real_t              sum;
  int                 ii;
  int                 ll;

  ii = 0;
  for (int i = 1; i <= n; i++) { // loop 12
    ll = indx(i);
    sum = b(ll);
    b(ll) = b(i);
    if (ii != 0) {
      for (int j  = ii; j <= i-1; j++) { // loop 11
        sum += - a(i,j)*b(j);
      } // end loop 11
    } else if (sum != 0.0) {
      ii = i;
    } // end else if (sum != 0.0)
    b(i) = sum;
  } // end loop 12

  for (int i = n; i >= 1; i--) { // loop 14
    sum = b(i);
    for (int j = i+1; j <= n; j++) { // loop 13
      sum += - a(i,j)*b(j);
    } // end loop 13
    b(i) = sum / a(i,i);
  } // end loop 14
}

//
//  inverse
//
void inverse(real_t *a_, int n, int np, real_t *y_)
{
  int                 indx_[n];
  ViewMatrixTypeReal  y(y_,n,n);
  ViewMatrixTypeInt   indx(indx_,n);
  ViewMatrixTypeReal  a(a_,n,n);
  real_t              d = 0.0;

  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      y(i,j) = 0.0;
    }
    y(i,i) = 1.0;
  }

  ludcmpc(a_,n,n,indx_,d);

  for (int j = 1; j <= n; j++) {
    lubksbc(a_,n,n,indx_,&y(1,j));
  }

  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      a(i,j) = y(i,j);
    }
  }
}

//
//  ludcmpc
//
void ludcmpc(real_t *a_, int n, int np, int *indx_, real_t d)
{
  ViewMatrixTypeReal  a(a_,n,n);
  ViewMatrixTypeInt   indx(indx_,n);
  int                 imax;
  const int           nmax = 500;
  real_t              vv_[nmax];
  ViewMatrixTypeReal  vv(vv_,nmax);
  real_t              aamax;
  real_t              dum;
  real_t              sum;
  const real_t        tiny = 1.0e-20;

  d = 1.0;
  for (int i = 1; i <= n; i++) { // loop 12
    aamax = 0.0;
    for (int j = 1; j <= n; j++) { // loop 11
      if (abs(a(i,j)) > aamax) aamax = abs(a(i,j));
    } // end loop 11
    if (aamax == 0.0) {
      printf("singular matrix in ludcmp\n");
      PAUSE;
    } // end if (aamax == 0.0)
    vv(i) = 1.0 / aamax;
  } // end loop 12

  for (int j = 1; j <= n; j++) { // loop 19

    for (int i = 1; i <= j-1; i++) { // loop 14
      sum = a(i,j);
      for (int k = 1; k <= i-1; k++) { // loop 13
        sum += - a(i,k)*a(k,j);
      } // end loop 13
      a(i,j) = sum;
    } // end loop 14
    aamax = 0.0;

    for (int i = j; i <= n; i++) { // loop 16
      sum = a(i,j);
      for (int k = 1; k <= j-1; k++) { // loop 15
        sum += - a(i,k)*a(k,j);
      } // end loop 15
      a(i,j) = sum;
      dum = vv(i)*abs(sum);
      if (dum >= aamax) {
        imax = i;
        aamax = dum;
      } // end if (dum >= aamax)
    } // end loop 16

    if (j != imax) {
      for (int k = 1; k <= n; k++) { // loop 17
        dum = a(imax,k);
        a(imax,k) = a(j,k);
        a(j,k) = dum;
      } // end loop 17
      d = -d;
      vv(imax) = vv(j);
    } // end if (j != imax)
    indx(j) = imax;

    if (a(j,j) == 0.0) a(j,j) = tiny;

    if (j != n) {
      dum = 1.0 / a(j,j);
      for (int i = j+1; i <= n; i++) { // loop 18
        a(i,j) = a(i,j)*dum;
      } // end loop 18
    } // end if (j != n)

  } // end loop 19

}

//
//  lubksb
//
void lubksbc(real_t *a_, int n, int np, int *indx_, real_t *b_)
{
  // No need to retype the coding in lubksb for lubksbc
  lubksb(a_, n, np, indx_, b_);
}
