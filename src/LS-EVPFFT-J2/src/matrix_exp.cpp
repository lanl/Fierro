#include "matrix_exp.h"

KOKKOS_FUNCTION
void matrix_exp(real_t *mat_, real_t *mat_exp_)
{

  ViewMatrixTypeReal mat(mat_, 3,3);
  ViewMatrixTypeReal mat_exp(mat_exp_, 3,3);

  real_t mat2_[3*3];
  real_t mat3_[3*3];
  ViewMatrixTypeReal mat2(mat2_, 3,3);
  ViewMatrixTypeReal mat3(mat3_, 3,3);

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      mat2(i,j) = 0.0;
      for (int k = 1; k <=3; k++){
        mat2(i,j) += mat(i,k)*mat(k,j);
      }
    }
  }

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      mat3(i,j) = 0.0;
      for (int k = 1; k <=3; k++){
        mat3(i,j) += mat2(i,k)*mat(k,j);
      }
    }
  }

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      mat_exp(i,j) = mat(i,j) + mat2(i,j)/2.0 + mat3(i,j)/6.0;
    }
    mat_exp(i,i) = mat_exp(i,i) + 1.0;
  }
  
}
