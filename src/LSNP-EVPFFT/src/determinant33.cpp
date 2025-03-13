#include "determinant33.h"

KOKKOS_FUNCTION
void determinant33(real_t *mat_, real_t &det_mat)
{

  ViewMatrixTypeReal mat(mat_, 3,3);

  det_mat = 
     mat(1,1)*mat(2,2)*mat(3,3) 
    +mat(1,2)*mat(2,3)*mat(3,1) 
    +mat(1,3)*mat(2,1)*mat(3,2) 
    -mat(3,1)*mat(2,2)*mat(1,3) 
    -mat(3,2)*mat(2,3)*mat(1,1) 
    -mat(3,3)*mat(2,1)*mat(1,2);

  
}
