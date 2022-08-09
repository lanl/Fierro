#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"
#include "lin_alg.h"
#include "bernstein_polynomials.cpp"

using namespace utils;

void bernstein_vandermonde(ViewCArray <real_t> &B){
       
  for(int j = 0; j < elem.num_basis(); j++){
    for (int i = 0; i < elem.num_basis(); i++){
      int ref_position = elem.vert_node_map(i);
      B(i,j) = ref_elem.ref_nodal_basis(ref_position,j);
    }// end loop over i
  } // end loop over j
}// end B-V matrix


void BV_inv(ViewCArray <real_t> &B, ViewCArray <real_t> &B_inv){
  
  int index_a[elem.num_basis()];
  auto index = ViewCArray <int> (&index_a[0], elem.num_basis());
  for(int i=0; i < elem.num_basis(); i++) index(i) = 0;


  real_t col_a[elem.num_basis()];
  auto col = ViewCArray <real_t> (&col_a[0], elem.num_basis());
  for(int i=0; i < elem.num_basis(); i++) col(i) = 0;

  LU_invert(B, index, B_inv, col, elem.num_basis());
 
}// end BV inverse
