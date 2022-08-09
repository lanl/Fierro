#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"
#include "lin_alg.h"
#include "bernstein_polynomials.cpp"

using namespace utils;

void bernstein_vandermonde(ViewCArray <real_t> &B){
       // --- evaluate the basis at the nodal positions
  for(int j = 0; j < elem.num_basis(); j++){
    for(int i = 0; i < elem.num_basis(); i++){

      auto point = CArray <real_t> (elem.num_basis());

      // Get the nodal coordinates
      for(int dim = 0; dim < mesh.num_dim(); dim++){
        point(dim) = ref_elem.ref_node_positions(i, dim);
        B(i,j,dim) = bernstein::eval(elem.num_basis(), j, point(dim));
      } // end loop over dim
    } // end loop over i
  }// end loop over j
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
