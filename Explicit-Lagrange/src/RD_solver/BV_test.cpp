/* test LU BV inv */

#include <iostream>
#include <cmath>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"
#include "lin_alg.h"
#include "bernstein_polynomials.cpp"


using namespace utils;

bool Bern_test(ViewCArray <real_t> &B, ViewCArray <real_t> &B_inv){
  
  bool pass =false;

  real_t eps = 1e-15;

  auto Id_left = CArray <real_t> (ref_elem.num_basis(), ref_elem.num_basis(), mesh.num_dim());
  auto Id_right = CArray <real_t> (ref_elem.num_basis(), ref_elem.num_basis(), mesh.num_dim());

  for (int dim = 0; dim < mesh.num_dim(); dim++){
    for (int j = 0; j < ref_elem.num_basis(); j++){
      for( int i = 0; i < ref_elem.num_basis(); i++){
        Id_left(i,j,dim) = B_inv(j,i,dim)*B(i,j,dim);
	Id_right(i,j,dim) = B(j,i,dim)*B_inv(i,j,dim);
	if ( i==j){
	   real_t one = Id_left(i,j)*Id_right(i,j);
	   if ( abs( one - 1.0 ) > eps ){
             pass = false;
	     std::cout<< "not one" << std::endl;
	     std::cout<< "value is " << one << std::endl;
	   } // end if id false
	   else{
	     pass = true;
	   }// end else
	}// end if i=j
	if ( i != j){
	  real_t zero = Id_left(i,j)*Id_right(i,j);
	  if (abs(zero) > eps){
            pass = false;
	    std::cout << " not zero " << std::endl;
	    std::cout << " value is " << zero << std::endl;
          }//end if not zero 
	  else{
	    pass = true;
	  }// end else
	}// end if i!=j
      }// end loop over i
    } // end loop over j
  }// end loop over dim 
  
  return pass;


}// end test


int main(){
 
  std::cout.precision(15);
  std::cout << "test bernstein-vandermonde inverse" <<std::endl;
  
  bool pass = test();

  std::cout << int(pass) << std::endl;

  return 0;

}//end main

