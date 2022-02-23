// build and invert the mass matrix for each element in the mesh


#include <math.h>
#include <iostream>

#include "utilities.h"
#include "state.h"
#include "geometry/geometry.h"
#include "slam/slam.h"
#include "variables.h"

using namespace utils;


void calc_average_density(){
    
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        
        // initialize to 0
        elem_state.avg_density(elem_gid) = 0.0;
        real_t elem_vol = 0.0;
        
        // tally up all quadrature point values for element
        for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
            
            int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
            
            real_t vol_gauss = ref_elem.ref_node_g_weights(gauss_lid)*mesh.gauss_pt_det_j(gauss_gid);
            elem_vol += vol_gauss;
            elem_state.avg_density(elem_gid) += mat_pt.mass(gauss_gid);
            
        } // end loop gauss
        
        // the average element value
        elem_state.avg_density(elem_gid) = elem_state.avg_density(elem_gid)/elem_vol;
        
    } // end elem_gid loop
    
} // end calc avg specific total energy


// calculate the density from strong mass
void strong_mass_dg(){
    
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        
        // tally up all quadrature point values for element
        for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
            
            int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
            
            real_t vol_gauss = ref_elem.ref_node_g_weights(gauss_lid)*mesh.gauss_pt_det_j(gauss_gid);
            mat_pt.density(gauss_gid) = mat_pt.mass(gauss_gid)/vol_gauss;
            
        } // end loop gauss
        
    } // end elem_gid loop
    
} // end strong_mass_dg


void mass_mat_inverse(){

	int num_basis = ref_elem.num_basis();

	// used for LU problem
    int singular = 0;
    int permutation;   

	for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){


		real_t src_a[num_basis*num_basis];
		auto mass_mat = ViewCArray <real_t>(&src_a[0], num_basis, num_basis);

		// Initialize mass matrix to zero 
		for(int i=0; i<num_basis; i++){
			for(int j=0; j<num_basis; j++){
				mass_mat(i,j) = 0.0;
			}
		}

		for(int basis_m = 0; basis_m < num_basis; basis_m++){

			for(int basis_n = 0; basis_n < num_basis; basis_n++){

				elem_state.mass_mat_inv(elem_gid, basis_m, basis_n) = 0.0;

				for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

					int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);

					mass_mat(basis_m, basis_n) +=
						  mat_pt.density(gauss_gid)
						* ref_elem.ref_nodal_basis(gauss_lid, basis_m)
						* ref_elem.ref_nodal_basis(gauss_lid, basis_n)
						* mesh.gauss_pt_det_j(gauss_gid)
						* ref_elem.ref_node_g_weights(gauss_lid);

				} // end loop over gauss in element
			} // end loop over basis_n
		} // end loop over basis_m

		// Print out mass matrix
		// std::cout.precision(16);
  //       std::cout <<"mass mat for elem_gid = "<< elem_gid << std::endl;

  //       for(int basis_m = 0; basis_m < ref_elem.num_basis(); basis_m++){
  //           for(int basis_n = 0; basis_n < ref_elem.num_basis(); basis_n++){

  //               real_t val = mass_mat(basis_m, basis_n);

  //               std::cout << val <<", "; 
  //           }
  //           std::cout << std::endl;
  //       }

		auto mass_mat_lu = ViewCArray <real_t> (&mass_mat(0,0), num_basis, num_basis);

		int index_a[num_basis];
		auto index = ViewCArray <int>(&index_a[0], num_basis);
		
		for(int i=0; i<num_basis; i++) index(i) = 0;
		int parity = 0;

		// Get the LU decomposition of the Mass Matrix
    	singular = LU_decompos(mass_mat_lu, index, parity, num_basis);    // mass_mat is returned as the LU matrix  


    	real_t col_a[num_basis];
    	auto col = ViewCArray <real_t> (&col_a[0], num_basis);
    	for(int i=0; i<num_basis; i++) col(i) = 0;

    	auto mass_inv = ViewCArray <real_t> (&elem_state.mass_mat_inv(elem_gid,0,0), num_basis, num_basis);
    
    	LU_invert(mass_mat, index, mass_inv, col, num_basis);

    	// std::cout <<"mass mat inverse for elem_gid = "<< elem_gid << std::endl;

     //    for(int basis_m = 0; basis_m < ref_elem.num_basis(); basis_m++){
     //        for(int basis_n = 0; basis_n < ref_elem.num_basis(); basis_n++){

     //            real_t val = elem_state.mass_mat_inv(elem_gid,basis_m,basis_n);

     //            std::cout << val <<", "; 
     //        }
     //        std::cout << std::endl;
     //    }


	} // end of loop over elements
} // end of routine
