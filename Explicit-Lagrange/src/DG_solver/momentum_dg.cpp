#include <math.h>
#include <iostream>
#include <cmath> //abs
#include <algorithm> //min
#include <string.h>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"
#include "slam.h"

#define PI 3.14159265

using namespace utils;

void momentum_dg(real_t rk_alpha, int cycle){

	// Build the RHS vector for calculating the new vertex velocity


	/*
	std::string t1 = "limit_velocity_";
	std::string t2 = std::to_string(cycle);
	std::string combo_title = t1+t2;

	const char *title = combo_title.c_str();
	outfile = fopen(title, "w");

	fprintf(outfile, "elem_gid, r,  mega_alpha\n");
	*/

	//print file at the last cycle to show ke, mega alpha

	// FILE *out_alpha;
	// out_alpha = fopen("vel_alpha_end", "w");
	// fprintf(out_alpha, "elem_gid, gauss_gid, r, alpha, r, ke, r, avg_ke\n");

    // conservation check
    real_t RHS_momentum_check[mesh.num_dim()];
    for(int dim=0; dim<mesh.num_dim(); dim++){
        RHS_momentum_check[dim] = 0.0;
    }
    
	// Loop over the elements in the mesh.num_elems()
	for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){

		// Left and right parts of the RHS
		auto surface_int = CArray <real_t>(ref_elem.num_basis(), mesh.num_dim());
		auto volume_int = CArray <real_t>(ref_elem.num_basis(), mesh.num_dim());

		// RHS vector
	    auto RHS = CArray <real_t> (ref_elem.num_basis(), mesh.num_dim());

		// Initialize to zero
		for(int basis_id = 0; basis_id < ref_elem.num_basis(); basis_id++){
			for(int dim=0; dim<mesh.num_dim(); dim++){
				surface_int(basis_id, dim) = 0.0;
				volume_int(basis_id, dim) = 0.0;
				RHS(basis_id, dim) = 0.0;
			}
		}

		// loop over the basis associated with this element
		for(int basis_id = 0; basis_id < ref_elem.num_basis(); basis_id++){
	        
			// loop over the cells in the element and corners/nodes in the cell
			// to calculate \Phi a_{i} (\sigma^* n_{i})^T, which is the same 
			// as the basis at the node times the corner force 
	        

	        // Calculate the surface integral
	        for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){ 
	            
	            int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);

	            // Loop over the nodes/corners in the cell
	            for(int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){

	                // get global ids for nodes and corners
	                int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);

	                int corner_lid = node_lid;
	                int corner_gid = mesh.corners_in_cell(cell_gid, corner_lid);  // node_lid = corner_lid

	                // Get the id for the node in the reference element
	                int node_rid = ref_elem.cell_nodes_in_elem(cell_lid, node_lid); 

	                // Sum over all corners
	                for(int dim=0; dim<mesh.num_dim(); dim++){
	                	surface_int(basis_id, dim) += ref_elem.ref_nodal_basis(node_rid, basis_id) * corner.force(corner_gid, dim);
	                }

	            }// end loop over nodes/corners in a cell
	        } // end loop over cells in an element

	        // Calculate volume integral
	        for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

	            int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
	            int node_rid = gauss_lid; // Leveraging internal structured geometry

	            // Hope you like indicial notation :)
                for(int k=0; k<mesh.num_dim(); k++){
                    for(int j=0; j<mesh.num_dim(); j++){
                        for(int i=0; i<mesh.num_dim(); i++){

                            volume_int(basis_id, k) +=
                                ref_elem.ref_nodal_gradient(node_rid, basis_id, i)
                                * mesh.gauss_pt_jacobian_inverse(gauss_gid, i, j)
                                * mat_pt.stress(gauss_gid, j, k)
                                * mesh.gauss_pt_det_j(gauss_gid)
                                * ref_elem.ref_node_g_weights(gauss_lid);
                        } // end loop over i
                    } // end loop over j
	            } // end loop over k

	        } // end loop over the gauss points

	    } // end loop over the basis

	    // Zero if small, limits error from LU inverse for large spectral radius
	    for(int basis_id=0; basis_id < ref_elem.num_basis(); basis_id++){
	        for(int dim=0; dim<mesh.num_dim(); dim++){
                
                RHS(basis_id, dim) = surface_int(basis_id, dim) - volume_int(basis_id, dim);
                
                if(fabs(RHS(basis_id, dim)) < 5E-15) RHS(basis_id, dim) = 0.0;
                
                RHS_momentum_check[dim] += RHS(basis_id,dim);  // checking conservation
	    	}
	    }

	    // Calcuate (dU/dt) = M^{-1} * RHS
	    auto dudt = CArray <real_t> (ref_elem.num_basis(), mesh.num_dim());
		
		// Initialize to zero
		for(int basis_id=0; basis_id < ref_elem.num_basis(); basis_id++){
	    	for(int dim=0; dim<mesh.num_dim(); dim++){
	    		dudt(basis_id, dim) = 0.0;
	    	}
	    }
		
	    // Calculate the velocity per dimension
		for(int dim=0; dim<mesh.num_dim(); dim++){
			for(int basis_i=0; basis_i < ref_elem.num_basis(); basis_i++){
				for(int basis_j=0; basis_j < ref_elem.num_basis(); basis_j++){

					dudt(basis_i, dim) += 
						elem_state.mass_mat_inv(elem_gid, basis_i, basis_j)
						* RHS(basis_j, dim);

					
				} // end loop over basis j
		    } // end loop over basis i
		} // end loop over dim

		// update the momentum at the kinematic degrees of freedom
		for(int basis_id = 0; basis_id < elem.num_basis(); basis_id++){

			// get global and local index of things
			int node_lid = elem.vert_node_map(basis_id);
			int gauss_gid = mesh.gauss_in_elem(elem_gid, node_lid);
			int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

	        for (int dim = 0; dim < mesh.num_dim(); dim++){
	            
	            mat_pt.velocity(1, gauss_gid, dim) = 
	                 	mat_pt.velocity(0, gauss_gid, dim) + rk_alpha * dudt(basis_id, dim)*dt;

	        } // end loop over dim
	    } // end loop over the basis functions in an element
		
	    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

        	int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);

            // get the global id of the gauss point
            int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);

            real_t interp_vel[mesh.num_dim()]; 
            for(int i=0; i<mesh.num_dim(); i++) interp_vel[i] = 0.0;

            // Sum over the basis times the velocity defined at the basis vertex position
            for (int dim = 0; dim < mesh.num_dim(); dim++){
                for(int basis_id = 0; basis_id < elem.num_basis(); basis_id++){

                    int node_basis_id = elem.vert_node_map(basis_id);
                    int interp_gid = mesh.gauss_in_elem(elem_gid, node_basis_id);

                    interp_vel[dim] += mat_pt.velocity(1, interp_gid, dim) * ref_elem.ref_nodal_basis(gauss_lid, basis_id);
                }
            }
            
            // Save interpolated velocity back to gauss point
            for (int dim = 0; dim < mesh.num_dim(); dim++){   
                mat_pt.velocity(1, gauss_gid, dim) =  interp_vel[dim];
            }

        } // end loop over gauss_lid
        
    } // end loop over elements
    
    
    //std::cout << "######### RHS Momentum Error = " <<
    //RHS_momentum_check[0] << " , " <<
    //RHS_momentum_check[1] << " , " <<
    //RHS_momentum_check[2] << "   " << std::endl;
    
} // end of routine



// this function calculates the conservative average velocity in an element
void calc_average_velocity(){

    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        
        // initialize to 0
        for (int dim = 0; dim < mesh.num_dim(); dim++){
            elem_state.avg_velocity(elem_gid, dim) = 0.0;
        }
        real_t elem_mass = 0.0;
        
        // tally up all quadrature point values for element
        for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
            
            int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
            
            elem_mass += mat_pt.mass(gauss_gid);
            for(int dim = 0; dim < mesh.num_dim(); dim++) {
                elem_state.avg_velocity(elem_gid, dim) += mat_pt.velocity(1, gauss_gid, dim)*mat_pt.mass(gauss_gid);
            }
            
        } // end loop gauss
        
        // the average element value
        for(int dim = 0; dim < mesh.num_dim(); dim++) {
            elem_state.avg_velocity(elem_gid, dim) = elem_state.avg_velocity(elem_gid, dim)/elem_mass;
        }
        
    } // end elem_gid loop
    
} // end calc avg velocity



void gradvel_dg(){
    
    for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        
        // volume matrix
        real_t src_a[ref_elem.num_basis()* ref_elem.num_basis()]; // a 1D array that will be mapped to a 2D array
        auto elem_vol_matrx = ViewCArray <real_t> (&src_a[0], ref_elem.num_basis(), ref_elem.num_basis());
        
        // Left and right parts of the RHS
        auto surface_int = CArray <real_t>(ref_elem.num_basis(), mesh.num_dim(), mesh.num_dim());
        auto volume_int  = CArray <real_t>(ref_elem.num_basis(), mesh.num_dim(), mesh.num_dim());
        
        // RHS vector
        auto RHS = CArray <real_t> (ref_elem.num_basis(), mesh.num_dim(), mesh.num_dim());
        
        // Initialize to zero
        for(int basis_id = 0; basis_id < ref_elem.num_basis(); basis_id++){
            for (int i = 0; i < mesh.num_dim(); i++){
                for (int j = 0; j < mesh.num_dim(); j++){
                    surface_int(basis_id, i, j) = 0.0;
                    volume_int(basis_id, i, j) = 0.0;
                    RHS(basis_id, i, j) = 0.0;
                } // end for i
            } // end for j
        } // end for basis
        
        // Initialize to zero
        for(int basis_n = 0; basis_n < ref_elem.num_basis(); basis_n++){
            for(int basis_m = 0; basis_m < ref_elem.num_basis(); basis_m++){
                elem_vol_matrx(basis_n, basis_m) = 0.0;
            } // end basis_j
        } // end basis_i
        
        
        // --- Calculate the inverse of the volume matrix ---
        
        // used for LU problem
        int singular = 0;
        int permutation;
        
        int num_basis = ref_elem.num_basis();
        
        for(int basis_m = 0; basis_m < num_basis; basis_m++){
            
            for(int basis_n = 0; basis_n < num_basis; basis_n++){
                
                for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
                    
                    int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
                    
                    elem_vol_matrx(basis_m, basis_n) +=
                          ref_elem.ref_nodal_basis(gauss_lid, basis_m)
                        * ref_elem.ref_nodal_basis(gauss_lid, basis_n)
                        * mesh.gauss_pt_det_j(gauss_gid)
                        * ref_elem.ref_node_g_weights(gauss_lid);
                    
                } // end loop over gauss in element
            } // end loop over basis_n
        } // end loop over basis_m
        
        
        auto vol_matrx_lu = ViewCArray <real_t> (&elem_vol_matrx(0,0), num_basis, num_basis);
        
        int index_a[num_basis];
        auto index = ViewCArray <int>(&index_a[0], num_basis);
        
        for(int i=0; i<num_basis; i++) index(i) = 0;
        int parity = 0;
        
        // Get the LU decomposition of the Volume Matrix
        singular = LU_decompos(vol_matrx_lu, index, parity, num_basis);    // vol_matrix_lu is returned as the LU matrix
        
        real_t col_a[num_basis];
        auto col = ViewCArray <real_t> (&col_a[0], num_basis);
        for(int i=0; i<num_basis; i++) col(i) = 0;
        
        auto elem_vol_matrx_inv = CArray <real_t> (num_basis, num_basis);
        auto vol_inv = ViewCArray <real_t> (&elem_vol_matrx_inv(0,0), num_basis, num_basis);
        
        LU_invert(elem_vol_matrx, index, vol_inv, col, num_basis);

        // --- end of calculating inverse of the volume matrix ---
        
        /*
        std::cout << " elem = " << elem_gid << "\n ";
        for(int basis_m = 0; basis_m < num_basis; basis_m++){
            for(int basis_n = 0; basis_n < num_basis; basis_n++){
                //if ( fabs(elem_vol_matrx(basis_m, basis_n) - 0.000125) > 1e-12 )
                std::cout  <<  elem_vol_matrx_inv(basis_m, basis_n)  << " , " ;
            }
            std::cout << " \n ";
        }
        std::cout << " \n " << std::endl;
         */
        
        
        
        // loop over the basis associated with this element
        for(int basis_id = 0; basis_id < ref_elem.num_basis(); basis_id++){
            
            
            // Surface integral
            for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
                
                int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);
                
                // Loop over the nodes/corners in the cell
                for(int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){
                    
                    // get global ids for nodes and corners
                    int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);
                    
                    int corner_lid = node_lid;
                    int corner_gid = mesh.corners_in_cell(cell_gid, corner_lid);  // node_lid = corner_lid
                    
                    // Get the id for the node in the reference element
                    int node_rid = ref_elem.cell_nodes_in_elem(cell_lid, node_lid);
                    
                    // Sum over all corners
                    for (int i = 0; i < mesh.num_dim(); i++){
                        for (int j = 0; j < mesh.num_dim(); j++){
                        
                            surface_int(basis_id, i, j) +=
                                ref_elem.ref_nodal_basis(node_rid, basis_id)
                                * corner.normal(corner_gid, i)
                                * node.vel(1, node_gid, j); // Riemann velocity is stored at the node
                        } // end for j
                    } // end for i
                    
                } // end loop over nodes/corners in a cell
            } // end loop over cells in an element
            
            // Volume integral
            for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
                
                int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
                int node_rid = gauss_lid; // Leveraging internal structured geometry
                
                
                for(int i=0; i<mesh.num_dim(); i++){
                    for(int j=0; j<mesh.num_dim(); j++){
                        for(int k=0; k<mesh.num_dim(); k++){
                            
                            volume_int(basis_id, j, k) +=
                                ref_elem.ref_nodal_gradient(node_rid, basis_id, i)
                                * mesh.gauss_pt_jacobian_inverse(gauss_gid, i, j)
                                * mat_pt.velocity(1, gauss_gid, k)
                                * mesh.gauss_pt_det_j(gauss_gid)
                                * ref_elem.ref_node_g_weights(gauss_lid);
                        } // end for k
                    } // end for j
                } // end for i
                
            } // end loop over the gauss points
            
        } // end loop over the basis
        
        
        // save the RHS tensor
        for(int basis_id=0; basis_id < ref_elem.num_basis(); basis_id++){
            for(int i=0; i<mesh.num_dim(); i++){
                for(int j=0; j<mesh.num_dim(); j++){
                    RHS(basis_id,i,j) = surface_int(basis_id,i,j) - volume_int(basis_id,i,j);
                    
                    if (fabs( RHS(basis_id,i,j) ) <= 5E-15 ) RHS(basis_id,i,j) = 0;
                    
                } // end for j
            } // end for i
        } // end loop over the basis
        
        
        
        // Calcuate L = M^{-1} * RHS
        auto L = CArray <real_t> (ref_elem.num_basis(),mesh.num_dim(),mesh.num_dim());
        
        // Initialize to zero
        for(int basis_id=0; basis_id < ref_elem.num_basis(); basis_id++){
            for(int i=0; i<mesh.num_dim(); i++){
                for(int j=0; j<mesh.num_dim(); j++){
                    L(basis_id, i, j) = 0.0;
                } // end for j
            } // end for i
        } // end loop over basis
        
        // Calculate the time rate of change of total energy
        for(int i=0; i<mesh.num_dim(); i++){
            for(int j=0; j<mesh.num_dim(); j++){
                for(int basis_m=0; basis_m < ref_elem.num_basis(); basis_m++){
                    for(int basis_n=0; basis_n < ref_elem.num_basis(); basis_n++){
                
                        L(basis_m, i, j) +=
                            elem_vol_matrx_inv(basis_m, basis_n)
                            * RHS(basis_n, i, j);
                        
                    } // end basis n
                    
                    if (fabs(  L(basis_m,i,j) ) <= 5E-15 ) L(basis_m,i,j) = 0;
                    
                } // end basis m
            } // end for j
        } // end for i

        
        // update the grad vel at the kinematic degrees of freedom
        for(int basis_id = 0; basis_id < elem.num_basis(); basis_id++){
        
            int node_basis_id = elem.vert_node_map(basis_id);
            int gauss_gid = mesh.gauss_in_elem(elem_gid, node_basis_id);
            
            for (int i = 0; i < mesh.num_dim(); i++){
                for (int j = 0; j < mesh.num_dim(); j++){
                    mat_pt.grad_vel(gauss_gid, i, j) = L(basis_id, i, j);
                } // end for j
            } // end for i
            
        } // end loop over the basis and associated gauss points
        
        
        // Evaluate the specific volume polynomal at the nodes
        for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
            
            // get the global id of the quadrature point
            int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
            
            int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);
            
            
            auto interp_L = CArray <real_t> (mesh.num_dim(),mesh.num_dim());
            
            // initialize to 0
            for (int i = 0; i < mesh.num_dim(); i++){
                for (int j = 0; j < mesh.num_dim(); j++){
                    interp_L(i,j) = 0;
                } // end for j
            }// end for i
            
            
            // Sum over the basis times the L defined at the basis vertex position
            for(int basis_id = 0; basis_id < elem.num_basis(); basis_id++){
                
                int node_basis_id = elem.vert_node_map(basis_id);
                int interp_gid = mesh.gauss_in_elem(elem_gid, node_basis_id);
                
                for (int i = 0; i < mesh.num_dim(); i++){
                    for (int j = 0; j < mesh.num_dim(); j++){
                        
                        interp_L(i,j) += mat_pt.grad_vel(interp_gid, i, j)
                            * ref_elem.ref_nodal_basis(gauss_lid, basis_id);
                        
                    } // end for j
                }// end for i
                
            } // end for basis_id
            
            
            
            // Save interpolated L back to gauss point grad_vel and div_vel
            mat_pt.div_vel(gauss_gid) = 0.0;
            for (int i = 0; i < mesh.num_dim(); i++){
                
                for (int j = 0; j < mesh.num_dim(); j++){
                    mat_pt.grad_vel(gauss_gid, i, j)  = interp_L(i,j);
                } // end for j
                
                mat_pt.div_vel(gauss_gid) += interp_L(i,i);
            } // end for i
            
            //std::cout << mat_pt.div_vel(gauss_gid) << std::endl;
            
        } // end loop over gauss in element
        
        
    } // end of loop over elements
    
} // end of grad_vel_dg routine


// calculate the velocity gradient directly, no L2 projection
void gradvel_dg_direct(){
    
    for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        
        // initialize veloicty gradient and divergence to zero
        for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
            
            int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
            
            for (int dim_i = 0; dim_i < mesh.num_dim(); dim_i++){
                for (int dim_j = 0; dim_j < mesh.num_dim(); dim_j++){
                    mat_pt.grad_vel(gauss_gid, dim_i, dim_j) = 0.0;
                } // end for j
            } // end for i
            
            mat_pt.div_vel(gauss_gid) = 0.0;
            
        } // end of gauss_lid
    
        
        for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
        
            int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
            
            for (int dim_i = 0; dim_i < mesh.num_dim(); dim_i++){
                for (int dim_j = 0; dim_j < mesh.num_dim(); dim_j++){
                    
                    for(int basis_id = 0; basis_id < elem.num_basis(); basis_id++){
                        
                        int node_basis_lid = elem.vert_node_map(basis_id);
                        int node_gid = mesh.nodes_in_elem(elem_gid, node_basis_lid);
                        
                        mat_pt.grad_vel(gauss_gid, dim_i, dim_j) += ref_elem.ref_nodal_gradient(gauss_lid, basis_id, dim_i)
                                * node.vel(1, node_gid, dim_j);
                        
                    } // end loop over basis
                    
                } // end loop over dim_j
            } // end loop over dim_i

        } // end loop over gauss in element
        
        
        // loop over the guass points in the element
        for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
            
            int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
        
            for(int dim_i = 0; dim_i < mesh.num_dim(); dim_i++) {
                mat_pt.div_vel(gauss_gid) += mat_pt.grad_vel(gauss_gid, dim_i, dim_i);
            }

        } // for gauss_lid

    } // end of loop over elements
    
} // end of grad_vel_dg routine



void art_viscosity(){

    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
         for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
             
             int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
         
             // Calculate the q artifical viscosity q = \mu \l^{hat} div(velocity)
             
             real_t sspd = fmax(mat_pt.sspd(gauss_gid), 1.0e-3);
             //real_t char_length = pow(mesh.elem_vol(elem_gid)/mesh.num_cells_in_elem(), 0.333333);
             real_t char_length = pow(mesh.elem_vol(elem_gid), 0.333333);
             
             // Calculate the symmetric part of the velocity gradient
             auto sym_Velgrad = CArray <real_t> (mesh.num_dim(), mesh.num_dim());
             
             for (int dim_i = 0; dim_i < mesh.num_dim(); dim_i++){
                 for (int dim_j = 0; dim_j < mesh.num_dim(); dim_j++){
                     
                     sym_Velgrad(dim_i, dim_j) =  mat_pt.grad_vel(gauss_gid, dim_i, dim_j); //0.5*(mat_pt.grad_vel(gauss_gid, dim_i, dim_j) + mat_pt.grad_vel(gauss_gid, dim_j, dim_i));
                     
                 } // end loop over dim_j
             } // end loop over dim_i
             
             
             auto Q_visc = CArray <real_t> (mesh.num_dim(), mesh.num_dim());
             
             for (int dim_i = 0; dim_i < mesh.num_dim(); dim_i++) {
                 for (int dim_j = 0; dim_j < mesh.num_dim(); dim_j++) {
                     Q_visc(dim_i, dim_j) = 0.0;
                 }
             }
             
             real_t QC1 = 0.333;
             real_t QC2 = 1.333;
             real_t QC3 = 0.0;
             
             real_t shock_impeadance = mat_pt.density(gauss_gid) *
                                        (QC1*sspd +
                                         QC2*char_length*fabs(mat_pt.div_vel(gauss_gid)) +
                                         QC3*pow(char_length*mat_pt.div_vel(gauss_gid),2.0));
             //shock_impeadance = 5.0;
             
             for (int dim_i = 0; dim_i < mesh.num_dim(); dim_i++) {
                 for (int dim_j = 0; dim_j < mesh.num_dim(); dim_j++) {
                     Q_visc(dim_i, dim_j) =
                         (-shock_impeadance) * char_length * sym_Velgrad(dim_i, dim_j);
                 }
             }
             
             mat_pt.q_visc(gauss_gid) = 0.0;
             for (int dim_i = 0; dim_i < mesh.num_dim(); dim_i++) {
                 mat_pt.q_visc(gauss_gid) -= Q_visc(dim_i, dim_i);
             }
             
             
            // if ( fabs(sym_Velgrad(0, 0)) > 0.0001) std::cout << sym_Velgrad(0, 0) << " ratio = " << mat_pt.q_visc(gauss_gid)/mat_pt.pressure(gauss_gid) << std::endl;
             
             for (int dim_i = 0; dim_i < mesh.num_dim(); dim_i++) {
                 for (int dim_j = 0; dim_j < mesh.num_dim(); dim_j++) {
                     //corner.stress(corner_gid, dim_i, dim_j) -= Q_visc(dim_i, dim_j); //cell_state.stress(1, cell_gid, dim_i, dim_j);
                     mat_pt.stress(gauss_gid, dim_i, dim_j) -= Q_visc(dim_i, dim_j);
                 }
             }
 
         } // end of loop over guass points
    } // end of loop over elements

 
}// end of rountine


