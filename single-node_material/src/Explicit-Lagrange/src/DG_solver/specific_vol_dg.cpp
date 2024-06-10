#include <math.h>
#include <iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"


#define PI 3.14159265


using namespace utils;

void specific_vol_dg(real_t rk_alpha){

	// Build the RHS vector for calculating the specific volume evolution

	// Loop over the elements in the mesh
	for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){

		// Left and right parts of the RHS
		auto surface_int = CArray <real_t>(ref_elem.num_basis());
		auto volume_int = CArray <real_t>(ref_elem.num_basis());
		auto RHS = CArray <real_t> (ref_elem.num_basis());

		// Initialize to zero
		for(int basis_id = 0; basis_id < ref_elem.num_basis(); basis_id++){
			surface_int(basis_id) = 0.0;
			volume_int(basis_id) = 0.0;
		}

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
	                for(int dim=0; dim<mesh.num_dim(); dim++){
	                	
	                	surface_int(basis_id) += 
	                		ref_elem.ref_nodal_basis(node_rid, basis_id)
	                		* corner.normal(corner_gid, dim)
	                		* node.vel(1, node_gid, dim); // Riemann velocity is stored at the node
	                }

	            }// end loop over nodes/corners in a cell
	        } // end loop over cells in an element

	        // Volume integral
	        for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

	        	int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
	        	int node_rid = gauss_lid; // Leveraging internal structured geometry


	        	for(int i=0; i<mesh.num_dim(); i++){
	        		for(int j=0; j<mesh.num_dim(); j++){
        				
	        			volume_int(basis_id) += 
	        								ref_elem.ref_nodal_gradient(node_rid, basis_id, i)
	        								* mesh.gauss_pt_jacobian_inverse(gauss_gid, i, j)
	        								* mat_pt.velocity(1, gauss_gid, j)
	        								* mesh.gauss_pt_det_j(gauss_gid)
	        								* ref_elem.ref_node_g_weights(gauss_lid);
	        		}
	        	}

	        } // end loop over the gauss points
	    } // end loop over the basis

	    // save the RHS vector
	    for(int basis_id=0; basis_id < ref_elem.num_basis(); basis_id++){
	    	RHS(basis_id) = surface_int(basis_id) - volume_int(basis_id); 
	    }


	    // Calcuate (d spec_vol/dt) = M^{-1} * RHS
	    auto dgamma_dt = CArray <real_t> (ref_elem.num_basis());
		
		// Initialize to zero
		for(int basis_id=0; basis_id < ref_elem.num_basis(); basis_id++){
	    	dgamma_dt(basis_id) = 0.0;
	    }

	    // Calculate the time rate of change of total energy
		for(int basis_i=0; basis_i < ref_elem.num_basis(); basis_i++){
			for(int basis_j=0; basis_j < ref_elem.num_basis(); basis_j++){

				dgamma_dt(basis_i) += 
					elem_state.mass_mat_inv(elem_gid, basis_i, basis_j)
					* RHS(basis_j);
			}
	    }

		// update the spec vol at the kinematic degrees of freedom
		for(int basis_id = 0; basis_id < elem.num_basis(); basis_id++){

			int node_basis_id = elem.vert_node_map(basis_id);
			int gauss_gid = mesh.gauss_in_elem(elem_gid, node_basis_id);

            mat_pt.specific_volume(1, gauss_gid) = 
                mat_pt.specific_volume(0, gauss_gid) + rk_alpha * dgamma_dt(basis_id) * (dt) ; // /(mesh.gauss_pt_det_j(gauss_gid));
	    
	    } // end loop over the basis and associated gauss points




	    // Evaluate the specific volume polynomal at the nodes
	    real_t avg = 0.0;
	    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

            // get the global id of the quadrature point
            int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);

            int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);

            real_t interp_spec_vol = 0.0;

            // Sum over the basis times the spec vol defined at the basis vertex position
            for(int basis_id = 0; basis_id < elem.num_basis(); basis_id++){

                int node_basis_id = elem.vert_node_map(basis_id);
                int interp_gid = mesh.gauss_in_elem(elem_gid, node_basis_id);

                interp_spec_vol += mat_pt.specific_volume(1, interp_gid) 
                			   * ref_elem.ref_nodal_basis(gauss_lid, basis_id);
            }

            // Save interpolated spec vol back to gauss point
            mat_pt.specific_volume(1, gauss_gid) = interp_spec_vol;
        } // end loop over gauss in element


        // Evaluate the density
        //for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

            // get the global id of the quadrature point
           // int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);


           // mat_pt.density(gauss_gid) = 1.0/mat_pt.specific_volume(1, gauss_gid);

       // }


    } // end loop over elements

} // end of energy evolution routine



void calc_average_specific_vol(){
    
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        
        // initialize to 0
        elem_state.avg_specific_volume(elem_gid) = 0.0;
        real_t elem_mass = 0.0;
        
        elem_state.avg_density(elem_gid) = 0.0;
        real_t elem_vol = 0.0;
        
        // tally up all quadrature point values for element
        for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
            
            int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
            
            real_t vol_gauss = ref_elem.ref_node_g_weights(gauss_lid)*mesh.gauss_pt_det_j(gauss_gid);
            

            elem_vol += vol_gauss;
            elem_mass += mat_pt.mass(gauss_gid);
            
            // element_vol/element_mass
            elem_state.avg_specific_volume(elem_gid) += vol_gauss;
            
            //
            //elem_state.avg_density(elem_gid) += (1.0/mat_pt.specific_volume(1, gauss_gid))*vol_gauss;
            
        } // end loop gauss
        
        // the average element value
        elem_state.avg_specific_volume(elem_gid) = elem_state.avg_specific_volume(elem_gid)/elem_mass;

        // the average element value
        elem_state.avg_density(elem_gid) = 1.0/elem_state.avg_specific_volume(elem_gid);
        
    } // end elem_gid loop
    
} // end calc avg specific total energy
