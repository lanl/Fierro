#include <math.h>
#include <iostream>
#include <string.h>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"


#define PI 3.14159265


using namespace utils;

void energy_dg(real_t rk_alpha, int cyle){

	// Build the RHS vector for calculating the new vertex velocity

	/*
	std::string t1 = "limit_energy_";
	std::string t2 = std::to_string(cycle);
	std::string combo_title = t1+t2;

	const char *title = combo_title.c_str();

	outfile = fopen(title, "w");

	fprintf(outfile, "elem_gid, gauss_gid, r, mega_alpha, r, e_lim\n");

	//create file that keeps track of the negative energy values 


	FILE *out_negative;

	std::string n1 = "negative_energy_";
	std::string n2 = std::to_string(cycle);
	std::string combo_title2 = n1+n2;

	const char *title2 = combo_title2.c_str();
	
	out_negative = fopen(title2, "w");

	fprintf(out_negative, " elem_gid, gauss_gid, r, energy, avg, r, mega_alpha\n");

	*/
	
	// FILE *out_alpha;
	// out_alpha = fopen("energy_limiter_end", "w");
	// fprintf(out_alpha, " elem_gid, gauss_gid, r, alpha, r, energy, r, avg\n");

    // conservation check
    real_t RHS_energy_check = 0.0;
    
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
	        
			// loop over the cells in the element and corners/nodes in the cell
			// to calculate \Phi a_{i}n_{i} \cdot (\sigma^* u*), which is the same 
			// as the basis at the node times the corner force * riemann velocity at 
			// this node
	        
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
	                		* corner.force(corner_gid, dim)
	                		* node.vel(1, node_gid, dim); // Riemann velocity is stored at the node

	                		// check that force times velocity is positive
	                }

	            }// end loop over nodes/corners in a cell
	        } // end loop over cells in an element

	        // Volume integral
	        for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){


	                /*	
			//--- source term for taylor green ----
	        	int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);

	        	real_t x = mesh.node_coords(node_gid, 0);
	        	real_t y = mesh.node_coords(node_gid, 1);

	        	real_t front = 3.14159265/(4.0*((7.0/5.0) - 1.0));

	        	// // Source term for taylor green problem
	                real_t source = front * (cos(3.0*PI*x)*cos(1.0*PI*y) - cos(3.0*PI*y)*(cos(1.0*PI*x)));
                        //----- source term ends here----
			*/

	        	int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
	        	int node_rid = gauss_lid; // Leveraging internal structured geometry
		        // vel_pt = series expansion (over verts)  ---> mat-pt.vel

	        	for(int i=0; i<mesh.num_dim(); i++){
	        		for(int j=0; j<mesh.num_dim(); j++){
	        			for(int k=0; k<mesh.num_dim(); k++){
        				
	        				volume_int(basis_id) += 
	        								ref_elem.ref_nodal_gradient(node_rid, basis_id, i)
	        								* mesh.gauss_pt_jacobian_inverse(gauss_gid, i, j)
	        								* mat_pt.stress(gauss_gid, j, k)
	        								* mat_pt.velocity(1, gauss_gid, k)
	        								* mesh.gauss_pt_det_j(gauss_gid)
	        								* ref_elem.ref_node_g_weights(gauss_lid);
	        			}
	        		}
	        	}

	        } // end loop over the gauss points
	    } // end loop over the basis

	    // save the RHS vector
	    for(int basis_id=0; basis_id < ref_elem.num_basis(); basis_id++){
	    	RHS(basis_id) = surface_int(basis_id) - volume_int(basis_id);
            
            RHS_energy_check += RHS(basis_id);  // checking conservation
	    }


	    // Calcuate (dtau/dt) = M^{-1} * RHS
	    auto dtau_dt = CArray <real_t> (ref_elem.num_basis());
		
		// Initialize to zero
		for(int basis_id=0; basis_id < ref_elem.num_basis(); basis_id++){
	    	dtau_dt(basis_id) = 0.0;
	    }

	    // Calculate the time rate of change of total energy
		for(int basis_i=0; basis_i < ref_elem.num_basis(); basis_i++){
			for(int basis_j=0; basis_j < ref_elem.num_basis(); basis_j++){

				dtau_dt(basis_i) += 
					elem_state.mass_mat_inv(elem_gid, basis_i, basis_j)
					* RHS(basis_j);
			}
	    }

		// update the total energy at the kinematic degrees of freedom
		for(int basis_id = 0; basis_id < elem.num_basis(); basis_id++){

			int node_basis_id = elem.vert_node_map(basis_id);
			int gauss_gid = mesh.gauss_in_elem(elem_gid, node_basis_id);

            mat_pt.specific_total_energy(1, gauss_gid) = 
                mat_pt.specific_total_energy(0, gauss_gid) + rk_alpha * dtau_dt(basis_id) * (dt);
	    
	    } // end loop over the basis and associated gauss points


	    // Evaluate the specific total energy polynomal at the nodes

	    real_t elem_tot_energy_pre = 0.0;
	    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

            // get the global id of the quadrature point
            int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);

            int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);

            real_t interp_energy = 0.0;

            // Sum over the basis times the velocity defined at the basis vertex position
            for(int basis_id = 0; basis_id < elem.num_basis(); basis_id++){

                int node_basis_id = elem.vert_node_map(basis_id);
                int interp_gid = mesh.gauss_in_elem(elem_gid, node_basis_id);

                interp_energy += mat_pt.specific_total_energy(1, interp_gid) 
                			   * ref_elem.ref_nodal_basis(gauss_lid, basis_id);
            }

            // Save interpolated energy back to gauss point
            mat_pt.specific_total_energy(1, gauss_gid) = interp_energy;


        } // end loop over gauss in element
        
    } // end for elements
    
    
    //if (fabs(RHS_energy_check) >= 1.0e-14 ){
    //std::cout << " #### RHS Energy Error = " << RHS_energy_check << std::endl;
    //}
    
} // end of energy evolution routine



void calc_average_specific_total_energy(){
    
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        
        // initialize to 0
        elem_state.avg_specific_total_energy(elem_gid) = 0.0;
        real_t elem_mass = 0.0;
        
        // tally up all quadrature point values for element
        for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
            
            int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
            
            elem_mass += mat_pt.mass(gauss_gid);
            elem_state.avg_specific_total_energy(elem_gid) += mat_pt.specific_total_energy(1, gauss_gid)*mat_pt.mass(gauss_gid);
            
        } // end loop gauss
        
        // the average element value
        elem_state.avg_specific_total_energy(elem_gid) = elem_state.avg_specific_total_energy(elem_gid)/elem_mass;
        
    } // end elem_gid loop
    
} // end calc avg specific total energy
