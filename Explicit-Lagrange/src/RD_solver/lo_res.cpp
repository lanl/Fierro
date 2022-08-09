/* low order residual */
/* R^{r,m}_{p(h)} = \sum_{q}M_{qp}(v^{r,m}_p - v^n_p)
 *                                      + \int^{t^m}_{t^n}(Q^{r,m}_p + \int_{V_h}(\grad\varphi\cdot\sigma)dV)dt  */

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"
#include "jacobi_polynomials.h"

using namespace utils;

void get_lo_res(real_t sub_dt, int t_step, real_t sub_time){
  
  num_dim = mesh.num_dim();

  // Loop over elements //
  for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    // Loop over nodes in element //
    for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
      // Get node global id //
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid); 
      
      // Create a view of vel and vel_n //
      auto vel = ViewCArray <real_t> (&node.vel(t_step, node_gid, 0), num_dim);// Should correspond to current "corrected" velocity at subtime m (i.e. v^{r,m}) //
      auto vel_n = ViewCArray <real_t> (&node.vel(0,node_gid,0), num_dim);// Velocity at previous time (v^n) //
      
      // Create CArray for vel_bar used in artificial viscosity //
      auto vel_bar = CArray <real_t> (num_dim, t_step);

      // Create CArray for Q //
      auto Q = CArray <real_t> (num_dim, t_step);

      // Create CArray for sigma //
      auto sigma = CArray <real_t> (num_dim, num_dim, t_step);

      
      
      // Loop over cells in node //
      for (int cell_lid = 0; cell_lid < mesh.num_cells_in_node(node_gid); cell_lid++){
        // get cell_gid //
        int cell_gid = mesh.cells_in_node(node_gid, cell_lid);

        // Loop over node_lid in cells to compute vel_bar //
        for (int node_lid_in_cell = 0; node_lid<mesh.num_nodes_in_cell(); node_lid_in_cell++){
          // Get node_gid for each node_lid in cell  //
          int node_gid_from_cell = mesh.nodes_in_cell(cell_gid,node_lid);
          
	  //  replaced with line 71
	  // Create view of each vel and vel_n in cell //
          //auto vel_in_cell = ViewCArray <real_t> (&node.vel(sub_time_steps, node_gid_from_cell, 0), num_dim);
          
	  // Initialize sigma //
	  auto sigma = CArray <real_t> (num_dim, num_dim, t_step);

	   
          // Fill sigma //
	  for (int prev_times = 0; prev_times <= t_step; prev_times++){
	    for (int j=0; j < num_dim; j++){
	      for (int i = 0; i < num_dim; i++){
                 sigma(i, j, prev_times) = cell_state.stress(prev_times, cell_gid, i, j);
	      } // end loop over i
	    }// end loop over j
	  }//end loop over previous sub-times
	  

          // Loop over dimension to assign vel_bar values //
	  for (int prev_times = 0 ; prev_times <= t_step; prev_times++){
            for (int dim = 0; dim < num_dim; dim++){
              vel_bar(dim, prev_times) += node.vel(prev_times, node_gid, dim);
            }// end loop over dim for vel_bar and vel_bar_n
	  }// end loop over previous sub times
        }// end loop over node_lid in cell_gid
         

        // Assign values to Q //
        real_t alpha_k = 1.0; // set alpha_k to 1 for now. <--- (CHANGE THIS)
	for (int prev_times = 0; prev_times <= t_step; prev_times++){
          for (int dim = 0; dim < num_dim; dim++){
            Q(dim, prev_times) = alpha_k
		                     *(node.vel(prev_times,node_gid,dim) - vel_bar(dim,prev_times))
				     *jacobi::eval(t_step,-0.5, -0.5, (1-sub_time)/(1+sub_time) );
          }// end loop over dim for Q
        }// end loop over previous sub times
        
         
        // Create CArray to store volume integral over cell of force at each sub time //
        auto force_cell_volume = CArray <real_t> (num_dim, t_step);

        // Begin volume integration of "force" over cell //
        // \int_{V_h}(\grad\varphi\cdot\sigma)dV //
	// and rational chebyshev interpolation //
	for (int prev_times = 0 ; prev_times <= t_step; prev_times++){
          for (int dim_i=0; dim_i < num_dim; dim_i ++){
            for (int dim_j =0; dim_j < num_dim; dim_j++){
              for (int basis_id = 0; basis_id < elem.num_basis(); basis_id++){
              
                force_cell_volume(dim_i, prev_times) += mesh.gauss_cell_pt_jacobian_inverse(cell_gid, dim_i, dim_j)
                                          *ref_elem.ref_cell_gradient(cell_lid, basis_id, dim_j)
                                          *sigma(dim_i, dim_i, prev_times) /// --- everythong around this line needs to be fixed --- ///
                                          *ref_elem.ref_cell_g_weights(cell_lid)
                                          *mesh.gauss_cell_pt_det_j(cell_gid)
					  *jacobi::eval(t_step, -0.5, -0.5, (1-sub_time)/(1+sub_time)); // rational chebyshev approximation

              }// end loop over basis_id
            }// end loop over dim_j        
          }// end loop over dim_i
        }// end loop over previous sub_times
	

        // Create CArray to store time integral of force and artificial viscosity
        auto time_integral = CArray <real_t> (num_dim);
        
	// Begin time integration of force_cell_volume and Q //
        // int^{t^m}_{t^n}(Q^{r,m}_p + \int_{V_h}(\grad\varphi\cdot\sigma)dV)dt /
        
	// compute chebyshev nodes for tie integration
	
	auto cheb_nodes = CArray <real_t> (t_step);
	
	elements::chebyshev_nodes_1D( cheb_nodes, t_step);

        for (int dim = 0; dim < num_dim; dim++){
          for (int prev_times = 0; prev_times <= t_step; prev_times++){

            time_integral(dim) += cheb_nodes(prev_times)*(force_cell_volume(dim, prev_times)+Q(dim,prev_times));
	  
	  }// end loop over previous sub times
        }// end loop over dim for time_int_term
	

        // Store lo_res with node_gid and cell_gid (how to best do this?...)
        
        for (int dim = 0; dim < num_dim; dim++){
          node.lo_res(node_gid, cell_gid, dim) = vel(dim) - vel_n(dim) + time_integral(dim)/cell_state.lumped_mass(cell_gid,node_gid);  // mass_p undefined, node.lo_res(node_gid, cell_gid, dim) undefined //
        }// end loop over dim 

      }// end loop over cells in node 
            
    }// end loop over node_lid
  }// end loop over elements



}// end get_lo_res
