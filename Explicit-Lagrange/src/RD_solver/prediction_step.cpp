/* Initialize correction step */

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void prediction_step(real_t sub_dt, int prediction_step){

  int num_dim = mesh.num_dim();

  for (int elem_gid = 0 ; elem_gid < mesh.num_elems(); elem_gid++){
    for (int node_lid = 0 ; node_lid < mesh.num_nodes_in_elem(); node_lid++){

      // Get node_gid //
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
      
      // Create view of vel_n //
      auto vel_n = ViewCArray <real_t> (&node.vel(0, node_gid, 0), num_dim);
    
      // Loop over dims //
      for (int dim = 0; dim < num_dim; dim++){
        real_t sum_lo_res = 0.0;

        // Loop over cells to sum residual //
        for (int cell_lid = 0; cell_lid < mesh.num_cells_in_node(node_gid); cell_lid++){
          int cell_gid = mesh.cells_in_node(node_gid, cell_lid);
          sum_lo_res += node.lo_res(node_gid, cell_gid, dim);
        }// end loop over cell_lid
        
	// push initial correction to nodes
        node.vel(prediction_step, node_gid, dim) = vel_n(dim) - sub_dt*sum_lo_res/node.mass(node_gid);

      } // end loop over dims

    }// end loop over node_lid
  }// end loop over elements


}// end init_correction_step
