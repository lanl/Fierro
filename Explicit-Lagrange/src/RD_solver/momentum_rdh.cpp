/* momentum_rd.cpp*/

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_momentum_rd(int pred_step, int correction_step){

  num_dim = mesh.num_dim();

#pragma omp simd    
  for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++){

    // Create views of vel_r //
    auto vel_r = ViewCArray <real_t> (&node.vel(pred_step, node_gid, 0), num_dim);

    // Create CArray  to store summed residual //
    auto sum_res = CArray <real_t> (num_dim);

    // Initialize to zero //
    for (int dim = 0; dim < num_dim; dim++){

      sum_res(dim) = 0.0;

    }// end sum_res initialization
    
    // Update each dim of vel_{r+1} //
    for (int dim = 0; dim < num_dim; dim++){

      // Sum res in cells around node //
      for (int cell_lid = 0; cell_lid < mesh.num_cells_in_node(node_gid); cell_lid++){

        // Get cell_gid //
        int cell_gid = mesh.cells_in_node(node_gid, cell_gid);

        // Perform summation //
        sum_res(dim) += node.lo_res(node_gid, cell_gid, dim);

      }// end loop over cell_lid

      // Update momentum //
      node.vel(correction_step,node_gid,dim) = vel_r(dim) - sum_res(dim)/node.mass(node_gid);

    }// end loop over dim
  }//end loop over nodes
}// end get_momentum_rd()
