// -----------------------------------------------------------------------------
// This code handles the geometric information for the mesh
//------------------------------------------------------------------------------
#include <iostream>

#include "utilities.h"
#include "state.h"
#include "geometry/geometry.h"
#include "variables.h"

using namespace utils;



void update_position_cch(real_t rk_alpha){
    
    // walk over the points to evolve position
#pragma omp simd
    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
        
        // create view of the vertex velocities
        auto vel   = ViewCArray <real_t> (&node.vel(1, node_gid, 0), num_dim);

        for (int dim = 0; dim < 3; dim++){
            node.coords(1, node_gid, dim) = 
                node.coords(0, node_gid, dim) + rk_alpha * dt*(vel(dim));

            mesh.node_coords(node_gid, dim) = node.coords(1, node_gid, dim);
        }
    } // end for loop over points
} // end subroutine
