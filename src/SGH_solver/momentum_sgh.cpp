/* momentum.cpp */

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

// -----------------------------------------------------------------------------
// This function evolves the velocity at the vertices of the mesh
//------------------------------------------------------------------------------
void update_velocity_sgh(real_t rk_alpha){

    num_dim = mesh.num_dim();

    
    // walk over points to update velocity
#pragma omp simd
    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
    	
        // create view into vertex velocities and forces
        auto vel   = ViewCArray <real_t> (&node.vel(1, node_gid, 0), num_dim);
        auto vel_n = ViewCArray <real_t> (&node.vel(0, node_gid, 0), num_dim);

        auto force = ViewCArray <real_t> (&node.force(node_gid, 0), num_dim);
        
        // loop over dimension
        for (int dim = 0; dim < num_dim; dim++){
        	vel(dim) = vel_n(dim) + rk_alpha * dt*force(dim)/node.mass(node_gid);
        }

        // std::cout<<" node "<< node_gid<<" = "<<
        //     node.vel(1, node_gid, 0)<<", "<<
        //     node.vel(1, node_gid, 1)<<", "<<
        //     node.vel(1, node_gid, 2)<< std::endl;

        
    }
    
}
