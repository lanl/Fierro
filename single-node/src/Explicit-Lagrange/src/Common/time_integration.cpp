/* time_integration.cpp */

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"


using namespace utils;

// -----------------------------------------------------------------------------
// This function saves the variables at rk_stage = 0, which is t_n
//------------------------------------------------------------------------------
void rk_init(){

    // walk over points, save position and velocity
#pragma omp simd
    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {

        // Create view of vertex velocity
        auto vel     = ViewCArray <real_t> (&node.vel(1, node_gid, 0), mesh.num_dim());
        auto vel_n   = ViewCArray <real_t> (&node.vel(0, node_gid, 0), mesh.num_dim());
        
        // Save velocity at time=n
        for (int dim = 0; dim < mesh.num_dim(); dim++){
            vel_n(dim) = vel(dim);
        }

        // Save position at time=n
        for (int dim = 0; dim < mesh.num_dim(); dim++){
            node.coords(0, node_gid, dim) = node.coords(1, node_gid, dim);
        }

    }
    
    // Save internal energy at time=n
#pragma omp simd
    for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){
        
        cell_state.ie(0, cell_gid) = cell_state.ie(1, cell_gid);
        
        cell_state.total_energy(0, cell_gid) = cell_state.total_energy(1, cell_gid);
        
        auto vel  = ViewCArray <real_t> (&cell_state.velocity(1 , cell_gid, 0), mesh.num_dim());
        auto vel_n  = ViewCArray <real_t> (&cell_state.velocity(0 , cell_gid, 0), mesh.num_dim());
        
        // Save element velocity at time=n
        for (int dim = 0; dim < mesh.num_dim(); dim++){
            vel_n(dim) = vel(dim);
        }
    }
    
} // end rk_init


