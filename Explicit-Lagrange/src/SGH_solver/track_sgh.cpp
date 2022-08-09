#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

// -----------------------------------------------------------------------------
// This function tracks the energy for the SGH solver
//------------------------------------------------------------------------------

void track_sgh(real_t &x, real_t &y){

    /* track the total energy in the problem */
    real_t mass = 0.0;
    real_t mass1 = 0.0;
    real_t ke = 0.0;
    real_t ie = 0.0;


    // loop over all elements to find total internal energy 
#pragma omp simd
    for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++) {
        
        ie   += cell_state.mass(cell_gid) * cell_state.ie(1, cell_gid);
        mass += cell_state.mass(cell_gid);
    
    }

    // loop over all points to find total kinetic energy 
    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
        
        // create view into vertex velocity
        auto vel = ViewCArray <real_t> (&node.vel(1, node_gid, 0), num_dim);

        ke += 0.5 * node.mass(node_gid) * 
            (vel(0)*vel(0) + vel(1)*vel(1) + vel(2)*vel(2));

        mass1 += node.mass(node_gid);
    }

    x = ke;
    y = ie;
  
} // end of track
