#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

// -----------------------------------------------------------------------------
// This function tracks the energy for the CCH solver
//------------------------------------------------------------------------------

void track_cch(real_t &x, real_t &y){

    /* track the total energy in the problem */
    real_t ke = 0.0;
    real_t ie = 0.0;


    // loop over all cells to find total internal energy 
#pragma omp simd
    for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++) {
        
        ie += cell_state.mass(cell_gid) * cell_state.ie(1, cell_gid);
    
        ke += 0.5 * cell_state.mass(cell_gid) * 
              (cell_state.velocity(1, cell_gid, 0)*cell_state.velocity(1, cell_gid, 0) +
               cell_state.velocity(1, cell_gid, 1)*cell_state.velocity(1, cell_gid, 1) +
               cell_state.velocity(1, cell_gid, 2)*cell_state.velocity(1, cell_gid, 2));

    }

    x = ke;
    y = ie;
  
} // end of track