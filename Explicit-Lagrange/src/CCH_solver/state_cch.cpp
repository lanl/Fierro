/* state.cpp */

#include "utilities.h"
#include "state.h"
#include "geometry/geometry.h"
#include "variables.h"

using namespace utils;

// -----------------------------------------------------------------------------
// This function that calculates the cell pressure, density, sound speed
//------------------------------------------------------------------------------
void update_state_cch(){

#pragma omp simd
    for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++) {

        cell_properties(cell_gid);
        
    }
}


