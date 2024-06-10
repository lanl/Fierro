#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

// -----------------------------------------------------------------------------
// This function tracks the energy for the DGH solver
//------------------------------------------------------------------------------

void track_dgh(real_t &x, real_t &y){

    /* track the total energy in the problem */
    real_t ke = 0.0;
    real_t ie = 0.0;


    // loop over all elements to find total internal energy 
#pragma omp simd

    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        
        for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
            
            int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
            
            real_t mass = mat_pt.mass(gauss_gid); //= mat_pt.density(gauss_gid)*ref_elem.ref_node_g_weights(gauss_lid)*mesh.gauss_pt_det_j(gauss_gid);
            real_t sie = mat_pt.ie(gauss_gid);
            
            real_t vel_dot_vel = 0.0;
            for(int dim = 0; dim < mesh.num_dim(); dim++) {
                vel_dot_vel += mat_pt.velocity(1, gauss_gid, dim)*mat_pt.velocity(1, gauss_gid, dim);
            }
            ke += mass*vel_dot_vel/2.0;
            ie += mass*sie;
            
        }
    } // end elem_gid loop


    x = ke;
    y = ie;
  
} // end of track
