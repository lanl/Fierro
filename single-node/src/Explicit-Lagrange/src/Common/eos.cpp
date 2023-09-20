// -----------------------------------------------------------------------------
// This code contains the equation of state information (constitutive relations)
//------------------------------------------------------------------------------
#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"


using namespace utils;

// -----------------------------------------------------------------------------
// This function evaluates the equation of state (eos) at the cells
//------------------------------------------------------------------------------
void cell_properties(int cell_gid){
    
    // eos_function can return different variables
    eos_return_t eos_return;
    
    // calculate the density  
    cell_state.density(cell_gid) = cell_state.mass(cell_gid)/
                                mesh.cell_vol(cell_gid);


    // calculate the pressure
    cell_state.pressure(cell_gid) = 
        material[cell_state.mat_id(cell_gid)].eos_func(
            p_of_de,
            cell_state.mat_id(cell_gid), 
            cell_state.density(cell_gid), 
            cell_state.ie(1, cell_gid)
            );
    
    // calculate the sound speed
    cell_state.cs(cell_gid) = 
        material[cell_state.mat_id(cell_gid)].eos_func(
            sspd, 
            cell_state.mat_id(cell_gid), 
            cell_state.density(cell_gid), 
            cell_state.ie(1, cell_gid)
            );
        
    

} // end of function


// -----------------------------------------------------------------------------
// This function evaluates the equation of state (eos) at the material points
// (same as gauss points)
//------------------------------------------------------------------------------
void gauss_properties(int mat_pt_gid){
    
    // eos_function can return different variables
    eos_return_t eos_return;
    
    // calculate the density  
    // mat_pt.density(mat_pt_gid) = mat_pt.mass(mat_pt_gid)/
    //                             (mesh.gauss_pt_det_j(mat_pt_gid) * mat_pt.weight(mat_pt_gid));

    // mat_pt.density(mat_pt_gid) = 1.0;
    // Note: Limit density using the same trick. 

    // create a view of the nodal velocity, which is the Riemann velocity
    int node_gid = mesh.node_in_gauss(mat_pt_gid);


    // calculate the pressure
    mat_pt.pressure(mat_pt_gid) = 
        material[mat_pt.mat_id(mat_pt_gid)].eos_func(
            p_of_de,
            mat_pt.mat_id(mat_pt_gid), 
            mat_pt.density(mat_pt_gid), 
            mat_pt.ie(mat_pt_gid)
            );
    
    // calculate the sound speed
    mat_pt.sspd(mat_pt_gid) = 
        material[mat_pt.mat_id(mat_pt_gid)].eos_func(
            sspd, 
            mat_pt.mat_id(mat_pt_gid), 
            mat_pt.density(mat_pt_gid), 
            mat_pt.ie(mat_pt_gid)
            );
        
} // end of function



// -----------------------------------------------------------------------------
// This is the gamma-law gas eos
//------------------------------------------------------------------------------
real_t ideal_gas(int kind, int k, real_t den, real_t ie)
{

    switch (kind) {
            
        // return pressure
        case 1 :
        {
          return (material[k].g - 1.0)*ie*den;
        }
        // return specific internal energy
        case 2 :
        {
          // add method to calculate e
          return ie;
        }
        // return specific heat
        case 3 :
        {
            // add method to calculate cv
            return 1.0;
        }
        // return sound speed
        case 4 :
        {
            real_t sspd = material[k].g*(material[k].g - 1.0)*ie;

            if (sspd < material[k].csmin) sspd = material[k].csmin;

            return sqrt(sspd);
        }
            
        // return temperature
        case 5 :
            return ie/material[k].cv;

        default :
        {
            printf("eos kind is not understood\n");
            return 0.0;
        }

    } // end of switch of return

} // end of ideal_gas

