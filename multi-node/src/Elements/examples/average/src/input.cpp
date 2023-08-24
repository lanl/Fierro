/* input.c */                                                             

#include <string.h>
#include <sys/stat.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>


#include "utilities.h"
#include "state.h"
#include "header.h"
#include "elements.h"
#include "swage.h"

using namespace utils;

// -----------------------------------------------------------------------------
// This is the input function
//------------------------------------------------------------------------------

void input(){

    // ---- time varaibles and cycle info ---- //
    // time step
    TFINAL = 0.5;
    
    dt_min = 1.e-8;
    dt_max = 1.e-2;
    dt_start = 1.e-5;
    
    cycle_stop = 100000000;

    rk_num_stages = 1;

    rk_storage = 2;

    dt_cfl = 0.5;


    p_order = 3;

    // ---- graphics information ---- //
    graphics_cyc_ival = 10;
    graphics_dt_ival  = 0.1;    // every 0.1 time units


    // ---- fill instructions and intial conditions ---- //

    NF = 2; // number of fills
    
    mat_fill = (mat_fill_t *) malloc((size_t)(NF*sizeof(mat_fill_t)));
    
    // Global instructions
    mat_fill[0].volume = region::global;    // fill everywhere
    mat_fill[0].mat_id = 0;                 // material id
    mat_fill[0].field1 = 1.0;               // some field
    mat_fill[0].field2 = 0.0;               // some other field

    
    // Specific instructions
    mat_fill[1].volume = region::box;   // fill a sphere
    mat_fill[1].mat_id = 1;             // material id
    mat_fill[1].x1 = 0.0;
    mat_fill[1].x2 = 0.7;
    mat_fill[1].y1 = 0.0;
    mat_fill[1].y2 = 2.0;
    mat_fill[1].z1 = 0.0;
    mat_fill[1].z2 = 2.0;
    mat_fill[1].field1 = 10.0;  // some field
    mat_fill[1].field2 = 0.0;   // some other field

    // ---- boundary conditions ---- //
    NB = 6; // number of boundaries
    
    // allocate boundary memory
    boundary = (boundary_t *) malloc((size_t)(NB*sizeof(boundary_t)));
    
    // Tag X=0 plane
    boundary[0].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
    boundary[0].value = 0.0;
    boundary[0].thermal_bc = bdy::isothermal;
    
    // Tag Y=0 plane
    boundary[1].surface = bdy::y_plane;
    boundary[1].value = 0.0;
    boundary[1].thermal_bc = bdy::isothermal;
    
    // Tag Z=0 plane
    boundary[2].surface = bdy::z_plane;
    boundary[2].value = 0.0;
    boundary[2].thermal_bc = bdy::isothermal;
    
    // Tag X=0 plane
    boundary[3].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
    boundary[3].value = 2.0;
    boundary[3].thermal_bc = bdy::isothermal;
    
    // Tag Y=0 plane
    boundary[4].surface = bdy::y_plane;
    boundary[4].value = 2.0;
    boundary[4].thermal_bc = bdy::isothermal;
    
    // Tag Z=0 plane
    boundary[5].surface = bdy::z_plane;
    boundary[5].value = 2.0;
    boundary[5].thermal_bc = bdy::isothermal;


} // end of input


