/* input.c */                                                             

#include <string.h>
#include <sys/stat.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>


#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"


using namespace utils;

#define PI 3.14159265
// -----------------------------------------------------------------------------
// This is the input function
//------------------------------------------------------------------------------
void input(){


    // ---- Method Choice ---- //
    CCH = false;
    SGH = false ;
    DGH = true;


    // p_order =
    //  = 0 uses 2x2x2 (trapezodal rule) quadrature points
    //  = 1 uese 3x3x3 (Simpson's rule) quadrature points
    //  = 2 uses 5x5x5 quadrature points
    //  = N uses (2N+1)(2N+1)(2N+1) quadrature points
    p_order = 1;  // DG will use = 1,2,3,..,N
    
    if(SGH == true) p_order = 0;
    if(CCH == true) p_order = 0;

    // ---- time varaibles and cycle info ---- //
    // time step
    TFINAL = 1.0;  //1.0;
    
    //C1 = 1.0;
    //C2 = 1.333;

    dt_min = 1.e-8;
    dt_max = 1.e-2;
    dt_start = 1.e-5;
    
    cycle_stop = 100000000;

    rk_num_stages = 2;

    rk_storage = 2;

    dt_cfl = 0.4;

    // ---- graphics information ---- //
    graphics_cyc_ival = 1000000;
    graphics_dt_ival  = 0.05;    // every 0.1 time units




    // ---- EOS parameters and material model ---- //
    NR = 1;
    material = (material_t *) malloc((size_t)(NR*sizeof(material_t)));
    
    material[0].eos_func = ideal_gas; // EOS model
    material[0].cv       = 1.0;       // specific heat
    material[0].csmin    = 1.0E-14;       // minimum sound speed
    material[0].g        = 1.4;       // gamma value
    material[0].b1       = 1.3333;    // linear slope of UsUp for Riemann solver


    // ---- fill instructions and intial conditions ---- //
    
    // Problem choice //
    // 1 = sedov (set up for 3D 30x30x30 mesh)
    // 2 = sod x  (set up for 3D 200x1x1 mesh)
    // 3 = sod y  (set up for 3D 1x200x1 mesh)
    // 4 = sod z  (set up for 3D 1x1x200 mesh)
    // 5 = double rarefaction (set up for 3D 30x1x1 mesh)
    // 6 = Noh 3D
    // 7 = Triple point

    // 8 = Taylor Green Vortex (add source)

    // 9 = Shockless Noh (fix boundary conditions)

    
    int test_problem = 1;

    // SEDOV on a 30x30x30 mesh
    if (test_problem == 1){

        // currently set up for sedov blast wave 
        material[0].g        = 5.0/3.0;   // gamma value

        NF = 2; // number of fills
        
        mat_fill = (mat_fill_t *) malloc((size_t)(NF*sizeof(mat_fill_t)));
        
        // Global instructions
        mat_fill[0].volume = region::global;   // fill everywhere
        mat_fill[0].mat_id = 0;                // material id
        mat_fill[0].r = 1.0;                   // intial density
        mat_fill[0].ie = 1.e-10;               // intial specific internal energy
        // mat_fill[0].ie = 10;               // intial specific internal energy
        mat_fill[0].u = 0.0;                   // initial x-dir velocity
        mat_fill[0].v = 0.0;
        mat_fill[0].w = 0.0;
        
        // Specific instructions
        mat_fill[1].volume = region::sphere;// fill a sphere
        mat_fill[1].mat_id = 0;             // material id
        mat_fill[1].radius1 = 0.0;          // inner radius of fill region
        
        // mat_fill[1].radius2 = 1.2/30.0;     // outer radius of fill region
        //mat_fill[1].radius2 = 1.0;     // outer radius of fill region

        // mat_fill[1].radius2 = 1.2/6.0;     // outer radius of fill region

        mat_fill[1].radius2 = 1.2/8.0;     // outer radius of fill region

        //mat_fill[1].radius2 = 1.2/16.0;     // outer radius of fill region

        //mat_fill[1].radius2 = 1.2/32.0;     // outer radius of fill region
        
        
        mat_fill[1].r = 1.0;                // initial density


        // 8x8x8 cells
        mat_fill[1].ie = (963.652344*pow((1.2/30.0),3))/pow((1.2/8.0),3);

        // // 16x16x16 cells
        //mat_fill[1].ie = (963.652344*pow((1.2/30.0),3))/pow((1.2/16.0),3);


        // // 32x32x32cells
        //mat_fill[1].ie = (963.652344*pow((1.2/30.0),3))/pow((1.2/32.0),3);
        
        // 30x30x30 cells
        // mat_fill[1].ie = 963.652344;
        
        // 2x2x2 cells
        // mat_fill[1].ie = (963.652344*pow((1.2/30.0),3))/pow((1.2/2.0),3);

        // 6x6x6 cells
        // mat_fill[1].ie = (963.652344*pow((1.2/30.0),3))/pow((1.2/6.0),3);

        // 4x4x4 cells
        // mat_fill[1].ie = (963.652344*pow((1.2/30.0),3))/pow((1.2/4.0),3);
        
        // EXTENSIVE SOURCE TERM = 0.493390

        // 12x12x12 cells
        
        // WARNING: WEIRD  FIX
        // mat_fill[1].ie = (963.652344*pow((1.2/30.0),3))/pow((1.2/12.0),3);

        // mat_fill[1].ie = 61.67376/4.0;

        // mat_fill[1].ie = 18.27;

        // mat_fill[1].ie = 493.39;

        mat_fill[1].velocity = init_conds::cartesian;
        mat_fill[1].u = 0.0;              // initial x-dir velocity
        mat_fill[1].v = 0.0;
        mat_fill[1].w = 0.0;
        



        // ---- boundary conditions ---- //
        NB = 6; // number of boundaries
        // // allocate boundary memory
        boundary = (boundary_t *) malloc((size_t)(NB*sizeof(boundary_t)));
        
        // Tag X plane
        boundary[0].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
        boundary[0].value = 0.0;
        boundary[0].hydro_bc = bdy::reflected;
        
        // Tag Y plane
        boundary[1].surface = bdy::y_plane;
        boundary[1].value = 0.0;
        boundary[1].hydro_bc = bdy::reflected;
        
        // Tag Z plane
        boundary[2].surface = bdy::z_plane;
        boundary[2].value = 0.0;
        boundary[2].hydro_bc = bdy::reflected;



        // Tag X plane
        boundary[3].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
        boundary[3].value = 1.2;
        boundary[3].hydro_bc = bdy::reflected;
        
        // Tag Y plane
        boundary[4].surface = bdy::y_plane;
        boundary[4].value = 1.2;
        boundary[4].hydro_bc = bdy::reflected;
        
        // Tag Z plane
        boundary[5].surface = bdy::z_plane;
        boundary[5].value = 1.2;
        boundary[5].hydro_bc = bdy::reflected;


    }
    
    // Sod in X on a 200x1x1 mesh
    if (test_problem == 2){



        material[0].g        = 1.4;       // gamma value 

        NF = 2; // number of fills
        
        mat_fill = (mat_fill_t *) malloc((size_t)(NF*sizeof(mat_fill_t)));
        
        // Global instructions
        mat_fill[0].volume = region::global;   // fill everywhere
        mat_fill[0].mat_id = 0;                // material id
        mat_fill[0].r = 1.0;                   // intial density
        mat_fill[0].ie = 2.5;                  // intial specific internal energy
        mat_fill[0].velocity = init_conds::cartesian;
        mat_fill[0].u = 0.0;                   // initial x-dir velocity
        mat_fill[0].v = 0.0;
        mat_fill[0].w = 0.0;
        
        // Specific instructions
        mat_fill[1].volume = region::box; // fill a box
        mat_fill[1].mat_id = 0;             // material id
        
        mat_fill[1].x1 = 0.5;          // 
        mat_fill[1].x2 = 1.0;           // 


        // mat_fill[1].x1 = 0.0;          // 
        // mat_fill[1].x2 = 8.0;           // 
        
        mat_fill[1].y1 = 0.0;          // 
        mat_fill[1].y2 = 0.25;          // 
        
        mat_fill[1].z1 = 0.0;          // 
        mat_fill[1].z2 = 0.25;          // 
        
        mat_fill[1].r = 0.125;            
        
        // 30x1x1 cells
        mat_fill[1].ie = 2.5;
        
        mat_fill[1].velocity = init_conds::cartesian;
        mat_fill[1].u = 0.0;              // initial x-dir velocity
        mat_fill[1].v = 0.0;
        mat_fill[1].w = 0.0;
            

        // ---- boundary conditions ---- //
        NB = 6; // number of boundaries
        // allocate boundary memory
        boundary = (boundary_t *) malloc((size_t)(NB*sizeof(boundary_t)));
        
        // Tag X = 0 plane
        boundary[0].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
        boundary[0].value = 0.0;
        boundary[0].hydro_bc = bdy::reflected;
        
        // Tag Y = 0 plane
        boundary[1].surface = bdy::y_plane;
        boundary[1].value = 0.0;
        boundary[1].hydro_bc = bdy::reflected;
        
        // Tag Z = 0 plane
        boundary[2].surface = bdy::z_plane;
        boundary[2].value = 0.0;
        boundary[2].hydro_bc = bdy::reflected;
        
        // Tag X = 100 plane
        boundary[3].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
        // boundary[3].value = 100.0;
        boundary[3].value = 1.0;
        boundary[3].hydro_bc = bdy::reflected;
        
        // Tag Y = 1 plane
        boundary[4].surface = bdy::y_plane;
        boundary[4].value = 0.1;
        boundary[4].hydro_bc = bdy::reflected;
        
        // Tag Z = 1 plane
        boundary[5].surface = bdy::z_plane;
        boundary[5].value = 0.1;
        boundary[5].hydro_bc = bdy::reflected;
    
    }


    // Sod in y on a 1x200x1 mesh
    if (test_problem == 3){

        material[0].g = 1.4;   // gamma value

        NF = 2; // number of fills
        
        mat_fill = (mat_fill_t *) malloc((size_t)(NF*sizeof(mat_fill_t)));
        
        // Global instructions
        mat_fill[0].volume = region::global;   // fill everywhere
        mat_fill[0].mat_id = 0;                // material id
        mat_fill[0].r = 1.0;                   // intial density
        mat_fill[0].ie = 2.5;               // intial specific internal energy
        mat_fill[0].velocity = init_conds::cartesian;
        mat_fill[0].u = 0.0;                   // initial x-dir velocity
        mat_fill[0].v = 0.0;
        mat_fill[0].w = 0.0;
        
        // Specific instructions
        mat_fill[1].volume = region::box; // fill a box
        mat_fill[1].mat_id = 0;             // material id
        
        mat_fill[1].x1 = 0;          // 
        mat_fill[1].x2 = 1;           // 
        
        mat_fill[1].y1 = 50;          // 
        mat_fill[1].y2 = 100;          // 
        
        mat_fill[1].z1 = 0;          // 
        mat_fill[1].z2 = 1;          // 
        
        mat_fill[1].r = 0.1;            
        
        // 30x1x1 cells
        mat_fill[1].ie = 2.5;
        
        mat_fill[1].velocity = init_conds::cartesian;
        mat_fill[1].u = 0.0;              // initial x-dir velocity
        mat_fill[1].v = 0.0;
        mat_fill[1].w = 0.0;
            

        // ---- boundary conditions ---- //
        NB = 6; // number of boundaries
        // allocate boundary memory
        boundary = (boundary_t *) malloc((size_t)(NB*sizeof(boundary_t)));
        
        // Tag X = 0 plane
        boundary[0].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
        boundary[0].value = 0.0;
        boundary[0].hydro_bc = bdy::reflected;
        
        // Tag Y = 0 plane
        boundary[1].surface = bdy::y_plane;
        boundary[1].value = 0.0;
        boundary[1].hydro_bc = bdy::reflected;
        
        // Tag Z = 0 plane
        boundary[2].surface = bdy::z_plane;
        boundary[2].value = 0.0;
        boundary[2].hydro_bc = bdy::reflected;
        
        // Tag X = 1 plane
        boundary[3].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
        boundary[3].value = 1.0;
        boundary[3].hydro_bc = bdy::reflected;
        
        // Tag Y = 100 plane
        boundary[4].surface = bdy::y_plane;
        boundary[4].value = 100.0;
        boundary[4].hydro_bc = bdy::reflected;
        
        // Tag Z = 1 plane
        boundary[5].surface = bdy::z_plane;
        boundary[5].value = 1.0;
        boundary[5].hydro_bc = bdy::reflected;
    
    }


    // Sod in Z on a 1x1x200 mesh
    if (test_problem == 4){

        // currently set up for sedov blast wave 

        NF = 2; // number of fills
        
        mat_fill = (mat_fill_t *) malloc((size_t)(NF*sizeof(mat_fill_t)));
        
        // Global instructions
        mat_fill[0].volume = region::global;   // fill everywhere
        mat_fill[0].mat_id = 0;                // material id
        mat_fill[0].r = 1.0;                   // intial density
        mat_fill[0].ie = 2.5;               // intial specific internal energy
        mat_fill[0].velocity = init_conds::cartesian;
        mat_fill[0].u = 0.0;                   // initial x-dir velocity
        mat_fill[0].v = 0.0;
        mat_fill[0].w = 0.0;
        
        // Specific instructions
        mat_fill[1].volume = region::box; // fill a box
        mat_fill[1].mat_id = 0;             // material id
        
        mat_fill[1].x1 = 0;          // 
        mat_fill[1].x2 = 1;           // 
        
        mat_fill[1].y1 = 0;          // 
        mat_fill[1].y2 = 1;          // 
        
        mat_fill[1].z1 = 50;          // 
        mat_fill[1].z2 = 100;          // 
        
        mat_fill[1].r = 0.1;            
        
        // 30x1x1 cells
        mat_fill[1].ie = 2.5;
        
        mat_fill[1].velocity = init_conds::cartesian;
        mat_fill[1].u = 0.0;              // initial x-dir velocity
        mat_fill[1].v = 0.0;
        mat_fill[1].w = 0.0;
            

        // ---- boundary conditions ---- //
        NB = 6; // number of boundaries
        // allocate boundary memory
        boundary = (boundary_t *) malloc((size_t)(NB*sizeof(boundary_t)));
        
        // Tag X = 0 plane
        boundary[0].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
        boundary[0].value = 0.0;
        boundary[0].hydro_bc = bdy::reflected;
        
        // Tag Y = 0 plane
        boundary[1].surface = bdy::y_plane;
        boundary[1].value = 0.0;
        boundary[1].hydro_bc = bdy::reflected;
        
        // Tag Z = 0 plane
        boundary[2].surface = bdy::z_plane;
        boundary[2].value = 0.0;
        boundary[2].hydro_bc = bdy::reflected;
        
        // Tag X = 1 plane
        boundary[3].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
        boundary[3].value = 1.0;
        boundary[3].hydro_bc = bdy::reflected;
        
        // Tag Y = 1 plane
        boundary[4].surface = bdy::y_plane;
        boundary[4].value = 1.0;
        boundary[4].hydro_bc = bdy::reflected;
        
        // Tag Z = 100 plane
        boundary[5].surface = bdy::z_plane;
        boundary[5].value = 100.0;
        boundary[5].hydro_bc = bdy::reflected;
    
    }


    // Double rarefaction in X
    if (test_problem == 5){

        material[0].g = 1.4;   // gamma value

        NF = 3; // number of fills
        
        mat_fill = (mat_fill_t *) malloc((size_t)(NF*sizeof(mat_fill_t)));
        
        // Global instructions
        mat_fill[0].volume = region::global;   // fill everywhere
        mat_fill[0].mat_id = 0;                // material id
        mat_fill[0].r = 1.0;                   // intial density
        mat_fill[0].ie = 1.0;               // intial specific internal energy
        
        mat_fill[0].velocity = init_conds::cartesian;
        mat_fill[0].u = 0.0;                   // initial x-dir velocity
        mat_fill[0].v = 0.0;
        mat_fill[0].w = 0.0;
        
        // Specific instructions

        mat_fill[1].volume = region::box; // fill a box
        mat_fill[1].mat_id = 0;             // material id
        
        mat_fill[1].x1 = 0.0;          // 
        mat_fill[1].x2 = 0.5;           // 
        
        mat_fill[1].y1 = 0.0;          // 
        mat_fill[1].y2 = 1.0;          // 
        
        mat_fill[1].z1 = 0.0;          // 
        mat_fill[1].z2 = 1.0;          // 
        
        mat_fill[1].r = 1.0;            
        

        mat_fill[1].ie = 1.0;
        
        mat_fill[1].velocity = init_conds::cartesian;
        mat_fill[1].u = -2.0;              // initial x-dir velocity
        mat_fill[1].v = 0.0;
        mat_fill[1].w = 0.0;



        mat_fill[2].volume = region::box; // fill a box
        mat_fill[2].mat_id = 0;             // material id
        
        mat_fill[2].x1 = 0.5;          // 
        mat_fill[2].x2 = 1.0;           // 
        
        mat_fill[2].y1 = 0.0;          // 
        mat_fill[2].y2 = 1.0;          // 
        
        mat_fill[2].z1 = 0.0;          // 
        mat_fill[2].z2 = 1.0;          // 
        
        mat_fill[2].r = 1.0;            
        

        mat_fill[2].ie = 1.0;
        
        mat_fill[2].velocity = init_conds::cartesian;
        mat_fill[2].u = 2.0;              // initial x-dir velocity
        mat_fill[2].v = 0.0;
        mat_fill[2].w = 0.0;
            

        // ---- boundary conditions ---- //
        NB = 4; // number of boundaries
        // allocate boundary memory
        boundary = (boundary_t *) malloc((size_t)(NB*sizeof(boundary_t)));
        


        // Tag x = 0 plane
        // boundary[0].surface = bdy::x_plane;
        // boundary[0].value = 0.0;
        // boundary[0].hydro_bc = bdy::reflected;

        // Tag Y = 0 plane
        boundary[0].surface = bdy::y_plane;
        boundary[0].value = 0.0;
        boundary[0].hydro_bc = bdy::reflected;
        
        // Tag Z = 0 plane
        boundary[1].surface = bdy::z_plane;
        boundary[1].value = 0.0;
        boundary[1].hydro_bc = bdy::reflected;
        

         // Tag x = 1 plane
        // boundary[3].surface = bdy::x_plane;
        // boundary[3].value = 1.0;
        // boundary[3].hydro_bc = bdy::reflected;

        // Tag Y =  0.1 plane
        boundary[2].surface = bdy::y_plane;
        boundary[2].value = 0.1;
        boundary[2].hydro_bc = bdy::reflected;
        
        // Tag Z =  0.1 plane
        boundary[3].surface = bdy::z_plane;
        boundary[3].value = 0.1;
        boundary[3].hydro_bc = bdy::reflected;
    

    }

    // Noh 3D
    if (test_problem == 6){

        TFINAL = 0.6;
        
        material[0].g        = 5.0/3.0;   // gamma value

        NF = 1; // number of fills
        
        mat_fill = (mat_fill_t *) malloc((size_t)(NF*sizeof(mat_fill_t)));
        
        // Global instructions
        mat_fill[0].volume = region::global;   // fill everywhere
        mat_fill[0].mat_id = 0;                // material id
        mat_fill[0].r = 1.0;                   // intial density
        mat_fill[0].ie = 1e-9;             // intial specific internal energy

        mat_fill[0].velocity = init_conds::spherical;
        mat_fill[0].speed = -1.0;

        // ---- boundary conditions ---- //
        NB = 3; // number of boundaries
        
        // allocate boundary memory
        boundary = (boundary_t *) malloc((size_t)(NB*sizeof(boundary_t)));
        
        // Tag X = 0 plane
        boundary[0].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
        boundary[0].value = 0.0;
        boundary[0].hydro_bc = bdy::reflected;
        
        // Tag Y = 0 plane
        boundary[1].surface = bdy::y_plane;
        boundary[1].value = 0.0;
        boundary[1].hydro_bc = bdy::reflected;
        
        // Tag Z = 0 plane
        boundary[2].surface = bdy::z_plane;
        boundary[2].value = 0.0;
        boundary[2].hydro_bc = bdy::reflected;
        

        

    }

    // triple point problem
    if (test_problem == 7){


        material[0].g        = 5.0/3.0;   // gamma value

        NF = 3; // number of fills
        
        mat_fill = (mat_fill_t *) malloc((size_t)(NF*sizeof(mat_fill_t)));
        
        // Global instructions
        mat_fill[0].volume = region::global;   // fill everywhere
        mat_fill[0].mat_id = 0;                // material id
        mat_fill[0].r = 1.0;                   // intial density
        mat_fill[0].ie = 2.5;             // intial specific internal energy
        
        mat_fill[0].velocity = init_conds::cartesian;
        mat_fill[0].u = 0.0;                   // initial x-dir velocity
        mat_fill[0].v = 0.0;
        mat_fill[0].w = 0.0;
        

        // Specific instructions
        mat_fill[1].volume = region::box; // fill a box
        mat_fill[1].mat_id = 0;             // material id
        
        mat_fill[1].x1 = 1.0;          // 
        mat_fill[1].x2 = 7.0;           // 
        
        mat_fill[1].y1 = 0.0;          // 
        mat_fill[1].y2 = 1.5;          // 
        
        mat_fill[1].z1 = 0.0;          // 
        mat_fill[1].z2 = 1.0;          // 
        
        mat_fill[1].r = 1.0;            
        

        mat_fill[1].ie = 0.25;
        
        mat_fill[1].velocity = init_conds::cartesian;
        mat_fill[1].u = 0.00;              // initial x-dir velocity
        mat_fill[1].v = 0.0;
        mat_fill[1].w = 0.0;


        mat_fill[2].volume = region::box; // fill a box
        mat_fill[2].mat_id = 0;             // material id
        
        mat_fill[2].x1 = 1.0;          // 
        mat_fill[2].x2 = 7.0;           // 
        
        mat_fill[2].y1 = 1.5;          // 
        mat_fill[2].y2 = 3.0;          // 
        
        mat_fill[2].z1 = 0.0;          // 
        mat_fill[2].z2 = 1.0;          // 
        
        mat_fill[2].r = 0.1;            
        

        mat_fill[2].ie = 2.5;
        
        mat_fill[2].velocity = init_conds::cartesian;
        mat_fill[2].u = 0.0;              // initial x-dir velocity
        mat_fill[2].v = 0.0;
        mat_fill[2].w = 0.0;





        // ---- boundary conditions ---- //
        NB = 6; // number of boundaries
        
        // allocate boundary memory
        boundary = (boundary_t *) malloc((size_t)(NB*sizeof(boundary_t)));
        
        // Tag X = 0 plane
        boundary[0].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
        boundary[0].value = 0.0;
        boundary[0].hydro_bc = bdy::reflected;
        
        // Tag Y = 0 plane
        boundary[1].surface = bdy::y_plane;
        boundary[1].value = 0.0;
        boundary[1].hydro_bc = bdy::reflected;
        
        // Tag Z = 0 plane
        boundary[2].surface = bdy::z_plane;
        boundary[2].value = 0.0;
        boundary[2].hydro_bc = bdy::reflected;


        // Tag X = 7 plane
        boundary[3].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
        boundary[3].value = 6.0;  // some meshes are 7 and others are 6
        boundary[3].hydro_bc = bdy::reflected;
        
        // Tag Y = 3 plane
        boundary[4].surface = bdy::y_plane;
        boundary[4].value = 3.0;
        boundary[4].hydro_bc = bdy::reflected;
        
        // Tag Z = 1 plane
        boundary[5].surface = bdy::z_plane;
        boundary[5].value = 1.0;
        boundary[5].hydro_bc = bdy::reflected;
        

    }

    
    // TG vortex
    if (test_problem == 8){


        material[0].g        = 7.0/5.0;   // gamma value

        NF = 1; // number of fills
        
        mat_fill = (mat_fill_t *) malloc((size_t)(NF*sizeof(mat_fill_t)));
        
        // Global instructions
        mat_fill[0].volume = region::global;   // fill everywhere
        mat_fill[0].mat_id = 0;                // material id
        mat_fill[0].r = 1.0;                   // intial density
        mat_fill[0].ie = 1.0;             // intial specific internal energy
        
        mat_fill[0].velocity = init_conds::tg_vortex;  // note: pressure and energy are setup with this vel ic


        // ---- boundary conditions ---- //
        NB = 6; // number of boundaries
        
        // allocate boundary memory
        boundary = (boundary_t *) malloc((size_t)(NB*sizeof(boundary_t)));
        
        // Tag X = 0 plane
        boundary[0].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
        boundary[0].value = 0.0;
        boundary[0].hydro_bc = bdy::reflected;
        
        // Tag Y = 0 plane
        boundary[1].surface = bdy::y_plane;
        boundary[1].value = 0.0;
        boundary[1].hydro_bc = bdy::reflected;
        
        // Tag Z = 0 plane
        boundary[2].surface = bdy::z_plane;
        boundary[2].value = 0.0;
        boundary[2].hydro_bc = bdy::reflected;


        // Tag X = 7 plane
        boundary[3].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
        boundary[3].value = 1.0;
        boundary[3].hydro_bc = bdy::reflected;
        
        // Tag Y = 3 plane
        boundary[4].surface = bdy::y_plane;
        boundary[4].value = 1.0;
        boundary[4].hydro_bc = bdy::reflected;
        
        // Tag Z = 1 plane
        boundary[5].surface = bdy::z_plane;
        boundary[5].value = 0.1;
        boundary[5].hydro_bc = bdy::reflected;
        

    }


    // Shockless Noh
    if (test_problem == 9){

        TFINAL = 0.6;
        
        material[0].g        = 7.0/5.0;   // gamma value

        NF = 1; // number of fills
        
        mat_fill = (mat_fill_t *) malloc((size_t)(NF*sizeof(mat_fill_t)));
        
        // Global instructions
        mat_fill[0].volume = region::global;   // fill everywhere
        mat_fill[0].mat_id = 0;                // material id
        mat_fill[0].r = 1.0;                   // intial density
        mat_fill[0].ie = 1.0;             // intial specific internal energy
        
        mat_fill[0].velocity = init_conds::spherical_linear;
        



        // ---- boundary conditions ---- //
        NB = 6; // number of boundaries
        
        // allocate boundary memory
        boundary = (boundary_t *) malloc((size_t)(NB*sizeof(boundary_t)));
        
        // Tag X = 0 plane
        boundary[0].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
        boundary[0].value = 0.0;
        boundary[0].hydro_bc = bdy::reflected;
        
        // Tag Y = 0 plane
        boundary[1].surface = bdy::y_plane;
        boundary[1].value = 0.0;
        boundary[1].hydro_bc = bdy::reflected;
        
        // Tag Z = 0 plane
        boundary[2].surface = bdy::z_plane;
        boundary[2].value = 0.0;
        boundary[2].hydro_bc = bdy::reflected;


        // Add Noh BC's here
        // Tag X = 7 plane
        boundary[3].surface = bdy::x_plane; // planes, cylinder, spheres, or a files
        boundary[3].value = 1.0;
        boundary[3].hydro_bc = bdy::velocity;
        
        // Tag Y = 3 plane
        boundary[4].surface = bdy::y_plane;
        boundary[4].value = 1.0;
        boundary[4].hydro_bc = bdy::velocity;
        
        // Tag Z = 1 plane
        boundary[5].surface = bdy::z_plane;
        boundary[5].value = 1.0;
        boundary[5].hydro_bc = bdy::velocity;
        

    }



} // end of input


