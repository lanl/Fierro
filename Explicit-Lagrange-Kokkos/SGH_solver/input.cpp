                                                           
// -----------------------------------------------------------------------------
// This code is the input file
//------------------------------------------------------------------------------
#include "matar.h"
#include "state.h"
#include "mesh.h"


// -----------------------------------------------------------------------------
// The user must input their parameters inside this function
//------------------------------------------------------------------------------
void input(CArrayKokkos <material_t> &material,
           CArrayKokkos <mat_fill_t> &mat_fill,
           CArrayKokkos <boundary_t> &boundary,
           CArrayKokkos <double> &state_vars,
           size_t &num_materials,
           size_t &num_fills,
           size_t &num_bcs,
           size_t &num_dims,
           size_t &num_state_vars,
           double &dt_start,
           double &time_final,
           double &dt_max,
           double &dt_min,
           double &dt_cfl,
           double &graphics_dt_ival,
           size_t &graphics_cyc_ival,
           size_t &cycle_stop,
           size_t &rk_num_stages){
    
    // Dimensions
    num_dims = 3;
    
    // ---- time varaibles and cycle info ----
    time_final = 1.0;
    dt_min = 1.e-8;
    dt_max = 1.e-2;
    dt_start = 1.e-5;
    cycle_stop = 100000000;


    // ---- graphics information ----
    graphics_cyc_ival = 1000000;
    graphics_dt_ival  = 0.025;

    
    // --- number of material regions ---
    num_materials = 1;
    material = CArrayKokkos <material_t> (num_materials); // create material
    
    
    // --- declare model state variable array size ---
    num_state_vars = 6;  // 6 values for gamma_law gas, it is a memory block
    state_vars = CArrayKokkos <double> (num_materials, num_state_vars); // init values
    
    
    // --- number of fill regions ---
    num_fills = 2;
    mat_fill = CArrayKokkos <mat_fill_t> (num_fills); // create fills
    
    
    // --- number of boundary conditions ---
    num_bcs= 6;
    boundary = CArrayKokkos <boundary_t> (num_bcs);  // create boundaries
    
    // --- test problems ---
    int test_problem = 1;
    
    
    // ---- fill instructions and intial conditions ---- //
    
    // Problem choice:
    // 1 = sedov (set up for 3D 30x30x30 mesh)
    // 2 = sod x  (set up for 3D 200x1x1 mesh)
    // 3 = sod y  (set up for 3D 1x200x1 mesh)
    // 4 = sod z  (set up for 3D 1x1x200 mesh)
    // 5 = double rarefaction (set up for 3D 30x1x1 mesh)
    // 6 = Noh 3D
    // 7 = Triple point
    // 8 = Taylor Green Vortex (add source)
    // 9 = Shockless Noh (fix boundary conditions)
    
    
    // Sedov blast wave test case
    if (test_problem == 1){
        
        RUN({
            // gamma law model
            // state_vars(0) = specific heat
            // state_vars(1) = ref temperature
            // state_vars(2) = ref density
            // state_vars(3) = ref specific internal energy
            // state_vars(4) = gamma
            // state_vars(5) = minimum sound speed
            
            material(0).mat_model = ideal_gas; // EOS model
            material(0).b1        = 1.3333;    // linear slope of UsUp for Riemann solver
            
            material(0).num_state_vars = 6;  // actual num_state_vars
            state_vars(0,0) = 1.0;     // specific heat
            state_vars(0,4) = 5.0/3.0; // gamma value
            state_vars(0,5) = 1.0E-14; // minimum sound speed
            
            
            // global initial conditions
            mat_fill(0).volume = region::global; // fill everywhere
            mat_fill(0).mat_id = 0;              // material id
            mat_fill(0).den = 1.0;               // intial density
            mat_fill(0).sie = 1.e-10;            // intial specific internal energy
            
            mat_fill(0).velocity = init_conds::cartesian;
            mat_fill(0).u = 0.0;   // initial x-dir velocity
            mat_fill(0).v = 0.0;   // initial y-dir velocity
            mat_fill(0).w = 0.0;   // initial z-dir velocity
            
            // energy source initial conditions
            mat_fill(1).volume = region::sphere; // fill a sphere
            mat_fill(1).mat_id = 0;              // material id
            mat_fill(1).radius1 = 0.0;           // inner radius of fill region
            mat_fill(1).radius2 = 1.2/8.0;       // outer radius of fill region
            mat_fill(1).den = 1.0;               // initial density
            mat_fill(1).sie = (963.652344*
                               pow((1.2/30.0),3))/pow((mat_fill(1).radius2),3);
            
            mat_fill(1).velocity = init_conds::cartesian;
            mat_fill(1).u = 0.0;   // initial x-dir velocity
            mat_fill(1).v = 0.0;   // initial y-dir velocity
            mat_fill(1).w = 0.0;   // initial z-dir velocity



            // ---- boundary conditions ---- //
            
            // Tag X plane
            boundary(0).surface = bdy::x_plane; // planes, cylinder, spheres, or a files
            boundary(0).value = 0.0;
            boundary(0).hydro_bc = bdy::reflected;
            
            // Tag Y plane
            boundary(1).surface = bdy::y_plane;
            boundary(1).value = 0.0;
            boundary(1).hydro_bc = bdy::reflected;
            
            // Tag Z plane
            boundary(2).surface = bdy::z_plane;
            boundary(2).value = 0.0;
            boundary(2).hydro_bc = bdy::reflected;
            
            
            // Tag X plane
            boundary(3).surface = bdy::x_plane; // planes, cylinder, spheres, or a files
            boundary(3).value = 1.2;
            boundary(3).hydro_bc = bdy::reflected;
            
            // Tag Y plane
            boundary(4).surface = bdy::y_plane;
            boundary(4).value = 1.2;
            boundary(4).hydro_bc = bdy::reflected;
            
            // Tag Z plane
            boundary(5).surface = bdy::z_plane;
            boundary(5).value = 1.2;
            boundary(5).hydro_bc = bdy::reflected;
            
        });  // end RUN

    } // end if Sedov
    

} // end of input


