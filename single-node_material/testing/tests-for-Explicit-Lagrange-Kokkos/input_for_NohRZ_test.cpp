                                                           
// -----------------------------------------------------------------------------
// This code is the input file
//------------------------------------------------------------------------------
#include "matar.h"
#include "state.h"
#include "mesh.h"


// Problem choice:
//     sedov 3D and RZ
//     sod3D x  (set up for 3D 200x1x1 mesh)
//     sod3D y  (set up for 3D 1x200x1 mesh)
//     sod3D z  (set up for 3D 1x1x200 mesh)
//     double rarefaction (set up for 3D 30x1x1 mesh)
//     Noh 3D and RZ
//     Triple point
//     Taylor Green Vortex (add source)
//     Shockless Noh (fix boundary conditions)
//     Taylor Anvil 3D and RZ
namespace test
{
    
    // applying initial conditions
    enum setup
    {
        none = 0,
        Sedov3D = 1,
        SedovRZ = 2,
        
        Noh3D = 3,
        NohRZ = 4,
        
        SodZ = 5,
        Sod3DX = 6,
        Sod3DY = 7,
        Sod3DZ = 8,
        
        TriplePoint = 9,
        TaylorAnvil3D = 10,
        TaylorAnvilRZ = 11
    };
    
} // end of initial conditions namespace

test::setup test_problem;

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
           size_t &max_num_state_vars,
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
    num_dims = 2;
    
    // ---- time varaibles and cycle info ----
    time_final = 1.0;  // 1.0 for Sedov
    dt_min = 1.e-8;
    dt_max = 1.e-2;
    dt_start = 1.e-5;
    cycle_stop = 100000;


    // ---- graphics information ----
    graphics_cyc_ival = 1000000;
    graphics_dt_ival  = 0.25;

    
    // --- number of material regions ---
    num_materials = 1;
    material = CArrayKokkos <material_t> (num_materials); // create material
    
    
    // --- declare model state variable array size ---
    max_num_state_vars = 6;  // it is a memory block
    state_vars = CArrayKokkos <double> (num_materials, max_num_state_vars); // init values
    
    
    // --- number of fill regions ---
    num_fills = 1;  
    mat_fill = CArrayKokkos <mat_fill_t> (num_fills); // create fills
    
    
    // --- number of boundary conditions ---
    num_bcs=2;  
    boundary = CArrayKokkos <boundary_t> (num_bcs);  // create boundaries
    
    // --- test problems ---
    test_problem = test::NohRZ;
    
    // ---- fill instructions and intial conditions ---- //

    // use mesh: mesh_NohRZ_64.geo
    // Noh 2D for 64 x 64 mesh
    if (test_problem == test::NohRZ){

        time_final = 0.6;
        
        RUN({
            
            material(0).eos_model = ideal_gas; // EOS model
            material(0).q1        = 1.0;       // accoustic coefficient
            material(0).q2        = 1.3333;    // linear slope of UsUp for Riemann solver
            material(0).q1ex      = 1.0;       // accoustic coefficient in expansion
            material(0).q2ex      = 0.0;       // linear slope of UsUp in expansion
            
            material(0).num_state_vars = 3;  // actual num_state_vars
            state_vars(0,0) = 5.0/3.0; // gamma value
            state_vars(0,1) = 1.0E-14; // minimum sound speed
            state_vars(0,2) = 1.0;     // specific heat c_v
            
            // Global instructions
            mat_fill(0).volume = region::global;   // fill everywhere
            mat_fill(0).mat_id = 0;                // material id
            mat_fill(0).den = 1.0;                   // intial density
            mat_fill(0).sie = 1e-9;             // intial specific internal energy
            
            mat_fill(0).velocity = init_conds::radial;
            mat_fill(0).speed = -1.0;
            
            // ---- boundary conditions ---- //
            
            // Tag X plane
            boundary(0).surface = bdy::x_plane; // planes, cylinder, spheres, or a files
            boundary(0).value = 0.0;
            boundary(0).hydro_bc = bdy::reflected;
            
            // Tag Y plane
            boundary(1).surface = bdy::y_plane;
            boundary(1).value = 0.0;
            boundary(1).hydro_bc = bdy::reflected;
            
        });  // end RUN
            
    } // end if Noh

    return;

} // end of input