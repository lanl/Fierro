/**********************************************************************************************
 © 2020. Triad National Security, LLC. All rights reserved.
 This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
 National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
 Department of Energy/National Nuclear Security Administration. All rights in the program are
 reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
 Security Administration. The Government is granted for itself and others acting on its behalf a
 nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
 derivative works, distribute copies to the public, perform publicly and display publicly, and
 to permit others to do so.
 This program is open source under the BSD-3 License.
 Redistribution and use in source and binary forms, with or without modification, are permitted
 provided that the following conditions are met:
 
 1.  Redistributions of source code must retain the above copyright notice, this list of
 conditions and the following disclaimer.
 
 2.  Redistributions in binary form must reproduce the above copyright notice, this list of
 conditions and the following disclaimer in the documentation and/or other materials
 provided with the distribution.
 
 3.  Neither the name of the copyright holder nor the names of its contributors may be used
 to endorse or promote products derived from this software without specific prior
 written permission.
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **********************************************************************************************/

#include "utilities.h"
#include "matar.h"
#include "Simulation_Parameters_SGH.h"

using namespace utils;

Simulation_Parameters_SGH::Simulation_Parameters_SGH() : Simulation_Parameters(){
  
  //equate parent class multi maps for parser code
  possible_options = sgh_possible_options;
  possible_options_nested2 = sgh_possible_options_nested2;
  possible_options_nested3 = sgh_possible_options_nested3;

  //initialize data and flags to defaults
  output_strain_flag = false;
  output_stress_flag = false;
  displaced_mesh_flag = false;
  report_runtime_flag = false;
  unit_scaling = 1;
  strain_max_flag = false;
  gravity_flag = false;
  // ---- boundary conditions ---- //
  NB = 0; 
  NBSF = 0; 
  NBV = 0;

  // --- Graphics output variables ---
  graphics_id = 0;
  graphics_cyc_ival = 50;

  graphics_times = CArray<double>(2000);
  graphics_dt_ival = 1.0e8;
  graphics_time = graphics_dt_ival;  // the times for writing graphics dump


  // --- Time and cycling variables ---
  time_value = 0.0;
  time_final = 1.e16;
  dt = 1.e-8;
  dt_max = 1.0e-2;
  dt_min = 1.0e-8;
  dt_cfl = 0.4;
  dt_start = 1.0e-8;

  rk_num_stages = 2;
  rk_num_bins = 2;

  cycle = 0;
  cycle_stop = 1000000000;


  // --- Precision variables ---
  fuzz = 1.0e-16;  // machine precision
  tiny = 1.0e-12;  // very very small (between real_t and single)
  small= 1.0e-8;   // single precision
}

Simulation_Parameters_SGH::~Simulation_Parameters_SGH(){
}

void Simulation_Parameters_SGH::input(){
  
  Simulation_Parameters::input();
  //output settings
  output_velocity_flag = true;
  //requires displacement flag to be true
  displaced_mesh_flag = true;
  
  output_strain_flag = true;
  output_stress_flag = false;

  //simulation spatial dimension
  num_dim = 3;
  unit_scaling = 1;

  //polynomial interpolation order
  p_order = 0;
  
  //Gauss-Legendre integration order
  num_gauss_points = 2;

  //debug and performance report flags
  report_runtime_flag = true;

  // ---- boundary conditions ---- //
  NB = 6; // number of boundaries
  NBSF = 4; //number of surface density force conditions
  NBV = 2; //number of surface sets used to specify a fixed displacement on nodes belonging to respective surfaces

  //apply body forces
  gravity_flag = false;
  gravity_vector[0] = 9.81;
  gravity_vector[1] = 0;
  gravity_vector[2] = 0;
    
  // ---- time varaibles and cycle info ----
  time_final = 1.0;  // 1.0 for Sedov
  dt_min = 1.e-8;
  dt_max = 1.e-2;
  dt_start = 1.e-5;
  cycle_stop = 2000000;

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
  num_fills = 2;  // =2 for Sedov
  mat_fill = CArrayKokkos <mat_fill_t> (num_fills); // create fills
    
  // --- number of boundary conditions ---
  num_bcs=6;  // =6 for Sedov
  boundary = CArrayKokkos <boundary_t> (num_bcs);  // create boundaries
    
  // --- test problems ---
  test_problem = Sedov3D;
    
  // ---- fill instructions and intial conditions ---- //
    
  time_value = 0.0;
  dt = dt_start;
  graphics_id = 0;
  graphics_times(0) = 0.0;
  graphics_time = graphics_dt_ival;  // the times for writing graphics dump
    
    // Sedov blast wave test case
    if (test_problem == Sedov3D){
        time_final = 1.0;  // 1.0 for Sedov
        
        RUN_CLASS({
            // gamma law model
            // statev(0) = gamma
            // statev(1) = minimum sound speed
            // statev(2) = specific heat
            // statev(3) = ref temperature
            // statev(4) = ref density
            // statev(5) = ref specific internal energy
            
            material(0).eos_model = ideal_gas; // EOS model is required
            
            material(0).strength_type = model::none;
            material(0).strength_model = NULL;  // not needed, but illistrates the syntax
            
            material(0).q1        = 1.0;       // accoustic coefficient
            material(0).q2        = 1.3333;    // linear slope of UsUp for Riemann solver
            material(0).q1ex      = 1.0;       // accoustic coefficient in expansion
            material(0).q2ex      = 0.0;       // linear slope of UsUp in expansion
            
            material(0).num_state_vars = 3;  // actual num_state_vars
            material(0).read_state_vars = 0; // no, state_vars declared here
            state_vars(0,0) = 5.0/3.0; // gamma value
            state_vars(0,1) = 1.0E-14; // minimum sound speed
            state_vars(0,2) = 1.0;     // specific heat
            
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
            mat_fill(1).radius2 = 1.2/32.0;       // outer radius of fill region
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
            
        });  // end RUN_CLASS

    } // end if Sedov
    
    // 2D RZ Sedov blast wave test case
    if (test_problem == SedovRZ){
        time_final = 1.0;  // 1.0 for Sedov
        RUN_CLASS({
            // gamma law model
            // statev(0) = gamma
            // statev(1) = minimum sound speed
            // statev(2) = specific heat
            // statev(3) = ref temperature
            // statev(4) = ref density
            // statev(5) = ref specific internal energy
            
            material(0).eos_model = ideal_gas; // EOS model is required
            
            material(0).strength_type = model::none;
            material(0).strength_model = NULL;  // not needed, but illistrates the syntax
            
            material(0).q1        = 1.0;       // accoustic coefficient
            material(0).q2        = 1.3333;    // linear slope of UsUp for Riemann solver
            material(0).q1ex      = 1.0;       // accoustic coefficient in expansion
            material(0).q2ex      = 0.0;       // linear slope of UsUp in expansion
            
            material(0).num_state_vars = 3;  // actual num_state_vars
            material(0).read_state_vars = 0; // no, state_vars declared here
            state_vars(0,0) = 5.0/3.0; // gamma value
            state_vars(0,1) = 1.0E-14; // minimum sound speed
            state_vars(0,2) = 1.0;     // specific heat
            
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
            mat_fill(1).radius1 = 0.01;           // inner radius of fill region
            mat_fill(1).radius2 = (1.2-mat_fill(1).radius1)/50 + mat_fill(1).radius1;       // outer radius of fill region
            mat_fill(1).den = 1.0;               // initial density
            double vol = PI*( pow((mat_fill(1).radius2),3)
                            - pow((mat_fill(1).radius1),3) );
            //vol = 4./3.* PI * ( pow((mat_fill(1).radius2),3) - pow((mat_fill(1).radius1),3) )/2.0;
            mat_fill(1).sie = (0.5*0.49339/vol);
            
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
            
            
            // Tag X plane
            //boundary(2).surface = bdy::x_plane; // planes, cylinder, spheres, or a files
            //boundary(2).value = 1.2;
            //boundary(2).hydro_bc = bdy::reflected;
            //
            //// Tag Y plane
            //boundary(3).surface = bdy::y_plane;
            //boundary(3).value = 1.2;
            //boundary(3).hydro_bc = bdy::reflected;
            
            
            // Tag inner cylinder
            boundary(2).surface = bdy::cylinder;
            boundary(2).value = 0.01;
            boundary(2).hydro_bc = bdy::fixed;
            
        });  // end RUN_CLASS

    } // end if Sedov
    
    
    // Noh 3D
    if (test_problem == Noh3D){

        time_final = 0.6;
        
        RUN_CLASS({
            
            material(0).eos_model = ideal_gas; // EOS model
            material(0).q1        = 1.0;       // accoustic coefficient
            material(0).q2        = 1.3333;    // linear slope of UsUp for Riemann solver
            material(0).q1ex      = 1.0;       // accoustic coefficient in expansion
            material(0).q2ex      = 0.0;       // linear slope of UsUp in expansion
            
            material(0).num_state_vars = 6;  // actual num_state_vars
            state_vars(0,0) = 5.0/3.0; // gamma value
            state_vars(0,1) = 1.0E-14; // minimum sound speed
            state_vars(0,2) = 1.0;     // specific heat c_v
            
            // Global instructions
            mat_fill(0).volume = region::global;   // fill everywhere
            mat_fill(0).mat_id = 0;                // material id
            mat_fill(0).den = 1.0;                   // intial density
            mat_fill(0).sie = 1e-9;             // intial specific internal energy
            
            mat_fill(0).velocity = init_conds::spherical;
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
            
            // Tag Z plane
            boundary(2).surface = bdy::z_plane;
            boundary(2).value = 0.0;
            boundary(2).hydro_bc = bdy::reflected;
            
        });  // end RUN_CLASS
            
    } // end if Noh
    
    
    // Noh 2D
    if (test_problem == NohRZ){

        time_final = 0.6;
        
        RUN_CLASS({
            
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
            
            // Tag Y plane
            boundary(2).surface = bdy::cylinder;
            boundary(2).value = 0.01;
            boundary(2).hydro_bc = bdy::fixed;
            
            
        });  // end RUN_CLASS
            
    } // end if Noh
    
    
    // Sod in Z direction in RZ coordinates (eq. to x-dir)
    if (test_problem == SodZ){
        
        time_final = 0.2;  // 1.0 for Sedov
        
        RUN_CLASS({
            // gamma law model
            // statev(0) = gamma
            // statev(1) = minimum sound speed
            // statev(2) = specific heat
            // statev(3) = ref temperature
            // statev(4) = ref density
            // statev(5) = ref specific internal energy
            
            material(0).eos_model = ideal_gas; // EOS model is required
            
            material(0).strength_type = model::none;
            material(0).strength_model = NULL;  // not needed, but illistrates the syntax
            
            material(0).q1        = 1.0;       // accoustic coefficient
            material(0).q2        = 1.3333;    // linear slope of UsUp for Riemann solver
            material(0).q1ex      = 1.0;       // accoustic coefficient in expansion
            material(0).q2ex      = 0.0;       // linear slope of UsUp in expansion
            
            material(0).num_state_vars = 3;  // actual num_state_vars
            material(0).read_state_vars = 0; // no, state_vars declared here
            state_vars(0,0) = 1.4; // gamma value
            state_vars(0,1) = 1.0E-14; // minimum sound speed
            state_vars(0,2) = 1.0;     // specific heat
            
            // global initial conditions
            mat_fill(0).volume = region::global; // fill everywhere
            mat_fill(0).mat_id = 0;              // material id
            mat_fill(0).den = 1.0;               // intial density
            mat_fill(0).sie = 2.5;            // intial specific internal energy
            
            mat_fill(0).velocity = init_conds::cartesian;
            mat_fill(0).u = 0.0;   // initial x-dir velocity
            mat_fill(0).v = 0.0;   // initial y-dir velocity
            mat_fill(0).w = 0.0;   // initial z-dir velocity
            
            // energy source initial conditions
            mat_fill(1).volume = region::box;    // fill a box
            mat_fill(1).mat_id = 0;              // material id
            mat_fill(1).x1 = 0.5;           //
            mat_fill(1).x2 = 1.0;           //
            mat_fill(1).y1 = 0.0;           //
            mat_fill(1).y2 = 1.0;           //
            mat_fill(1).z1 = 0.0;           //
            mat_fill(1).z2 = 1.0;           //
            mat_fill(1).den = 0.125;        // initial density
            mat_fill(1).sie = 2.5;          // initial specific internal energy
            
            
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
            
            
            // Tag X plane
            boundary(2).surface = bdy::x_plane; // planes, cylinder, spheres, or a files
            boundary(2).value = 1.0;
            boundary(2).hydro_bc = bdy::reflected;

            // Tag Y plane
            boundary(3).surface = bdy::y_plane;
            boundary(3).value = 0.1;
            boundary(3).hydro_bc = bdy::reflected;
            
        });  // end RUN_CLASS

    } // end if SodZ
    
    
    // Triple point
    if (test_problem == TriplePoint){
        
        time_final = 4.0; 
        
        RUN_CLASS({
            // gamma law model
            // statev(0) = gamma
            // statev(1) = minimum sound speed
            // statev(2) = specific heat
            // statev(3) = ref temperature
            // statev(4) = ref density
            // statev(5) = ref specific internal energy
            
            material(0).eos_model = ideal_gas; // EOS model is required
            
            material(0).strength_type = model::none;
            material(0).strength_model = NULL;  // not needed, but illistrates the syntax
            
            material(0).q1        = 1.0;       // accoustic coefficient
            material(0).q2        = 1.3333;    // linear slope of UsUp for Riemann solver
            material(0).q1ex      = 1.0;       // accoustic coefficient in expansion
            material(0).q2ex      = 0.0;       // linear slope of UsUp in expansion
            
            material(0).num_state_vars = 3;  // actual num_state_vars
            material(0).read_state_vars = 0; // no, state_vars declared here
            state_vars(0,0) = 5.0/3.0; // gamma value
            state_vars(0,1) = 1.0E-14; // minimum sound speed
            state_vars(0,2) = 1.0;     // specific heat
            
            // global initial conditions
            mat_fill(0).volume = region::global; // fill everywhere
            mat_fill(0).mat_id = 0;              // material id
            mat_fill(0).den = 1.0;               // intial density
            mat_fill(0).sie = 2.5;            // intial specific internal energy
            
            mat_fill(0).velocity = init_conds::cartesian;
            mat_fill(0).u = 0.0;   // initial x-dir velocity
            mat_fill(0).v = 0.0;   // initial y-dir velocity
            mat_fill(0).w = 0.0;   // initial z-dir velocity
            
            
            // initial conditions, region 1
            mat_fill(1).volume = region::box;    // fill a box
            mat_fill(1).mat_id = 0;              // material id
            mat_fill(1).x1 = 1.0;           //
            mat_fill(1).x2 = 7.0;           //
            mat_fill(1).y1 = 0.0;           //
            mat_fill(1).y2 = 1.5;           //
            mat_fill(1).z1 = 0.0;           //
            mat_fill(1).z2 = 1.0;           //
            mat_fill(1).den = 1.0;          // initial density
            mat_fill(1).sie = 0.25;         // initial specific internal energy
            
            mat_fill(1).velocity = init_conds::cartesian;
            mat_fill(1).u = 0.0;   // initial x-dir velocity
            mat_fill(1).v = 0.0;   // initial y-dir velocity
            mat_fill(1).w = 0.0;   // initial z-dir velocity
            
            // initial conditions, region 2
            mat_fill(2).volume = region::box;    // fill a box
            mat_fill(2).mat_id = 0;              // material id
            mat_fill(2).x1 = 1.0;           //
            mat_fill(2).x2 = 7.0;           //
            mat_fill(2).y1 = 1.5;           //
            mat_fill(2).y2 = 3.0;           //
            mat_fill(2).z1 = 0.0;           //
            mat_fill(2).z2 = 1.0;           //
            mat_fill(2).den = 0.1;        // initial density
            mat_fill(2).sie = 2.5;          // initial specific internal energy
            
            mat_fill(2).velocity = init_conds::cartesian;
            mat_fill(2).u = 0.0;   // initial x-dir velocity
            mat_fill(2).v = 0.0;   // initial y-dir velocity
            mat_fill(2).w = 0.0;   // initial z-dir velocity



            // ---- boundary conditions ---- //
            
            // Tag X = 0 plane
            boundary(0).surface = bdy::x_plane; // planes, cylinder, spheres, or a files
            boundary(0).value = 0.0;
            boundary(0).hydro_bc = bdy::reflected;
            
            // Tag Y = 0 plane
            boundary(1).surface = bdy::y_plane;
            boundary(1).value = 0.0;
            boundary(1).hydro_bc = bdy::reflected;
            
            // Tag Z = 0 plane
            boundary(2).surface = bdy::z_plane;
            boundary(2).value = 0.0;
            boundary(2).hydro_bc = bdy::reflected;


            // Tag X = 7 plane
            boundary(3).surface = bdy::x_plane; // planes, cylinder, spheres, or a files
            boundary(3).value = 6.0;  // some meshes are 7 and others are 6
            boundary(3).hydro_bc = bdy::reflected;
            
            // Tag Y = 3 plane
            boundary(4).surface = bdy::y_plane;
            boundary(4).value = 3.0;
            boundary(4).hydro_bc = bdy::reflected;
            
            // Tag Z = 1 plane
            boundary(5).surface = bdy::z_plane;
            boundary(5).value = 1.0;
            boundary(5).hydro_bc = bdy::reflected;
            
        });  // end RUN_CLASS

    } // end if SodZ
    
    
    // Taylor Anvil
    if (test_problem == TaylorAnvil){

        time_final = 25.0;
        
        RUN_CLASS({
            
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
            
            mat_fill(0).velocity = init_conds::cartesian;
            mat_fill(0).u = 0.0;
            mat_fill(0).v = 0.0;
            mat_fill(0).w = -1.0;
            
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
            
            
        });  // end RUN_CLASS
        
    } // end if Taylor Anvil

}

std::string Simulation_Parameters_SGH::yaml_input(std::string filename){
  Yaml::Node root;
  std::string current_option_outer, current_setting_outer;
  std::string current_option_inner, current_setting_inner;
  std::string current_option_center, current_setting_center;
  std::string error = "success";
  bool inner_found, inner_found_nest;
    try
    {
        Yaml::Parse(root, filename.c_str());
    }
    catch (const Yaml::Exception e)
    {
        std::cout << "Exception " << e.Type() << ": " << e.what() << std::endl;
    }
    
    
    std::cout << "print root Size = " << root.Size() << "\n";
    std::cout << "print root = \n";
    //for (size_t i=0; i<root.Size(); i++){
        
        // get the outer stuff
        Yaml::Node & outer_item = root;
        
        std::cout << "\n";
        std::cout << "size = " << outer_item.Size() << std::endl;
        if (outer_item.Size()!=0){
            for(auto outer_it = outer_item.Begin(); outer_it != outer_item.End(); outer_it++)
            {   
                current_option_outer = (*outer_it).first;
                current_setting_outer = (*outer_it).second.As<std::string>();
                
                //find the keyword for this option out of the three multimaps of possible options with different nesting structure
                multimap_iterator multi_outer_iterator = sgh_possible_options.find(current_option_outer);
                multimap_iterator_nested2 multi_outer_iterator_nested2 = sgh_possible_options_nested2.find(current_option_outer);
                multimap_iterator_nested3 multi_outer_iterator_nested3 = sgh_possible_options_nested3.find(current_option_outer);

                //check if this keyword is an option for this solver type
                if(multi_outer_iterator==sgh_possible_options.end()){
                    if(multi_outer_iterator_nested2==sgh_possible_options_nested2.end()){
                        if(multi_outer_iterator_nested3==sgh_possible_options_nested3.end()){
                            std::string message = "Unsupported option requested in YAML input file: ";
                            error = message + current_option_outer;
                        }
                    }
                  //return error;
                }
                
                std::cout << current_option_outer << " " << current_setting_outer << std::endl;
            
                Yaml::Node & inner_item = (*outer_it).second;
            
                // inner layer
                if (inner_item.Size()!=0){
                    for(auto inner_it = inner_item.Begin(); inner_it != inner_item.End();   inner_it++)
                    {
                        current_option_inner = (*inner_it).first;
                        current_setting_inner = (*inner_it).second.As<std::string>();
                        
                        //check if this option is supported
                        std::pair<multimap_iterator_nested2,multimap_iterator_nested2> iterator_range = sgh_possible_options_nested2.equal_range(current_option_outer);
                        std::pair<multimap_iterator_nested3,multimap_iterator_nested3> iterator_range_nested = sgh_possible_options_nested3.equal_range(current_option_outer);
                        //multimap_iterator multi_outer_iterator_nested2 = sgh_possible_options_nested2[current_option_outer].find(current_option_outer);
                        //multimap_iterator_nested2 multi_outer_iterator_nested3 = sgh_possible_options_nested3[current_option_outer].find(current_option_outer);
                        
                        inner_found = inner_found_nest = false;
                        for(auto temp_it = iterator_range.first; temp_it != iterator_range.second; temp_it++){
                            //test if map element corresponds to a possible option
                            multimap_iterator multi_inner_iterator = temp_it->second.find(current_option_inner);
                            //std::cout << "Test print " << temp_it->first << std::endl;
                            if(multi_inner_iterator!=temp_it->second.end()){
                                inner_found = true;
                            }
                        }

                        for(auto temp_it = iterator_range_nested.first; temp_it != iterator_range_nested.second; temp_it++){
                            //test if map element corresponds to a possible option
                            multimap_iterator_nested2 multi_inner_iterator = temp_it->second.find(current_option_inner);
                            //std::cout << "Test print " << temp_it->first << std::endl;
                            if(multi_inner_iterator!=temp_it->second.end()){
                                inner_found_nest = true;
                            }
                        }
                        
                        if(!inner_found&&!inner_found_nest){
                            std::string message = "Unsupported option requested in YAML input file: ";
                            error = message + current_option_inner;
                        }

                        //std::cout << "    " << current_option_inner << " " << current_setting_inner << " " << inner_found << " " << inner_found_nest << std::endl;
            
                        // inner_most layer
                        Yaml::Node & center_item = (*inner_it).second;
            
                        //std::cout << "  \n";
                        if (center_item.Size()!=0){
                            for(auto center_it = center_item.Begin(); center_it !=  center_item.End();   center_it++)
                            {
                                current_option_center = (*center_it).first;
                                current_setting_center = (*center_it).second.As<std::string>();

                                //check if this option is supported
                                for(auto temp_it = iterator_range_nested.first; temp_it != iterator_range_nested.second; temp_it++){
                                    //test if map element corresponds to a possible option
                                    std::pair<multimap_iterator_nested2,multimap_iterator_nested2> iterator_range_nested3 = temp_it->second.equal_range(current_option_inner);
                                    for(auto temp_it_inner = iterator_range_nested3.first; temp_it_inner != iterator_range_nested3.second; temp_it_inner++){
                                        //test if map element corresponds to a possible option
                                        multimap_iterator multi_center_iterator = temp_it_inner->second.find(current_option_center);
                                        if(multi_center_iterator!=temp_it_inner->second.end()){
                                            std::string message = "Unsupported option requested in YAML input file: ";
                                            error = message + current_option_center;
                                        }
                                    }
                                }

                                //std::cout << "        " << current_option_center << "   " << current_setting_center << std::endl;

                                set_options[current_option_center] = current_setting_center;
                            } // end for
                        }
                        else{
                            set_options[current_option_inner] = current_setting_inner;
                        }
            
                    } // end for
                } // end if inner_item.Size
                else{
                    //there were no inner items; add setting to map of user defined settings to query later
                    set_options[current_option_outer] = current_setting_outer;
                    //if(current_option_outer=="solver_type"){
                      //std::string test_string = "solver_type";
                      //std::cout << "solver type: " << set_options[test_string] << std::endl;
                    //}

                } //end else for outer item with no inner items
            
            } // end for outer_it
        } // end if outer_it
    //} // end for
    return error;
}

void Simulation_Parameters_SGH::FEA_module_setup(){
  
  //initial buffer size for FEA module list storage
  int buffer_size = 10 + nfea_modules;
  FEA_Module_List.resize(buffer_size);
  fea_module_must_read.resize(buffer_size);
  int start_module = nfea_modules;

  //decides which FEA modules to setup based on user decided implicit solves
  FEA_Module_List[nfea_modules] = "SGH";
  nfea_modules++;
  //example for later
  if(nfea_modules==buffer_size){
    buffer_size += 10;
    FEA_Module_List.resize(buffer_size);
    fea_module_must_read.resize(buffer_size);
  }

  //initialize
  for(int imodule = start_module; imodule < nfea_modules; imodule++){
    fea_module_must_read[imodule] = false;
  }
}
