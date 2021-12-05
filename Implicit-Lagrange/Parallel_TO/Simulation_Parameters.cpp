#include "utilities.h"
#include "Simulation_Parameters.h"
#include "elements.h"
#include "swage.h"

using namespace utils;

Simulation_Parameters::Simulation_Parameters(){

  //initialize data and flags to defaults
  output_strain_flag = 0;
  report_runtime_flag = 0;
  nodal_density_flag = 1;
  unit_scaling = 1;
  strain_max_flag = 0;
  optimization_output_freq = 100;
  direct_solver_flag = false;
  penalty_power = 3;
  restart_file = false;
  tecplot_input = false;
}

Simulation_Parameters::~Simulation_Parameters(){
}

void Simulation_Parameters::input(){
  //file input flags
  tecplot_input = false;
  restart_file = false;

  //simulation spatial dimension
  num_dim = 3;
  unit_scaling = 1;

  // time step
  TFINAL = 0.5;
  
  dt = 0.001;
  dt_min = 1.e-8;
  dt_max = 1.e-2;
  dt_start = 1.e-5;
    
  cycle_stop = 100000000;

  rk_num_stages = 1;

  rk_storage = 2;

  dt_cfl = 0.5;

  //polynomial interpolation order
  p_order = 0;

  // ---- graphics information ---- //
  graphics_cyc_ival = 10;
  graphics_dt_ival  = 0.1;    // every 0.1 time units
  output_strain_flag = 1;

  // ---- fill instructions and intial conditions ---- //

  NF = 2; // number of fills
  
  //Static isotropic parameters to move into a child class later
  Elastic_Modulus = 10000;
  Poisson_Ratio = 0.3;
  words_per_line = 1;
  tecplot_words_per_line = 3;
  elem_words_per_line = 8;

  //ensight file readin for Isotropical elasticity

  num_gauss_points = 2;

  //debug and performance report flags
  report_runtime_flag = 1;

  //Topology Optimization flags
  nodal_density_flag = 1;

  //Topology Optimization parameters
  penalty_power = 6;
  maximum_strain = 0.02;
  maximum_strain_energy = 1000;

  // ---- boundary conditions ---- //
  NB = 5; // number of boundaries
  NBSF = 4; //number of surface density force conditions
  NBD = 1; //number of surface sets used to specify a fixed displacement on nodes belonging to respective surfaces
           //note this only implies a fixed displacement on the surface if no other basis functions have support on the surface

}