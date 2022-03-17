#include "utilities.h"
#include "Simulation_Parameters_Elasticity.h"

using namespace utils;

Simulation_Parameters_Elasticity::Simulation_Parameters_Elasticity() : Simulation_Parameters(){

  //initialize data and flags to defaults
  output_displacement_flag = false;
  output_strain_flag = false;
  output_stress_flag = false;
  displaced_mesh_flag = false;
  report_runtime_flag = false;
  unit_scaling = 1;
  strain_max_flag = false;
  direct_solver_flag = false;
  gravity_flag = false;
  multigrid_timers = false;
  equilibrate_matrix_flag = false;
}

Simulation_Parameters_Elasticity::~Simulation_Parameters_Elasticity(){
}

void Simulation_Parameters_Elasticity::input(){

  //output settings
  output_displacement_flag = false;
  displaced_mesh_flag = true;
  output_strain_flag = true;
  output_stress_flag = false;

  //multigrid_timers = true;
  equilibrate_matrix_flag = false;

  //simulation spatial dimension
  num_dim = 3;
  unit_scaling = 1;

  //polynomial interpolation order
  p_order = 0;
  
  //Static isotropic parameters to move into a child class later
  Elastic_Modulus = 10000;
  Poisson_Ratio = 0.3;
  
  //Gauss-Legendre integration order
  num_gauss_points = 2;

  //debug and performance report flags
  report_runtime_flag = true;

  // ---- boundary conditions ---- //
  NB = 6; // number of boundaries
  NBSF = 4; //number of surface density force conditions
  NBD = 2; //number of surface sets used to specify a fixed displacement on nodes belonging to respective surfaces

  //apply body forces
  gravity_flag = false;
  gravity_vector[0] = 9.81;
  gravity_vector[1] = 0;
  gravity_vector[2] = 0;

}