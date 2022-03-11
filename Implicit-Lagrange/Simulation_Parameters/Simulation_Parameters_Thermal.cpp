#include "utilities.h"
#include "Simulation_Parameters_Thermal.h"

using namespace utils;

Simulation_Parameters_Thermal::Simulation_Parameters_Thermal() : Simulation_Parameters(){

  //initialize data and flags to defaults
  output_heat_flux_flag = false;
  report_runtime_flag = false;
  unit_scaling = 1;
  flux_max_flag = false;
  direct_solver_flag = false;
  thermal_flag = false;
  multigrid_timers = false;
  equilibrate_matrix_flag = false;
}

Simulation_Parameters_Thermal::~Simulation_Parameters_Thermal(){
}

void Simulation_Parameters_Thermal::input(){

  //multigrid_timers = true;
  equilibrate_matrix_flag = false;

  //simulation spatial dimension
  num_dim = 3;
  unit_scaling = 1;

  //polynomial interpolation order
  p_order = 0;

  // ---- graphics information ---- //
  output_heat_flux_flag = true;
  
  //Isotropic Conductivity parameters to move into a child class later
  Thermal_Conductivity = 17;

  //Gauss-Legendre parameters
  num_gauss_points = 2;

  //debug and performance report flags
  report_runtime_flag = true;

  // ---- boundary conditions ---- //
  NB = 6; // number of boundary conditions for this module
  NBSF = 4; //number of surface heat flux conditions
  NBT = 2; //number of surface sets used to specify a fixed temperature on nodes belonging to respective surfaces

  //apply body forces
  thermal_flag = false;
  specific_internal_energy_rate = 1;

}