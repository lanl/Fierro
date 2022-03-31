#include "utilities.h"
#include "Simulation_Parameters_Inertial.h"

using namespace utils;

Simulation_Parameters_Inertial::Simulation_Parameters_Inertial() : Simulation_Parameters(){

  //initialize data and flags to defaults
  report_runtime_flag = false;
  unit_scaling = 1;
}

Simulation_Parameters_Inertial::~Simulation_Parameters_Inertial(){
}

void Simulation_Parameters_Inertial::input(){
  
  Simulation_Parameters::input();

  //simulation spatial dimension
  num_dim = 3;
  unit_scaling = 1;

  //polynomial interpolation order
  p_order = 0;
  
  //Gauss-Legendre integration order
  num_gauss_points = 2;

  //debug and performance report flags
  report_runtime_flag = true;

}