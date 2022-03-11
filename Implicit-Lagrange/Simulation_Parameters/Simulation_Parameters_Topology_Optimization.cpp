#include "utilities.h"
#include "Simulation_Parameters_Topology_Optimization.h"
#include "elements.h"
#include "swage.h"

using namespace utils;

Simulation_Parameters_Topology_Optimization::Simulation_Parameters_Topology_Optimization(){

  //initialize data and flags to defaults
  report_runtime_flag = false;
  nodal_density_flag = true;
  penalty_power = 3;
}

Simulation_Parameters_Topology_Optimization::~Simulation_Parameters_Topology_Optimization(){
}

void Simulation_Parameters_Topology_Optimization::input(){

  //simulation spatial dimension
  num_dim = 3;

  //polynomial interpolation order
  p_order = 0;

  //debug and performance report flags
  report_runtime_flag = true;

  //Topology Optimization flags
  nodal_density_flag = true;

  //Topology Optimization parameters
  penalty_power = 3;

  // ---- boundary conditions ---- //
  NB = 6; // number of boundaries
  NBD = 2; //number of surface sets used to specify a fixed density

}