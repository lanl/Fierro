#ifndef SIMULATION_PARAMETERS_TOPOLOGY_OPTIMIZATION_H
#define SIMULATION_PARAMETERS_TOPOLOGY_OPTIMIZATION_H

#include "utilities.h"
#include "Simulation_Parameters.h"
using namespace utils;

class Simulation_Parameters_Topology_Optimization : public Simulation_Parameters
{
 public:
  Simulation_Parameters_Topology_Optimization();
  virtual ~Simulation_Parameters_Topology_Optimization();
  virtual void input();
  //==============================================================================
  //   Mesh Variables
  //==============================================================================

  // --- Mesh regions and material fills ---
  int NB; // number of boundary patch sets to tag
  int NBD; //number of density boundary conditions
  
  // --- Dimensional and mesh constants ---
  int num_dim;
  int p_order;

  //debug and performance reporting flags
  bool report_runtime_flag;

  //Topology Optimization flags
  bool nodal_density_flag;

  //Topology Optimization parameters
  real_t penalty_power;

  //volumes to hold density constant
};

#endif // end HEADER_H
