#ifndef SIMULATION_PARAMETERS_INERTIAL_H
#define SIMULATION_PARAMETERS_INERTIAL_H

#include "utilities.h"
#include "Simulation_Parameters.h"
using namespace utils;

class Simulation_Parameters_Inertial : public Simulation_Parameters
{
 public:
  Simulation_Parameters_Inertial();
  virtual ~Simulation_Parameters_Inertial();
  virtual void input();
  
  //==============================================================================
  //   Mesh Variables
  //==============================================================================

  // --- Graphics output variables ---
  bool output_flag;

  // -- Integration rule
  int num_gauss_points;

  // --- Dimensional and mesh constants ---
  int num_dim;
  int p_order;
  real_t unit_scaling;

  //debug and performance reporting flags
  bool report_runtime_flag;
};

#endif // end HEADER_H
