#ifndef SIMULATION_PARAMETERS_THERMAL_H
#define SIMULATION_PARAMETERS_THERMAL_H

#include "utilities.h"
#include "Simulation_Parameters.h"
using namespace utils;

class Simulation_Parameters_Thermal: public Simulation_Parameters
{
 public:
  Simulation_Parameters_Thermal();
  virtual ~Simulation_Parameters_Thermal();
  virtual void input();
  //==============================================================================
  //   Thermal FEA problem parameters
  //==============================================================================

  // --- Boundary Conditions ---
  int NB; // number of boundary patch sets to tag
  int NBSF; //number of surface heat flux boundary conditions
  int NBT; //number of temperature boundary conditions

  // --- Graphics output variables ---
  bool output_heat_flux_flag, flux_max_flag;

  // --- Constitutive Parameters ---
  real_t Thermal_Conductivity;

  // --- Integration Scheme
  int num_gauss_points;

  // --- Dimensional and mesh constants ---
  int num_dim;
  int p_order;
  real_t unit_scaling;

  //debug and performance reporting flags
  bool report_runtime_flag, multigrid_timers, direct_solver_flag;

  //Body force parameters
  bool thermal_flag, electric_flag;
  real_t specific_internal_energy_rate;

  //Linear Solver Flags
  bool equilibrate_matrix_flag;

  //Topology Optimization parameters
  real_t maximum_strain, maximum_strain_energy;
};

#endif // end HEADER_H
