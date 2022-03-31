#ifndef SIMULATION_PARAMETERS_ELASTICITY_H
#define SIMULATION_PARAMETERS_ELASTICITY_H

#include "utilities.h"
#include "Simulation_Parameters.h"
using namespace utils;

class Simulation_Parameters_Elasticity : public Simulation_Parameters
{
 public:
  Simulation_Parameters_Elasticity();
  virtual ~Simulation_Parameters_Elasticity();
  virtual void input();
  //==============================================================================
  //   Mesh Variables
  //==============================================================================

  // --- Mesh regions and material fills ---
  int NB; // number of boundary patch sets to tag
  int NBSF; //number of surface force density boundary conditions
  int NBD; //number of displacement boundary conditions


  // --- Graphics output variables ---
  bool output_displacement_flag, output_stress_flag, output_strain_flag, strain_max_flag, displaced_mesh_flag;

  // --- Isotropic Elastic Parameters
  real_t Elastic_Modulus, Poisson_Ratio;

  // -- Integration rule
  int num_gauss_points;

  // --- Dimensional and mesh constants ---
  int num_dim;
  int p_order;
  real_t unit_scaling;

  //debug and performance reporting flags
  bool report_runtime_flag;

  //Body force parameters
  bool gravity_flag;
  real_t gravity_vector[3];

  //Linear Solver Flags
  bool direct_solver_flag, multigrid_timers, equilibrate_matrix_flag;
};

#endif // end HEADER_H
