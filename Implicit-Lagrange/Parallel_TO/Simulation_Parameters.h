#ifndef ELEMENTS_SIMULATION_PARAMETERS_H
#define ELEMENTS_SIMULATION_PARAMETERS_H

#include "utilities.h"
using namespace utils;

class Simulation_Parameters
{
 public:
  Simulation_Parameters();
  virtual ~Simulation_Parameters();
  virtual void input();
  //==============================================================================
  //   Mesh Variables
  //==============================================================================

  // --- Mesh regions and material fills ---
  int NR; // number of Regions
  int NC; // number of contours
  int NF; // number of fill
  int NB; // number of boundary patch sets to tag
  int NBSF; //number of surface force density boundary conditions
  int NBD; //number of displacement boundary conditions


  // --- Graphics output variables ---
  int graphics_id;
  int graphics_cyc_ival;
  int output_strain_flag, strain_max_flag;

  real_t graphics_times[250];
  real_t graphics_dt_ival;
  real_t graphics_time;  // the times for writing graphics dump


  // --- Time and cycling variables ---
  real_t TIME;
  real_t TFINAL;
  real_t dt;
  real_t dt_max;
  real_t dt_min;
  real_t dt_cfl;
  real_t dt_start;
  real_t Elastic_Modulus, Poisson_Ratio;
  int num_gauss_points;

  int rk_num_stages;
  int rk_storage;
  int rk_stage;

  int cycle;
  int cycle_stop;
  int stop_calc;    // a flag to end the calculation when = 1

  real_t percent_comp;

  // --- Precision variables ---
  real_t fuzz;  // machine precision
  real_t tiny;  // very very small (between real_t and single)
  real_t small;   // single precision


  // --- Dimensional and mesh constants ---
  int num_dim;
  int p_order;

  //file input parameters 
  int words_per_line, elem_words_per_line, tecplot_words_per_line;
  char *format_specification;  //per line file format when reading dofs
  real_t unit_scaling;
  bool restart_file;
  bool tecplot_input;

  //file output parameters
  int optimization_output_freq;

  //debug and performance reporting flags
  int report_runtime_flag;

  //Body force parameters
  bool gravity_flag;
  real_t gravity_vector[3];

  //Linear Solver Flags
  bool direct_solver_flag, multigrid_timers, equilibrate_matrix_flag;

  //Topology Optimization flags
  bool nodal_density_flag;

  //Topology Optimization parameters
  real_t maximum_strain, maximum_strain_energy, penalty_power;
};

#endif // end HEADER_H
