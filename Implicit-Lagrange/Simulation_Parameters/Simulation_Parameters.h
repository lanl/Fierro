#ifndef IMPLICIT_SIMULATION_PARAMETERS_H
#define IMPLICIT_SIMULATION_PARAMETERS_H

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

  // --- Graphics output variables ---
  int output_strain_flag;

  // --- Dimensional and mesh constants ---
  int num_dim;

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
};

#endif // end HEADER_H
