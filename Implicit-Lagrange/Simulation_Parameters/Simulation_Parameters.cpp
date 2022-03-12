#include "utilities.h"
#include "Simulation_Parameters.h"
#include "elements.h"
#include "swage.h"

using namespace utils;

Simulation_Parameters::Simulation_Parameters(){

  //initialize data and flags to defaults
  report_runtime_flag = 0;
  unit_scaling = 1;
  optimization_output_freq = 100;
  restart_file = false;
  tecplot_input = false;
}

Simulation_Parameters::~Simulation_Parameters(){
}

void Simulation_Parameters::input(){
  //file input flags
  tecplot_input = false;
  restart_file = false;

  //simulation spatial dimension
  num_dim = 3;
  unit_scaling = 1;
    
  // ---- graphics information ---- //
  output_strain_flag = 1;
  
  //file readin parameters
  words_per_line = 1;
  tecplot_words_per_line = 3;
  elem_words_per_line = 8;

  //debug and performance report flags
  report_runtime_flag = 1;

}