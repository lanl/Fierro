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
  tecplot_input = ansys_dat_input = false;
  nfea_modules = 0;
  element_type = "Hex8";
}

Simulation_Parameters::~Simulation_Parameters(){
}

void Simulation_Parameters::input(){
  //file input flags
  tecplot_input = false;
  restart_file = false;
  ansys_dat_input = true;

  //simulation spatial dimension
  num_dim = 3;
  unit_scaling = 1;
  
  //file readin parameters
  words_per_line = 1;
  tecplot_words_per_line = 3;
  ansys_dat_node_words_per_line = 4;
  elem_words_per_line = 8;
  ansys_dat_elem_words_per_line = 8;
  element_type = "Hex8";

  //debug and performance report flags
  report_runtime_flag = 1;

}

void Simulation_Parameters::FEA_module_setup(){
  
  //initial buffer size for FEA module list storage
  int buffer_size = 10;
  FEA_Module_List = std::vector<std::string>(buffer_size);

  //decides which FEA modules to setup based on user decided implicit solves
  FEA_Module_List[0] = "Elasticity";
  nfea_modules++;
  FEA_Module_List[1] = "Heat_Conduction";
  nfea_modules++;
  //example for later
  if(nfea_modules==buffer_size){
    buffer_size += 10;
    FEA_Module_List.resize(buffer_size);
  }
}
