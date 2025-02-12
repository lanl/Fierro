#pragma once
#include "yaml-serializable.h"
#include "Fields.h"
#include <set>

SERIALIZABLE_ENUM(OUTPUT_FORMAT, vtk, pvtu, none)
SERIALIZABLE_ENUM(TIMER_VERBOSITY, standard, thorough)
struct Output_Options : Yaml::DerivedFields, Yaml::ValidatedYaml {
  TIMER_VERBOSITY timer_output_level = TIMER_VERBOSITY::standard;
  bool include_default_output_fields = true;
  std::set<FIELD> output_fields;

  OUTPUT_FORMAT output_file_format = OUTPUT_FORMAT::vtk;
  size_t max_num_user_output_vars=0;
  bool write_initial = true;
  bool write_final   = true;
  bool convert_to_vtk = false;
  bool convert_to_tecplot = false;
  bool restart_file = false;
  bool optimization_restart_file = false;
  size_t restart_step_interval = 100000;
  size_t optimization_restart_step_interval = 0;
    
  std::string output_file_location = "vtk/";
};

IMPL_YAML_SERIALIZABLE_FOR(Output_Options, 
  timer_output_level, output_fields, include_default_output_fields,
  output_file_format, write_initial, write_final, max_num_user_output_vars,
  convert_to_vtk, convert_to_tecplot, output_file_location, restart_file,
  optimization_restart_file, restart_step_interval, optimization_restart_step_interval
)
