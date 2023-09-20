#pragma once
#include "yaml-serializable.h"
#include <string>

SERIALIZABLE_ENUM(OUTPUT_FORMAT, vtk, vtu, none)

struct Output_Options : Yaml::DerivedFields {
  int graphics_step_frequency;
  double graphics_step;
  OUTPUT_FORMAT output_file_format;
  size_t max_num_user_output_vars=0;
  bool write_initial = true;
  bool write_final = true;
};

IMPL_YAML_SERIALIZABLE_FOR(Output_Options, 
    graphics_step_frequency, graphics_step, output_file_format, 
    write_initial, write_final, max_num_user_output_vars
)
