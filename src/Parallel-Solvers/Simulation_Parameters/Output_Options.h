#pragma once
#include "yaml-serializable.h"
#include <string>

SERIALIZABLE_ENUM(OUTPUT_FORMAT, vtk, vtu, none)

struct Output_Options : Yaml::DerivedFields, Yaml::ValidatedYaml {
  bool output_velocity = true;
  bool output_strain   = true;
  bool output_stress   = true;
  bool output_displacement = true;

  bool displaced_mesh      = true;

  bool output_temperature = false;
  bool output_temperature_gradient = false;
  bool output_heat_flux = false;
  
  bool strain_max       = false;

  OUTPUT_FORMAT output_file_format = OUTPUT_FORMAT::vtk;
  size_t max_num_user_output_vars=0;
  bool write_initial = true;
  bool write_final = true;

  void validate() {
    if (displaced_mesh && !output_displacement)
      throw Yaml::ConfigurationException("`displaced_mesh` cannot be set to true if `output_displacement` is false.");
  }
};

IMPL_YAML_SERIALIZABLE_FOR(Output_Options, 
    output_file_format, output_velocity, output_strain, output_displacement,
    output_stress, displaced_mesh, 
    write_initial, write_final, max_num_user_output_vars
)
