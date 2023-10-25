#pragma once
#include "yaml-serializable.h"

struct Graphics_Options : Yaml::DerivedFields {
  bool output_velocity_flag = true;
  bool output_strain_flag   = true;
  bool output_stress_flag   = false;

  bool displaced_mesh_flag  = true;
  bool strain_max_flag      = false;

  double graphics_dt_ival   = 0.25;

  // Non-serialized
  size_t graphics_id = 0;
  double graphics_time;  // the times for writing graphics dump
  CArray <double> graphics_times;

  void derive() {
    graphics_time     = graphics_dt_ival;
    graphics_times    = CArray<double>(2000);
    graphics_times(0) = 0.0;
  }
};
IMPL_YAML_SERIALIZABLE_FOR(Graphics_Options, 
  output_velocity_flag, output_stress_flag, output_strain_flag,
  strain_max_flag, displaced_mesh_flag, graphics_dt_ival
)
