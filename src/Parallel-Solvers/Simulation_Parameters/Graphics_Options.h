#pragma once
#include "yaml-serializable.h"

struct Graphics_Options : Yaml::DerivedFields {
  double graphics_step;
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
  graphics_step, graphics_dt_ival
)
