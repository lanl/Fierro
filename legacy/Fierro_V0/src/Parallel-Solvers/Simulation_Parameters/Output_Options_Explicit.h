#pragma once
#include "yaml-serializable.h"
#include "Output_Options.h"
#include "matar.h"

struct Output_Options_Explicit : Output_Options {
  double graphics_step    = 0.25;
  double graphics_dt_ival = 0.25;
  int graphics_cyc_ival   = 1000000;

  // Non-serialized
  size_t graphics_id = 0;
  double graphics_time = 0;  // the times for writing graphics dump
  mtr::CArray <double> graphics_times;

  void derive() {
    graphics_time     = graphics_dt_ival;
    graphics_times    = mtr::CArray<double>(2000);
    graphics_times(0) = 0.0;
  }
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Output_Options_Explicit, Output_Options,
  graphics_step, graphics_dt_ival, graphics_cyc_ival
)
