#pragma once
#include "yaml-serializable.h"

SERIALIZABLE_ENUM(TIME_OUTPUT_LEVEL,
    none,   // output no time sequence information of forward solve
    low,    // output just final timestep of forward solve
    high,   // output time sequence for forward solve
    extreme // output time sequence for all phases of solver
)

struct Dynamic_Options : Yaml::DerivedFields {
    unsigned long cycle_stop = 2000000;
    double time_initial = 0.0;
    double time_final = 1.0;
    double dt_min     = 1e-8;
    double dt_max     = 1e-2;
    double dt_start   = 1e-5;
    double dt_cfl     = 0.4;
    double fuzz       = 1e-16; // machine precision
    double tiny       = 1e-12; // very very small (between real_t and single)
    double small      = 1e-8;  // single precision
    TIME_OUTPUT_LEVEL output_time_sequence_level = TIME_OUTPUT_LEVEL::high;
    int rk_num_stages = 2;
    int rk_num_bins   = -1;
    double time_value = -1;

    // Non-serialized Fields
    double dt;
    void derive() {
      if (time_value < 0)
        time_value = time_initial;
      dt = dt_start;

      if (rk_num_bins < 0)
        rk_num_bins = rk_num_stages;
    }
};
IMPL_YAML_SERIALIZABLE_FOR(Dynamic_Options, output_time_sequence_level,
  time_initial, time_value, time_final, dt_min, dt_max, dt_start, dt_cfl,
  cycle_stop, fuzz, tiny, small, rk_num_stages, rk_num_bins
)