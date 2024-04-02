

#ifndef FIERRO_DYNAMIC_OPTIONS_H
#define FIERRO_DYNAMIC_OPTIONS_H
#include <stdio.h>
#include "matar.h"




// Time and cycle options
struct dynamic_options_t
{

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
    
    int rk_num_stages = 2;
    int rk_num_bins   = 2;

}; // output_options_t

// ----------------------------------
// valid inputs for dynamic options
// ----------------------------------
static std::vector <std::string> str_dyn_opts_inps
{
    "time_initial",
    "time_final",
    "dt_min",
    "dt_max",
    "dt_start",
    "dt_cfl",
    "cycle_stop",
    "fuzz",
    "tiny",
    "small",
    "rk_num_stages",
    "rk_num_bins"
};



#endif // end Header Guard