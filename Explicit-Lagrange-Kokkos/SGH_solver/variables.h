#ifndef VARIABLES_H
#define VARIABLES_H 

#include "state.h"
#include "mesh.h"




// --- Graphics output variables ---
extern int graphics_cyc_ival;
extern double graphics_dt_ival;
extern double graphics_time;


// --- Time and cycling variables ---
extern double time_final;  // the final time of the simulation
extern double dt;          // time step size
extern double dt_start;    // the initial time step size
extern double dt_min;
extern double dt_max;
extern double dt_cfl;

extern size_t cycle;          // the time cycle number
extern size_t cycle_stop;     // stop calculation at this cycle number
extern size_t stop_calc;      // a flag to stop the calculation


// --- Precision variables ---
extern double fuzz;  // machine precision
extern double tiny;  // very very small (between double and single)
extern double small; // single precision




#endif //VARIABLES_H
