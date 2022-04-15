#ifndef VARIABLES_H
#define VARIABLES_H 

#include "state.h"
#include "mesh.h"

// --- State variables ---
extern node_t  node;
extern elem_t  elem;
extern corner_t  corner;
extern material_t material;

extern size_t num_state_vars;  // number of variables in consitutive model


// --- Mesh ---
extern mesh_t mesh;

extern size_t num_dims;

extern size_t num_materials;
extern size_t num_fills;
extern size_t num_boundaries;



// --- Method variables ---


// --- Graphics output variables ---
extern int graphics_id;
extern int graphics_cyc_ival;

extern double graphics_times[];
extern double graphics_dt_ival;
extern double graphics_time;



// --- Time and cycling variables ---
extern double time_value;
extern double time_final;  // the final time of the simulation
extern double dt;          // time step size
extern double dt_start;    // the initial time step size
extern double dt_min;
extern double dt_max;
extern double dt_cfl;

extern size_t rk_num_stages;  	// number of rk stages, 1,2,3,4
extern size_t rk_num_bins;      // number of bins for rk storage


extern size_t cycle;          // the time cycle number
extern size_t cycle_stop;     // stop calculation at this cycle number
extern size_t stop_calc;      // a flag to stop the calculation



// --- Precision variables ---
extern double fuzz;  // machine precision
extern double tiny;  // very very small (between double and single)
extern double small; // single precision



// --- Dimensional and mesh constants ---

extern size_t num_bdy_sets;



#endif //VARIABLES_H
