#ifndef VARIABLES_H
#define VARIABLES_H 

extern size_t num_dims;


// --- Mesh state variables ---
extern node_t  node;


extern int num_regions;   // number or Regions
extern int num_contours;  // number of contours
extern int num_fills;     // number of fill regions
extern int num_boundarys; // number of boundary patch sets to tag

// --- Method variables ---


// --- Graphics output variables ---
extern int graphics_id;
extern int graphics_cyc_ival;

extern double graphics_times[];
extern double graphics_dt_ival;
extern double graphics_time;



// --- Time and cycling variables ---
extern double TIME;
extern double TFINAL;      // the final time of the simulation
extern double dt;          // time step size
extern double dt_start;    // the initial time step size
extern double dt_min;
extern double dt_max;
extern double dt_cfl;

extern size_t rk_num_stages;  	// number of rk stages, 1,2,3,4
extern size_t rk_storage;      // number of bins for rk storage


extern int cycle;          // the time cycle number
extern int cycle_stop;     // stop calculation at this cycle number
extern int stop_calc;      // a flag to stop the calculation


extern double te_0;


// --- Precision variables ---
extern double fuzz;  // machine precision
extern double tiny;  // very very small (between double and single)
extern double small; // single precision


// --- Artificial viscosity coefficients ---
extern double C1;
extern double C2;


// --- Dimensional and mesh constants ---

extern int num_bdy_sets;



#endif //VARIABLES_H
