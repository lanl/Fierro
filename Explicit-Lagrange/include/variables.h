/* variables.h   */
#ifndef VARIABLES_H
#define VARIABLES_H 

extern int num_dim;


// --- Reference element variables ---
extern elements::HexN      elem;
extern elements::ref_element  ref_elem;

// --- Mesh state variables ---
extern elem_state_t elem_state;
extern cell_state_t cell_state;
extern node_t      	node;
extern corner_t     corner;
extern mat_pt_t  	mat_pt;
// extern gauss_pt_t  	gauss_pt;


extern material_t  	*material;
extern mat_fill_t  	*mat_fill;
extern boundary_t  	*boundary;


extern int NR;     // number or Regions
extern int NC;     // number of contours
extern int NF;     // number of fill regions
extern int NB;     // number of boundary patch sets to tag

// --- Method variables ---
extern bool CCH;
extern bool SGH;
extern bool DGH;
extern bool RDH;

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

extern int rk_num_stages;  	// number of rk stages, 1,2,3,4
extern int rk_storage;      // number of bins for rk storage
//extern int rk_stage;		// current rk stage in calculation
extern int num_prediction_steps;  // number of sub time stages in RD
extern int num_correction_steps;  // number of correction stages for RD

extern int cycle;          // the time cycle number
extern int cycle_stop;     // stop calculation at this cycle number
extern int stop_calc;      // a flag to stop the calculation


extern double te_0;
extern double ke;
extern double ie;

// --- Precision variables ---
extern double fuzz;  // machine precision
extern double tiny;  // very very small (between double and single)
extern double small; // single precision


// --- Artificial viscosity coefficients ---
extern double C1;
extern double C2;


// --- Dimensional and mesh constants ---

extern int p_order;
extern int num_bdy_sets;



#endif //VARIABLES_H
