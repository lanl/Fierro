#include <stdio.h>
#include <iostream>
#include <cmath>
#include <sys/stat.h>
#include <ctime>

#include "state.h"
#include "variables.h"
#include "matar.h"



//==============================================================================
//   Variables
//==============================================================================

size_t num_dims = 3;

// --- state data declarations ---
node_t  node;


// --- Mesh regions and material fills ---
int num_regions;   // number or Regions
int num_contours;  // number of contours
int num_fills;     // number of fill regions
int num_boundarys; // number of boundary patch sets to tag


// --- Graphics output variables ---
int graphics_id = 0;
int graphics_cyc_ival = 50;

real_t graphics_times[2000];
real_t graphics_dt_ival = 1.0e8;
real_t graphics_time = graphics_dt_ival;  // the times for writing graphics dump


// --- Time and cycling variables ---
real_t TIME = 0.0;
real_t TFINAL = 1.e16;
real_t dt = 1.e-8;
real_t dt_max = 1.0e-2;
real_t dt_min = 1.0e-8;
real_t dt_cfl = 0.3;
real_t dt_start = 1.0e-8;

size_t rk_num_stages = 1;
size_t rk_storage = 2;


int cycle = 0;
int cycle_stop = 1000000000;
int stop_calc = 0;    // a flag to end the calculation when = 1


// --- Precision variables ---
real_t fuzz = 1.0e-16;  // machine precision
real_t tiny = 1.0e-12;  // very very small (between real_t and single)
real_t small= 1.0e-8;   // single precision


// --- Artificial viscosity coefficients ---
real_t C1 = 0.0;
real_t C2 = 0.0;


// --- Energy Metrics ---
real_t te_0;



//==============================================================================
//    Main
//==============================================================================
int main(int argc, char *argv[]){


    printf("Starting Lagrangian SGH code\n");

    
    // allocating node state
    
    rk_storage = 2;
    size_t num_nodes = 100;  // WARNING: a filler
    num_dims = 3;
    
    
    node.initialize(rk_storage, num_nodes, num_dims);
    
    
    // create Dual Views of the node struct variables
    auto node_coords =  DViewCArrayKokkos <double> (&node.coords(0,0,0),
                                           rk_storage, num_nodes, num_dims);
    //DViewCArrayKokkos <double> node_vel(&node.coords(0,0,0),
   //                                     rk_storage, num_nodes,num_dims);
   // DViewCArrayKokkos <double> node_mass(&node.coords(0,0,0), num_nodes);
    
    // WARNING: just filling
    for (size_t rk_lid=0; rk_lid<rk_storage;  rk_lid++){
        for (size_t node_gid=0; node_gid<num_nodes; node_gid++){
            for (size_t dim=0; dim<num_dims; dim++){
                node.coords(rk_lid, node_gid, dim) = 0.0;
                node.vel(rk_lid, node_gid, dim) = 0.0;
                node.mass(node_gid) = 0.0;
            } // end for
        } // end for
    } // end for
    
    

    
    
    // output files
	FILE *out_elem_state;  //element average state

    
    // output files
	out_elem_state  = fopen("elem_state_file", "w");


    

    
    // Save total energy at time=0
    double ke = 0.0;
    double ie = 0.0;
    te_0 = ke + ie;




    //----intialize time, time_step, and cycles----//
    TIME = 0.0;
    dt = dt_start;


    // write grophics outputs
    graphics_id = 0;
    graphics_times[0] = 0.0;
    graphics_time = graphics_dt_ival;  // the times for writing graphics dump
    

    


    
    fclose(out_elem_state);





    printf("Finished\n");


    return 0;
}// end of main function

