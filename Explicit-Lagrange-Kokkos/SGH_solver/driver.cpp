// -----------------------------------------------------------------------------
// This is the main function
//------------------------------------------------------------------------------
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <sys/stat.h>
#include <ctime>


#include "mesh.h"
#include "state.h"
#include "variables.h"
#include "matar.h"

size_t num_dims = 3;


//==============================================================================
//   Variables, setting default inputs
//==============================================================================

size_t num_materials;
size_t num_state_vars = 6;  // gamma_law has 6 parameters


size_t num_fills;
size_t num_boundaries;



// --- Graphics output variables ---
int graphics_id = 0;
int graphics_cyc_ival = 50;

double graphics_times[2000];
double graphics_dt_ival = 1.0e8;
double graphics_time = graphics_dt_ival;  // the times for writing graphics dump


// --- Time and cycling variables ---
double time_value = 0.0;
double time_final = 1.e16;
double dt = 1.e-8;
double dt_max = 1.0e-2;
double dt_min = 1.0e-8;
double dt_cfl = 0.4;
double dt_start = 1.0e-8;

size_t rk_num_stages = 2;
size_t rk_num_bins = 2;

size_t cycle = 0;
size_t cycle_stop = 1000000000;
size_t stop_calc = 0;    // a flag to end the calculation when = 1


// --- Precision variables ---
double fuzz = 1.0e-16;  // machine precision
double tiny = 1.0e-12;  // very very small (between real_t and single)
double small= 1.0e-8;   // single precision




//==============================================================================
//    main
//==============================================================================
int main(int argc, char *argv[]){
    
    
    // check to see of a mesh was supplied when running the code
    if (argc == 1) {
        std::cout << "\n\n**********************************\n\n";
        std::cout << " ERROR:\n";
        std::cout << " Please supply a mesh \n";
        std::cout << "   ./fierro my_mesh.geo \n\n";
        std::cout << "**********************************\n\n" << std::endl;
        std::exit(EXIT_FAILURE);
    } // end if


    printf("Starting Lagrangian SGH code\n");


    
    
    
    // The kokkos scope
    Kokkos::initialize();
    {
        // --- state data declarations ---
        node_t  node;
        elem_t  elem;
        corner_t  corner;
        CArrayKokkos <material_t> material;
        
        
        // --- mesh declaration ---
        mesh_t mesh;
        CArrayKokkos <mat_fill_t> mat_fill;
        CArrayKokkos <boundary_t> boundary;
        
        
        // --- mesh supplied mesh ---
        read_mesh_ensight(argv[1], mesh, node, elem, corner, num_dims);
        mesh.build_corner_connectivity();
        
        
        // --- read the input file ---
        CArrayKokkos <double> state_vars0; // array to hold init model variables
        input(material, mat_fill, boundary, state_vars0);
        
        
        // short hand names
        size_t num_nodes = mesh.num_nodes;
        size_t num_elems = mesh.num_elems;
        size_t num_corners = mesh.num_corners;
        
        // allocate elem_statev
        elem.statev = CArray <double> (num_elems, num_state_vars);
        
        
        // --- make dual views of data on CPU and GPU ---

        // create Dual Views of the node struct variables
        DViewCArrayKokkos <double> node_coords(&node.coords(0,0,0),
                                               rk_num_bins,
                                               num_nodes,
                                               num_dims);
        DViewCArrayKokkos <double> node_vel(&node.coords(0,0,0),
                                            rk_num_bins,
                                            num_nodes,
                                            num_dims);
        DViewCArrayKokkos <double> node_mass(&node.coords(0,0,0),
                                             num_nodes);
        
        // create Dual Views of the elem struct variables
        DViewCArrayKokkos <double> elem_den(&elem.den(0),
                                            num_elems);
        DViewCArrayKokkos <double> elem_pres(&elem.pres(0),
                                             num_elems);
        DViewCArrayKokkos <double> elem_stress(&elem.stress(0,0,0,0),
                                               rk_num_bins,
                                               num_elems,
                                               num_dims,
                                               num_dims);
        DViewCArrayKokkos <double> elem_sspd(&elem.sspd(0),
                                             num_elems);
        DViewCArrayKokkos <double> elem_sie(&elem.sie(0),
                                            rk_num_bins,
                                            num_elems);
        DViewCArrayKokkos <double> elem_vol(&elem.vol(0),
                                            num_elems);
        DViewCArrayKokkos <double> elem_mass(&elem.mass(0),
                                             num_elems);
        DViewCArrayKokkos <size_t> elem_mat_id(&elem.mat_id(0),
                                               num_elems);
        DViewCArrayKokkos <double> elem_statev(&elem.statev(0,0),
                                               num_elems,
                                               num_state_vars );
        
        
        // create Dual Views of the corner struct variables
        DViewCArrayKokkos <double> corner_force(&corner.force(0,0,0), num_corners, num_dims);
        DViewCArrayKokkos <double> corner_mass (&corner.mass(0), num_corners);
        
        
        // --- calculate geometry ---
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
            get_vol_hex(elem_vol, elem_gid, node_coords, mesh);
        });

        // --- setup the IC's and BC's ---
        // call setup here
        
        
        
        
        // output files
        FILE *out_elem_state;  //element average state
        
        
        // output files
        out_elem_state  = fopen("elem_state_file", "w");
        
        
        
        
        
        // Save total energy at time=0
        double ke = 0.0;
        double ie = 0.0;
        double te_0 = ke + ie;
        
        
        
        
        //----intialize time, time_step, and cycles----//
        time_value = 0.0;
        dt = dt_start;
        
        
        // write grophics outputs
        graphics_id = 0;
        graphics_times[0] = 0.0;
        graphics_time = graphics_dt_ival;  // the times for writing graphics dump
        
        fclose(out_elem_state);
        
        
        
    } // end of kokkos scope
    Kokkos::finalize();
    


    printf("Finished\n");


    return 0;
}// end of main function

