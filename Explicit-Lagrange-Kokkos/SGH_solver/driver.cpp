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




//==============================================================================
//   Variables, setting default inputs
//==============================================================================



// --- num vars ----
size_t num_dims = 3;

size_t num_materials;
size_t num_state_vars;

size_t num_fills;
size_t num_bcs;


// --- Graphics output variables ---
size_t graphics_id = 0;
int graphics_cyc_ival = 50;

CArray <double> graphics_times(2000);
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
     
        // ---------------------------------------------------------------------
        //    state data type declarations
        // ---------------------------------------------------------------------
        node_t  node;
        elem_t  elem;
        corner_t  corner;
        CArrayKokkos <material_t> material;
        CArrayKokkos <double> state_vars; // array to hold init model variables
        
        
        // ---------------------------------------------------------------------
        //    mesh data type declarations
        // ---------------------------------------------------------------------
        mesh_t mesh;
        CArrayKokkos <mat_fill_t> mat_fill;
        CArrayKokkos <boundary_t> boundary;
        

        // ---------------------------------------------------------------------
        //    read the input file
        // ---------------------------------------------------------------------  
        input(material, mat_fill, boundary, state_vars,
              num_materials, num_fills, num_bcs,
              num_dims, num_state_vars);


        // ---------------------------------------------------------------------
        //    read in supplied mesh
        // --------------------------------------------------------------------- 
        read_mesh_ensight(argv[1], mesh, node, elem, corner, num_dims, rk_num_bins);
        mesh.build_corner_connectivity();
        mesh.build_elem_elem_connectivity();
        mesh.build_patch_connectivity();
        
        
        // ---------------------------------------------------------------------
        //    allocate memory
        // ---------------------------------------------------------------------

        // shorthand names
        const size_t num_nodes = mesh.num_nodes;
        const size_t num_elems = mesh.num_elems;
        const size_t num_corners = mesh.num_corners;

        
        // allocate elem_statev
        elem.statev = CArray <double> (num_elems, num_state_vars);

        // --- make dual views of data on CPU and GPU ---
        //  Notes:
        //     Instead of using a struct of dual types like the mesh type, 
        //     individual dual views will be made for all the state 
        //     variables.  The motivation is to reduce memory movement 
        //     when passing state into a function.  Passing a struct by 
        //     reference will copy the meta data and pointers for the 
        //     variables held inside the struct.  Since all the mesh 
        //     variables are typically used by most functions, a single 
        //     mesh struct or passing the arrays will be roughly equivalent 
        //     for memory movement.

        
        // create Dual Views of the individual node struct variables
        DViewCArrayKokkos <double> node_coords(&node.coords(0,0,0),
                                               rk_num_bins,
                                               num_nodes,
                                               num_dims);

        DViewCArrayKokkos <double> node_vel(&node.vel(0,0,0),
                                            rk_num_bins,
                                            num_nodes,
                                            num_dims);

        DViewCArrayKokkos <double> node_mass(&node.mass(0),
                                             num_nodes);
        
        
        // create Dual Views of the individual elem struct variables
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
        DViewCArrayKokkos <double> corner_force(&corner.force(0,0,0), 
                                                num_corners, 
                                                num_dims);

        DViewCArrayKokkos <double> corner_mass (&corner.mass(0), 
                                                num_corners);
        
        
        // ---------------------------------------------------------------------
        //   calculate geometry
        // ---------------------------------------------------------------------
        node_coords.update_device();
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
            get_vol_hex(elem_vol, elem_gid, node_coords, mesh);
        });
        Kokkos::fence();


        // ---------------------------------------------------------------------
        //   setup the IC's and BC's
        // ---------------------------------------------------------------------
        setup( material,
               mat_fill,
               boundary,
               mesh,
               node_coords,
               node_vel,
               node_mass,      
               elem_den,
               elem_pres,
               elem_stress,
               elem_sspd,       
               elem_sie,
               elem_vol,
               elem_mass,
               elem_mat_id,
               elem_statev,
               state_vars,
               num_fills,
               rk_num_bins,
               num_bcs
           );
        
        // intialize time, time_step, and cycles
        time_value = 0.0;
        dt = dt_start;
        graphics_id = 0;
        graphics_times(0) = 0.0;
        graphics_time = graphics_dt_ival;  // the times for writing graphics dump

        
        // ---------------------------------------------------------------------
        //   t=0 ensight and state output
        // ---------------------------------------------------------------------
        elem_den.update_device();
        elem_pres.update_device();
        elem_stress.update_device();
        elem_sspd.update_device();
        elem_sie.update_device();
        elem_vol.update_device();
        elem_mass.update_device();
        elem_mat_id.update_device();

        // write out ensight file
        ensight( mesh,
                 node_coords,
                 node_vel,
                 node_mass,
                 elem_den,
                 elem_pres,
                 elem_stress,
                 elem_sspd, 
                 elem_sie,
                 elem_vol,
                 elem_mass,
                 elem_mat_id,
                 graphics_times,
                 graphics_id,
                 time_value);

        // output files
        FILE *out_elem_state;  //element average state

        // output files
        out_elem_state  = fopen("elem_state_t0", "w");

        fclose(out_elem_state);
        

        // Calculate total energy at time=0
        double ke = 0.0;
        double ie = 0.0;
        double te_0 = ke + ie;

        // get_timestep();

        // call hydro here

        // calculate total energy at time=t_end
        
        
    } // end of kokkos scope
    Kokkos::finalize();
    


    printf("Finished\n");


    return 0;
}// end of main function


// for easy copy and paste to functions
//DViewCArrayKokkos <double> node_coords,
//DViewCArrayKokkos <double> node_vel,
//DViewCArrayKokkos <double> node_mass,
//DViewCArrayKokkos <double> elem_den,
//DViewCArrayKokkos <double> elem_pres,
//DViewCArrayKokkos <double> elem_stress,
//DViewCArrayKokkos <double> elem_sspd, 
//DViewCArrayKokkos <double> elem_sie,
//DViewCArrayKokkos <double> elem_vol,
//DViewCArrayKokkos <double> elem_mass,
//DViewCArrayKokkos <size_t> elem_mat_id,
//DViewCArrayKokkos <double> elem_statev
