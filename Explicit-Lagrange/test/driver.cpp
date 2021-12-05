#include <stdio.h>
#include <iostream>  // std::cout etc.
#include <cmath>
#include <sys/stat.h>
#include <ctime>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;



//==============================================================================
//   Mesh Variables
//==============================================================================

int num_dim = 3;

// --- Mesh state data declarations ---
elem_state_t    elem_state;
cell_state_t    cell_state;

node_t          node;
corner_t        corner;
mat_pt_t        mat_pt;
// gauss_pt_t      gauss_pt;


material_t  * material;
mat_fill_t  * mat_fill;

boundary_t   * boundary;



// Reference element
//elements::ref_element  ref_elem;

// Everything is HexN for now
//elements::HexN      elem;


// --- Mesh regions and material fills ---
int NR; // number of Regions
int NC; // number of contours
int NF; // number of fill
int NB; // number of boundary patch sets to tag


// Global variables
int p_order = 0;


// --- Graphics output variables ---
int graphics_id = 0;
int graphics_cyc_ival = 50;

real_t graphics_times[250];
real_t graphics_dt_ival = 1.0e8;
real_t graphics_time = graphics_dt_ival;  // the times for writing graphics dump


// --- Choose force routine(s) ---
bool CCH = false;
bool SGH = false;
bool DGH = false;

// --- Time and cycling variables ---
real_t TIME = 0.0;
real_t TFINAL = 1.e16;
real_t dt = 1.e-8;
real_t dt_max = 1.0e-2;
real_t dt_min = 1.0e-8;
real_t dt_cfl = 0.3;
real_t dt_start = 1.0e-8;

int rk_num_stages = 1;
int rk_storage = 2;
int rk_stage = 0;

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


// --- Threading (e.g., OpenMP) coloring order ---
int *cell_order; 



// --- Energy Metrics ---
real_t te_0;
real_t ie;
real_t ke;

//==============================================================================
//    Main
//==============================================================================




// Main function for the testing program
int main(int argc, char *argv[]){

    // check to see of a mesh was supplied when running the code
    if (argc == 1) {
        std::cout << "\n\n**********************************\n\n";
        std::cout << " ERROR:\n";
        std::cout << " Please supply a mesh \n";
        std::cout << "   ./fierro my_mesh.geo \n\n";
        std::cout << "**********************************\n\n" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    
    std::cout<<"MAIN START"<<std::endl;

    // Initialize RK stage
    rk_stage = 0;
    
    // output files
    FILE *out_cell_state;  // cell state values
    FILE *out_energy;      // energy and partiion as a function of time
    FILE *out_matpt_state; // matpt state
	FILE *out_elem_state;  //element average state

    
    // output files
    out_energy      = fopen("energy_file","w");
    out_cell_state  = fopen("cell_state_file", "w");
    out_matpt_state = fopen("matpt_state_file","w");
	out_elem_state  = fopen("elem_state_file", "w");

    
    real_t te_0, ke, ie;
    real_t percent_comp = 0.0;

    
    //----read input----//
    input();

    std::cout << "CCH = "<< CCH << std::endl;
    std::cout << "SGH = "<< SGH << std::endl;
    std::cout << "DGH = "<< DGH << std::endl;

    //---- Read intial mesh and build connectivity ----//
    if(CCH == true) setup_cch(argv[1]);
    if(SGH == true) setup_sgh(argv[1]);
    if(DGH == true) setup_dgh(argv[1]);


    // calculate the total energy at the beginning of the calculation 
    if(CCH == true) track_cch(ke, ie);
    if(SGH == true) track_sgh(ke, ie);
    if(DGH == true) track_dgh(ke, ie);
    

    // Save total energy at time=0
    te_0 = ke + ie;


    std::cout<<"Num cells read in = "<<init_mesh.num_cells()<<std::endl;
    std::cout<<"Num elements in mesh = "<<mesh.num_elems()<<std::endl;
    std::cout<<"Num cells in mesh = "<<mesh.num_cells()<<std::endl; 

    //----intialize time, time_step, and cycles----//
    TIME = 0.0;
    dt = dt_start;



    get_timestep();
 
    std::cout<<"After first timestep call"<<std::endl;

    std::cout<<"First Time Step = "<<dt<<std::endl;
    std::cout<<std::endl;

    // write grophics outputs
    graphics_id = 0;
    graphics_times[0] = 0.0;
    graphics_time = graphics_dt_ival;  // the times for writing graphics dump
    

    


    std::cout<<"Before first ensight "<<std::endl;
    
    ensight();
    std::cout<<"After first ensight "<<std::endl;

    std::cout<<std::endl;
    std::cout<<"Calling Hydro Solver"<<std::endl;
    std::cout<<std::endl;
    std::cout<<std::endl;

    // get the wall clock time
    std::clock_t start_time = std::clock();

    if(CCH == true) cch_hydro();
    
    if(SGH == true) sgh_hydro();

    if(DGH == true) dg_hydro();


    // graphics output after end of hydro solve
    ensight();


    if(CCH == true){

        // print to file energy diagonostics
        fprintf(out_energy,
                "Cycle=%d Time=%e ie=%e ke=%e te_error=%e",
                cycle,
                TIME,
                ie,
                ke,
                te_0-(ie+ke));

        
        // state output for plotting
        fprintf(out_cell_state, "cell_gid x y z r den pres sie vol\n");

        real_t x, y, z;
        std::cout<<"Num Cells: " << mesh.num_cells() << std::endl; 
        for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){

            x = mesh.cell_coords(cell_gid, 0);
            y = mesh.cell_coords(cell_gid, 1);
            z = mesh.cell_coords(cell_gid, 2);


            
            fprintf(out_cell_state, "%i, %e, %e, %e, %e, %e, %e, %e, %e\n",
                    cell_gid,
                    x, y, z,
                    sqrt(x*x + y*y + z*z),
                    cell_state.density(cell_gid), 
                    cell_state.pressure(cell_gid), 
                    cell_state.ie(1, cell_gid),
                    mesh.cell_vol(cell_gid)
                    );
        }
    }

    if(SGH == true){

        // print to file energy diagonostics
        fprintf(out_energy,
                "Cycle=%d Time=%e ie=%e ke=%e te_error=%e",
                cycle,
                TIME,
                ie,
                ke,
                te_0-(ie+ke));

        
        // state output for plotting
        fprintf(out_cell_state, "cell_gid x y z r den pres sie vol\n");

        real_t x, y, z;
        
        for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){

            x = mesh.cell_coords(cell_gid, 0);
            y = mesh.cell_coords(cell_gid, 1);
            z = mesh.cell_coords(cell_gid, 2);


            
            fprintf(out_cell_state, "%i, %e, %e, %e, %e, %e, %e, %e, %e\n",
                    cell_gid,
                    x, y, z,
                    sqrt(x*x + y*y + z*z),
                    cell_state.density(cell_gid), 
                    cell_state.pressure(cell_gid), 
                    cell_state.ie(1, cell_gid),
                    mesh.cell_vol(cell_gid)
                    );
        }

    }

    if(DGH == true){

        // print to file energy diagonostics
        fprintf(out_energy,
                "Cycle=%d Time=%e ie=%e ke=%e te_error=%e",
                cycle,
                TIME,
                ie,
                ke,
                te_0-(ie+ke));


        fprintf(out_elem_state, "x, y, z, r, den, pres, sie, ske, ste \n");


        for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++) {
            
            real_t xe = 0.0, ye = 0.0, ze = 0.0;
            real_t sie = 0.0, ske = 0.0;  // partitions of specific total energy
            real_t pres = 0.0;
            
            real_t num = (real_t)mesh.num_nodes_in_elem();
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
                
                int gauss_lid = node_lid;  // they overlap inside an element
                int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
                int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

                xe += mesh.node_coords(node_gid, 0)/num;
                ye += mesh.node_coords(node_gid, 1)/num;
                ze += mesh.node_coords(node_gid, 2)/num;
                
                sie += mat_pt.ie(gauss_gid)/num;
                ske += mat_pt.ke(gauss_gid)/num;
                pres += mat_pt.pressure(gauss_gid)/num;
            }
            
            real_t r = sqrt( xe*xe + ye*ye + ze*ze);

            fprintf(out_elem_state, "%e, %e, %e, %e, %e, %e, %e, %e, %e \n",
                    xe, ye, ze, r,
                    elem_state.avg_density(elem_gid),
                    pres,
                    sie,
                    ske,
                    elem_state.avg_specific_total_energy(elem_gid)
                    );

		} //end loop over elements

        
        // state output for plotting
        fprintf(out_matpt_state, "x, y, z, r, den, pres, sie, vol \n");

        real_t x, y, z;
        for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){

            for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

                int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);

                int node_lid = gauss_lid;
                int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);


                x = mesh.node_coords(node_gid, 0);
                y = mesh.node_coords(node_gid, 1);
                z = mesh.node_coords(node_gid, 2);

				real_t r = sqrt(x*x + y*y + z*z);

                fprintf(out_matpt_state, "%e, %e, %e, %e, %e, %e, %e, %e\n",
                    x, y, z, r,
                    mat_pt.density(gauss_gid), 
                    mat_pt.pressure(gauss_gid), 
                    mat_pt.ie(gauss_gid),
                    mesh.gauss_pt_det_j(gauss_gid) * mat_pt.weight(gauss_gid)
                    );

            }
        }

    }
    
    
    fclose(out_energy);      // energy and partion as a function of time
    fclose(out_cell_state);  // cell state values
    fclose(out_matpt_state);
    fclose(out_elem_state);





    // calculate the total energy at the end of the calculation
    if(CCH == true) track_cch(ke, ie);
    if(SGH == true) track_sgh(ke, ie);
    if(DGH == true) track_dgh(ke, ie);
    

    std::cout<<"Kinetic Energy at time = " << TIME << " is = " << ke <<std::endl;
    std::cout<<"Internal Energy at time = " << TIME << " is = " << ie <<std::endl;
    std::cout<<"Total energy at time = 0 = "<< te_0 <<std::endl;
    std::cout<<"Total energy at time = " << TIME << " is = "<< ke+ie <<std::endl;

    std::cout<<"Energy Error = "<< (te_0) - (ke+ie) <<std::endl;

    // get the wall clock time
    std::clock_t end_time = std::clock();
    std::cout << "Program Run Time = " << (end_time-start_time)/CLOCKS_PER_SEC << " sec " << std::endl;


    std::cout<<"MAIN END"<<std::endl;



    return 0;
}// end  of main function

