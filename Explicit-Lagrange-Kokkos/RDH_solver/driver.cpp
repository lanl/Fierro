// -----------------------------------------------------------------------------
// This is the main function
//------------------------------------------------------------------------------
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <sys/stat.h>
#include <ctime>

#include "ref_elem.h"
#include "mesh.h"
#include "state.h"
#include "matar.h"

using namespace mtr;


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
size_t graphics_cyc_ival = 50;

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
        std::cout << "   ./FierroRDH my_mesh.geo \n\n";
        std::cout << "**********************************\n\n" << std::endl;
        std::exit(EXIT_FAILURE);
    } // end if


    printf("Starting Lagrangian RDH code\n");


    // The kokkos scope
    Kokkos::initialize();
    {
     
        // ---------------------------------------------------------------------
        //    state data type declarations
        // ---------------------------------------------------------------------
        node_t  node;
        elem_t  elem;
	corner_t corner;
        ref_elem_t ref_elem;
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
        input(material,
              mat_fill,
              boundary,
              state_vars,
              num_materials,
              num_fills,
              num_bcs,
              num_dims,
              num_state_vars,
              dt_start,
              time_final,
              dt_max,
              dt_min,
              dt_cfl,
              graphics_dt_ival,
              graphics_cyc_ival,
              cycle_stop,
              rk_num_stages
              );
        



        // ---------------------------------------------------------------------
        //    read in supplied mesh
        // ---------------------------------------------------------------------
        std::string str;
        std::string delimiter = ".";
        std::string arg1(argv[1]);
        std::vector<std::string> split_file_str = split (arg1, delimiter);
        size_t len = split_file_str.size();
        if(split_file_str[len-1]=="geo"){
            read_mesh_ensight(argv[1], mesh, node, elem, corner, num_dims, rk_num_bins);
        }
        else
        {
            // arbitrary order elements
            //printf("Calling readVTKPn \n");
            readVTKPn(argv[1], mesh, node, elem, corner, ref_elem, num_dims, rk_num_bins);
        }
        
        // write VTKPn
        printf("writing VTK file \n");
        VTKHexN(mesh, node);
        
	printf("building corners \n");
        mesh.build_corner_connectivity();
        printf("building elem_elem \n");
        mesh.build_elem_elem_connectivity();
        printf("building patches \n");
        mesh.build_patch_connectivity();
        printf("building node_node \n");
        mesh.build_node_node_connectivity();
        printf("done building connectivity \n");
        
        // ---------------------------------------------------------------------
        //    allocate memory
        // ---------------------------------------------------------------------

        // shorthand names
        const size_t num_nodes = mesh.num_nodes;
        const size_t num_elems = mesh.num_elems;
        const size_t num_corners = mesh.num_corners;
        const size_t num_zones = mesh.num_elems*mesh.num_zones_in_elem;
        const size_t num_lob_pts = elem.num_lob_pts*num_elems;
        const size_t num_leg_pts = elem.num_leg_pts*num_elems;
        
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
        
        DViewCArrayKokkos <double> node_div(&node.div(0),
                                            num_nodes);
        
        
        // create Dual Views of the individual elem struct variables
        DViewCArrayKokkos <double> elem_den(&elem.den(0),
                                            num_zones);

        DViewCArrayKokkos <double> elem_pres(&elem.pres(0),
                                             num_zones);

        DViewCArrayKokkos <double> elem_stress(&elem.stress(0,0,0,0),
                                               rk_num_bins,
                                               num_zones,
                                               3,
                                               3); // always 3D even in 2D-RZ

        DViewCArrayKokkos <double> elem_sspd(&elem.sspd(0),
                                             num_zones);

        DViewCArrayKokkos <double> elem_sie(&elem.sie(0,0),
                                            rk_num_bins,
                                            num_zones);

        DViewCArrayKokkos <double> elem_div(&elem.div(0),
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
        
	DViewCArrayKokkos <double> lobatto_jacobian(&elem.gauss_lobatto_jacobian(0,0,0),
                                                    num_lob_pts,
                                                    num_dims,
                                                    num_dims);
	
	DViewCArrayKokkos <double> legendre_jacobian(&elem.gauss_legendre_jacobian(0,0,0),
                                                    num_leg_pts,
                                                    num_dims,
                                                    num_dims);

	DViewCArrayKokkos <double> lobatto_jacobian_inverse(&elem.gauss_lobatto_jacobian_inverse(0,0,0),
                                                            num_lob_pts,
                                                            num_dims,
                                                            num_dims);
	
	DViewCArrayKokkos <double> legendre_jacobian_inverse(&elem.gauss_legendre_jacobian_inverse(0,0,0),
                                                            num_leg_pts,
                                                            num_dims,
                                                            num_dims);
	
	DViewCArrayKokkos <double> lobatto_det(&elem.gauss_lobatto_det_j(0),
                                               num_lob_pts);

	DViewCArrayKokkos <double> legendre_det(&elem.gauss_legendre_det_j(0),
                                               num_leg_pts);

        // create Dual Views of the corner struct variables
        DViewCArrayKokkos <double> corner_force(&corner.force(0,0),
                                                num_corners, 
                                                num_dims);

        DViewCArrayKokkos <double> corner_mass(&corner.mass(0),
                                               num_corners);
        
        
        // ---------------------------------------------------------------------
        //   calculate geometry
        // ---------------------------------------------------------------------
        node_coords.update_device();
        Kokkos::fence();

        get_gauss_leg_pt_jacobian(mesh,
                                  elem,
                                  ref_elem,
                                  node_coords,
                                  legendre_jacobian,
                                  legendre_det,
                                  legendre_jacobian_inverse);
        
        get_vol(elem_vol, node_coords, mesh, elem, ref_elem);
        
        /* 
        double vol_check = 0.0;
        for (int i = 0; i < mesh.num_elems; i++){
           vol_check += elem_vol(i);
        }
        printf("calculated volume is: %f \n", vol_check); 
        */
       /* 
        // check jacobian inverse works //
        double temp_left = 0.0;
        double temp_right = 0.0;
        
        std::cout << "left inverse " << std::endl;
        for (int i = 0; i < num_leg_pts; i++){
          std::cout << " At gauss pt " << i << std::endl;
          std::cout << " ######################## " << std::endl;
          for (int dim_1 = 0; dim_1 < mesh.num_dims; dim_1++){
            for (int dim_2 = 0; dim_2 < mesh.num_dims; dim_2++){
              for (int k = 0; k < mesh.num_dims; k++){
                temp_left += legendre_jacobian_inverse(i,dim_1,k)*legendre_jacobian(i, k, dim_2); 
              }
              std::cout<<  temp_left << ", ";
              temp_left = 0.0;
            }
            std::cout<< " "<< std::endl;
          }
          std::cout << " ######################## " << std::endl;
        }
        
        std::cout << "right inverse " << std::endl;
        for (int i = 0; i < num_leg_pts; i++){
          std::cout << " At gauss pt " << i << std::endl;
          std::cout << " ######################## " << std::endl;
          for (int dim_1 = 0; dim_1 < mesh.num_dims; dim_1++){
            for (int dim_2 = 0; dim_2 < mesh.num_dims; dim_2++){
              for (int k = 0; k < mesh.num_dims; k++){
                temp_right += legendre_jacobian(i,dim_1,k)*legendre_jacobian_inverse(i, k, dim_2); 
              }
              std::cout<< temp_right <<", ";
              temp_right = 0.0;
            }
            std::cout<< " "<< std::endl;
          }
          std::cout << " ######################## " << std::endl;
        }
        */
        // ---------------------------------------------------------------------
        //   setup the IC's and BC's
        // ---------------------------------------------------------------------
        setup(material,
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
              corner_mass,
              num_fills,
              rk_num_bins,
              num_bcs,
              num_materials,
              num_state_vars);
        
        // intialize time, time_step, and cycles
        time_value = 0.0;
        dt = dt_start;
        graphics_id = 0;
        graphics_times(0) = 0.0;
        graphics_time = graphics_dt_ival;  // the times for writing graphics dump
        
        
        // --- testing ---
        // testing high-order mesh initialization
        // --- testing ---
        VTKHexN(mesh,
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
        
        

        // ---------------------------------------------------------------------
        //   Calculate the RDH solution
        // ---------------------------------------------------------------------
        rdh_solve(material,
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
		  elem_div,
                  elem_mass,
                  elem_mat_id,
                  elem_statev,
                  corner_force,
                  corner_mass,
                  time_value,
                  time_final,
                  dt_max,
                  dt_min,
                  dt_cfl,
                  graphics_time,
                  graphics_cyc_ival,
                  graphics_dt_ival,
                  cycle_stop,
                  rk_num_stages,
                  dt,
                  fuzz,
                  tiny,
                  small,
                  graphics_times,
                  graphics_id);


        // calculate total energy at time=t_end
        
        
        VTKHexN(mesh,
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
