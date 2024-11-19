// -----------------------------------------------------------------------------
// This is the main function
//------------------------------------------------------------------------------
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <sys/stat.h>
#include <ctime>

#include "ref_elem.h"
#include "ref_surf_elem.h"
#include "mesh.h"
#include "rdh.h"
#include "state.h"
#include "matar.h"
#include "linear_algebra.h"

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
double dt = 1.e-5;
double dt_max = 1.0e-2;
double dt_min = 1.0e-8;
double fe_order = 1.0;
double dt_cfl = 1.0/(2.0*fe_order+1.0);
double dt_start = 1.0e-5;

size_t rk_num_stages = 2;
size_t rk_num_bins = 2;

size_t cycle = 0;
size_t cycle_stop = 1000000000;


// --- Precision variables ---
double fuzz = 1.0e-16;  // machine precision
double tiny = 1.0e-12;  // very very small (between real_t and single)
double small= 1.0e-8;   // single precision

bool viscosity_cond = false;
bool source_cond = false;


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
        zone_t zone;
        mat_pt_t mat_pt;
	    corner_t corner;
        fe_ref_elem_t ref_elem;
        fe_ref_surf_t ref_surf;
        CArrayKokkos <material_t> material;
        CArrayKokkos <double> state_vars; // array to hold init model variables
        
        
        // ---------------------------------------------------------------------
        //    mesh data type declarations
        // ---------------------------------------------------------------------
        mesh_t mesh;
        CArrayKokkos <mat_fill_t> mat_fill;
        CArrayKokkos <boundary_t> boundary;
        
        //printf("Before input \n");
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
              rk_num_stages,
              viscosity_cond,
              source_cond);
        

        printf("After input \n");

        // ---------------------------------------------------------------------
        //    read in supplied mesh
        // ---------------------------------------------------------------------
        std::string str;
        std::string delimiter = ".";
        std::string arg1(argv[1]);
        std::vector<std::string> split_file_str = split (arg1, delimiter);
        size_t len = split_file_str.size();
        if(split_file_str[len-1]=="geo"){
            //read_mesh_ensight(argv[1], mesh, node, elem, corner, num_dims, rk_num_bins);
        }
        else
        {   
            //printf("inside conditional to readVTKPn \n");
            // arbitrary order elements
            //printf("Calling readVTKPn \n");
            readVTKPn(argv[1], mesh, node, elem, zone, mat_pt, corner, ref_elem, num_dims, rk_num_bins);
            //printf("after readVTKPn \n");
        }

        
        // write VTKPn
        //printf("writing VTK file \n");
        //VTKHexN(mesh, node);
        
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
        const size_t num_lob_pts = mat_pt.num_lob_pts;//*num_elems;
        const size_t num_leg_pts = mat_pt.num_leg_pts;//*num_elems;
        const size_t num_lob_pts_per_elem = ref_elem.num_gauss_lob_in_elem;
        const size_t num_leg_pts_per_elem = ref_elem.num_gauss_leg_in_elem;
        // printf(" num_leg_pts_per_elem = %d\n", num_leg_pts_per_elem);
        //printf(" num_zones_in_elem = %zu and num_zones = %zu \n", mesh.num_zones_in_elem, num_zones);
        
        // allocate elem_statev
        elem.statev = CArray <double> (num_elems, num_state_vars);
        // mat_pt.statev = CArrayKokkos <double> (num_leg_pts, num_state_vars);
        
        
        // create Dual Views of the individual node struct variables
        DViewCArrayKokkos <double> node_coords(&node.coords(0,0,0),
                                               rk_num_bins,
                                               num_nodes,
                                               num_dims);

        // DViewCArrayKokkos <double> mat_pt_coords(&mat_pt.coords(0,0),
        //                                         num_leg_pts,
        //                                         num_dims);

        DViewCArrayKokkos <double> node_vel(&node.vel(0,0,0),
                                            rk_num_bins,
                                            num_nodes,
                                            num_dims);

        DViewCArrayKokkos <double> M_u(&elem.M(0,0,0),
                                            num_elems,
                                            mesh.num_nodes_in_elem,
                                            mesh.num_nodes_in_elem);
        
        DViewCArrayKokkos <double> F_u(&elem.F_u(0,0,0,0),
                                            rk_num_stages,
                                            num_elems,
                                            mesh.num_nodes_in_elem,
                                            mesh.num_dims);
        DViewCArrayKokkos <double> PHI(&elem.PHI(0,0,0,0),
                                            rk_num_stages,
                                            num_elems,
                                            mesh.num_nodes_in_elem,
                                            mesh.num_dims);

        DViewCArrayKokkos <double> F_e(&elem.F_e(0,0,0),
                                            rk_num_stages,
                                            num_elems,
                                            mesh.num_zones_in_elem);

        DViewCArrayKokkos <double> S(&zone.S(0,0,0),
                                            rk_num_stages,
                                            num_elems,
                                            mesh.num_zones_in_elem);

        DViewCArrayKokkos <double> PSI(&elem.PSI(0,0,0),
                                            rk_num_stages,
                                            num_elems,
                                            mesh.num_zones_in_elem);

        // DViewCArrayKokkos <double> mat_pt_vel(&mat_pt.vel(0,0),
        //                                     num_leg_pts,
        //                                     num_dims);

        DViewCArrayKokkos <double> node_mass(&node.mass(0),
                                             num_nodes);
        
        DViewCArrayKokkos <double> node_div(&node.div(0),
                                            num_nodes);
        
        
        // create Dual Views of the individual elem struct variables
        DViewCArrayKokkos <double> mat_pt_den(&mat_pt.den(0),
                                            num_leg_pts);

        DViewCArrayKokkos <double> mat_pt_pressure(&mat_pt.pres(0),
                                             num_leg_pts);

        DViewCArrayKokkos <double> mat_pt_stress(&mat_pt.stress(0,0,0,0),
                                               rk_num_bins,
                                               num_leg_pts,
                                               3,
                                               3); // always 3D even in 2D-RZ

        DViewCArrayKokkos <double> mat_pt_sspd(&mat_pt.sspd(0),
                                             num_leg_pts);

        DViewCArrayKokkos <double> zone_sie(&zone.sie(0,0),
                                            rk_num_bins,
                                            num_zones);
        
        DViewCArrayKokkos <double> M_e(&zone.M(0,0,0),
                                            num_elems,
                                            mesh.num_zones_in_elem,
                                            mesh.num_zones_in_elem);

        DViewCArrayKokkos <double> zone_mass(&zone.zonal_mass(0),
                                             num_zones);

        DViewCArrayKokkos <double> mat_pt_sie(&mat_pt.sie(0),
                                             num_leg_pts);

        // DViewCArrayKokkos <double> mat_pt_div(&mat_pt.div(0),
        //                                     num_leg_pts);

        DViewCArrayKokkos <double> elem_vol(&elem.vol(0),
                                            num_elems);
        

        // DViewCArrayKokkos <double> mat_pt_mass(&mat_pt.mass(0),
        //                                      num_leg_pts);

        DViewCArrayKokkos <int> elem_mat_id(&elem.mat_id(0),
                                               num_elems);

        DViewCArrayKokkos <double> elem_statev(&elem.statev(0,0),
                                               num_elems,
                                               num_state_vars);
        

        // create Dual Views of the corner struct variables
        // DViewCArrayKokkos <double> corner_force(&corner.force(0,0),
        //                                         num_corners, 
        //                                         num_dims);

        // DViewCArrayKokkos <double> corner_mass(&corner.mass(0),
        //                                        num_corners);

        DViewCArrayKokkos <double> mat_pt_h(&mat_pt.h(0),
                                                num_leg_pts);

        DViewCArrayKokkos <double> gauss_legendre_jacobian(&mat_pt.gauss_legendre_jacobian(0,0,0),
                                                            num_leg_pts, 3, 3);

        DViewCArrayKokkos <double> gauss_legendre_jacobian_inverse(&mat_pt.gauss_legendre_jacobian_inverse(0,0,0),
                                                            num_leg_pts, 3, 3);

        DViewCArrayKokkos <double> SigmaJacInv(&mat_pt.SigmaJacInv(0,0,0,0),
                                                rk_num_bins, num_leg_pts, 3, 3);

        DViewCArrayKokkos <double> gauss_legendre_det_j(&mat_pt.gauss_legendre_det_j(0), num_leg_pts);
        
        DViewCArrayKokkos <double> den0DetJac0(&mat_pt.den0DetJac0(0), num_leg_pts);

        DViewCArrayKokkos <double> h0(&mat_pt.h0(0), num_leg_pts);

        DViewCArrayKokkos <double> Jac0Inv(&mat_pt.Jac0Inv(0),
                                                num_leg_pts, 3, 3);
        
        
        // ---------------------------------------------------------------------
        //   calculate geometry
        // ---------------------------------------------------------------------
        node_coords.update_device();
        Kokkos::fence();

        get_gauss_leg_pt_jacobian(mesh,
                                  elem,
                                  ref_elem,
                                  node_coords,
                                  gauss_legendre_jacobian,
                                  gauss_legendre_det_j,
                                  gauss_legendre_jacobian_inverse,
                                  0);
        Kokkos::fence();

        gauss_legendre_jacobian.update_host();
        gauss_legendre_jacobian_inverse.update_host();
        gauss_legendre_det_j.update_host();

        get_J0Inv(gauss_legendre_jacobian, Jac0Inv, mat_pt);
        Kokkos::fence();
        
        Jac0Inv.update_host();

        get_vol(elem_vol, node_coords, gauss_legendre_det_j, mesh, elem, ref_elem);
        Kokkos::fence();

        get_h0(elem_vol, h0, mesh, ref_elem);
        Kokkos::fence();

        h0.update_host();

        // double vol_check = 0.0;
        // for (int i = 0; i < mesh.num_elems; i++){
        //    vol_check += elem_vol(i);
        // }
        // printf("calculated volume is: %f \n", vol_check); 
        
        
        // // check jacobian inverse works //
        // double temp_left = 0.0;
        // double temp_right = 0.0;
        
        // std::cout << "left inverse " << std::endl;
        // for (int i = 0; i < num_leg_pts; i++){
        //   std::cout << " At gauss pt " << i << std::endl;
        //   std::cout << " ######################## " << std::endl;
        //   for (int dim_1 = 0; dim_1 < mesh.num_dims; dim_1++){
        //     for (int dim_2 = 0; dim_2 < mesh.num_dims; dim_2++){
        //       for (int k = 0; k < mesh.num_dims; k++){
        //         temp_left += mat_pt.gauss_legendre_jacobian_inverse(i,dim_1,k)*mat_pt.gauss_legendre_jacobian(i, k, dim_2); 
        //       }
        //       std::cout<<  temp_left << ", ";
        //       temp_left = 0.0;
        //     }
        //     std::cout<< " "<< std::endl;
        //   }
        //   std::cout << " ######################## " << std::endl;
        // }
        
        // std::cout << "right inverse " << std::endl;
        // for (int i = 0; i < num_leg_pts; i++){
        //   std::cout << " At gauss pt " << i << std::endl;
        //   std::cout << " ######################## " << std::endl;
        //   for (int dim_1 = 0; dim_1 < mesh.num_dims; dim_1++){
        //     for (int dim_2 = 0; dim_2 < mesh.num_dims; dim_2++){
        //       for (int k = 0; k < mesh.num_dims; k++){
        //         temp_right += mat_pt.gauss_legendre_jacobian(i,dim_1,k)*mat_pt.gauss_legendre_jacobian_inverse(i, k, dim_2); 
        //       }
        //       std::cout<< temp_right <<", ";
        //       temp_right = 0.0;
        //     }
        //     std::cout<< " "<< std::endl;
        //   }
        //   std::cout << " ######################## " << std::endl;
        // }
        
        // ---------------------------------------------------------------------
        //   setup the IC's and BC's
        // ---------------------------------------------------------------------
        setup(material,
              mat_fill,
              boundary,
              mesh,
              elem,
              zone,
              mat_pt,
              ref_elem,
              node_coords,
              node_vel,
              node_mass,
              mat_pt_den,
              mat_pt_pressure,
              mat_pt_stress,
              mat_pt_sspd,
              zone_sie,
              mat_pt_sie,
              elem_vol,
              elem_mat_id,
              elem_statev,
              state_vars,
              num_fills,
              rk_num_bins,
              num_bcs,
              num_materials,
              num_state_vars);
        Kokkos::fence();
        printf("after setup \n");


        // intialize time, time_step, and cycles
        time_value = 0.0;
        dt = dt_start;
        graphics_id = 0;
        graphics_times(0) = 0.0;
        graphics_time = graphics_dt_ival;  // the times for writing graphics dump
        
        mat_pt_den.update_host();

        printf("before den0 deet(J_0) \n");
        get_den0DetJac0(den0DetJac0, mat_pt_den, gauss_legendre_det_j, mat_pt);
        Kokkos::fence();

        den0DetJac0.update_host();
        mat_pt_pressure.update_host();
        
        printf("before stress \n");
        // initalize stress
        for (int stage = 0; stage < rk_num_bins; stage++){
            get_stress(mesh, mat_pt, mat_pt_pressure, mat_pt_stress, stage);
        }
        Kokkos::fence();

        mat_pt_stress.update_host();
        // initialize sigma J^{-1}
        printf("before sigma J^{-1} \n");
        for (int stage = 0; stage < rk_num_bins; stage++){
            get_SigmaJacInv(mesh, mat_pt, mat_pt_stress, gauss_legendre_jacobian_inverse, SigmaJacInv, stage);
        }
        Kokkos::fence();

        SigmaJacInv.update_host();

        mat_pt_sspd.update_host();
        // mat_pt_vel.update_host();
        // mat_pt_coords.update_host();
        mat_pt_sie.update_host();
        mat_pt_h.update_host();
        zone_sie.update_host();
        elem_vol.update_host();
        // mat_pt_mass.update_host();
        elem_mat_id.update_host();
        
        node_coords.update_host();
        node_vel.update_host();
        node_div.update_host();

        Kokkos::fence();

        // printf("before vtk \n");
        VTKHexN(mesh,
                node_coords,
                node_vel,
                node_mass,
                mat_pt_den,
                mat_pt_pressure,
                mat_pt_stress,
                mat_pt_sspd,
                zone_sie,
                elem_vol,
                elem_mat_id,
                graphics_times,
                graphics_id,
                time_value);
        
        // assemble consistent mass matrices (per element) //
        printf("before mass matrices \n");
        assemble_mass_matrices(mesh, ref_elem, mat_pt_den, gauss_legendre_det_j,
                               M_u, M_e);

        M_u.update_host();
        M_e.update_host();

        // compute lumped masses //
        printf("before lumped mass \n");
        get_lumped_mass(mesh, ref_elem, gauss_legendre_det_j, mat_pt_den, M_u, M_e, node_mass, zone_mass);
        node_mass.update_host();
        zone_mass.update_host();
        


        //////////////////////////////////////////////////
        ////////////////// Run Solver ////////////////////
        //////////////////////////////////////////////////
        //////////////////////////////////////////////////
        printf("running solver\n");
        run(mesh,
            ref_elem,
            ref_surf,
            elem,
            mat_pt,
            boundary,
            node_vel,
            node_mass,
            M_u,
            zone_sie,
            zone_mass,
            M_e,
            PHI,
            PSI,
            F_u,
            F_e,
            S,
            node_coords,
            gauss_legendre_jacobian,
            gauss_legendre_jacobian_inverse,
            gauss_legendre_det_j,
            Jac0Inv,
            h0,
            mat_pt_stress,
            SigmaJacInv,
            mat_pt_den,
            den0DetJac0,
            mat_pt_pressure,
            mat_pt_sspd,
            elem_vol,
            elem_mat_id,
            elem_statev,
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
            graphics_id,
            viscosity_cond,
            source_cond);


        //////////////////////////////////////////////////
        ////////////////// End Solver ////////////////////
        //////////////////////////////////////////////////


        mat_pt_den.update_host();
        mat_pt_pressure.update_host();
        mat_pt_stress.update_host();
        mat_pt_sspd.update_host();
        // mat_pt_vel.update_host();
        // mat_pt_coords.update_host();
        mat_pt_sie.update_host();
        mat_pt_h.update_host();
        zone_sie.update_host();
        zone_mass.update_host();
        elem_vol.update_host();
        // mat_pt_mass.update_host();
        elem_mat_id.update_host();
        
        node_coords.update_host();
        node_vel.update_host();
        node_mass.update_host();
        node_div.update_host();

        gauss_legendre_jacobian.update_host();
        gauss_legendre_jacobian_inverse.update_host();
        gauss_legendre_det_j.update_host();

        Kokkos::fence();

        state_file( mesh, node_coords, node_vel,
                    mat_pt_h, node_mass, mat_pt_den, mat_pt_pressure, mat_pt_stress,
                    mat_pt_sspd, mat_pt_sie, elem_vol,
                    elem_mat_id, time_value );
            
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
//DViewCArrayKokkos <double> elem_pressure,
//DViewCArrayKokkos <double> elem_stress,
//DViewCArrayKokkos <double> elem_sspd, 
//DViewCArrayKokkos <double> elem_sie,
//DViewCArrayKokkos <double> elem_vol,
//DViewCArrayKokkos <double> elem_mass,
//DViewCArrayKokkos <size_t> elem_mat_id,
//DViewCArrayKokkos <double> elem_statev
