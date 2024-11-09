#include "matar.h"
#include "mesh.h"
#include "rdh.h"
#include "state.h"
#include "ref_elem.h"
#include "ref_surf_elem.h"
#include <cmath>
#include <chrono>

#define PI 3.141592653589793

using namespace mtr;

void run(const mesh_t &mesh,
         const fe_ref_elem_t &ref_elem,
         const fe_ref_surf_t &ref_surf,
         const elem_t &elem,
         const mat_pt_t &mat_pt,
         CArrayKokkos <boundary_t> &boundary,
         DViewCArrayKokkos <double> &node_vel,
         DViewCArrayKokkos <double> &node_mass,
         DViewCArrayKokkos <double> &M_u,
         DViewCArrayKokkos <double> &zone_sie,
         DViewCArrayKokkos <double> &zone_mass,
         DViewCArrayKokkos <double> &M_e,
         DViewCArrayKokkos <double> &PHI,
         DViewCArrayKokkos <double> &PSI,
         DViewCArrayKokkos <double> &F_u,
         DViewCArrayKokkos <double> &F_e,
         DViewCArrayKokkos <double> &S,
         DViewCArrayKokkos <double> &node_coords,
         DViewCArrayKokkos <double> &Jacobian,
         DViewCArrayKokkos <double> &JacInv,
         DViewCArrayKokkos <double> &DetJac,
         DViewCArrayKokkos <double> &stress,
         DViewCArrayKokkos <double> &SigmaJacInv,
         DViewCArrayKokkos <double> &den,
         DViewCArrayKokkos <double> &den0DetJac0,
         DViewCArrayKokkos <double> &pres,
         DViewCArrayKokkos <double> &sspd,
         DViewCArrayKokkos <double> &elem_vol,
         DViewCArrayKokkos <int> &elem_mat_id,
         const DViewCArrayKokkos <double> &elem_state_vars,
         double &time_value,
         const double time_final,
         const double dt_max,
         const double dt_min,
         const double dt_cfl,
         double &graphics_time,
         size_t graphics_cyc_ival,
         double graphics_dt_ival,
         const size_t cycle_stop,
         const size_t num_stages,
         double dt,
         const double fuzz,
         const double tiny,
         const double small,
         CArray <double> &graphics_times,
         size_t &graphics_id,
         bool &viscosity_cond,
         bool &source_cond){

    
    size_t stop_calc=0;
    
    auto time_1 = std::chrono::high_resolution_clock::now();

    CArrayKokkos <double> time_int_weights(2,2);


    // HARD CODED FOR NOW //
    RUN({
        time_int_weights(0,0) = 1.0;
        time_int_weights(1,0) = 0.5;
        time_int_weights(1,1) = 0.5;
    });
    
	// loop over the max number of time integration cycles
	for (size_t cycle = 0; cycle < cycle_stop; cycle++) {

        // stop calculation if flag
	    if (stop_calc == 1) break;

        if (cycle==0){
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        }
        // print time step every 10 cycles
        else if (cycle%20==0){
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        } // end if

        init_tn(mesh, node_coords, node_vel, zone_sie, stress, num_stages);

        get_timestep_HexN(mesh, node_coords, node_vel, sspd, elem_vol,
                         time_value, graphics_time, time_final, dt_max,
                         dt_min, dt_cfl, dt, fuzz);

        for (int stage = 0; stage < num_stages; stage++){

            get_gauss_leg_pt_jacobian(mesh,
                                        elem,
                                        ref_elem,
                                        node_coords,
                                        Jacobian,
                                        DetJac,
                                        JacInv,
                                        stage);
            Jacobian.update_host();
            DetJac.update_host();
            JacInv.update_host();

            pointwise_mass_conservation(den, DetJac, den0DetJac0, mat_pt);

            den.update_host();

            update_thermo(mesh, ref_elem, zone_sie, den, pres, stress, sspd, stage, elem_state_vars);

            pres.update_host();
            sspd.update_host();

            get_stress(mesh, mat_pt, pres, stress, stage);

            stress.update_host();

            get_SigmaJacInv(mesh, mat_pt, stress, JacInv, SigmaJacInv, stage);

            SigmaJacInv.update_host();

            get_momentum_rhs(mesh, ref_elem, ref_surf, SigmaJacInv, DetJac, F_u, stage, viscosity_cond);

            F_u.update_host();

            get_momentum_residual(mesh, M_u, node_vel, F_u, PHI, dt, stage, time_int_weights);

            PHI.update_host();

            update_momentum(node_vel, PHI, node_mass, mesh, stage);
            // FOR_ALL(node_gid, 0, mesh.num_nodes, {
            //     node_vel(stage, node_gid, 0) = sin( PI * node_coords(1, node_gid, 0) ) * cos(PI * node_coords(1, node_gid, 1)); 
            //     node_vel(stage, node_gid, 1) =  -1.0*cos(PI * node_coords(1, node_gid, 0)) * sin(PI * node_coords(1, node_gid, 1)); 
            //     node_vel(stage, node_gid, 2) = 0.0;
            // });// end fill vel with TG exact
            // Kokkos::fence();
        
            node_vel.update_host();

            boundary_velocity(mesh, boundary, node_vel, time_value);

            node_vel.update_host();

            get_energy_rhs(mesh, ref_elem, DetJac, SigmaJacInv, node_vel, stress, node_coords, F_e, S, stage, viscosity_cond, source_cond);

            F_e.update_host();

            get_energy_residual(mesh, M_e, zone_sie, F_e, PSI, dt,  stage, time_int_weights);

            PSI.update_host();

            update_energy(zone_sie, PSI, zone_mass, mesh, stage);
            // FOR_ALL(elem_gid, 0, mesh.num_elems, {

            //     for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
                        
            //         int zone_gid = mesh.zones_in_elem(elem_gid, zone_lid);

                    
            //         double zone_coords[3]; 
            //         zone_coords[0] = 0.0;
            //         zone_coords[1] = 0.0;
            //         zone_coords[2] = 0.0;

            //         // get the coordinates of the zone center
            //         for (int node_lid = 0; node_lid < mesh.num_nodes_in_zone; node_lid++){
                        
            //             zone_coords[0] += node_coords(stage, mesh.nodes_in_zone(zone_gid, node_lid), 0);
            //             zone_coords[1] += node_coords(stage, mesh.nodes_in_zone(zone_gid, node_lid), 1);
            //             zone_coords[2] += node_coords(stage, mesh.nodes_in_zone(zone_gid, node_lid), 2);
                        
            //         } // end loop over nodes in element
            //         zone_coords[0] = zone_coords[0]/mesh.num_nodes_in_zone;
            //         zone_coords[1] = zone_coords[1]/mesh.num_nodes_in_zone;
            //         zone_coords[2] = zone_coords[2]/mesh.num_nodes_in_zone;

            //         // p = rho*ie*(gamma - 1)
            //         double gamma = elem_state_vars(elem_gid,0); // gamma value

            //         double temp_pres = 0.0;
            //         temp_pres = 0.25*( cos(2.0*PI*zone_coords[0]) + cos(2.0*PI*zone_coords[1]) ) + 1.0;
                    
            //         zone_sie(stage, zone_gid) = temp_pres/((gamma - 1.0));
                    
            //     }// end loop over zones
            // });// end fill sie with TG exact
            // Kokkos::fence();

            zone_sie.update_host();

            update_position_rdh(stage, dt, mesh, node_coords, node_vel);

            node_coords.update_host();

        }// end stages

    
        // increment the time
	    // time_value+=dt;
        double min_detJ = 0.0;
        double local_min_detJ = 1.0e+12;
        REDUCE_MIN(gauss_gid, 0, mat_pt.num_leg_pts, local_min_detJ, {
            local_min_detJ = local_min_detJ < DetJac(gauss_gid) ? local_min_detJ : DetJac(gauss_gid); 
        }, min_detJ);

            
        //     double lcl_dt_vol = 1.0e+12;
        //     double dt_vol = 0.0;
        //     REDUCE_MIN(gauss_gid, 0, mat_pt.num_leg_pts, lcl_dt_vol, {
        //         lcl_dt_vol = lcl_dt_vol > 0.1/(abs(div_vel(gauss_gid)) + fuzz) ? 0.1/(abs(div_vel(gauss_gid)) + fuzz) : lcl_dt_vol ;
        //     }, dt_vol);

        //     double dt_vol_old = dt;
        //     if (dt_vol < dt){
        //         dt = dt_vol;
        //     }
        
        //     // make dt be exact for final time
        //     dt = fmin(dt, time_final-time_value);
        // // increment the time
        //     if (dt < 0.80*dt_vol_old){
        //         time_value = time_value;
        //     cycle -= 1;
        //         printf("dt_vol decreased too much; restarting step. \n");
        //     }
        if (min_detJ < 0.0){
            // repeat step and decrease dt
            time_value = time_value;	
            printf("min detJ < 0 ; restarting step. \n");
            cycle -= 1;
            dt *= 0.85;
            // make dt be exact for final time
            dt = fmin(dt, time_final-time_value);

            if (dt < dt_min){
                printf("time step crashed : dt = %10.16lf \n", dt);
            }
        }
        else{
            time_value += dt;
        }
        
        size_t write = 0;
        if ((cycle+1)%graphics_cyc_ival == 0 && cycle>0){
            write = 1;
        }
        else if (cycle == cycle_stop) {
            write = 1;
        }
        else if (time_value >= time_final){
            write = 1;
        }
        else if (time_value >= graphics_time){
            write = 1;
        }
            
        // write outputs
        if (write == 1){
            printf("Writing outputs to file at %f \n", graphics_time);
            VTKHexN(mesh,
                node_coords,
                node_vel,
                node_mass,
                den,
                pres,
                stress,
                sspd,
                zone_sie,
                elem_vol,
                elem_mat_id,
                graphics_times,
                graphics_id,
                time_value);
            graphics_time = time_value + graphics_dt_ival;
        } // end if
        
        
        // end of calculation
        if (time_value>=time_final) break;
    }// end cycle
    
    auto time_2 = std::chrono::high_resolution_clock::now();
    auto calc_time = std::chrono::duration_cast
                           <std::chrono::nanoseconds>(time_2 - time_1).count();
    
    printf("\nCalculation time in seconds: %f \n", calc_time * 1e-9);
    
    return;
}