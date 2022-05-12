// Call cycle loop for the SGH solver


#include "state.h"
#include "mesh.h"
#include <chrono>

void sgh_solve(CArrayKokkos <material_t> &material,
               CArrayKokkos <boundary_t> &boundary,
               mesh_t &mesh,
               DViewCArrayKokkos <double> &node_coords,
               DViewCArrayKokkos <double> &node_vel,
               DViewCArrayKokkos <double> &node_mass,
               DViewCArrayKokkos <double> &elem_den,
               DViewCArrayKokkos <double> &elem_pres,
               DViewCArrayKokkos <double> &elem_stress,
               DViewCArrayKokkos <double> &elem_sspd,
               DViewCArrayKokkos <double> &elem_sie,
               DViewCArrayKokkos <double> &elem_vol,
               DViewCArrayKokkos <double> &elem_div,
               DViewCArrayKokkos <double> &elem_mass,
               DViewCArrayKokkos <size_t> &elem_mat_id,
               DViewCArrayKokkos <double> &elem_statev,
               DViewCArrayKokkos <double> corner_force,
               double &time_value,
               const double time_final,
               const double dt_max,
               const double dt_min,
               const double dt_cfl,
               double &graphics_time,
               size_t graphics_cyc_ival,
               double graphics_dt_ival,
               const size_t cycle_stop,
               const size_t rk_num_stages,
               double dt,
               const double fuzz,
               const double tiny,
               const double small,
               CArray <double> &graphics_times,
               size_t &graphics_id){
    
    
    printf("Writing outputs to file at %f \n", time_value);
    write_outputs(mesh,
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
    
    // a flag to exit the calculation
    size_t stop_calc=0;
    
    auto time_1 = std::chrono::high_resolution_clock::now();
    
	// loop over the max number of time integration cycles
	for (size_t cycle = 0; cycle < cycle_stop; cycle++) {

	    // stop calculation if flag
	    if (stop_calc == 1) break;

	    // get the step
	    get_timestep(mesh,
                     node_coords,
                     node_vel,
                     elem_sspd,
                     elem_vol,
                     time_value,
                     graphics_time,
                     time_final,
                     dt_max,
                     dt_min,
                     dt_cfl,
                     dt,
                     fuzz);
        
        if (cycle==0){
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        }
        // print time step every 10 cycles
        else if (cycle%10==0){
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        } // end if
        
        
        // ---------------------------------------------------------------------
        //  integrate the solution forward to t(n+1) via Runge Kutta (RK) method
        // ---------------------------------------------------------------------
	     
        // save the values at t_n
        rk_init(node_coords,
                node_vel,
                elem_sie,
                elem_stress,
                mesh.num_dims,
                mesh.num_elems,
                mesh.num_nodes);
	        
        // integrate solution forward in time
        for (size_t rk_stage = 0; rk_stage < rk_num_stages; rk_stage++){

            
            // ---- RK coefficient ----
            double rk_alpha = 1.0/((double)rk_num_stages - (double)rk_stage);
            
            
            // ---- Calculate velocity diveregence for the element ----
            get_divergence(elem_div,
                           mesh,
                           node_coords,
                           node_vel,
                           elem_vol);
            
            
            // ---- calculate the forces on the vertices and evolve stress (hypo model) ----
            get_force_sgh(material,
                          mesh,
                          node_coords,
                          node_vel,
                          elem_den,
                          elem_sie,
                          elem_pres,
                          elem_stress,
                          elem_sspd,
                          elem_vol,
                          elem_div,
                          elem_mat_id,
                          corner_force,
                          fuzz,
                          small,
                          elem_statev,
                          dt,
                          rk_alpha);

            
            // ---- Update nodal velocities ---- //
            update_velocity_sgh(rk_alpha,
                                dt,
                                mesh,
                                node_vel,
                                node_mass,
                                corner_force);
            
            // ---- apply force boundary conditions to the boundary patches----
            boundary_velocity(mesh, boundary, node_vel);
            
            
            // ---- Update specific internal energy in the elements ----
            update_energy_sgh(rk_alpha,
                              dt,
                              mesh,
                              node_vel,
                              elem_sie,
                              elem_mass,
                              corner_force);
            
            
            // ---- Update nodal positions ----
            update_position_sgh(rk_alpha,
                                dt,
                                mesh.num_dims,
                                mesh.num_nodes,
                                node_coords,
                                node_vel);
            
            
            // ---- Calculate cell volume for next time step ----
            get_vol(elem_vol, node_coords, mesh);
            
            
            // ---- Calculate elem state (den, pres, sound speed, stress) for next time step ----
            update_state(material,
                         mesh,
                         node_coords,
                         node_vel,
                         elem_den,
                         elem_pres,
                         elem_stress,
                         elem_sspd,
                         elem_sie,
                         elem_vol,
                         elem_mass,
                         elem_mat_id,
                         elem_statev,
                         dt,
                         rk_alpha);
            
            // ----
            // Notes on strength:
            //    1) hyper-elastic strength models are called in update_state
            //    2) hypo-elastic strength models are called in get_force
            //    3) strength models must be added by the user in user_mat.cpp

            
        } // end of RK loop
	        

	    

	    // increment the time
	    time_value+=dt;
    	
        
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
            write_outputs(mesh,
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
            
            graphics_time = time_value + graphics_dt_ival;
        } // end if
        
        
        // end of calculation
        if (time_value>=time_final) break;
        
    } // end for cycle loop
    
    
    auto time_2 = std::chrono::high_resolution_clock::now();
    auto calc_time = std::chrono::duration_cast
                           <std::chrono::nanoseconds>(time_2 - time_1).count();
    
    printf("\nCalculation time in seconds: %f \n", calc_time * 1e-9);
    
}
