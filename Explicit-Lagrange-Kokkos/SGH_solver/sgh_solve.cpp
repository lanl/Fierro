// Call cycle loop for the SGH solver


#include "state.h"
#include "mesh.h"


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
               DViewCArrayKokkos <double> &elem_mass,
               DViewCArrayKokkos <size_t> &elem_mat_id,
               DViewCArrayKokkos <double> &elem_statev,
               double &time_value,
               const double time_final,
               const double dt_max,
               const double dt_min,
               const double dt_cfl,
               double &graphics_time,
               size_t graphics_cyc_ival,
               size_t graphics_dt_ival,
               const size_t cycle_stop,
               const size_t rk_num_stages,
               double dt,
               const double fuzz,
               const double tiny,
               const double small
               ){
    
    // a flag to exit the calculation
    size_t stop_calc = 0;
    
    
    printf("cycle stop = %lu \n", cycle_stop);
    
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
            printf("time step = %f \n", dt);
        }
        // print time step every 10 cycles
        else if (cycle%10==0){
            printf("time step = %f \n", dt);
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

            // ---- RK coefficient ---- //
            double rk_alpha = 1.0/((double)rk_num_stages - (double)rk_stage);
            
            
            // Calculate cell velocity gradient
            //get_velgrad();
            
            // ---- calculate the forces on the vertices ---- //
            //get_force_sgh();
            
            // ---- apply force boundary conditions to the boundary patches---- //
            //boundary_force();
            boundary_velocity(mesh, boundary, node_vel);
            
            
            // ---- Update nodal velocities ---- //
            //update_velocity_sgh(rk_alpha);
            
            // ---- Update specific internal energy in the elements ---- //
            //update_energy_sgh(rk_alpha);
            
            // ---- Update nodal positions ---- //
            //update_position_sgh(rk_alpha);
            
            // ---- Calculate cell volume for next time step ---- //
            //get_vol_hex(mesh, ref_elem);
            
            // ---- Calculate cell state (den,pres,sound speed) for next time step ---- //
            //update_state_sgh();
            
        } // end of RK loop
	        

	    

	    // increment the time
	    time_value+=dt;

	   	// end of calculation
    	if (time_value>=time_final) break;
    	
	    //run_info(cycle);
	    
	} // end for cycle loop
}
