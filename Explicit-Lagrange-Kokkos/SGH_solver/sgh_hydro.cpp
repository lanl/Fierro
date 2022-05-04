// Call cycle loop for the SGH solver


#include "state.h"
#include "mesh.h"


void sgh_hydro(CArrayKokkos <material_t> &material,
               CArrayKokkos <mat_fill_t> &mat_fill,
               CArrayKokkos <boundary_t> &boundary,
               mesh_t &mesh,
               DViewCArrayKokkos <double> &node_coords,
               DViewCArrayKokkos <double> &node_vel,
               DViewCArrayKokkos <double> &node_mass,
               const DViewCArrayKokkos <double> &elem_den,
               const DViewCArrayKokkos <double> &elem_pres,
               const DViewCArrayKokkos <double> &elem_stress,
               const DViewCArrayKokkos <double> &elem_sspd,
               const DViewCArrayKokkos <double> &elem_sie,
               const DViewCArrayKokkos <double> &elem_vol,
               const DViewCArrayKokkos <double> &elem_mass,
               const DViewCArrayKokkos <size_t> &elem_mat_id,
               const DViewCArrayKokkos <double> &elem_statev,
               const CArrayKokkos <double> &state_vars,
               double &time_value,
               const double time_final,
               const double dt_max,
               const double dt_min,
               const double dt_cfl,
               size_t &graphics_time,
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
    
    
    
	// loop over the max number of time integration cycles
	for (size_t cycle = 1; cycle <= cycle_stop; cycle++) {

	    // stop calculation if flag
	    if (stop_calc == 1) break;

	    // get the step
	    //get_timestep();
	    
	    // ensure time step hits the graphics time intervals
	    dt = fmin(dt, (graphics_time - time_value)+fuzz);
	    
	    // integrate the solution forward in time via Runge Kutta (RK) method
	    { // scope inside rk integration
	     
	        // save the values at t_n
	        //rk_init();

	        
	        // integrate solution forward in time
	        for (int rk_stage = 0; rk_stage < rk_num_stages; rk_stage++){

	            // k_init();

	            // ---- RK coefficient ---- //
	            double rk_alpha = 1.0/((double)rk_num_stages - (double)rk_stage);


	            // Calculate cell velocity gradient
	            // std::cout<<"Before velgrad"<<std::endl;
	            //get_velgrad();

	            // ---- calculate the forces on the vertices ---- //
	            //get_force_sgh();

	            // ---- apply force boundary conditions to the boundary patches---- //
	            //boundary_force();
	            
                boundary_velocity(mesh, boundary, node_vel);

	            // ---- RK update (Velocity, position, internal energy, stress) ---- //
	            
	            // ---- Update nodal velocities ---- //
	            //update_velocity_sgh(rk_alpha);

	            // ---- Update specific internal energy in the elements ---- //
	            //update_energy_sgh(rk_alpha);
	            
	            // ---- Update nodal positions ---- //
	            //update_position_sgh(rk_alpha);
	            
	            // ---- Calculate cell volume for next time step ---- //
	            //get_vol_hex(mesh, ref_elem);

	            // ---- Calculate cell B matrix for next time step ---- //
	            //get_bmatrix();

	            // ---- Calculate cell state (den,pres,sound speed) for next time step ---- //
	            //update_state_sgh();

	        } // end of RK loop
	        
	    } // end of scope for time integration
	    

	    // increment the time
	    time_value+=dt;

	   	// end of calculation
    	if (time_value>=time_final) break;
    	
	    //run_info(cycle);
	    
	    
	} // end of loop
}
