// Call cycle loop for the CCH P0 solver


#include "utilities.h"
#include "state.h"
#include "geometry/geometry.h"
#include "variables.h"

using namespace utils;

void cch_hydro(){
	
	// loop over the max number of time integration cycles
	for (cycle = 1; cycle <= cycle_stop; cycle++) {

	    // reset rk_stage to zero
	    //rk_stage = 0;
	    
	    // stop calculation if flag
	    if (stop_calc == 1) break;

	    // calculate the total energy
	    track_cch(ke, ie);

	    // get the step
	    get_timestep();
	    
	    // ensure time step hits the graphics time intervals
	    dt = fmin(dt, (graphics_time - TIME)+fuzz);
	    
	    // integrate the solution forward in time via RK
	    { // scope inside rk integration
	     
	        // save the values at t_n
	        rk_init();

	        
	        // integrate solution forward in time
	        for (int rk_stage = 0; rk_stage < rk_num_stages; rk_stage++){

	            // k_init();

	            // ---- RK coefficient ---- //
	            real_t rk_alpha = 1.0/((real_t)rk_num_stages - (real_t)rk_stage);



	            // ---- calculate the forces on the vertices ---- //
	            get_force_cch();
	            
	            // ---- RK update (Velocity, position, internal energy, stress) ---- //
	            
	            // ---- Update nodal velocities ---- //
	            // ---- Update specific internal energy in the elements ---- //
	            update_mom_energy_cch(rk_alpha);
	            
	            // ---- Update nodal positions ---- //
	            // std::cout<<"Before position"<<std::endl;
	            update_position_cch(rk_alpha);
	            
	            // ---- Calculate cell volume for next time step ---- //
	            // std::cout<<"Before volume"<<std::endl;

	            // Get jacobians at the gauss points
	            get_gauss_pt_jacobian(mesh, ref_elem);


	            // Get exact volume for a linear hexahedron
	            get_vol_hex(mesh, ref_elem);
	            
	            // ---- Calculate cell state (den,pres,sound speed) for next time step ---- //
	            // std::cout<<"Before state"<<std::endl;
	            update_state_cch();
	           	            
	            
	        } // end of RK loop
	        
	    } // end of scope for time integration
	    
	    

	    track_cch(ke, ie);



	    // increment the time
	    TIME+=dt;

	   	// end of calculation
    	if (TIME>=TFINAL) break;
    	
	    run_info(cycle);

	    
	    
	} // end of loop



	


	
}