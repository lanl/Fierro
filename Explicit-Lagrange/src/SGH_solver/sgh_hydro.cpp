// Call cycle loop for the SGH solver


#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"
using namespace utils;


void sgh_hydro(){
	// loop over the max number of time integration cycles
	for (cycle = 1; cycle <= cycle_stop; cycle++) {

	    // std::cout<<"********CYCLE =  "<< cycle<<std::endl;
	    
	    // std::cout<<"Cycle = "<< cycle <<std::endl;

	    // reset rk_stage to zero
	    //rk_stage = 0;
	    
	    // stop calculation if flag
	    if (stop_calc == 1) break;


	    // calculate the total energy
	    if(CCH == true) track_cch(ke, ie);
	    if(SGH == true) track_sgh(ke, ie);

	    if (cycle == 1) te_0 = ie + ke;

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


	            // Calculate cell velocity gradient
	            // std::cout<<"Before velgrad"<<std::endl;
	            get_velgrad();

	            // ---- calculate the forces on the vertices ---- //
	            // std::cout<<"Before force"<<std::endl;
	            get_force_sgh();

	            // ---- apply force boundary conditions to the boundary patches---- //
	            // std::cout<<"Before boundary"<<std::endl;
	            boundary_force();
	            
	            boundary_velocity();

	            // ---- RK update (Velocity, position, internal energy, stress) ---- //
	            
	            // ---- Update nodal velocities ---- //
	            update_velocity_sgh(rk_alpha);

	            // ---- Update specific internal energy in the elements ---- //
	            // std::cout<<"Before energy"<<std::endl;
	            update_energy_sgh(rk_alpha);
	            
	            // ---- Update nodal positions ---- //
	            // std::cout<<"Before position"<<std::endl;
	            update_position_sgh(rk_alpha);
	            
	            // ---- Calculate cell volume for next time step ---- //
	            // std::cout<<"Before volume"<<std::endl;
	            get_vol_hex(mesh, ref_elem);

	            // ---- Calculate cell B matrix for next time step ---- //
	            // std::cout<<"Before b matrix"<<std::endl;
	            get_bmatrix();

	            // ---- Calculate cell state (den,pres,sound speed) for next time step ---- //
	            // std::cout<<"Before state"<<std::endl;
	            update_state_sgh();

	        } // end of RK loop
	        
	    } // end of scope for time integration
	    
	    track_sgh(ke, ie);
	    

	    // increment the time
	    TIME+=dt;

	   	// end of calculation
    	if (TIME>=TFINAL) break;
    	
	    run_info(cycle);
	    
	    
	} // end of loop
}