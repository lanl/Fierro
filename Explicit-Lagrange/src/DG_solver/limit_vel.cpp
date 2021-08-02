#include <math.h>
#include <iostream>		//cout
#include <cmath>		//abs
#include <algorithm>	//min
#include <string>		//string -> char for file output

#include "state.h"
#include "geometry.h"
#include "variables.h"
#include "utilities.h"

using namespace utils;

//limit fields function:
//Inputs:
//	1. swage mesh
//	2. ref element
//	3. string: type of limiting function
//		- "BJ" = Barth-Jesperson = min(1, r)
//		- "V" = Venkatakrishan = (x^2 + 2x) / ( x^2 + x + 2)
//	4. int elem_gid: element global id, needed since the limiting is done 
//	   element by element
//	NOTE 1: whats the best way to switch between limiting? I don't want an if statement 
//	at every gauss point, its very expensive.
//
//	NOTE 2: should add a different parameter for what to limit, ie ke
//	or each component of the velocity
//
//

real_t epsilon = 1e-10;



void limit_vel(swage::mesh_t& mesh, elements::ref_element& ref_elem, std::string lim_type, int elem_gid) {

    // the element average for elem_gid
    real_t avg_vel[mesh.num_dim()];
    for(int dim = 0; dim < mesh.num_dim(); dim++){
        avg_vel[dim] = elem_state.avg_velocity(elem_gid, dim);
        //std::cout << "limiting routine, " << avg << std::endl;
    }
    
    // kinetic energy using averages
    real_t ke_avg = 0.0;
    real_t ke_max = 0.0;
    real_t ke_min = 0.0;
    

	// 1. find max and min bounds using elem_gid and neighboring elements
    real_t max_avg_vel[mesh.num_dim()];
    real_t min_avg_vel[mesh.num_dim()];
    
    for(int dim = 0; dim < mesh.num_dim(); dim++){
	    max_avg_vel[dim] = avg_vel[dim];
	    min_avg_vel[dim] = avg_vel[dim];
        
        ke_avg += avg_vel[dim]*avg_vel[dim];
        ke_max += avg_vel[dim]*avg_vel[dim];
        ke_min += avg_vel[dim]*avg_vel[dim];
    }
    ke_avg = ke_avg/2.0;
    ke_max = ke_max/2.0;
    ke_min = ke_min/2.0;
    
    
    for (int neigh_elem = 0; neigh_elem < mesh.num_elems_in_elem(elem_gid); neigh_elem++) {
        
        //a. get neighboor gid
        int neighboor_gid = mesh.elems_in_elem(elem_gid, neigh_elem);
        
        //b. check if avgs are new max/mins
        for(int dim = 0; dim < mesh.num_dim(); dim++){
            real_t avg_vel_neigh = elem_state.avg_velocity(neighboor_gid, dim);

            if ( avg_vel_neigh > max_avg_vel[dim]) {
                max_avg_vel[dim] = avg_vel_neigh;
            }
            if ( avg_vel_neigh < min_avg_vel[dim]) {
                min_avg_vel[dim] = avg_vel_neigh;
            }
            
        } // end dim
        
        // check if kinetic energy avgs are new max/mins
        real_t avg_ke_neigh = 0.0;
        for(int dim = 0; dim < mesh.num_dim(); dim++){
            avg_ke_neigh += elem_state.avg_velocity(neighboor_gid, dim)*elem_state.avg_velocity(neighboor_gid, dim);
        }
        avg_ke_neigh = avg_ke_neigh/2.0;
        
        if ( avg_ke_neigh > ke_max) {
            ke_max = avg_ke_neigh;
        }
        if ( avg_ke_neigh < ke_min) {
            ke_min = avg_ke_neigh;
        }
        
    } //end loop over neighboors to find max/min avg values

    ke_avg = std::max(epsilon, ke_avg);
    ke_max = std::max(epsilon, ke_max);
    ke_min = std::max(epsilon, ke_min);



    

	// 2. Now compute the limiting scalar value using "type" limiter
	// "BJ" = Barth-Jesperson
	// "V" = Venkatrashina
	real_t elem_alpha_vel = 1.0;  // element limiting coefficient

    
    // for shock detector
    real_t char_length = pow(mesh.elem_vol(elem_gid)/mesh.num_cells_in_elem(), 0.333333);
    real_t alpha_shock = 1.0;  // 1=no shock, 0=definately a shock

    
	for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++) {

        int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
        
        real_t diff[mesh.num_dim()];
        real_t alpha_vel[mesh.num_dim()];
        real_t alpha_ke;
        
        // limit velocity
        int lim_approach = 1;  // 0 = vel, 1 = ke
        
        
        // limit based on the velocity
        if (lim_approach == 0){
            
            // a. initialize alpha and calculate diff
            for(int dim = 0; dim < mesh.num_dim(); dim++) {
                diff[dim] = mat_pt.velocity(1, gauss_gid, dim) - avg_vel[dim];
            }
            
            
            // b. find ratio and alpha for each component of velocity
            for(int dim = 0; dim < mesh.num_dim(); dim++) {
                real_t ratio = 0.0;
                
                if (diff[dim] > epsilon) {
                    ratio = (max_avg_vel[dim] - avg_vel[dim]) / diff[dim];
                }
                else if ( diff[dim] < -epsilon) {
                    ratio = (min_avg_vel[dim] - avg_vel[dim]) / diff[dim];
                }
                else {
                    ratio = 1.0;
                }
                
                ratio = std::max(ratio, 0.0);
                if (lim_type == "BJ") {
                    alpha_vel[dim] = std::min(ratio, 1.0);
                }
                else if ( lim_type == "V") {
                    alpha_vel[dim] = (ratio*ratio + 2.0*ratio) / (ratio*ratio + ratio + 2.0);
                }
                else if ( lim_type == "P0") {
                    alpha_vel[dim] = 0.0;
                }
                else{
                    alpha_vel[dim] = 1.0; // no limiting at all
                }
                
            } //end finding alpha for each velocity component
            
            
            // c. find the smallest limiting value for all velocity components and use the smallest for the element
            real_t smallest = alpha_vel[0];
            for(int dim = 1; dim < mesh.num_dim(); dim++) {
                smallest = std::min(smallest, alpha_vel[dim]);
            }
            elem_alpha_vel = std::min(elem_alpha_vel, smallest);
            
        } // end of limiting approach based on velocity
        else{
            
            real_t diff_ke = 0.0;
            real_t ke = 0.0;  // the value at a quadrature point
            
            // a. calculate ke at this gauss point and then the diff
            for(int dim = 0; dim < mesh.num_dim(); dim++) {
                ke += mat_pt.velocity(1, gauss_gid, dim)*mat_pt.velocity(1, gauss_gid, dim);
            }
            ke = ke/2.0;
            

            diff_ke = ke - ke_avg;  // the difference
            
            
            // b. find ratio and alpha for ke
            real_t ratio = 0.0;
                
            if (diff_ke > epsilon) {
                ratio = (ke_max - ke_avg)/diff_ke;
            }
            else if (diff_ke < -epsilon) {
                ratio = (ke_min - ke_avg)/diff_ke;
            }
            else {
                ratio = 1.0;
            }
            

            
            ratio = std::max(ratio, 0.0);
            if (lim_type == "BJ") {
                alpha_ke = std::min(ratio, 1.0);
            }
            else if ( lim_type == "V") {
                alpha_ke = (ratio*ratio + 2.0*ratio) / (ratio*ratio + ratio + 2.0);
            }
            else if ( lim_type == "P0") {
                alpha_ke = 0.0;
            }
            else{
                alpha_ke = 1.0; // no limiting at all
            }
            
            elem_alpha_vel = std::min(elem_alpha_vel, alpha_ke);
            
        }  // end if on limiting using vel versus ke
            
        
        elem_alpha_vel = std::min(mat_pt.te_phi(gauss_gid), elem_alpha_vel);  // limiter must be smaller than energy limiter
        
        // calculate shock detector
        real_t ssp = std::max(mat_pt.sspd(gauss_gid), 1.0E-14);
        real_t ratio = ssp/( fabs(char_length*mat_pt.div_vel(gauss_gid)) + 1.0E-16 );
        ratio = std::min(1.0, ratio/200.0);
        alpha_shock = std::min(alpha_shock, ratio);
        
        // expansion
        //if (mat_pt.div_vel(gauss_gid) >= 0.0 && 2.0*char_length*mat_pt.div_vel(gauss_gid) <= ssp){
        //    alpha_shock = 1.0;  // don't limit if in expansion and if the exapansion is less than the speed of sound
        //}
        
    } //end loop over gauss points in element to find elem_alpha
    
    // d. use shock detector
    //elem_alpha_vel = std::max(elem_alpha_vel, alpha_shock);

    
    // verifying that the internal energy is positive
    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++) {
        
        int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
        
        real_t ke = 0;
        
        // a. calculate the ke
        for(int dim = 0; dim < mesh.num_dim(); dim++) {
            ke += mat_pt.velocity(1, gauss_gid, dim)*mat_pt.velocity(1, gauss_gid, dim);
        }
        ke = ke/2.0;
        
        
        real_t alpha_ke = std::max(0.0, mat_pt.specific_total_energy(1, gauss_gid)/(ke+1.0e-14));
        alpha_ke = std::min( 1.0, alpha_ke);
        
        elem_alpha_vel = std::min(elem_alpha_vel, alpha_ke);
    }
    
    
    // 3. Limit velocity at every gaussian point in element
    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++) {
    
        int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
    
        for(int dim = 0; dim < mesh.num_dim(); dim++) {
            real_t current_vel = mat_pt.velocity(1, gauss_gid, dim);
            mat_pt.velocity(1, gauss_gid, dim) = avg_vel[dim] + elem_alpha_vel*(current_vel - avg_vel[dim]);
        }
        
        mat_pt.vel_phi(gauss_gid) = elem_alpha_vel;
        
    
    } //end loop over gauss in element


} // end limiting function



