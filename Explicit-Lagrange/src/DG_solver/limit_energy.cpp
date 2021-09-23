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
//	NOTE: FOR NOW ONLY BJ IS IMPLEMENTED, SO THE STRING VARIABLE IS USELESS

void limit_energy(swage::mesh_t& mesh, elements::ref_element& ref_elem, std::string lim_type, int elem_gid) {

    // the element average for elem_gid
    real_t avg_te = elem_state.avg_specific_total_energy(elem_gid);

    // 1. find max and min bounds using elem_gid and neighboring elements
    real_t max_avg_te = avg_te;
    real_t min_avg_te = avg_te;
    
    for (int neigh_elem = 0; neigh_elem < mesh.num_elems_in_elem(elem_gid); neigh_elem++) {
        
        //a. get neighboor gid
        int neighboor_gid = mesh.elems_in_elem(elem_gid, neigh_elem);
        
        //b. check if avgs are new max/mins
        real_t avg_te_neigh = elem_state.avg_specific_total_energy(neighboor_gid);
        
        if ( avg_te_neigh > max_avg_te) {
            max_avg_te = avg_te_neigh;
        }
        if ( avg_te_neigh < min_avg_te) {
            min_avg_te = avg_te_neigh;
        }
        
    } //end loop over neighboors to find max/min avg energy


    //2. Find the scalar limiting value alpha using "type" limiter
    real_t epsilon = 1e-10;
    real_t elem_alpha_te = 1.0;  // element limiter coefficient is the largest one in the element

    // for shock detector
    real_t char_length = pow(mesh.elem_vol(elem_gid)/mesh.num_cells_in_elem(), 0.333333);
    real_t alpha_shock = 1.0;  // 1=no shock, 0=definately a shock
    
    
    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++) {

        int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);

        // a. initialize alpha and calculate diff
        real_t diff  = mat_pt.specific_total_energy(1, gauss_gid) - avg_te;
        real_t ratio = 0.0;
        real_t alpha_te = 1.0; // limiter value at a guass point

        // b. find ratio and alpha
        if ( diff > epsilon) {
            ratio = (max_avg_te - avg_te) / diff;
        }
        else if (diff < -epsilon) {
            ratio = (min_avg_te - avg_te) / diff;
        }
        else {
            ratio = 1.0;
        }
        
        ratio = std::max(ratio, 0.0);
        if (lim_type == "BJ") {
            alpha_te = std::min(ratio, 1.0);
        }
        else if ( lim_type == "V") {
            alpha_te = (ratio*ratio + 2.0*ratio) / (ratio*ratio + ratio + 2.0);
        }
        else if ( lim_type == "P0") {
            alpha_te = 0.0;
        }
        else{
            alpha_te = 1.0; // no limiting
        }
        
        // c. find the smallest limiting value for the entire element
        elem_alpha_te = std::min(elem_alpha_te, alpha_te);
        
        
        // calculate shock detector
        real_t ssp = std::max(mat_pt.sspd(gauss_gid), 1.0E-14);
        ratio = ssp/( fabs(char_length*mat_pt.div_vel(gauss_gid)) + 1.0E-16 );
        ratio = std::min(1.0, ratio/20.0);
        alpha_shock = std::min(alpha_shock, ratio);
        
        // expansion
        if (mat_pt.div_vel(gauss_gid) >= 0.0){ // && 2.0*char_length*mat_pt.div_vel(gauss_gid) <= ssp){
            alpha_shock = std::max(0.0, std::min( 1.0, 1.5-char_length*mat_pt.div_vel(gauss_gid)/ssp ));  // don't limit if in expansion and if the exapansion is less than the speed of sound
        }

    } //end finding alpha for limiting

    
    // d. use shock detector
    //elem_alpha_te = std::max(elem_alpha_te, alpha_shock);
    
    
    // verifying that the internal energy is positive
    /*
    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++) {
        
        int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
        
        real_t ke = 0;
        
        // a. calculate the ke
        for(int dim = 0; dim < mesh.num_dim(); dim++) {
            ke += mat_pt.velocity(1, gauss_gid, dim)*mat_pt.velocity(1, gauss_gid, dim);
        }
        ke = ke/2.0;
        
        
        real_t alpha_ke = std::max(0.0, mat_pt.specific_total_energy(1, gauss_gid)/(ke+1.0e-14));
        alpha_ke = std::min( 1.0, alpha_ke );
        
        elem_alpha_te = std::min(elem_alpha_te, alpha_ke);
    }
     */

    // 3. Limit specific total energy at every gaussian point in element
    // energy_limited = avg + elem_alpha*(energy - avg)
    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++) {

        int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
        
        real_t curent_tot_energy = mat_pt.specific_total_energy(1, gauss_gid); // unlimited energy
        mat_pt.specific_total_energy(1, gauss_gid) = avg_te + elem_alpha_te*(curent_tot_energy - avg_te); // limited energy

        mat_pt.te_phi(gauss_gid) = elem_alpha_te;
        
    } //end loop over gauss in element

} // end limiting function
