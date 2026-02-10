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

void limit_specific_volume(swage::mesh_t& mesh, elements::ref_element& ref_elem, std::string lim_type, int elem_gid) {

    // the element average for elem_gid
    real_t avg_spec_vol = elem_state.avg_specific_volume(elem_gid);

    // 1. find max and min bounds using elem_gid and neighboring elements
    real_t max_avg_spec_vol = avg_spec_vol;
    real_t min_avg_spec_vol = avg_spec_vol;
    
    for (int neigh_elem = 0; neigh_elem < mesh.num_elems_in_elem(elem_gid); neigh_elem++) {
        
        //a. get neighboor gid
        int neighboor_gid = mesh.elems_in_elem(elem_gid, neigh_elem);
        
        //b. check if avgs are new max/mins
        real_t avg_spec_vol_neigh = elem_state.avg_specific_volume(neighboor_gid);
        
        if ( avg_spec_vol_neigh > max_avg_spec_vol) {
            max_avg_spec_vol = avg_spec_vol_neigh;
        }
        if ( avg_spec_vol_neigh < min_avg_spec_vol) {
            min_avg_spec_vol = avg_spec_vol_neigh;
        }
        
    } //end loop over neighboors to find max/min avg specific_volume


    //2. Find the scalar limiting value alpha using "type" limiter
    real_t epsilon = 1e-10;
    real_t elem_alpha_spec_vol = 1.0;  // element limiter coefficient is the largest one in the element

    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++) {

        int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);

        real_t diff = mat_pt.specific_volume(1, gauss_gid) - avg_spec_vol;
        real_t ratio = 0.0;
        real_t alpha_spec_vol = 0.0; // limiter value at a guass point

        if ( diff > epsilon) {
            ratio = (max_avg_spec_vol - avg_spec_vol) / diff;
        }
        else if (diff < -epsilon) {
            ratio = (min_avg_spec_vol - avg_spec_vol) / diff;
        }
        else {
            ratio = 1.0;
        }
        
        ratio = std::max(ratio, 0.0);
        if (lim_type == "BJ") {
            alpha_spec_vol = std::min(ratio, 1.0);
        }
        else if ( lim_type == "V") {
            alpha_spec_vol = (ratio*ratio + 2.0*ratio) / (ratio*ratio + ratio + 2.0);
        }
        else if ( lim_type == "P0") {
            alpha_spec_vol = 0.0;
        }
        else{
            alpha_spec_vol = 1.0; // no limiting
        }
        elem_alpha_spec_vol = std::min(elem_alpha_spec_vol, alpha_spec_vol);

    } //end finding alpha for limiting


    // 3. Limit specific_volume at every gaussian point in element
    // specific_volume_limited = avg + elem_alpha*(specific_volume - avg)
    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++) {

        int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
        
        real_t curent_specific_volume = mat_pt.specific_volume(1, gauss_gid); // unlimited specific_volume
        mat_pt.specific_volume(1, gauss_gid) = avg_spec_vol + elem_alpha_spec_vol*(curent_specific_volume - avg_spec_vol); // limited specific_volume
        
        mat_pt.den_phi(gauss_gid) = elem_alpha_spec_vol;
        
        
        // update the density
        mat_pt.density(gauss_gid) = 1.0/mat_pt.specific_volume(1, gauss_gid);
        
    } //end loop over gauss in element



} // end limiting function
