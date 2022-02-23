#include <math.h>
#include <iostream>		//cout
#include <cmath>		//abs
#include <algorithm>	//min
#include <string>		//string -> char for file output

#include "state.h"
#include "geometry/geometry.h"
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

void limit_density(swage::mesh_t& mesh, elements::ref_element& ref_elem, std::string lim_type, int elem_gid) {

    // the element average for elem_gid
    real_t avg_den = elem_state.avg_density(elem_gid);

    // 1. find max and min bounds using elem_gid and neighboring elements
    real_t max_avg_den = avg_den;
    real_t min_avg_den = avg_den;
    
    for (int neigh_elem = 0; neigh_elem < mesh.num_elems_in_elem(elem_gid); neigh_elem++) {
        
        //a. get neighboor gid
        int neighboor_gid = mesh.elems_in_elem(elem_gid, neigh_elem);
        
        //b. check if avgs are new max/mins
        real_t avg_den_neigh = elem_state.avg_density(neighboor_gid);
        
        if ( avg_den_neigh > max_avg_den) {
            max_avg_den = avg_den_neigh;
        }
        if ( avg_den_neigh < min_avg_den) {
            min_avg_den = avg_den_neigh;
        }
        
    } //end loop over neighboors to find max/min avg density


    //2. Find the scalar limiting value alpha using "type" limiter
    real_t epsilon = 1e-10;
    real_t elem_alpha_den = 1.0;  // element limiter coefficient is the largest one in the element

    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++) {

        int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);

        real_t diff = mat_pt.density(gauss_gid) - avg_den;
        real_t ratio = 0.0;
        real_t alpha_den = 1.0; // limiter value at a guass point

        if ( diff > epsilon) {
            ratio = (max_avg_den - avg_den) / diff;
        }
        else if (diff < -epsilon) {
            ratio = (min_avg_den - avg_den) / diff;
        }
        else {
            ratio = 1.0;
        }
        
        ratio = std::max(ratio, 0.0);
        if (lim_type == "BJ") {
            alpha_den = std::min(ratio, 1.0);
        }
        else if ( lim_type == "V") {
            alpha_den = (ratio*ratio + 2.0*ratio) / (ratio*ratio + ratio + 2.0);
        }
        else if ( lim_type == "P0") {
            alpha_den = 0.0;
        }
        else{
            alpha_den = 1.0; // no limiting
        }
        elem_alpha_den = std::min(elem_alpha_den, alpha_den);

    } //end finding alpha for limiting


    // 3. Limit density at every gaussian point in element
    // density_limited = avg + elem_alpha*(density - avg)
    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++) {

        int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
        
        real_t curent_density = mat_pt.density(gauss_gid); // unlimited density
        mat_pt.density(gauss_gid) = avg_den + elem_alpha_den*(curent_density - avg_den); // limited density
        
        
        mat_pt.den_phi(gauss_gid) = elem_alpha_den;

        // recalculate the mass at the material point
        //mat_pt.mass(gauss_gid) = mat_pt.density(gauss_gid)*ref_elem.ref_node_g_weights(gauss_lid)*mesh.gauss_pt_det_j(gauss_gid);

    } //end loop over gauss in element



} // end limiting function
