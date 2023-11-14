
#include "mesh.h"
#include "state.h"
#include "FEA_Module_Dynamic_Elasticity.h"

// -----------------------------------------------------------------------------
// This function saves the variables at rk_stage = 0, which is t_n
//------------------------------------------------------------------------------
void FEA_Module_Dynamic_Elasticity::rk_init(DViewCArrayKokkos <double> &node_coords,
             DViewCArrayKokkos <double> &node_vel,
             DViewCArrayKokkos <double> &elem_sie,
             DViewCArrayKokkos <double> &elem_stress,
             const size_t num_elems,
             const size_t num_nodes){

    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;
    int num_dims = simparam->num_dims;
    // save elem quantities
    FOR_ALL_CLASS(elem_gid, 0, num_elems, {

        // stress is always 3D even with 2D-RZ
        for(size_t i=0; i<3; i++){
            for(size_t j=0; j<3; j++){
                elem_stress(0,elem_gid,i,j) = elem_stress(rk_level,elem_gid,i,j);
            }
        }  // end for

        elem_sie(0,elem_gid) = elem_sie(rk_level,elem_gid);

    }); // end parallel for
    
    
    // save nodal quantities
    FOR_ALL_CLASS(node_gid, 0, num_nodes, {
        
        for(size_t i=0; i<num_dims; i++){
            node_coords(0,node_gid,i) = node_coords(rk_level,node_gid,i);
            node_vel(0,node_gid,i) = node_vel(rk_level,node_gid,i);
        }
    }); // end parallel for
    Kokkos::fence();
    
    return;
    
} // end rk_init




// -----------------------------------------------------------------------------
// This function calculates the time step by finding the shortest distance
// between any two nodes in the mesh
//------------------------------------------------------------------------------
// WARNING WARNING :  Only works for 3D, 8 node elements
void FEA_Module_Dynamic_Elasticity::get_timestep(mesh_t &mesh,
                  DViewCArrayKokkos <double> &node_coords,
                  DViewCArrayKokkos <double> &node_vel,
                  DViewCArrayKokkos <double> &elem_sspd,
                  DViewCArrayKokkos <double> &elem_vol){

    
    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;

    // increase dt by 10%, that is the largest dt value
    dt = dt*1.1;
    int num_dims = simparam->num_dims;
    double dt_lcl;
    double min_dt_calc;
    REDUCE_MIN_CLASS(elem_gid, 0, rnum_elem, dt_lcl, {
        
        double coords0[24];  // element coords
        ViewCArrayKokkos <double> coords(coords0, 8, 3);

        
        double distance0[28];  // array for holding distances between each node
        ViewCArrayKokkos <double> dist(distance0, 28);
        
        // Getting the coordinates of the element
        for(size_t node_lid = 0; node_lid < 8; node_lid++){

            for (size_t dim = 0; dim < num_dims; dim++){
                coords(node_lid, dim) = node_coords(rk_level,  nodes_in_elem(elem_gid, node_lid), dim);
            } // end for dim
            
        } // end for loop over node_lid

        // loop conditions needed for distance calculation
        size_t countA = 0;
        size_t countB = 1;
        size_t a;
        size_t b;
        size_t loop = 0;
        
        // Only works for 3D
        // Solving for the magnitude of distance between each node
        for (size_t i = 0; i < 28; i++){
            
            a = countA;
            b = countB;
            
            // returns magnitude of distance between each node, 28 total options
            dist(i) = fabs(sqrt(( pow((coords(b, 0) - coords(a, 0)), 2.0)
                                + pow((coords(b, 1) - coords(a, 1)), 2.0)
                                + pow((coords(b, 2) - coords(a, 2)), 2.0))));

            countB++;
            countA++;
            
            // tricky indexing
            if (countB > 7) {
                loop++;
                countB = 1 + loop;
                countA = 0;
            }
            
        } // endo for i


        double dist_min = dist(0);
        
        for(int i = 0; i < 28; ++i){
            dist_min = fmin(dist(i), dist_min);
        }
        
        // local dt calc based on CFL
        double dt_lcl_ = dt_cfl*dist_min/(elem_sspd(elem_gid) + fuzz);
        
        // make dt be in bounds
        dt_lcl_ = fmin(dt_lcl_, dt_max);               // make dt small than dt_max
        dt_lcl_ = fmax(dt_lcl_, dt_min);               // make dt larger than dt_min
        
       
        if (dt_lcl_ < dt_lcl) dt_lcl = dt_lcl_;
        
    }, min_dt_calc);  // end parallel reduction
    Kokkos::fence();
    
    // save the min dt
    if(min_dt_calc < dt) dt = min_dt_calc;
    
    // ensure time step hits the graphics time intervals
    dt = fmin(dt, (graphics_time - time_value)+fuzz);
    
    // make dt be exact for final time
    dt = fmin(dt, time_final-time_value);
    
    return;
    
} // end get_timestep

// -----------------------------------------------------------------------------
// This function calculates the time step by finding the shortest distance
// between any two nodes in the mesh
//------------------------------------------------------------------------------
// WARNING WARNING :  Only works for 3D, 8 node elements
void FEA_Module_Dynamic_Elasticity::get_timestep2D(mesh_t &mesh,
                    DViewCArrayKokkos <double> &node_coords,
                    DViewCArrayKokkos <double> &node_vel,
                    DViewCArrayKokkos <double> &elem_sspd,
                    DViewCArrayKokkos <double> &elem_vol){

    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;

    // increase dt by 10%, that is the largest dt value
    dt = dt*1.1;
    int num_dims = simparam->num_dims;
    double dt_lcl;
    double min_dt_calc;

    REDUCE_MIN_CLASS(elem_gid, 0, rnum_elem, dt_lcl, {
        
        double coords0[8];  // element coords
        ViewCArrayKokkos <double> coords(coords0, 4, 2);

        
        double distance0[6];  // array for holding distances between each node
        ViewCArrayKokkos <double> dist(distance0, 6);
        
        // Getting the coordinates of the nodes of the element
        for(size_t node_lid = 0; node_lid < 4; node_lid++){

            for (size_t dim = 0; dim < num_dims; dim++){
                coords(node_lid, dim) = node_coords(rk_level,  nodes_in_elem(elem_gid, node_lid), dim);
            } // end for dim
            
        } // end for loop over node_lid

        
        // Only works for 2D
        // Solving for the magnitude of distance between each node
        size_t count = 0;
        for (size_t i = 0; i < 3; i++){
            for (size_t j=i+1; j<=3; j++){
                
                // returns magnitude of distance between each node, 6 total options
                dist(count) = fabs(
                                 sqrt( pow((coords(i, 0) - coords(j, 0)), 2.0)
                                     + pow((coords(i, 1) - coords(j, 1)), 2.0) )
                                  );
                count ++;
            } // end for j
        } // end for i


        double dist_min = dist(0);
        
        for(int i = 0; i < 6; ++i){
            dist_min = fmin(dist(i), dist_min);
        }
        
        // local dt calc based on CFL
        double dt_lcl_ = dt_cfl*dist_min/(elem_sspd(elem_gid) + fuzz);
        
        // make dt be in bounds
        dt_lcl_ = fmin(dt_lcl_, dt_max);    // make dt small than dt_max
        dt_lcl_ = fmax(dt_lcl_, dt_min);    // make dt larger than dt_min
        
       
        if (dt_lcl_ < dt_lcl) dt_lcl = dt_lcl_;
        
    }, min_dt_calc);  // end parallel reduction
    Kokkos::fence();
    
    // save the min dt
    if(min_dt_calc < dt) dt = min_dt_calc;
    
    // ensure time step hits the graphics time intervals
    dt = fmin(dt, (graphics_time - time_value)+fuzz);
    
    // make dt be exact for final time
    dt = fmin(dt, time_final-time_value);
    
    return;
    
} // end get_timestep2D


