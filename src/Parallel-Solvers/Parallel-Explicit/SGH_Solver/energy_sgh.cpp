
#include "mesh.h"
#include "state.h"
#include "FEA_Module_SGH.h"
#include "Simulation_Parameters_SGH.h"

void FEA_Module_SGH::update_energy_sgh(double rk_alpha,
                       const mesh_t &mesh,
                       const DViewCArrayKokkos <double> &node_vel,
                       const DViewCArrayKokkos <double> &node_coords,
                       DViewCArrayKokkos <double> &elem_sie,
                       const DViewCArrayKokkos <double> &elem_mass,
                       const DViewCArrayKokkos <double> &corner_force){
   
    const size_t rk_level = simparam.rk_num_bins - 1; 
    int num_dims = simparam.num_dims;

    // loop over all the elements in the mesh
    FOR_ALL_CLASS (elem_gid, 0, rnum_elem, {

        double elem_power = 0.0;

        // --- tally the contribution from each corner to the element ---

        // Loop over the nodes in the element
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
            
            size_t corner_lid = node_lid;
            
            // Get node global id for the local node id
            size_t node_gid = nodes_in_elem(elem_gid, node_lid);
            
            // Get the corner global id for the local corner id
            size_t corner_gid = corners_in_elem(elem_gid, corner_lid);
            
            double node_radius = 1;
            if(num_dims==2){
                node_radius = node_coords(rk_level,node_gid,1);
            }

            // calculate the Power=F dot V for this corner
            for (size_t dim=0; dim<num_dims; dim++){
                
                double half_vel = (node_vel(rk_level, node_gid, dim) + node_vel(0, node_gid, dim))*0.5;
                elem_power += corner_force(corner_gid, dim)*node_radius*half_vel;
                
            } // end for dim
            
        } // end for node_lid

        // update the specific energy
        elem_sie(rk_level, elem_gid) = elem_sie(0, elem_gid) -
                                rk_alpha*dt/elem_mass(elem_gid) * elem_power;

    }); // end parallel loop over the elements
    
    return;
} // end subroutine
