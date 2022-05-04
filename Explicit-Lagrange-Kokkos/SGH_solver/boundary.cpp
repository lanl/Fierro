// -------------------------------------------------------
// This function applys the boundary condition
// to points on a list of patches created at setup
//--------------------------------------------------------

#include "mesh.h"
#include "state.h"


void boundary_velocity(const mesh_t &mesh,
                       const CArrayKokkos <boundary_t> &boundary,
                       DViewCArrayKokkos <double> &node_vel){

    
    // Loop over boundary sets
    FOR_ALL(bdy_set,0, mesh.num_bdy_sets,{
        
	// reflected (boundary array is on the device)
        if(boundary(bdy_set).hydro_bc == 1){
	        
            // Loop over boundary nodes in a boundary set
            for (size_t bdy_node_lid=0; bdy_node_lid<mesh.num_bdy_nodes_in_set(bdy_set); bdy_node_lid++) {
            
                // directions with hydro_bc:
                // x_plane  = 0,
                // y_plane  = 1,
                // z_plane  = 2,
                size_t direction = boundary(bdy_set).surface;
                
                size_t bdy_node_gid = mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid);
                    
                // Set velocity to zero in that directdion
                node_vel(1, bdy_node_gid, direction) = 0.0;
                
            } // end for bdy_node_lid
	    
	} // end if
       
    }); // end FOR_ALL bdy_set
    
} // end boundary_velocity function
