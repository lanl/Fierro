// -------------------------------------------------------
// This function applys the boundary condition
// to points on a list of patches created at setup
//--------------------------------------------------------

#include "mesh.h"
#include "state.h"


void boundary_velocity(const mesh_t &mesh,
                       const CArrayKokkos <boundary_t> &boundary,
                       DViewCArrayKokkos <double> &node_vel,
                       const double time_value){

    
    // Loop over boundary sets
    for (size_t bdy_set=0; bdy_set<mesh.num_bdy_sets; bdy_set++){
        
        // Loop over boundary nodes in a boundary set
        FOR_ALL(bdy_node_lid, 0, mesh.num_bdy_nodes_in_set.host(bdy_set), {
                
            // reflected (boundary array is on the device)
            if (boundary(bdy_set).hydro_bc == bdy::reflected){
            
                // directions with hydro_bc:
                // x_plane  = 0,
                // y_plane  = 1,
                // z_plane  = 2,
                size_t direction = boundary(bdy_set).surface;
                
                size_t bdy_node_gid = mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid);
                    
                // Set velocity to zero in that directdion
                node_vel(1, bdy_node_gid, direction) = 0.0;
                        
            }
            else if (boundary(bdy_set).hydro_bc == bdy::fixed){
                
                
                size_t bdy_node_gid = mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid);
                
                for(size_t dim=0; dim<mesh.num_dims; dim++){
                    // Set velocity to zero
                    node_vel(1, bdy_node_gid, dim) = 0.0;
                }
                
            }// end if
            else if (boundary(bdy_set).hydro_bc == bdy::velocity){
                
                
                size_t bdy_node_gid = mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid);
                
                
                // directions with hydro_bc:
                // x_plane  = 0,
                // y_plane  = 1,
                // z_plane  = 2,
                size_t direction = boundary(bdy_set).surface;
                
                // Set velocity to that directdion to specified value
                // if t_end > time > t_start
                // v(t) = v0 exp(-v1*(time - time_start) )
                if (time_value >= boundary(bdy_set).hydro_bc_vel_t_start &&
                    time_value <= boundary(bdy_set).hydro_bc_vel_t_end){
                    
                    double time_delta = time_value - boundary(bdy_set).hydro_bc_vel_t_start;
                    
                    node_vel(1, bdy_node_gid, direction) =
                        boundary(bdy_set).hydro_bc_vel_0 *
                        exp(-boundary(bdy_set).hydro_bc_vel_1 * time_delta );
                    
                } // end if on time
                
            }// end if
            
                
        }); // end for bdy_node_lid
	    
    } // end for bdy_set
    
    return;
} // end boundary_velocity function
