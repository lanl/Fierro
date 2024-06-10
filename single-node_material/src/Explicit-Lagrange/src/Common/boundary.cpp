// -------------------------------------------------------
// This function applys the boundary condition
// to points on a list of patches created at setup
//--------------------------------------------------------

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void boundary_force(){

    // Loop over boundary sets
    for(int bdy_set = 0; bdy_set < mesh.num_bdy_sets(); bdy_set++){

        int direction = boundary[bdy_set].surface;

        // Loop over boundary patches in boundary set
        for (int bdy_patch_gid = 0; bdy_patch_gid < mesh.num_bdy_patches_in_set(bdy_set); bdy_patch_gid++){

            // get the global id for this boundary patch
            int patch_gid = mesh.bdy_patches_in_set(bdy_set, bdy_patch_gid);

            // apply boundary condition at nodes on boundary
            for(int node_lid = 0; node_lid < 4; node_lid++){

                int node_gid = mesh.node_in_patch(patch_gid, node_lid);

                // Set nodal force to zero
                node.force(node_gid, direction) = 0.0;
            }
        }
    }
}

void boundary_velocity(){

    // Loop over boundary sets
    for(int bdy_set = 0; bdy_set < mesh.num_bdy_sets(); bdy_set++){

        int direction = boundary[bdy_set].surface;

        // Loop over boundary patches in boundary set
        for (int bdy_patch_gid = 0; bdy_patch_gid < mesh.num_bdy_patches_in_set(bdy_set); bdy_patch_gid++){
                
            // get the global id for this boundary patch
            int patch_gid = mesh.bdy_patches_in_set(bdy_set, bdy_patch_gid);

            // apply boundary condition at nodes on boundary
            for(int node_lid = 0; node_lid < 4; node_lid++){
                
                int node_gid = mesh.node_in_patch(patch_gid, node_lid);

                // Set nodal force to zero
                node.vel(1, node_gid, direction) = 0.0;
                // node.force(node_gid, direction) = 0.0;

            }
        }
    }
}