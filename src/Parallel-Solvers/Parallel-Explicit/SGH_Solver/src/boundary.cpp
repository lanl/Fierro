#include "mesh.h"
#include "state.h"
#include "FEA_Module_SGH.h"
#include "Simulation_Parameters/FEA_Module/Boundary_Conditions.h"

void FEA_Module_SGH::boundary_velocity(const mesh_t&                    mesh,
                                       const DCArrayKokkos<boundary_t>& boundary,
                                       DViewCArrayKokkos<double>&       node_vel)
{
    // error and debug flag
    // DCArrayKokkos<bool> print_flag(1, "print_flag");
    // print_flag.host(0) = false;
    // print_flag.update_device();

    const size_t rk_level = rk_num_bins - 1;
    int          num_dims = num_dim;
    // Loop over boundary sets
    for (size_t bdy_set = 0; bdy_set < num_bdy_sets; bdy_set++)
    {
        // Loop over boundary nodes in a boundary set
        FOR_ALL_CLASS(bdy_node_lid, 0, num_bdy_nodes_in_set.host(bdy_set), {
            // reflected (boundary array is on the device)
            if (boundary(bdy_set).type == BOUNDARY_CONDITION_TYPE::reflected)
            {
                // directions with hydro_bc:
                // x_plane  = 0,
                // y_plane  = 1,
                // z_plane  = 2,
                size_t direction = boundary(bdy_set).surface.planar_surface_index();

                size_t bdy_node_gid = bdy_nodes_in_set(bdy_set, bdy_node_lid);

                // Set velocity to zero in that directdion
                node_vel(rk_level, bdy_node_gid, direction) = 0.0;
            }
            else if (boundary(bdy_set).type == BOUNDARY_CONDITION_TYPE::fixed_position)
            {
                size_t bdy_node_gid = bdy_nodes_in_set(bdy_set, bdy_node_lid);

                // debug clause
                // if(bdy_node_gid==549412) print_flag(0) = true;

                for (size_t dim = 0; dim < num_dims; dim++)
                {
                    // Set velocity to zero
                    node_vel(rk_level, bdy_node_gid, dim) = 0.0;
                }
            }
            else if (boundary(bdy_set).type == BOUNDARY_CONDITION_TYPE::velocity)
            {
                size_t bdy_node_gid = mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid);

                node_vel(rk_level, bdy_node_gid, 0) = boundary(bdy_set).u;
                node_vel(rk_level, bdy_node_gid, 1) = boundary(bdy_set).v;
                if (mesh.num_dims == 3)
                {
                    node_vel(rk_level, bdy_node_gid, 2) = boundary(bdy_set).w;
                }
            } // end if
        }); // end for bdy_node_lid
    } // end for bdy_set

    // debug check
    // print_flag.update_host();
    // if(print_flag.host(0)) std::cout << "found boundary node with id 549412" << std::endl;
    return;
} // end boundary_velocity function