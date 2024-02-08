
#include "mesh.h"
#include "state.h"

void update_energy_sgh(double                           rk_alpha,
                       double                           dt,
                       const mesh_t&                    mesh,
                       const DViewCArrayKokkos<double>& node_vel,
                       const DViewCArrayKokkos<double>& node_coords,
                       DViewCArrayKokkos<double>&       elem_sie,
                       const DViewCArrayKokkos<double>& elem_mass,
                       const DViewCArrayKokkos<double>& corner_force)
{
    // loop over all the elements in the mesh
    FOR_ALL(elem_gid, 0, mesh.num_elems, {
        double elem_power = 0.0;

        // --- tally the contribution from each corner to the element ---

        // Loop over the nodes in the element
        for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++)
        {
            size_t corner_lid = node_lid;

            // Get node global id for the local node id
            size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

            // Get the corner global id for the local corner id
            size_t corner_gid = mesh.corners_in_elem(elem_gid, corner_lid);

            double node_radius = 1;
            if (mesh.num_dims == 2)
            {
                node_radius = node_coords(1, node_gid, 1);
            }

            // calculate the Power=F dot V for this corner
            for (size_t dim = 0; dim < mesh.num_dims; dim++)
            {
                double half_vel = (node_vel(1, node_gid, dim) + node_vel(0, node_gid, dim)) * 0.5;
                elem_power     += corner_force(corner_gid, dim) * node_radius * half_vel;
            } // end for dim
        } // end for node_lid

        // update the specific energy
        elem_sie(1, elem_gid) = elem_sie(0, elem_gid) -
                                rk_alpha * dt / elem_mass(elem_gid) * elem_power;
    }); // end parallel loop over the elements

    return;
} // end subroutine
