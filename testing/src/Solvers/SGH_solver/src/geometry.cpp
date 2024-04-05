
#include "sgh_solver.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn update_position
///
/// \brief Updates the nodal positions based on the nodal velocity
///
/// \param Runge Kutta time integration alpha value
/// \param Time step size
/// \param Number of dimensions in the mesh (REMOVE)
/// \param Number of nodes in the mesh
/// \param View of nodal position data
/// \param View of nodal velocity data
///
/////////////////////////////////////////////////////////////////////////////
void SGH::update_position(double rk_alpha,
    double dt,
    const size_t num_dims,
    const size_t num_nodes,
    DCArrayKokkos<double>& node_coords,
    const DCArrayKokkos<double>& node_vel)
{
    // loop over all the nodes in the mesh
    FOR_ALL(node_gid, 0, num_nodes, {
        for (int dim = 0; dim < num_dims; dim++) {
            double half_vel = (node_vel(1, node_gid, dim) + node_vel(0, node_gid, dim)) * 0.5;
            node_coords(1, node_gid, dim) = node_coords(0, node_gid, dim) + rk_alpha * dt * half_vel;
        }
    }); // end parallel for over nodes
} // end subroutine