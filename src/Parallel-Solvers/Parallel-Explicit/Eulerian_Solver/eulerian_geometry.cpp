// -----------------------------------------------------------------------------
// This code handles the geometric information for the mesh for the SHG solver
//------------------------------------------------------------------------------
#include "matar.h"
#include "mesh.h"
#include "state.h"
#include "FEA_Module_Eulerian.h"

void FEA_Module_Eulerian::example_function(double rk_alpha,
                         const size_t num_nodes,
                         DViewCArrayKokkos <double> &node_coords,
                         const DViewCArrayKokkos <double> &node_vel){

    const size_t rk_level = rk_num_bins - 1;
    int num_dims = num_dim;

    // example device loop
    FOR_ALL_CLASS(node_gid, 0, num_nodes, {

        for (int dim = 0; dim < num_dims; dim++){
            double half_vel = (node_vel(rk_level, node_gid, dim) + node_vel(0, node_gid, dim))*0.5;
            node_coords(rk_level, node_gid, dim) = node_coords(0, node_gid, dim) + rk_alpha*dt*half_vel;
        }
        
    }); // end parallel for over nodes
    
} // end subroutine


// -----------------------------------------------------------------------------
//  This function claculates
//    B_p =  J^{-T} \cdot (\nabla_{xi} \phi_p w
//  where
//    \phi_p is the basis function for vertex p
//    w is the 1 gauss point for the cell (everything is evaluted at this point)
//    J^{-T} is the inverse transpose of the Jacobi matrix
//    \nabla_{xi} is the gradient opperator in the reference coordinates
//
//   B_p is the OUTWARD corner area normal at node p
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
void FEA_Module_Eulerian::example_device_function(const ViewCArrayKokkos <double> &B_matrix,
                 const size_t elem_gid,
                 const DViewCArrayKokkos <double> &node_coords,
                 const ViewCArrayKokkos <size_t>  &elem_node_gids,
                 const size_t rk_level) const {


} // end subroutine


