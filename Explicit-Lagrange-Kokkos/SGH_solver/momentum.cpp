

#include "mesh.h"
#include "state.h"



// -----------------------------------------------------------------------------
// This function evolves the velocity at the nodes of the mesh
//------------------------------------------------------------------------------
void update_velocity_sgh(double rk_alpha,
                         double dt,
                         const mesh_t &mesh,
                         DViewCArrayKokkos <double> &node_vel,
                         const DViewCArrayKokkos <double> &node_mass,
                         const DViewCArrayKokkos <double> &corner_force
                         ){

    
    const size_t num_dims = mesh.num_dims;
    
    // walk over the nodes to update the velocity
    FOR_ALL(node_gid, 0, mesh.num_nodes, {

        double node_force[3];
        for (size_t dim = 0; dim < num_dims; dim++){
            node_force[dim] = 0.0;
        } // end for dim
        
        // loop over all corners around the node and calculate the nodal force
        for (size_t corner_lid=0; corner_lid<mesh.num_corners_in_node(node_gid); corner_lid++){
        
            // Get corner gid
            size_t corner_gid = mesh.corners_in_node(node_gid, corner_lid);
            
            // loop over dimension
            for (size_t dim = 0; dim < num_dims; dim++){
                node_force[dim] += corner_force(corner_gid, dim);
            } // end for dim
            
        } // end for corner_lid
        
        // update the velocity
        for (int dim = 0; dim < num_dims; dim++){
            node_vel(1, node_gid, dim) = node_vel(0, node_gid, dim) +
                                         rk_alpha * dt*node_force[dim]/node_mass(node_gid);
        } // end for dim
        
    }); // end for parallel for over nodes
    
} // end subroutine update_velocity




// -----------------------------------------------------------------------------
// This function calculates the velocity gradient
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
void get_velgrad(ViewCArrayKokkos <double> &vel_grad,
                 const mesh_t &mesh,
                 const DViewCArrayKokkos <double> &node_vel,
                 const ViewCArrayKokkos <double> &b_matrix,
                 const double elem_vol,
                 const size_t elem_gid){
    

    const size_t num_nodes_in_elem = 8;
    
    double u_array[num_nodes_in_elem];
    double v_array[num_nodes_in_elem];
    double w_array[num_nodes_in_elem];
    ViewCArrayKokkos <double> u(u_array, num_nodes_in_elem); // x-dir vel component
    ViewCArrayKokkos <double> v(v_array, num_nodes_in_elem); // y-dir vel component
    ViewCArrayKokkos <double> w(w_array, num_nodes_in_elem); // z-dir vel component
    
    // get the vertex velocities for the cell
    for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
        
        // Get node gid
        size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

        u(node_lid) = node_vel(1, node_gid, 0);
        v(node_lid) = node_vel(1, node_gid, 1);
        w(node_lid) = node_vel(1, node_gid, 2);
        
    } // end for

    
    // --- calculate the velocity gradient terms ---
    double inverse_vol = 1.0/elem_vol;
    
    // x-dir
    vel_grad(0,0) = (u(0)*b_matrix(0,0) + u(1)*b_matrix(1,0)
                   + u(2)*b_matrix(2,0) + u(3)*b_matrix(3,0)
                   + u(4)*b_matrix(4,0) + u(5)*b_matrix(5,0)
                   + u(6)*b_matrix(6,0) + u(7)*b_matrix(7,0))*inverse_vol;
    
    vel_grad(0,1) = (u(0)*b_matrix(0,1) + u(1)*b_matrix(1,1)
                   + u(2)*b_matrix(2,1) + u(3)*b_matrix(3,1)
                   + u(4)*b_matrix(4,1) + u(5)*b_matrix(5,1)
                   + u(6)*b_matrix(6,1) + u(7)*b_matrix(7,1))*inverse_vol;
    
    vel_grad(0,2) = (u(0)*b_matrix(0,2) + u(1)*b_matrix(1,2)
                   + u(2)*b_matrix(2,2) + u(3)*b_matrix(3,2)
                   + u(4)*b_matrix(4,2) + u(5)*b_matrix(5,2)
                   + u(6)*b_matrix(6,2) + u(7)*b_matrix(7,2))*inverse_vol;
    
    // y-dir
    vel_grad(1,0) = (v(0)*b_matrix(0,0) + v(1)*b_matrix(1,0)
                   + v(2)*b_matrix(2,0) + v(3)*b_matrix(3,0)
                   + v(4)*b_matrix(4,0) + v(5)*b_matrix(5,0)
                   + v(6)*b_matrix(6,0) + v(7)*b_matrix(7,0))*inverse_vol;
    
    vel_grad(1,1) = (v(0)*b_matrix(0,1) + v(1)*b_matrix(1,1)
                   + v(2)*b_matrix(2,1) + v(3)*b_matrix(3,1)
                   + v(4)*b_matrix(4,1) + v(5)*b_matrix(5,1)
                   + v(6)*b_matrix(6,1) + v(7)*b_matrix(7,1))*inverse_vol;
    vel_grad(1,2) = (v(0)*b_matrix(0,2) + v(1)*b_matrix(1,2)
                   + v(2)*b_matrix(2,2) + v(3)*b_matrix(3,2)
                   + v(4)*b_matrix(4,2) + v(5)*b_matrix(5,2)
                   + v(6)*b_matrix(6,2) + v(7)*b_matrix(7,2))*inverse_vol;
    
    // z-dir
    vel_grad(2,0) = (w(0)*b_matrix(0,0) + w(1)*b_matrix(1,0)
                   + w(2)*b_matrix(2,0) + w(3)*b_matrix(3,0)
                   + w(4)*b_matrix(4,0) + w(5)*b_matrix(5,0)
                   + w(6)*b_matrix(6,0) + w(7)*b_matrix(7,0))*inverse_vol;
    
    vel_grad(2,1) = (w(0)*b_matrix(0,1) + w(1)*b_matrix(1,1)
                   + w(2)*b_matrix(2,1) + w(3)*b_matrix(3,1)
                   + w(4)*b_matrix(4,1) + w(5)*b_matrix(5,1)
                   + w(6)*b_matrix(6,1) + w(7)*b_matrix(7,1))*inverse_vol;
    
    vel_grad(2,2) = (w(0)*b_matrix(0,2) + w(1)*b_matrix(1,2)
                   + w(2)*b_matrix(2,2) + w(3)*b_matrix(3,2)
                   + w(4)*b_matrix(4,2) + w(5)*b_matrix(5,2)
                   + w(6)*b_matrix(6,2) + w(7)*b_matrix(7,2))*inverse_vol;

} // end function



// -----------------------------------------------------------------------------
// This subroutine to calculate the velocity divergence in all elements
//------------------------------------------------------------------------------
void get_divergence(DViewCArrayKokkos <double> &elem_div,
                    const mesh_t &mesh,
                    const DViewCArrayKokkos <double> &node_coords,
                    const DViewCArrayKokkos <double> &node_vel,
                    const DViewCArrayKokkos <double> &elem_vol){
    
    // --- calculate the forces acting on the nodes from the element ---
    FOR_ALL (elem_gid, 0, mesh.num_elems, {
    
        const size_t num_nodes_in_elem = 8;
        const size_t num_dims = 3;
        
        double u_array[num_nodes_in_elem];
        double v_array[num_nodes_in_elem];
        double w_array[num_nodes_in_elem];
        ViewCArrayKokkos <double> u(u_array, num_nodes_in_elem); // x-dir vel component
        ViewCArrayKokkos <double> v(v_array, num_nodes_in_elem); // y-dir vel component
        ViewCArrayKokkos <double> w(w_array, num_nodes_in_elem); // z-dir vel component
        
        
        // The b_matrix are the outward corner area normals
        double b_matrix_array[24];
        ViewCArrayKokkos <double> b_matrix(b_matrix_array, num_nodes_in_elem, num_dims);
        get_bmatrix(b_matrix,
                    elem_gid,
                    node_coords,
                    mesh);
        
        // get the vertex velocities for the cell
        for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
        
            // Get node gid
            size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
    
            u(node_lid) = node_vel(1, node_gid, 0);
            v(node_lid) = node_vel(1, node_gid, 1);
            w(node_lid) = node_vel(1, node_gid, 2);
        
        } // end for
    
        
        // --- calculate the velocity divergence terms ---
        double inverse_vol = 1.0/elem_vol(elem_gid);
        
        elem_div(elem_gid) = 0.0;
        
        // x-dir
        elem_div(elem_gid) += (u(0)*b_matrix(0,0) + u(1)*b_matrix(1,0)
                             + u(2)*b_matrix(2,0) + u(3)*b_matrix(3,0)
                             + u(4)*b_matrix(4,0) + u(5)*b_matrix(5,0)
                             + u(6)*b_matrix(6,0) + u(7)*b_matrix(7,0))*inverse_vol;
        
        // y-dir
        elem_div(elem_gid) += (v(0)*b_matrix(0,1) + v(1)*b_matrix(1,1)
                             + v(2)*b_matrix(2,1) + v(3)*b_matrix(3,1)
                             + v(4)*b_matrix(4,1) + v(5)*b_matrix(5,1)
                             + v(6)*b_matrix(6,1) + v(7)*b_matrix(7,1))*inverse_vol;
        
        // z-dir
        elem_div(elem_gid) += (w(0)*b_matrix(0,2) + w(1)*b_matrix(1,2)
                             + w(2)*b_matrix(2,2) + w(3)*b_matrix(3,2)
                             + w(4)*b_matrix(4,2) + w(5)*b_matrix(5,2)
                             + w(6)*b_matrix(6,2) + w(7)*b_matrix(7,2))*inverse_vol;
        
    });  // end parallel for over elem_gid
    
} // end subroutine

