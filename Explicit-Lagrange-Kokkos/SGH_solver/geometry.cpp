// -----------------------------------------------------------------------------
// This code handles the geometric information for the mesh for the SHG solver
//------------------------------------------------------------------------------
#include "matar.h"
#include "mesh.h"
#include "state.h"


void update_position_sgh(double rk_alpha,
                         double dt,
                         const mesh_t &mesh,
                         DViewCArrayKokkos <double> &node_coords,
                         const DViewCArrayKokkos <double> &node_vel){
    
    // loop over all the nodes in the mesh
    FOR_ALL(node_gid, 0, mesh.num_nodes, {

        for (int dim = 0; dim < mesh.num_dims; dim++){
            double half_vel = (node_vel(1, node_gid, dim) + node_vel(0, node_gid, dim))*0.5;
            node_coords(1, node_gid, dim) = node_coords(0, node_gid, dim) + rk_alpha*dt*half_vel;
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
void get_bmatrix(const ViewCArrayKokkos <double> &B_matrix,
                 const size_t elem_gid,
                 const DViewCArrayKokkos <double> &node_coords,
                 const ViewCArrayKokkos <size_t>  &elem_node_gids){

    const size_t num_nodes = 8;

    double x_array[8];
    double y_array[8];
    double z_array[8];
    
    // x, y, z coordinates of elem vertices
    auto x  = ViewCArrayKokkos <double> (x_array, num_nodes);
    auto y  = ViewCArrayKokkos <double> (y_array, num_nodes);
    auto z  = ViewCArrayKokkos <double> (z_array, num_nodes);

    // get the coordinates of the nodes(rk,elem,node) in this element
    for (int node_lid = 0; node_lid < num_nodes; node_lid++){
        x(node_lid) = node_coords(1, elem_node_gids(node_lid), 0);
        y(node_lid) = node_coords(1, elem_node_gids(node_lid), 1);
        z(node_lid) = node_coords(1, elem_node_gids(node_lid), 2);
    } // end for

    double twelth = 1./12.;
        
    B_matrix(0,0) = (-y(1)*z(2) -y(1)*z(3) +y(1)*z(4) +y(1)*z(5) +y(2)*z(1) -y(2)*z(3)
                  +y(3)*z(1) +y(3)*z(2) -y(3)*z(4) -y(3)*z(7) -y(4)*z(1) +y(4)*z(3)
                  -y(4)*z(5) +y(4)*z(7) -y(5)*z(1) +y(5)*z(4) +y(7)*z(3) -y(7)*z(4) )*twelth;
        
    B_matrix(1,0) = (+y(0)*z(2) +y(0)*z(3) -y(0)*z(4) -y(0)*z(5) -y(2)*z(0) -y(2)*z(3)
                  +y(2)*z(5) +y(2)*z(6) -y(3)*z(0) +y(3)*z(2) +y(4)*z(0) -y(4)*z(5)
                  +y(5)*z(0) -y(5)*z(2) +y(5)*z(4) -y(5)*z(6) -y(6)*z(2) +y(6)*z(5) )*twelth;

    B_matrix(2,0) = (-y(0)*z(1) +y(0)*z(3) +y(1)*z(0) +y(1)*z(3) -y(1)*z(5) -y(1)*z(6)
                  -y(3)*z(0) -y(3)*z(1) +y(3)*z(6) +y(3)*z(7) +y(5)*z(1) -y(5)*z(6)
                  +y(6)*z(1) -y(6)*z(3) +y(6)*z(5) -y(6)*z(7) -y(7)*z(3) +y(7)*z(6) )*twelth;

    B_matrix(3,0) = (-y(0)*z(1) -y(0)*z(2) +y(0)*z(4) +y(0)*z(7) +y(1)*z(0) -y(1)*z(2)
                  +y(2)*z(0) +y(2)*z(1) -y(2)*z(6) -y(2)*z(7) -y(4)*z(0) +y(4)*z(7)
                  +y(6)*z(2) -y(6)*z(7) -y(7)*z(0) +y(7)*z(2) -y(7)*z(4) +y(7)*z(6) )*twelth;

    B_matrix(4,0) = (+y(0)*z(1) -y(0)*z(3) +y(0)*z(5) -y(0)*z(7) -y(1)*z(0) +y(1)*z(5)
                  +y(3)*z(0) -y(3)*z(7) -y(5)*z(0) -y(5)*z(1) +y(5)*z(6) +y(5)*z(7)
                  -y(6)*z(5) +y(6)*z(7) +y(7)*z(0) +y(7)*z(3) -y(7)*z(5) -y(7)*z(6) )*twelth;

    B_matrix(5,0) = (+y(0)*z(1) -y(0)*z(4) -y(1)*z(0) +y(1)*z(2) -y(1)*z(4) +y(1)*z(6)
                  -y(2)*z(1) +y(2)*z(6) +y(4)*z(0) +y(4)*z(1) -y(4)*z(6) -y(4)*z(7)
                  -y(6)*z(1) -y(6)*z(2) +y(6)*z(4) +y(6)*z(7) +y(7)*z(4) -y(7)*z(6))*twelth;

    B_matrix(6,0) = (+y(1)*z(2) -y(1)*z(5) -y(2)*z(1) +y(2)*z(3) -y(2)*z(5) +y(2)*z(7)
                  -y(3)*z(2) +y(3)*z(7) +y(4)*z(5) -y(4)*z(7) +y(5)*z(1) +y(5)*z(2)
                  -y(5)*z(4) -y(5)*z(7) -y(7)*z(2) -y(7)*z(3) +y(7)*z(4) +y(7)*z(5) )*twelth;

    B_matrix(7,0) = (-y(0)*z(3) +y(0)*z(4) +y(2)*z(3) -y(2)*z(6) +y(3)*z(0) -y(3)*z(2)
                  +y(3)*z(4) -y(3)*z(6) -y(4)*z(0) -y(4)*z(3) +y(4)*z(5) +y(4)*z(6)
                  -y(5)*z(4) +y(5)*z(6) +y(6)*z(2) +y(6)*z(3) -y(6)*z(4) -y(6)*z(5) )*twelth;

    B_matrix(0,1) = (-z(1)*x(2) -z(1)*x(3) +z(1)*x(4) +z(1)*x(5) +z(2)*x(1) -z(2)*x(3)
                  +z(3)*x(1) +z(3)*x(2) -z(3)*x(4) -z(3)*x(7) -z(4)*x(1) +z(4)*x(3)
                  -z(4)*x(5) +z(4)*x(7) -z(5)*x(1) +z(5)*x(4) +z(7)*x(3) -z(7)*x(4) )*twelth;

    B_matrix(1,1) = (+z(0)*x(2) +z(0)*x(3) -z(0)*x(4) -z(0)*x(5) -z(2)*x(0) -z(2)*x(3)
                  +z(2)*x(5) +z(2)*x(6) -z(3)*x(0) +z(3)*x(2) +z(4)*x(0) -z(4)*x(5)
                  +z(5)*x(0) -z(5)*x(2) +z(5)*x(4) -z(5)*x(6) -z(6)*x(2) +z(6)*x(5) )*twelth;

    B_matrix(2,1) = (-z(0)*x(1) +z(0)*x(3) +z(1)*x(0) +z(1)*x(3) -z(1)*x(5) -z(1)*x(6)
                  -z(3)*x(0) -z(3)*x(1) +z(3)*x(6) +z(3)*x(7) +z(5)*x(1) -z(5)*x(6)
                  +z(6)*x(1) -z(6)*x(3) +z(6)*x(5) -z(6)*x(7) -z(7)*x(3) +z(7)*x(6) )*twelth;

    B_matrix(3,1) = (-z(0)*x(1) -z(0)*x(2) +z(0)*x(4) +z(0)*x(7) +z(1)*x(0) -z(1)*x(2)
                  +z(2)*x(0) +z(2)*x(1) -z(2)*x(6) -z(2)*x(7) -z(4)*x(0) +z(4)*x(7)
                  +z(6)*x(2) -z(6)*x(7) -z(7)*x(0) +z(7)*x(2) -z(7)*x(4) +z(7)*x(6) )*twelth;

    B_matrix(4,1) = (+z(0)*x(1) -z(0)*x(3) +z(0)*x(5) -z(0)*x(7) -z(1)*x(0) +z(1)*x(5)
                  +z(3)*x(0) -z(3)*x(7) -z(5)*x(0) -z(5)*x(1) +z(5)*x(6) +z(5)*x(7)
                  -z(6)*x(5) +z(6)*x(7) +z(7)*x(0) +z(7)*x(3) -z(7)*x(5) -z(7)*x(6) )*twelth;

    B_matrix(5,1) = (+z(0)*x(1) -z(0)*x(4) -z(1)*x(0) +z(1)*x(2) -z(1)*x(4) +z(1)*x(6)
                  -z(2)*x(1) +z(2)*x(6) +z(4)*x(0) +z(4)*x(1) -z(4)*x(6) -z(4)*x(7)
                  -z(6)*x(1) -z(6)*x(2) +z(6)*x(4) +z(6)*x(7) +z(7)*x(4) -z(7)*x(6) )*twelth;

    B_matrix(6,1) = (+z(1)*x(2) -z(1)*x(5) -z(2)*x(1) +z(2)*x(3) -z(2)*x(5) +z(2)*x(7)
                  -z(3)*x(2) +z(3)*x(7) +z(4)*x(5) -z(4)*x(7) +z(5)*x(1) +z(5)*x(2)
                  -z(5)*x(4) -z(5)*x(7) -z(7)*x(2) -z(7)*x(3) +z(7)*x(4) +z(7)*x(5) )*twelth;

    B_matrix(7,1) = (-z(0)*x(3) +z(0)*x(4) +z(2)*x(3) -z(2)*x(6) +z(3)*x(0) -z(3)*x(2)
                  +z(3)*x(4) -z(3)*x(6) -z(4)*x(0) -z(4)*x(3) +z(4)*x(5) +z(4)*x(6)
                  -z(5)*x(4) +z(5)*x(6) +z(6)*x(2) +z(6)*x(3) -z(6)*x(4) -z(6)*x(5) )*twelth;

    B_matrix(0,2) = (-x(1)*y(2) -x(1)*y(3) +x(1)*y(4) +x(1)*y(5) +x(2)*y(1) -x(2)*y(3)
                  +x(3)*y(1) +x(3)*y(2) -x(3)*y(4) -x(3)*y(7) -x(4)*y(1) +x(4)*y(3)
                  -x(4)*y(5) +x(4)*y(7) -x(5)*y(1) +x(5)*y(4) +x(7)*y(3) -x(7)*y(4) )*twelth;

    B_matrix(1,2) = (+x(0)*y(2) +x(0)*y(3) -x(0)*y(4) -x(0)*y(5) -x(2)*y(0) -x(2)*y(3)
                  +x(2)*y(5) +x(2)*y(6) -x(3)*y(0) +x(3)*y(2) +x(4)*y(0) -x(4)*y(5)
                  +x(5)*y(0) -x(5)*y(2) +x(5)*y(4) -x(5)*y(6) -x(6)*y(2) +x(6)*y(5) )*twelth;

    B_matrix(2,2) = (-x(0)*y(1) +x(0)*y(3) +x(1)*y(0) +x(1)*y(3) -x(1)*y(5) -x(1)*y(6)
                  -x(3)*y(0) -x(3)*y(1) +x(3)*y(6) +x(3)*y(7) +x(5)*y(1) -x(5)*y(6)
                  +x(6)*y(1) -x(6)*y(3) +x(6)*y(5) -x(6)*y(7) -x(7)*y(3) +x(7)*y(6) )*twelth;

    B_matrix(3,2) = (-x(0)*y(1) -x(0)*y(2) +x(0)*y(4) +x(0)*y(7) +x(1)*y(0) -x(1)*y(2)
                  +x(2)*y(0) +x(2)*y(1) -x(2)*y(6) -x(2)*y(7) -x(4)*y(0) +x(4)*y(7)
                  +x(6)*y(2) -x(6)*y(7) -x(7)*y(0) +x(7)*y(2) -x(7)*y(4) +x(7)*y(6) )*twelth;

    B_matrix(4,2) = (+x(0)*y(1) -x(0)*y(3) +x(0)*y(5) -x(0)*y(7) -x(1)*y(0) +x(1)*y(5)
                  +x(3)*y(0) -x(3)*y(7) -x(5)*y(0) -x(5)*y(1) +x(5)*y(6) +x(5)*y(7)
                  -x(6)*y(5) +x(6)*y(7) +x(7)*y(0) +x(7)*y(3) -x(7)*y(5) -x(7)*y(6) )*twelth;

    B_matrix(5,2) = (+x(0)*y(1) -x(0)*y(4) -x(1)*y(0) +x(1)*y(2) -x(1)*y(4) +x(1)*y(6)
                  -x(2)*y(1) +x(2)*y(6) +x(4)*y(0) +x(4)*y(1) -x(4)*y(6) -x(4)*y(7)
                  -x(6)*y(1) -x(6)*y(2) +x(6)*y(4) +x(6)*y(7) +x(7)*y(4) -x(7)*y(6) )*twelth;

    B_matrix(6,2) = (+x(1)*y(2) -x(1)*y(5) -x(2)*y(1) +x(2)*y(3) -x(2)*y(5) +x(2)*y(7)
                  -x(3)*y(2) +x(3)*y(7) +x(4)*y(5) -x(4)*y(7) +x(5)*y(1) +x(5)*y(2)
                  -x(5)*y(4) -x(5)*y(7) -x(7)*y(2) -x(7)*y(3) +x(7)*y(4) +x(7)*y(5) )*twelth;

    B_matrix(7,2) = (-x(0)*y(3) +x(0)*y(4) +x(2)*y(3) -x(2)*y(6) +x(3)*y(0) -x(3)*y(2)
                  +x(3)*y(4) -x(3)*y(6) -x(4)*y(0) -x(4)*y(3) +x(4)*y(5) +x(4)*y(6)
                  -x(5)*y(4) +x(5)*y(6) +x(6)*y(2) +x(6)*y(3) -x(6)*y(4) -x(6)*y(5) )*twelth;

} // end subroutine



void get_vol(const DViewCArrayKokkos <double> &elem_vol,
             const DViewCArrayKokkos <double> &node_coords,
             const mesh_t &mesh){
    
    const size_t num_dims = mesh.num_dims;
    
    if (num_dims == 2){
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
            
            // cut out the node_gids for this element
            ViewCArrayKokkos <size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), 2);
            get_vol_quad(elem_vol, elem_gid, node_coords, elem_node_gids);
            
        });
        Kokkos::fence();
    }
    else {
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
            
            // cut out the node_gids for this element
            ViewCArrayKokkos <size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), 8);
            get_vol_hex(elem_vol, elem_gid, node_coords, elem_node_gids);
            
        });
        Kokkos::fence();
    } // end if
    
} // end subroutine


// Exact volume for a hex element
KOKKOS_FUNCTION
void get_vol_hex(const DViewCArrayKokkos <double> &elem_vol,
                 const size_t elem_gid,
                 const DViewCArrayKokkos <double> &node_coords,
                 const ViewCArrayKokkos <size_t>  &elem_node_gids){

    const size_t num_nodes = 8;

    double x_array[8];
    double y_array[8];
    double z_array[8];
    
    // x, y, z coordinates of elem vertices
    auto x  = ViewCArrayKokkos <double> (x_array, num_nodes);
    auto y  = ViewCArrayKokkos <double> (y_array, num_nodes);
    auto z  = ViewCArrayKokkos <double> (z_array, num_nodes);
    
    // get the coordinates of the nodes(rk,elem,node) in this element
    for (int node_lid = 0; node_lid < num_nodes; node_lid++){
        x(node_lid) = node_coords(1, elem_node_gids(node_lid), 0);
        y(node_lid) = node_coords(1, elem_node_gids(node_lid), 1);
        z(node_lid) = node_coords(1, elem_node_gids(node_lid), 2);
    } // end for

    double twelth = 1./12.;
        
    // element volume
    elem_vol(elem_gid) =
        (x(1)*(y(3)*(-z(0) + z(2)) + y(4)*( z(0) - z(5)) + y(0)*(z(2) + z(3) - z(4) - z(5)) + y(6)*(-z(2) + z(5)) + y(5)*(z(0) - z(2) + z(4) - z(6)) + y(2)*(-z(0) - z(3) + z(5) + z(6))) +
         x(7)*(y(0)*(-z(3) + z(4)) + y(6)*( z(2) + z(3)  - z(4) - z(5)) + y(2)*(z(3) - z(6)) + y(3)*(z(0) - z(2) + z(4) - z(6)) + y(5)*(-z(4) + z(6)) + y(4)*(-z(0) - z(3) + z(5) + z(6))) +
         x(3)*(y(1)*( z(0) - z(2)) + y(7)*(-z(0) + z(2)  - z(4) + z(6)) + y(6)*(z(2) - z(7)) + y(2)*(z(0) + z(1) - z(6) - z(7)) + y(4)*(-z(0) + z(7)) + y(0)*(-z(1) - z(2) + z(4) + z(7))) +
         x(5)*(y(0)*( z(1) - z(4)) + y(7)*( z(4) - z(6)) + y(2)*(-z(1) + z(6)) + y(1)*(-z(0) + z(2) - z(4) + z(6)) + y(4)*(z(0) + z(1) - z(6) - z(7)) + y(6)*(-z(1) - z(2) + z(4) + z(7))) +
         x(6)*(y(1)*( z(2) - z(5)) + y(7)*(-z(2) - z(3)  + z(4) + z(5)) + y(5)*(z(1) + z(2) - z(4) - z(7)) + y(4)*(z(5) - z(7)) + y(3)*(-z(2) + z(7)) + y(2)*(-z(1) + z(3) - z(5) + z(7))) +
         x(0)*(y(2)*( z(1) - z(3)) + y(7)*( z(3) - z(4)) + y(5)*(-z(1) + z(4)) + y(1)*(-z(2) - z(3) + z(4) + z(5)) + y(3)*(z(1) + z(2) - z(4) - z(7)) + y(4)*(-z(1) + z(3) - z(5) + z(7))) +
         x(2)*(y(0)*(-z(1) + z(3)) + y(5)*( z(1) - z(6)) + y(1)*(z(0) + z(3) - z(5) - z(6)) + y(7)*(-z(3) + z(6)) + y(6)*(z(1) - z(3) + z(5) - z(7)) + y(3)*(-z(0) - z(1) + z(6) + z(7))) +
         x(4)*(y(1)*(-z(0) + z(5)) + y(7)*( z(0) + z(3)  - z(5) - z(6)) + y(3)*(z(0) - z(7)) + y(0)*(z(1) - z(3) + z(5) - z(7)) + y(6)*(-z(5) + z(7)) + y(5)*(-z(0) - z(1) + z(6) + z(7))))*twelth;

} // end subroutine



KOKKOS_FUNCTION
void get_bmatrix2D(const ViewCArrayKokkos <double> &B_matrix,
                   const size_t elem_gid,
                   const DViewCArrayKokkos <double> &node_coords,
                   const ViewCArrayKokkos <size_t>  &elem_node_gids){

    const size_t num_nodes = 4;

    double x_array[4];
    double y_array[4];
    
    // x, y coordinates of elem vertices
    auto x  = ViewCArrayKokkos <double> (x_array, num_nodes);
    auto y  = ViewCArrayKokkos <double> (y_array, num_nodes);

    // get the coordinates of the nodes(rk,elem,node) in this element
    for (int node_lid = 0; node_lid < num_nodes; node_lid++){
        x(node_lid) = node_coords(1, elem_node_gids(node_lid), 0);
        y(node_lid) = node_coords(1, elem_node_gids(node_lid), 1);
    } // end for

    /* ensight node order   0 1 2 3
       Flanaghan node order 3 4 1 2
    */
    
    B_matrix(0,0) = 0.5*(y(3)-y(1));
                  
    B_matrix(1,0) = 0.5*(y(0)-y(2));
                  
    B_matrix(2,0) = 0.5*(y(1)-y(3));
                  
    B_matrix(3,0) = 0.5*(y(2)-y(0));
                 

    B_matrix(0,1) = 0.5*(x(1)-x(3));
                 
    B_matrix(1,1) = 0.5*(x(2)-x(0));
                  
    B_matrix(2,1) = 0.5*(x(3)-x(1));
                  
    B_matrix(3,1) = 0.5*(x(0)-x(2));


} // end subroutine


KOKKOS_FUNCTION
void get_vol_quad(const DViewCArrayKokkos <double> &elem_vol,
                  const size_t elem_gid,
                  const DViewCArrayKokkos <double> &node_coords,
                  const ViewCArrayKokkos <size_t>  &elem_node_gids){

    const size_t num_nodes = 4;

    double x_array[4];
    double y_array[4];
    
    // x, y coordinates of elem vertices
    auto x  = ViewCArrayKokkos <double> (x_array, num_nodes);
    auto y  = ViewCArrayKokkos <double> (y_array, num_nodes);
     
    // get the coordinates of the nodes(rk,elem,node) in this element
    for (int node_lid = 0; node_lid < num_nodes; node_lid++){
        x(node_lid) = node_coords(1, elem_node_gids(node_lid), 0);
        y(node_lid) = node_coords(1, elem_node_gids(node_lid), 1);
    } // end for

    /* ensight node order   0 1 2 3
       Flanaghan node order 3 4 1 2
    */
    
    // element volume
    elem_vol(elem_gid) = 0.5*((x(0)-x(2))*(y(1)-y(3))+(x(3)-x(1))*(y(0)-y(2)));

} // end subroutine


