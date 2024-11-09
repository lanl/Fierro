// -----------------------------------------------------------------------------
// This code handles the geometric information for the mesh for the hydro solver
//------------------------------------------------------------------------------
#include "matar.h"
#include "mesh.h"
#include "state.h"
#include "ref_elem.h"
#include "rdh.h"


void update_position_rdh(const int stage,
                         double dt,
                         const mesh_t &mesh,
                         DViewCArrayKokkos <double> &node_coords,
                         const DViewCArrayKokkos <double> &node_vel){
    
    // loop over all the nodes in the mesh
    FOR_ALL(node_gid, 0, mesh.num_nodes, {

        for (int dim = 0; dim < mesh.num_dims; dim++){
            
            double half_vel = 0.5*( node_vel(1, node_gid, dim) + node_vel(0, node_gid, dim) );
            node_coords(1, node_gid, dim) = node_coords(0, node_gid, dim) + dt*half_vel;
        }
        
    }); // end parallel for over nodes
    
} // end subroutine



void get_gauss_leg_pt_jacobian(const mesh_t &mesh,
                               const elem_t &elem,
                               const fe_ref_elem_t &ref_elem,
                               const DViewCArrayKokkos <double> &node_coords, 
                               DViewCArrayKokkos <double> &gauss_legendre_jacobian,
                               DViewCArrayKokkos <double> &gauss_legendre_det_j,
                               DViewCArrayKokkos <double> &gauss_legendre_jacobian_inverse, 
                               const int stage){

    // loop over the mesh
    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        
        for(int gauss_lid = 0; gauss_lid < ref_elem.num_gauss_leg_in_elem; gauss_lid++){

            int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);

            for(int dim_i = 0; dim_i < mesh.num_dims; dim_i++){
                for(int dim_j = 0; dim_j < mesh.num_dims; dim_j++){

                    gauss_legendre_jacobian(gauss_gid, dim_i, dim_j) = 0.0;
                    
                    //printf("num_basis = %d \n", ref_elem.num_basis);
                    // Sum over the basis functions and vertices where they are defined
                    for(int node_lid = 0; node_lid < ref_elem.num_basis; node_lid++){
                        
                        int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
                        //printf(" legendre grad basis =  : %f \n", ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, dim_j));
                        //printf(" node_coords at node %d and dim %d is %f \n", node_gid, dim_i, node_coords(1, node_gid , dim_i));

                        gauss_legendre_jacobian(gauss_gid, dim_i, dim_j) += ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, dim_j) 
                                                                            * node_coords(stage, node_gid , dim_i);
                        

                    }// end loop node_id
                } // end dim_j
            } // end dim_i

            gauss_legendre_det_j(gauss_gid) =
                gauss_legendre_jacobian(gauss_gid, 0, 0) * (gauss_legendre_jacobian(gauss_gid, 1, 1) * gauss_legendre_jacobian(gauss_gid, 2, 2) - gauss_legendre_jacobian(gauss_gid, 2, 1) * gauss_legendre_jacobian(gauss_gid, 1, 2)) -
                gauss_legendre_jacobian(gauss_gid, 0, 1) * (gauss_legendre_jacobian(gauss_gid, 1, 0) * gauss_legendre_jacobian(gauss_gid, 2, 2) - gauss_legendre_jacobian(gauss_gid, 1, 2) * gauss_legendre_jacobian(gauss_gid, 2, 0)) +
                gauss_legendre_jacobian(gauss_gid, 0, 2) * (gauss_legendre_jacobian(gauss_gid, 1, 0) * gauss_legendre_jacobian(gauss_gid, 2, 1) - gauss_legendre_jacobian(gauss_gid, 1, 1) * gauss_legendre_jacobian(gauss_gid, 2, 0)); 
            
//            printf(" leg determinant : %f \n", gauss_legendre_det_j(gauss_gid));
            
           real_t invdet = 1.0 / gauss_legendre_det_j(gauss_gid);
 
           gauss_legendre_jacobian_inverse(gauss_gid, 0, 0) = ( gauss_legendre_jacobian(gauss_gid, 1, 1) * gauss_legendre_jacobian(gauss_gid, 2, 2) - gauss_legendre_jacobian(gauss_gid, 2, 1) * gauss_legendre_jacobian(gauss_gid, 1, 2)) * invdet;
           gauss_legendre_jacobian_inverse(gauss_gid, 0, 1) = ( gauss_legendre_jacobian(gauss_gid, 0, 2) * gauss_legendre_jacobian(gauss_gid, 2, 1) - gauss_legendre_jacobian(gauss_gid, 0, 1) * gauss_legendre_jacobian(gauss_gid, 2, 2)) * invdet;
           gauss_legendre_jacobian_inverse(gauss_gid, 0, 2) = ( gauss_legendre_jacobian(gauss_gid, 0, 1) * gauss_legendre_jacobian(gauss_gid, 1, 2) - gauss_legendre_jacobian(gauss_gid, 0, 2) * gauss_legendre_jacobian(gauss_gid, 1, 1)) * invdet;
           
           gauss_legendre_jacobian_inverse(gauss_gid, 1, 0) = ( gauss_legendre_jacobian(gauss_gid, 1, 2) * gauss_legendre_jacobian(gauss_gid, 2, 0) - gauss_legendre_jacobian(gauss_gid, 1, 0) * gauss_legendre_jacobian(gauss_gid, 2, 2)) * invdet;
           gauss_legendre_jacobian_inverse(gauss_gid, 1, 1) = ( gauss_legendre_jacobian(gauss_gid, 0, 0) * gauss_legendre_jacobian(gauss_gid, 2, 2) - gauss_legendre_jacobian(gauss_gid, 0, 2) * gauss_legendre_jacobian(gauss_gid, 2, 0)) * invdet;
           gauss_legendre_jacobian_inverse(gauss_gid, 1, 2) = ( gauss_legendre_jacobian(gauss_gid, 1, 0) * gauss_legendre_jacobian(gauss_gid, 0, 2) - gauss_legendre_jacobian(gauss_gid, 0, 0) * gauss_legendre_jacobian(gauss_gid, 1, 2)) * invdet;
           gauss_legendre_jacobian_inverse(gauss_gid, 2, 0) = ( gauss_legendre_jacobian(gauss_gid, 1, 0) * gauss_legendre_jacobian(gauss_gid, 2, 1) - gauss_legendre_jacobian(gauss_gid, 2, 0) * gauss_legendre_jacobian(gauss_gid, 1, 1)) * invdet;
           gauss_legendre_jacobian_inverse(gauss_gid, 2, 1) = ( gauss_legendre_jacobian(gauss_gid, 2, 0) * gauss_legendre_jacobian(gauss_gid, 0, 1) - gauss_legendre_jacobian(gauss_gid, 0, 0) * gauss_legendre_jacobian(gauss_gid, 2, 1)) * invdet;
           gauss_legendre_jacobian_inverse(gauss_gid, 2, 2) = ( gauss_legendre_jacobian(gauss_gid, 0, 0) * gauss_legendre_jacobian(gauss_gid, 1, 1) - gauss_legendre_jacobian(gauss_gid, 1, 0) * gauss_legendre_jacobian(gauss_gid, 0, 1)) * invdet;

        } // end loop over gauss in element
    }); // end FOR_ALL over elements
} // end subroutine



void get_vol(DViewCArrayKokkos <double> &elem_vol,
             const DViewCArrayKokkos <double> &node_coords,
             const DViewCArrayKokkos <double> &legendre_jacobian_det,
             const mesh_t &mesh,
             const elem_t &elem,
             const fe_ref_elem_t &ref_elem) {
    
    const int num_dims = mesh.num_dims;
    
    
    FOR_ALL(elem_gid, 0, mesh.num_elems, {
    
        // // cut out the node_gids for this element
        // // ViewCArrayKokkos <size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), 8);
        // // printf("number of legendre points : %d \n", elem.num_leg_pts); 
        elem_vol(elem_gid) = 0.0;
        for (int gauss_lid = 0; gauss_lid < ref_elem.num_gauss_leg_in_elem; gauss_lid++){
            int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
            
            //printf("legendre weight : %f\n", legendre_weights(gauss_lid) );//ref_elem.gauss_leg_weights(gauss_lid));
            //printf("legendre J det : %f\n", legendre_jacobian_det(gauss_gid) );//elem.gauss_legendre_det_j(gauss_gid));
            //printf("legendre gid : %d\n", gauss_gid);

            elem_vol(elem_gid) += ref_elem.gauss_leg_weights(gauss_lid)*legendre_jacobian_det(gauss_gid);
            // Kokkos::atomic_add(&elem_vol(elem_gid), ref_elem.gauss_leg_weights(gauss_lid)*legendre_jacobian_det(gauss_gid));
        }
        //printf("elem_vol for elem %d inside loop = %f \n", elem_gid, elem_vol(elem_gid));
        //get_vol_jacobi(elem_vol, elem_gid, node_coords, elem_node_gids);//  
    
    });
    Kokkos::fence();
    
    return;
    
} // end subroutine
