#include "ref_elem.h"
#include "mesh.h"
#include "state.h"

void assemble_A1(   CArrayKokkos <double> &A1,
                    CArrayKokkos <double> &residual_in_elem,
                    const size_t stage,
                    const double dt,
                    const mesh_t &mesh,
                    CArrayKokkos <double> &M_dot_u,
                    CArrayKokkos <double> &F_dot_ones,
                    const CArrayKokkos <double> &force_tensor,
                    const CArrayKokkos <double> &mass_matrix,
                    const DViewCArrayKokkos <double> &node_vel ){


    // Assemble residual in each element //
    FOR_ALL( node_gid_1, 0, mesh.num_nodes,{
    //FOR_ALL(elem_gid, 0, mesh.num_elems, {

        // Computes M.\delta^k u
        //for (int node_lid_1 = 0; node_lid_1 < mesh.num_nodes_in_elem; node_lid_1++){
            //int node_gid_1 = mesh.nodes_in_elem(elem_gid, node_lid_1);

            for (int dim = 0; dim < mesh.num_dims; dim++){
                

                for (int node_gid_2 = 0; node_gid_2 < mesh.num_nodes; node_gid_2++){
                // for (int node_lid_2 = 0; node_lid_2 < mesh.num_nodes_in_elem; node_lid_2++){
                //     int node_gid_2 = mesh.nodes_in_elem(elem_gid, node_lid_2);

                    M_dot_u(node_gid_1, dim ) += mass_matrix(node_gid_1, node_gid_2)*( node_vel(stage, node_gid_2, dim) - node_vel(0, node_gid_2, dim) );

                }// end loop over node_lid_2
                
                // Compute \int 1.F.sigma dt
                for (int zone_gid = 0; zone_gid < mesh.num_zones; zone_gid++){
                // for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
                //     int zone_gid = mesh.zones_in_elem(elem_gid, zone_lid);
                
                    F_dot_ones(node_gid_1, dim) += 0.5*( force_tensor(stage, node_gid_1, zone_gid, dim) + force_tensor(0, node_gid_1, zone_gid, dim) );

                }// end loop over zone_lid

                //residual_in_elem(elem_gid, node_gid_1, dim) = M_dot_u(node_gid_1, dim) + dt*F_dot_ones(node_gid_1, dim);
                A1(stage, node_gid_1, dim) = M_dot_u(node_gid_1, dim)/dt + F_dot_ones(node_gid_1, dim);
                
            }// end loop over dim
            
        //}// end loop over node_lid_1
    } );// end FOR_ALL
    Kokkos::fence();

    // for (int elem_gid =0; elem_gid < mesh.num_elems; elem_gid++){
    //     for (int node_gid = 0; node_gid < mesh.num_nodes; node_gid++){
    //         for (int dim = 0; dim < mesh.num_dims; dim++){
    //             printf("residual_in_elem(%d, %d, %d) = %f \n", elem_gid, node_gid, dim, residual_in_elem(elem_gid, node_gid, dim));
    //         }// end loop over dim
    //     }// end loop over node_gid
    // }// end loop over elem_gid

    // Assemble A1 operator, A1 = sum_{E \ni i} \Phi^E_i(u^k)
    // FOR_ALL(node_gid, 0, mesh.num_nodes, {

    //     int num_elems_attached_to_node = mesh.num_corners_in_node(node_gid);
    //     //printf("num_elems attached to node %d is %d \n", node_gid, num_elems_attached_to_node);

    //     for (int dim = 0; dim < mesh.num_dims; dim++){
    
    //         for (int elem_in_node = 0; elem_in_node < num_elems_attached_to_node; elem_in_node++){
    //             int elem_in_node_gid = mesh.elems_in_node(node_gid, elem_in_node);
    //             //printf("elem_in_node_gid = %d for node  %d with %d elems in node\n", elem_in_node_gid, node_gid, num_elems_attached_to_node);
    //             A1(node_gid, dim) += residual_in_elem(elem_in_node_gid, node_gid, dim);

    //         }// end loop over elem_in_node
    //     }// end loop over dim
    // });// end FOR_ALL
    // Kokkos::fence();


}// end assemble residual