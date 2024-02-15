#include "ref_elem.h"
#include "mesh.h"
#include "state.h"

void assemble_A1(   CArrayKokkos <double> &A1,
                    CArrayKokkos <double> &residual_in_elem,
                    const size_t stage,
                    const double factor,
                    const double dt,
                    const mesh_t &mesh,
                    CArrayKokkos <double> &M_dot_u,
                    CArrayKokkos <double> &ones_dot_F_dot_sigma,
                    const CArrayKokkos <double> &force_tensor,
                    const CArrayKokkos <double> &mass_matrix,
                    const DViewCArrayKokkos <double> &stress_tensor,
                    const DViewCArrayKokkos <double> &node_vel ){
    // "Factor" is the weight for the time integration.  
    // Since we only consider 2nd order, it is 1 in the first stage and 1/2 in the second stage,
    // i.e. 1/(1+stage) since stage \in \{ 0, 1 \}.


    // Assemble residual in each element //
    FOR_ALL(elem_gid, 0, mesh.num_elems, {

        // Computes M.\delta^k u
        for (int node_lid_1 = 0; node_lid_1 < mesh.num_nodes_in_elem; node_lid_1++){
            int node_gid_1 = mesh.nodes_in_elem(elem_gid, node_lid_1);

            for (int node_lid_2 = 0; node_lid_2 < mesh.num_nodes_in_elem; node_lid_2++){
                int node_gid_2 = mesh.nodes_in_elem(elem_gid, node_lid_2);
                
                for (int dim = 0; dim < mesh.num_dims; dim++){

                    M_dot_u( node_gid_1, dim ) += mass_matrix(node_gid_1, node_gid_2)*( node_vel(stage, node_gid_2, dim) - node_vel(0, node_gid_2, dim) );

                }// end loop over dim
            }// end loop over node_lid_2
            

            // Compute \int 1.F.sigma dt
            for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
                int zone_gid = mesh.zones_in_elem(elem_gid, zone_lid);
                for (int contraction_zone_lid = 0; contraction_zone_lid < mesh.num_zones_in_elem; contraction_zone_lid++){
                    int contraction_zone_gid = mesh.zones_in_elem(elem_gid, contraction_zone_lid);
                    for (int dim_1 = 0; dim_1 < mesh.num_dims; dim_1++){
                        for (int dim_2 = 0; dim_2 < mesh.num_dims; dim_2++){
                    

                            ones_dot_F_dot_sigma(node_gid_1, dim_1) += factor*( force_tensor(stage, node_gid_1, zone_gid, contraction_zone_gid, dim_2)
                                                                                *stress_tensor(stage, zone_gid, dim_2, dim_1) 
                                                                                + force_tensor(0, node_gid_1, zone_gid, contraction_zone_gid, dim_2)
                                                                                *stress_tensor(0, zone_gid, dim_2, dim_1) );

                        }// end loop over dim_2
                    }// end loop over dim_1
                }// end loop over contraction_zone_lid
            }// end loop over zone_lid
            

            for (int dim = 0; dim < mesh.num_dims; dim++){

                residual_in_elem(elem_gid, node_gid_1, dim) = M_dot_u(node_gid_1, dim) + dt*ones_dot_F_dot_sigma(node_gid_1, dim);

            }
        }// end loop over node_lid_1
    } );// end FOR_ALL
    Kokkos::fence();

    // Assemble A1 operator, A1 = sum_{E \ni i} \Phi^E_i(u^k)
    FOR_ALL(node_gid, 0, mesh.num_nodes, {

        int num_elems_attached_to_node = mesh.num_corners_in_node(node_gid);

        for (int elem_in_node = 0; elem_in_node < num_elems_attached_to_node; elem_in_node++){
            int elem_in_node_gid = mesh.elems_in_node(node_gid, elem_in_node);
            
            for (int dim = 0; dim < mesh.num_dims; dim++){
    
                A1(node_gid, dim) += residual_in_elem(elem_in_node_gid, node_gid, dim);

            }// end loop over dim
        }// end loop over elem_in_node
    });// end FOR_ALL
    Kokkos::fence();


}// end assemble residual