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
                    CArrayKokkos <double> &F_dot_ones,
                    const CArrayKokkos <double> &force_tensor,
                    const CArrayKokkos <double> &mass_matrix,
                    const DViewCArrayKokkos <double> &node_vel ){
    
    // Assemble residual in each element //
    FOR_ALL(elem_gid, 0, mesh.num_elems, {


        for (int node_lid1 = 0; node_lid1 < mesh.num_nodes_in_elem; node_lid1++){
            int node_gid1 = mesh.nodes_in_elem(elem_gid, node_lid1);

            for (int node_lid2 = 0; node_lid2 < mesh.num_nodes_in_elem; node_lid2++){
                int node_gid2 = mesh.nodes_in_elem(elem_gid, node_lid2);
                
                for (int dim = 0; dim < mesh.num_dims; dim++){

                    M_dot_u( node_gid1, dim ) += mass_matrix(node_gid1, node_gid2)*( node_vel(stage, node_gid2, dim) - node_vel(0, node_gid2, dim) );

                }// end loop over dim
            }// end loop over node_lid2
            

            for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
                int zone_gid = mesh.zones_in_elem(elem_gid, zone_lid);

                for (int dim = 0; dim < mesh.num_dims; dim++){

                    F_dot_ones(node_gid1, dim) += factor*( force_tensor(stage, node_gid1, zone_gid, dim) + force_tensor(0, node_gid1, zone_gid, dim) );

                }// end loop over dim
            }// end loop over zone_lid

            for (int dim = 0; dim < mesh.num_dims; dim++){

                residual_in_elem(elem_gid, node_gid1, dim) = M_dot_u(elem_gid, node_gid1, dim) + dt*F_dot_ones(elem_gid, node_gid1, dim);

            }
        }// end loop over node_lid1
    } );// end FOR_ALL
    Kokkos::fence();

    // Assemble A1 operator, A1 = sum_{E \ni i} \Phi^E_i(u^k)
    FOR_ALL(node_gid, 0, mesh.num_nodes, {

        int num_elems_attached_to_node = 1;//mesh.num_elems_in_node(node_gid);

        for (int elem_in_node = 0; elem_in_node < num_elems_attached_to_node; elem_in_node++){
            int elem_in_node_gid = mesh.elems_in_node(node_gid, elem_in_node);
            
            for (int dim = 0; dim < mesh.num_dims; dim++){
    
                A1(node_gid, dim) += residual_in_elem(elem_in_node_gid, node_gid, dim);

            }// end loop over dim
        }// end loop over elem_in_node
    });// end FOR_ALL
    Kokkos::fence();


}// end assemble residual