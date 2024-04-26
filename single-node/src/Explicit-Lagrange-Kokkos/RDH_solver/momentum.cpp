#include "ref_elem.h"
#include "mesh.h"
#include "state.h"


void update_momentum(DViewCArrayKokkos <double> &node_vel,
                     const size_t stage,
                     const mesh_t &mesh,
                     const double dt,
                     const CArrayKokkos <double> &A1,
                     const CArrayKokkos <double> &lumped_mass){
    
    // FOR_ALL(elem_gid, 0, mesh.num_elems,{
    //     for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
    //         int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
    //         for (int dim = 0; dim < mesh.num_dims; dim++){
    //             node_vel( stage, node_gid, dim ) = node_vel(stage, node_gid, dim) - A1(node_gid, dim)/lumped_mass(elem_gid, node_lid);
    //         }
    //     }
    // });//end for all
    FOR_ALL(node_gid, 0, mesh.num_nodes,{

        double factor = 1.0/(2.0 - double(stage));
        for (int dim = 0; dim < mesh.num_dims; dim++){
            
            // printf(" lumped mass = %f \n", lumped_mass(node_gid));
            // printf(" A1/ lumped mass = %f \n", A1(node_gid, dim)/lumped_mass(node_gid));

            node_vel( 1, node_gid, dim ) = node_vel(stage, node_gid, dim) - dt*A1(stage, node_gid, dim)/lumped_mass(node_gid);

            // node_vel( 1, node_gid, dim ) = node_vel(0, node_gid, dim) 
            //                             - dt*factor*( A1(stage, node_gid, dim) + A1(0, node_gid, dim) )/lumped_mass(node_gid);

        }
        
    });//end for all
    Kokkos::fence();
    

}// end update momentum