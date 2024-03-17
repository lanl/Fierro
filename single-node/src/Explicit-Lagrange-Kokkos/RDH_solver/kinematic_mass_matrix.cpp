#include "ref_elem.h"
#include "mesh.h"
#include "state.h"

void assemble_kinematic_mass_matrix(CArrayKokkos <double> &M_V,
                                    CArrayKokkos <double> &lumped_mass,
                                    const mesh_t &mesh,
                                    const CArrayKokkos <double> &basis,
                                    const CArrayKokkos <double> &legendre_weights,
                                    const CArrayKokkos <double> &legendre_jacobian_det,
                                    const DViewCArrayKokkos <double> &density){

    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for (int i = 0; i < mesh.num_nodes_in_elem; i++){
            int global_i = mesh.nodes_in_elem(elem_gid, i);
            
            for (int j = 0; j < mesh.num_nodes_in_elem; j++){
                int global_j = mesh.nodes_in_elem(elem_gid, j);
                
                for (int legendre_lid = 0; legendre_lid < mesh.num_leg_gauss_in_elem; legendre_lid++){
                    int legendre_gid = mesh.legendre_in_elem(elem_gid, legendre_lid);

                    M_V( global_i, global_j) += density(legendre_gid)*legendre_weights(legendre_lid)
                                               *legendre_jacobian_det(legendre_gid)
                                               *basis(legendre_lid, i)
                                               *basis(legendre_lid, j); 

                }// end loop over legendre_lid

                // compute lumped mass as row sum of mass matrix
                lumped_mass(global_i) += M_V(global_i, global_j);

            }// end loop over j
        }// end loop over i
    });// end FOR_ALL
     Kokkos::fence();
}// end assemble kinematic mass matrix


void compute_lumped_mass(CArrayKokkos <double> &lumped_mass,
                         const mesh_t &mesh,
                         const CArrayKokkos <double> &basis,
                         const CArrayKokkos <double> &legendre_weights,
                         const CArrayKokkos <double> &legendre_jacobian_det,
                         const DViewCArrayKokkos <double> &density){
    
    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
            int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
            int num_elems_around_node = mesh.num_corners_in_node(node_gid);
            for (int elem_in_node_lid = 0; elem_in_node_lid < num_elems_around_node; elem_in_node_lid++){
                int elem_in_node_gid = mesh.corners_in_node(node_gid, elem_in_node_lid);
                for (int legendre_lid = 0; legendre_lid < mesh.num_leg_gauss_in_elem; legendre_lid++){
                    int legendre_gid = mesh.legendre_in_elem(elem_in_node_gid, legendre_lid);
                    lumped_mass(node_gid) += density(legendre_gid)*legendre_weights(legendre_lid)
                                                *legendre_jacobian_det(legendre_gid)
                                                *basis(legendre_lid, node_lid);
                }// end loop over legendre_lid
            }// end loop over elem_in_node_lid
        }// end loop over node_lid
    });// end FOR_ALL
    Kokkos::fence(); 

}// end compute lumped mass
