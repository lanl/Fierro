#include "ref_elem.h"
#include "mesh.h"
#include "state.h"


void assemble_thermodynamic_mass_matrix( CArrayKokkos <double> &M,
                                    const mesh_t &mesh,
                                    const CArrayKokkos <double> &basis,
                                    const CArrayKokkos <double> &legendre_weights,
                                    const DViewCArrayKokkos <double> &legendre_jacobian_det,
                                    const DViewCArrayKokkos <double> &density ){

    FOR_ALL( elem_gid, 0, mesh.num_elems, {
        for (int i = 0; i < mesh.num_zones_in_elem; i++){
            int global_i = mesh.zones_in_elem(elem_gid, i);
            
            for (int j = 0; j < mesh.num_zones_in_elem; j++){
                int global_j = mesh.zones_in_elem(elem_gid, j);
                
                for (int legendre_lid = 0; legendre_lid < mesh.num_leg_gauss_in_elem; legendre_lid++){
                    int legendre_gid = mesh.legendre_in_elem(elem_gid, legendre_lid);

                    M(global_i, global_j) += density(legendre_gid)
                                               *legendre_weights(legendre_lid)
                                               *legendre_jacobian_det(legendre_gid)
                                               *basis(i, legendre_lid)
                                               *basis(j, legendre_lid); 
                                               
                }// end loop over legendre_lid
            }// end loop over j
        }// end loop over i
    }); // end FOR_ALL
}// end assemble kinematic mass matrix