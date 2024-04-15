#include "ref_elem.h"
#include "mesh.h"
#include "state.h"


void assemble_thermodynamic_mass_matrix( CArrayKokkos <double> &M,
                                        CArrayKokkos <double> &m,
                                        const mesh_t &mesh,
                                        const CArrayKokkos <double> &basis,
                                        const CArrayKokkos <double> &legendre_weights,
                                        const CArrayKokkos <double> &legendre_jacobian_det,
                                        const DViewCArrayKokkos <double> &density ){

    FOR_ALL(i, 0, mesh.num_zones, 
            j, 0, mesh.num_zones,{
                M(i,j) = 0.0;
    });
    Kokkos::fence();

    FOR_ALL(i, 0, mesh.num_zones, {
        m(i) = 0.0;
    });
    Kokkos::fence();

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
                                               *basis(legendre_lid, i)
                                               *basis(legendre_lid, j); 
                                               
                }// end loop over legendre_lid
                //printf("thermo mass = %f at Tdofs %d, %d \n", M(global_i, global_j), global_i, global_j);

                m(global_i) += M(global_i, global_j);

            }// end loop over j

        }// end loop over i
    }); // end FOR_ALL
    Kokkos::fence();
    
    // RUN({
    //     for (int i =  0; i <  mesh.num_zones; i++){ 
    //         for( int j = 0; j < mesh.num_zones; j++) {
            
    //             printf("%f", M(i, j));
    //         }
    //         printf("\n");
    //     }
    // });

    FOR_ALL( i, 0, mesh.num_zones, {
        if (m(i) <= 0.0){
            printf("NEGATIVE thermo lumped mass in zone %d and val = %f \n", i, m(i));
            //stop_calc = 1;
        }
    });
    printf("\n");
    Kokkos::fence();
}// end assemble kinematic mass matrix