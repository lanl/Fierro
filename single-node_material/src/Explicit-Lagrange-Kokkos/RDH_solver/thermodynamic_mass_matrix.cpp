#include "ref_elem.h"
#include "mesh.h"
#include "state.h"


void assemble_thermodynamic_mass_matrix( CArrayKokkos <double> &M,
                                        CArrayKokkos <double> &m,
                                        CArrayKokkos <double> &M_inv,
                                        const mesh_t &mesh,
                                        const CArrayKokkos <double> &basis,
                                        const CArrayKokkos <double> &legendre_weights,
                                        const CArrayKokkos <double> &legendre_jacobian_det,
                                        const DViewCArrayKokkos <double> &density){

    FOR_ALL(i, 0, mesh.num_zones, 
            j, 0, mesh.num_zones,{
                M(i,j) = 0.0;
                M_inv(i,j) = 0.0;
    });
    Kokkos::fence();

    FOR_ALL(i, 0, mesh.num_zones, {
        m(i) = 0.0;
    });
    Kokkos::fence();
    
    
    FOR_ALL( elem_gid, 0, mesh.num_elems, {

        CArrayKokkos <double> temp(mesh.num_zones_in_elem, mesh.num_zones_in_elem);
        CArrayKokkos <double> temp_inv(mesh.num_zones_in_elem, mesh.num_zones_in_elem);
        for (int i = 0; i < mesh.num_zones_in_elem; i++){
            for (int j = 0; j < mesh.num_zones_in_elem; j++){
                temp(i, j) = 0.0;
                temp_inv(i, j) = 0.0;
            }
        }

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
                
                temp(i,j) = M(global_i, global_j);

                m(global_i) += M(global_i, global_j);

                //M_inv(global_i, global_i) = 1.0/m(global_i);

            }// end loop over j

        }// end loop over i

        invert_matrix(temp_inv, temp, mesh, mesh.num_zones_in_elem);

        for (int i = 0; i < mesh.num_zones_in_elem; i++){
            int global_i = mesh.zones_in_elem(elem_gid, i);
            
            for (int j = 0; j < mesh.num_zones_in_elem; j++){
                int global_j = mesh.zones_in_elem(elem_gid, j);
                M_inv(global_i, global_j) += temp_inv(i,j);
            }
        }


        
       

    }); // end FOR_ALL
    Kokkos::fence();

    
}// end assemble kinematic mass matrix

// CArrayKokkos <double> block_mass(mesh.num_zones_in_elem, mesh.num_zones_in_elem);
        // CArrayKokkos <double> block_mass_inv(mesh.num_zones_in_elem, mesh.num_zones_in_elem);
    
        // for (int i = 0; i < mesh.num_zones_in_elem; i++){
            
        //     for (int j = 0; j < mesh.num_zones_in_elem; j++){
        //         block_mass( i, j) = 0.0;
        //         block_mass_inv( i, j) = 0.0;
        //     }
        // }




        // block_mass( i, j ) += density(legendre_gid)
                    //                 *legendre_weights(legendre_lid)
                    //                 *legendre_jacobian_det(legendre_gid)
                    //                 *basis(legendre_lid, i)
                    //                 *basis(legendre_lid, j);





        // invert_matrix(block_mass_inv, block_mass, mesh, mesh.num_zones_in_elem);

        // for (int i = 0; i < mesh.num_zones_in_elem; i++){
        //     int global_i = mesh.zones_in_elem(elem_gid, i);
            
        //     for (int j = 0; j < mesh.num_zones_in_elem; j++){
        //         int global_j = mesh.zones_in_elem(elem_gid, j);

        //         M_inv(global_i, global_j) = block_mass_inv(i,j);
        //     }
        // }

        
        // CArrayKokkos <double> temp_left( mesh.num_zones_in_elem, mesh.num_zones_in_elem );
        // CArrayKokkos <double> temp_right( mesh.num_zones_in_elem, mesh.num_zones_in_elem );

        // for (int i = 0; i < mesh.num_zones_in_elem; i++){
        //     for (int j = 0; j < mesh.num_zones_in_elem; j++){
        //         temp_left(i,j) = 0.0;
        //         temp_right(i,j) = 0.0;
        //     }//
        // }//
        
        
        // for (int i = 0; i < mesh.num_zones_in_elem; i++){
        //     for (int j = 0; j < mesh.num_zones_in_elem; j++){
        //         for (int k = 0; k < mesh.num_zones_in_elem; k++){
        //             temp_left(i,j) += block_mass_inv(i,k)*block_mass(k,j); 
        //         } 
        //     }
        // }
          
        // for (int i = 0; i < mesh.num_zones_in_elem; i++){
        //     for (int j = 0; j < mesh.num_zones_in_elem; j++){
        //         for (int k = 0; k < mesh.num_zones_in_elem; k++){
        //             temp_right(i,j) += block_mass(i,k)*block_mass_inv(k,j); 
        //         }   
        //     }
        // }

        // printf(" ######################## \n");
        // printf(" ###### Left ##### \n");
        // printf(" ######################## \n");
        // for (int i = 0; i < mesh.num_zones_in_elem; i++){
        //     for (int j = 0; j < mesh.num_zones_in_elem; j++){
        //         printf(" %f ,", temp_left(i,j));
        //     }
        //     printf(" \n");
        // }
        // printf(" ######################## \n");
        // printf(" \n");

        // printf(" ######################## \n");
        // printf(" ###### Right ##### \n");
        // printf(" ######################## \n");
        // for (int i = 0; i < mesh.num_zones_in_elem; i++){
        //     for (int j = 0; j < mesh.num_zones_in_elem; j++){
        //         printf(" %f ,", temp_right(i,j));
        //     }
        //     printf(" \n");
        // }
        // printf(" ######################## \n");
        // printf(" \n");


        // FOR_ALL(i, 0, mesh.num_zones,
    //         j, 0, mesh.num_zones, {

    //     CArrayKokkos <double> temp_left(mesh.num_zones, mesh.num_zones);
    //     CArrayKokkos <double> temp_right(mesh.num_zones, mesh.num_zones);

    //     for (int i = 0; i < mesh.num_zones; i++){
    //         for (int j = 0; j < mesh.num_zones; j++){
    //             temp_left(i,j) = 0.0;
    //             temp_right(i,j) = 0.0;
    //         }
    //     }

    //     for (int i = 0; i < mesh.num_zones; i++){
    //         for (int j = 0; j < mesh.num_zones; j++){
    //             for (int k = 0; k < mesh.num_zones; k++){
    //             temp_left(i,j) += M_inv(i,k)*M(k, j); 
    //             } 
    //         }
    //     }
          
    //     for (int i = 0; i < mesh.num_zones; i++){
    //         for (int j = 0; j < mesh.num_zones; j++){
    //             for (int k = 0; k < mesh.num_zones; k++){
    //             temp_right(i,j) += M(i,k)*M_inv( k, j); 
    //             }   
    //         }
    //     }
        
    //     printf(" ######################## \n");
    //     printf(" ###### Left ##### \n");
    //     printf(" ######################## \n");
    //     for (int i = 0; i < mesh.num_zones_in_elem; i++){
    //         for (int j = 0; j < mesh.num_zones_in_elem; j++){
    //             printf(" %f ,", temp_left(i,j));
    //         }
    //         printf(" \n");
    //     }
    //     printf(" ######################## \n");

    //     printf(" ######################## \n");
    //     printf(" ###### Right ##### \n");
    //     printf(" ######################## \n");
    //     for (int i = 0; i < mesh.num_zones_in_elem; i++){
    //         for (int j = 0; j < mesh.num_zones_in_elem; j++){
    //             printf(" %f ,", temp_right(i,j));
    //         }
    //         printf(" \n");
    //     }
    //     printf(" ######################## \n");

    // });
    
    // RUN({
    //     for (int i =  0; i <  mesh.num_zones; i++){ 
    //         for( int j = 0; j < mesh.num_zones; j++) {
            
    //             printf("%f", M(i, j));
    //         }
    //         printf("\n");
    //     }
    // });

    // FOR_ALL( i, 0, mesh.num_zones, {
    //     if (m(i) <= 0.0){
    //         printf("NEGATIVE thermo lumped mass in zone %d and val = %f \n", i, m(i));
    //         //stop_calc = 1;
    //     }
    // });
    // printf("\n");
    // Kokkos::fence();