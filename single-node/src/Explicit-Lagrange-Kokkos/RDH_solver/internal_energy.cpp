#include "ref_elem.h"
#include "mesh.h"
#include "state.h"


void update_internal_energy(DViewCArrayKokkos <double> &zone_sie,
                     const size_t stage,
                     const mesh_t &mesh,
                     const CArrayKokkos <double> &M_e_inv,
                     const CArrayKokkos <double> &force_tensor,
                     CArrayKokkos <double> &F_dot_u,
                     const CArrayKokkos <double> &Fc,
                     CArrayKokkos <double> &Fc_dot_u,
                     const CArrayKokkos <double> &source,
                     const DViewCArrayKokkos <double> &node_vel,
                     const CArrayKokkos <double> &lumped_mass,
                     const double dt){
                    //  const CArrayKokkos <double> &A1,
                    //  const CArrayKokkos <double> &lumped_mass){
    
    
    FOR_ALL(zone_gid_1, 0, mesh.num_zones,{

        for (int node_gid = 0; node_gid < mesh.num_nodes; node_gid++){
            for (int dim = 0; dim < mesh.num_dims; dim++){
                double vel = 0.0;
                if (stage == 0){
                    vel = node_vel(1, node_gid, dim);
                }
                if (stage == 1){
                    vel = 0.5*( node_vel(1, node_gid, dim) + node_vel(0, node_gid, dim) );
                }
                F_dot_u(zone_gid_1) += 0.5*( force_tensor(stage, node_gid, zone_gid_1, dim)
                                              + force_tensor(0, node_gid, zone_gid_1, dim) )
                                       *vel;//0.5*( node_vel(1, node_gid, dim) + node_vel(0, node_gid, dim) );
                // Kokkos::atomic_add(&F_dot_u(zone_gid_1), 0.5*( force_tensor(stage, node_gid, zone_gid_1, dim)
                //                               + force_tensor(0, node_gid, zone_gid_1, dim))*vel);

            }
        }// end loop over zone_lid

        for (int node_gid = 0; node_gid < mesh.num_nodes; node_gid++){
            for (int dim = 0; dim < mesh.num_dims; dim++){
                double vel = 0.0;
                if (stage == 0){
                    vel = node_vel(1, node_gid, dim);
                }
                if (stage == 1){
                    vel = 0.5*( node_vel(1, node_gid, dim) + node_vel(0, node_gid, dim) );
                }
                Fc_dot_u(zone_gid_1) += (1.0/(1.0 + double(stage)))*( Fc(node_gid, zone_gid_1, dim))
                                        *vel;//0.5*( node_vel(1, node_gid, dim) + node_vel(0, node_gid, dim) );
                // Kokkos::atomic_add(&Fc_dot_u(zone_gid_1), (1.0/(1.0 + double(stage)))*( Fc(node_gid, zone_gid_1, dim))*vel);

            }
        }// end loop over zone_lid
    });
    Kokkos::fence();

    FOR_ALL(zone_gid_1, 0, mesh.num_zones,{

        double RHS1 = 0.0;
        double RHS2 = 0.0;
        double RHS = 0.0;

        for (int zone_gid_2 = 0; zone_gid_2 < mesh.num_zones; zone_gid_2++){
            RHS1 += M_e_inv(zone_gid_1, zone_gid_2)*(F_dot_u(zone_gid_2) + 0.0*Fc_dot_u(zone_gid_2));
            // Kokkos::atomic_add(&RHS1, M_e_inv(zone_gid_1, zone_gid_2)*(F_dot_u(zone_gid_2) + 0.0*Fc_dot_u(zone_gid_2)) );
        
            RHS2 += 0.5*M_e_inv(zone_gid_1, zone_gid_2)*
                    (source(stage, zone_gid_2) + source(0, zone_gid_2));
            // Kokkos::atomic_add(&RHS2, 0.5*M_e_inv(zone_gid_1, zone_gid_2)*(source(stage, zone_gid_2) + source(0, zone_gid_2)) );
        }

        RHS = RHS1 + RHS2;
        // printf("RHS %f \n", RHS);

        zone_sie( 1, zone_gid_1 ) = zone_sie( 0, zone_gid_1 ) + dt*RHS;


        if (zone_sie( 1, zone_gid_1 ) <= 0.0){
            // printf("NEGATIVE INTERNAL ENERGY %f \n", zone_sie( 1, zone_gid_1 ));
            printf("Negative sie. Switching to lumped mass \n");// at thermo dof %d \n", zone_gid_1);

            zone_sie(1, zone_gid_1) = 0.0;
            zone_sie(1, zone_gid_1) = zone_sie( 0, zone_gid_1 ) + dt*( F_dot_u(zone_gid_1) + 0.0*Fc_dot_u(zone_gid_1) )/lumped_mass(zone_gid_1);
            zone_sie(1, zone_gid_1) += dt*0.5*(source(stage, zone_gid_1) + source(0, zone_gid_1))/lumped_mass(zone_gid_1);
            // Kokkos::atomic_add(&zone_sie(1, zone_gid_1),dt*0.5*(source(stage, zone_gid_1) + source(0, zone_gid_1))/lumped_mass(zone_gid_1) );


            if (zone_sie( 1, zone_gid_1 ) <= 0.0){

                printf("INTERNAL ENERGY IS NEGATIVE AFTER MASS LUMPING \n");// AND REMOVING CORRECTION \n");
                
                // zone_sie(1, zone_gid_1) = 0.0;
                // zone_sie(1, zone_gid_1) = zone_sie( 0, zone_gid_1 ) + dt*( F_dot_u(zone_gid_1) )/lumped_mass(zone_gid_1);
                // //zone_sie(1, zone_gid_1) += dt*0.5*(source(stage, zone_gid_1) + source(0, zone_gid_1))/lumped_mass(zone_gid_1);
                // if (zone_sie( 1, zone_gid_1 ) <= 0.0){
                //     printf("INTERNAL ENERGY IS NEGATIVE AFTER MASS LUMPING \n");// AND REMOVING CORRECTION \n");
                // }
            }
        }

    });//end for all
    Kokkos::fence();
    

}// end update internal energy
