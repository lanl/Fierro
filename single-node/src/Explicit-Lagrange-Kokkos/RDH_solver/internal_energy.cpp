#include "ref_elem.h"
#include "mesh.h"
#include "state.h"


void update_internal_energy(DViewCArrayKokkos <double> &zone_sie,
                     const size_t stage,
                     const mesh_t &mesh,
                     const CArrayKokkos <double> &M_e_inv,
                     const CArrayKokkos <double> &force_tensor,
                     CArrayKokkos <double> &F_dot_u,
                     const CArrayKokkos <double> &source,
                     const DViewCArrayKokkos <double> &node_vel,
                     const double dt){
                    //  const CArrayKokkos <double> &A1,
                    //  const CArrayKokkos <double> &lumped_mass){
    
    
    FOR_ALL(zone_gid_1, 0, mesh.num_zones,{

        for (int node_gid = 0; node_gid < mesh.num_nodes; node_gid++){
            for (int dim = 0; dim < mesh.num_dims; dim++){
                F_dot_u(zone_gid_1) += 0.5*( force_tensor(stage, node_gid, zone_gid_1, dim)
                                              + force_tensor(0, node_gid, zone_gid_1, dim) )
                                       *0.5*( node_vel(1, node_gid, dim) + node_vel(0, node_gid, dim) );
            }
        }// end loop over zone_lid
    });
    Kokkos::fence();

    FOR_ALL(zone_gid_1, 0, mesh.num_zones,{

        double RHS1 = 0.0;
        double RHS2 = 0.0;
        double RHS = 0.0;

        for (int zone_gid_2 = 0; zone_gid_2 < mesh.num_zones; zone_gid_2++){
            RHS1 += M_e_inv(zone_gid_1, zone_gid_2)*F_dot_u(zone_gid_2);
        
            RHS2 += 0.5*M_e_inv(zone_gid_1, zone_gid_2)*
                    (source(stage, zone_gid_2) + source(0, zone_gid_2));
        }

        RHS = RHS1 + RHS2;

        zone_sie( 1, zone_gid_1 ) = zone_sie( 0, zone_gid_1 ) + dt*RHS;

    });//end for all
    Kokkos::fence();
    

}// end update internal energy


void update_internal_energy_DeC(DViewCArrayKokkos <double> &zone_sie,
                     const size_t stage,
                     const mesh_t &mesh,
                     const CArrayKokkos <double> &m,
                     const CArrayKokkos <double> &L_2,
                     const CArrayKokkos <double> &correction,
                     const double dt){
                    //  const CArrayKokkos <double> &A1,
                    //  const CArrayKokkos <double> &lumped_mass){
    
    
    
    FOR_ALL(zone_gid, 0, mesh.num_zones,{

        zone_sie( 1, zone_gid ) = zone_sie( stage, zone_gid ) - dt*L_2(zone_gid)/m(zone_gid) + dt*double(stage)*correction(zone_gid)/m(zone_gid);

        if (zone_sie( 1, zone_gid ) <= 0.0){
            printf("NEGATIVE INTERNAL ENERGY %f \n", zone_sie( 1, zone_gid ));
        }
    });//end for all
    Kokkos::fence();
    

}// end update internal energy