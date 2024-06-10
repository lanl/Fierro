#include "ref_elem.h"
#include "mesh.h"
#include "state.h"

void assemble_Thermo_L2( CArrayKokkos <double> &L2,
                    const size_t stage,
                    const double dt,
                    const mesh_t &mesh,
                    CArrayKokkos <double> &M_dot_e,
                    CArrayKokkos <double> &F_dot_u,
                    const CArrayKokkos <double> &force_tensor,
                    CArrayKokkos <double> &source,
                    const CArrayKokkos <double> &mass_matrix,
                    const DViewCArrayKokkos <double> &zone_sie,
                    const DViewCArrayKokkos <double> &node_vel ){


    // Assemble residual in each element //
    FOR_ALL( zone_gid_1, 0, mesh.num_zones,{
 

        for (int zone_gid_2 = 0; zone_gid_2 < mesh.num_zones; zone_gid_2++){
            M_dot_e(zone_gid_1) += mass_matrix(zone_gid_1, zone_gid_2)*zone_sie(stage, zone_gid_2) 
                                    - mass_matrix(zone_gid_1, zone_gid_2)*zone_sie(0, zone_gid_2);
        }// end loop over node_lid_2
        
        // Compute \int 1.F.sigma dt
        double tempF = 0.0;
        for (int node_gid = 0; node_gid < mesh.num_nodes; node_gid++){
            for (int dim = 0; dim < mesh.num_dims; dim++){
                
                tempF += 0.5*( force_tensor(stage, node_gid, zone_gid_1, dim)
                                + force_tensor(0, node_gid, zone_gid_1, dim) )
                               *0.5*( node_vel(stage, node_gid, dim)  + node_vel(0, node_gid, dim) );
            }
        }// end loop over node_gid
        F_dot_u(zone_gid_1) = tempF;

        double source_term = 0.0;

        source_term = 0.5*( source(stage, zone_gid_1) + source(0, zone_gid_1) );

        L2(zone_gid_1) = M_dot_e(zone_gid_1)/dt - F_dot_u(zone_gid_1) - source_term;
        
        //}// end loop over node_lid_1
    } );// end FOR_ALL
    Kokkos::fence();
}// end assemble T_A1