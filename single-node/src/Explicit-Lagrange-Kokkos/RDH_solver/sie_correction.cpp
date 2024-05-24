#include "ref_elem.h"
#include "mesh.h"
#include "state.h"
#include <math.h>


void get_sie_correction(CArrayKokkos <double> correction,
                        const CArrayKokkos <double> &ke_correction,
                        const CArrayKokkos <double> &M_dot_e,
                        const mesh_t &mesh){
                    // const DViewCArrayKokkos <double> &node_vel,
                    // const CArrayKokkos <double> &F,
                    // const CArrayKokkos <double> &vel_tilde,
                    // const CArrayKokkos <double> &ke_correction,
                    // const CArrayKokkos <double> &M_dot_e,
                    // const mesh_t &mesh){
                        
    FOR_ALL(zone_gid, 0 , mesh.num_zones, {
        correction(zone_gid) = 0.0;
    });// node_gid

    FOR_ALL(zone_gid, 0, mesh.num_zones,{


        // for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
        //     int zone_gid = mesh.zones_in_elem(elem_gid, zone_lid);

            // double temp = 0.0;
            
            // for (int node_gid = 0; node_gid < mesh.num_nodes; node_gid++){
            //     for (int dim = 0; dim < mesh.num_dims; dim++){
                
            //         temp += 0.5*F(0, node_gid, zone_gid, dim)*(node_vel(1, node_gid, dim));// - vel_tilde(node_gid, dim));

            //     }//dim
            // }// node_lid 

            // // if( 1.0e-06 < fabs(temp)){
                // printf("correction %f \n", temp);
            // // }// if

            // correction(zone_gid) = temp;

            correction(zone_gid) += ke_correction(zone_gid)/mesh.num_zones;
            
            correction(zone_gid) -= M_dot_e(zone_gid);

            // if ( 1.0e-07 < fabs(correction(zone_gid)) ){
            //     printf("correction %f \n", correction(zone_gid));
            // }//if
            
        //}//zone_lid
            
    });// end elem_gid
    Kokkos::fence();

    
}// end get_sie_source