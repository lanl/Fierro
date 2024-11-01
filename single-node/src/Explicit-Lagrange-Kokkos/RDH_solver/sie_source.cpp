#include "ref_elem.h"
#include "mesh.h"
#include "state.h"
#include <math.h>


void get_sie_source(CArrayKokkos <double> source,
                    const DViewCArrayKokkos <double> &node_coords,
                    const mat_pt_t &mat_pt,
                    const mesh_t &mesh,
                    const zone_t &zone,
                    const fe_ref_elem_t &ref_elem,
                    const size_t stage
                    ){

    CArrayKokkos <double> temp(mesh.num_elems, mesh.num_zones_in_elem, mesh.num_zones_in_elem);
    FOR_ALL(elem_gid, 0, mesh.num_elems,{

        for (int z_lid = 0; z_lid < mesh.num_zones_in_elem; z_lid++){
            

            for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
                temp(elem_gid, z_lid, zone_lid) = 0.0;
            }
        }
    });
    Kokkos::fence();

    CArrayKokkos <double> zone_coords(mesh.num_elems, mesh.num_zones_in_elem, mesh.num_dims);
    FOR_ALL(elem_gid, 0, mesh.num_elems,{

        for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
            
            for (int dim = 0; dim < mesh.num_dims; dim++){
                zone_coords(elem_gid, zone_lid, dim) = 0.0;
            }
        }
    });
    Kokkos::fence();

    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        // get the coordinates of the zone center
        for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
            int zone_gid = mesh.zones_in_elem(elem_gid, zone_lid);

            for (int dim = 0; dim < mesh.num_dims; dim++){
                for (int node_lid = 0; node_lid < mesh.num_nodes_in_zone; node_lid++){
                    
                    zone_coords(elem_gid, zone_lid, dim) += node_coords(stage, mesh.nodes_in_zone(zone_gid, node_lid), dim);
                    // Kokkos::atomic_add(&zone_coords(elem_gid, zone_lid, dim), node_coords(stage, mesh.nodes_in_zone(zone_gid, node_lid), dim));
                    
                } // end loop over nodes in element


                zone_coords(elem_gid, zone_lid, dim) = zone_coords(elem_gid, zone_lid, dim)/mesh.num_nodes_in_zone;
            
            }
        }
    });
    Kokkos::fence();

    FOR_ALL(elem_gid, 0, mesh.num_elems,{

        for (int z_lid = 0; z_lid < mesh.num_zones_in_elem; z_lid++){
            

            for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){

                for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                    int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);

                    temp(elem_gid, z_lid, zone_lid) += ( (3.0*PI)/8.0 )*mat_pt.gauss_legendre_det_j(gauss_gid)
                                            *ref_elem.gauss_leg_weights(gauss_lid)
                                            *ref_elem.gauss_leg_elem_basis(gauss_lid, z_lid)
                                            *ref_elem.gauss_leg_elem_basis(gauss_lid, zone_lid);
                    
                    // Kokkos::atomic_add(&temp(elem_gid, z_lid, zone_lid),( (3.0*PI)/8.0 )*mat_pt.gauss_legendre_det_j(gauss_gid)
                    //                         *ref_elem.gauss_leg_weights(gauss_lid)
                    //                         *ref_elem.gauss_leg_elem_basis(gauss_lid, z_lid)
                    //                         *ref_elem.gauss_leg_elem_basis(gauss_lid, zone_lid));
                }// gauss_lid 
            }//zone_lid
            
        }// z_lid
    });// end elem_gid
    Kokkos::fence();

    FOR_ALL(elem_gid, 0, mesh.num_elems,{

        for (int z_lid = 0; z_lid < mesh.num_zones_in_elem; z_lid++){
            int zone_gid = mesh.zones_in_elem(elem_gid, z_lid);


            for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
                

                source(stage, zone_gid) += temp(elem_gid, z_lid, zone_lid)*
                                ( cos( 3.0*PI*zone_coords(elem_gid, zone_lid, 0) )*cos( PI*zone_coords(elem_gid, zone_lid, 1) ) 
                                - cos( PI*zone_coords(elem_gid, zone_lid, 0) )*cos( 3.0*PI*zone_coords(elem_gid, zone_lid, 1) ) );
                // Kokkos::atomic_add(&source(stage, zone_gid), temp(elem_gid, z_lid, zone_lid)*
                //                 ( cos( 3.0*PI*zone_coords(elem_gid, zone_lid, 0) )*cos( PI*zone_coords(elem_gid, zone_lid, 1) ) 
                //                 - cos( PI*zone_coords(elem_gid, zone_lid, 0) )*cos( 3.0*PI*zone_coords(elem_gid, zone_lid, 1) ) ));
            }// zone_lid
        }// z_lid

    });// elem_gid
    Kokkos::fence();


}// end get_sie_source