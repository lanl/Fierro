#include "matar.h"
#include "mesh.h"
#include "rdh.h"
#include "state.h"
#include "ref_elem.h"
#include <cmath>

#define PI 3.141592653589793

using namespace mtr;


// compute lumped masses for kinematic and thermodynamic vars
void get_lumped_mass(mesh_t &mesh,
                     fe_ref_elem_t & ref_elem,
                     const DViewCArrayKokkos <double> &DetJac,
                     const DViewCArrayKokkos <double> &den,
                     const DViewCArrayKokkos <double> &M_u,
                     const DViewCArrayKokkos <double> &M_e,
                     DViewCArrayKokkos <double> & nodal_mass,
                     DViewCArrayKokkos <double> & zonal_mass){
    
    
    // compute \sum_{K \ni i} \int_K \rho \varphi_i dx
    FOR_ALL( node_gid, 0, mesh.num_nodes, {

        int num_elems_in_node = mesh.num_corners_in_node(node_gid);

        for (int elem_lid = 0; elem_lid < num_elems_in_node; elem_lid++){
            int elem_gid = mesh.elems_in_node(node_gid, elem_lid);
            
            int node_lid = mesh.local_node_id_in_elem(node_gid, elem_gid);
                    
            for (int gauss_lid = 0; gauss_lid < ref_elem.num_gauss_leg_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                Kokkos::atomic_add(&nodal_mass(node_gid), 
                                    den(gauss_gid)*DetJac(gauss_gid)
                                    *ref_elem.gauss_leg_basis(gauss_lid, node_lid)
                                    *ref_elem.gauss_leg_weights(gauss_lid) );
            }// gauss_lid

        }// elem_lid 
    });

    Kokkos::fence();

    //FOR_ALL(elem_gid, 0, mesh.num_elems,{
    // 
    //    for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
    //        int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

    //        for (int i = 0; i < mesh.num_nodes_in_elem; i++){
    //            nodal_mass(node_gid) += M_u(elem_gid, node_lid, i);
    //        }
    //    }
    //});

    // Kokkos::fence();

    // compute \sum_{K \ni i} \int_K \rho \psi_i dx
    // unlike the above, there is only one K per i as the i's are completely internal to the K's
    // i.e. no thermo degree of freedom is shared by two or more elements
    
	//FOR_ALL(zone_gid, 0, mesh.num_zones,{
    //    int elem_gid = mesh.elems_in_zone(zone_gid);
    //    int zone_lid = mesh.local_zone_id_in_elem(zone_gid, elem_gid);

    //    for (int gauss_lid = 0; gauss_lid < ref_elem.num_gauss_leg_in_elem; gauss_lid++){
    //        int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);

    //        zonal_mass(zone_gid) += den(gauss_gid)*DetJac(gauss_gid)
    //                                        *ref_elem.gauss_leg_elem_basis(gauss_lid, zone_lid)
    //                                        *ref_elem.gauss_leg_weights(gauss_lid);

    //    }
    //});

    //Kokkos::fence();

    FOR_ALL(elem_gid, 0, mesh.num_elems,{
     
        for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
            int zone_gid = mesh.zones_in_elem(elem_gid, zone_lid);

            for (int i = 0; i < mesh.num_zones_in_elem; i++){
                zonal_mass(zone_gid) += M_e(elem_gid, zone_lid, i);
            }
        }
    });

    Kokkos::fence();


    // Check positivity of lumped masses
    FOR_ALL( node_gid, 0, mesh.num_nodes, {
        if (nodal_mass(node_gid) < 0.0){
            printf(" negative lumped kinematic mass at node id %d with val: %f \n", node_gid, nodal_mass(node_gid));
        }
    });

    Kokkos::fence();

    FOR_ALL(zone_gid, 0, mesh.num_zones,{
        if (zonal_mass(zone_gid) < 0.0){
            printf(" negative lumped thermo mass at zone id %d with val: %f \n", zone_gid, zonal_mass(zone_gid));
        }
    });

    Kokkos::fence();


}// end get_lumped_mass
