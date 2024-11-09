#include "matar.h"
#include "mesh.h"
#include "rdh.h"
#include "state.h"
#include "ref_elem.h"
#include "ref_surf_elem.h"
#include <cmath>

#define PI 3.141592653589793

using namespace mtr;

void update_momentum(DViewCArrayKokkos <double> &node_vel,
                     const DViewCArrayKokkos <double> &PHI,
                     const DViewCArrayKokkos <double> &m,
                     const mesh_t &mesh,
                     const int stage){

    FOR_ALL( node_gid, 0, mesh.num_nodes, {

        int num_elems_in_node = mesh.num_corners_in_node(node_gid);

        double L2x = 0.;
        double L2y = 0.;
        double L2z = 0.;

        // Compute \sum_{elem_gid \in node_gid } PHI^{elem_gid}_{node_lid} (u(stage)) //
        for (int elem_lid = 0; elem_lid < num_elems_in_node; elem_lid++){
            int elem_gid = mesh.elems_in_node(node_gid, elem_lid);
            
            int node_lid = mesh.local_node_id_in_elem(node_gid, elem_gid);
                    
            Kokkos::atomic_add(&L2x, PHI(stage, elem_gid, node_lid, 0));
            Kokkos::atomic_add(&L2y, PHI(stage, elem_gid, node_lid, 1));
            Kokkos::atomic_add(&L2z, PHI(stage, elem_gid, node_lid, 2));

        }// elem_lid 

        node_vel(1, node_gid, 0 ) = node_vel(stage, node_gid, 0) - L2x/m(node_gid);
        // printf("vel at dim 0 is : %f \n", node_vel(1, node_gid, 0 ));
        node_vel(1, node_gid, 1 ) = node_vel(stage, node_gid, 1) - L2y/m(node_gid);
        // printf("vel at dim 1 is : %f \n", node_vel(1, node_gid, 1 ));
        node_vel(1, node_gid, 2 ) = node_vel(stage, node_gid, 2) - L2z/m(node_gid);
        // printf("vel at dim 2 is : %f \n", node_vel(1, node_gid, 2 ));


    });

    Kokkos::fence();

}// end update energy