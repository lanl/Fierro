// computes stress tensor at each time stage //
#include "ref_elem.h"
#include "mesh.h"
#include "state.h"


void get_stress_tensor(CArrayKokkos <double> &elem_stress,
                       const size_t stage,
                       const mesh_t &mesh,
                       const DViewCArrayKokkos <double> &elem_pressure){
    
    // FOR_ALL(elem_gid, 0, mesh.num_elems,{
    //     for (int legendre_node_lid = 0; legendre_node_lid < mesh.num_leg_gauss_in_elem; legendre_node_lid++){
    //         int legendre_node = mesh.legendre_in_elem(elem_gid, legendre_node_lid);
    //         for (int dim = 0; dim < mesh.num_dims; dim++){
    //             elem_stress(stage, legendre_node, dim, dim) = -1.0*elem_pressure(legendre_node);
    //         }// end loop over dim
    //     }
    // });// end FOR_ALL over legendre_nodes
    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){
            int zone_gid = mesh.zones_in_elem(elem_gid, zone_lid);
            for (int dim = 0; dim < mesh.num_dims; dim++){
                elem_stress(stage, zone_gid, dim, dim) = -1.0*elem_pressure(zone_gid);
            }// end loop over dim
        }
    });// end FOR_ALL over legendre_nodes
    Kokkos::fence();
}// end stress tensor computation


