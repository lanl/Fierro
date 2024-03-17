// computes stress tensor at each time stage //
#include "ref_elem.h"
#include "mesh.h"
#include "state.h"


void get_stress_tensor( DViewCArrayKokkos <double> &mat_pt_stress,
                       const size_t stage,
                       const mesh_t &mesh,
                       const DViewCArrayKokkos <double> &mat_pt_pressure){
    
    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for (int legendre_node_lid = 0; legendre_node_lid < mesh.num_leg_gauss_in_elem; legendre_node_lid++){
            int legendre_node = mesh.legendre_in_elem(elem_gid, legendre_node_lid);
            for (int dim = 0; dim < mesh.num_dims; dim++){
                mat_pt_stress(stage, legendre_node, dim, dim) = -1.0*mat_pt_pressure(legendre_node);
            }// end loop over dim
        }
    });// end FOR_ALL over legendre_nodes
    Kokkos::fence();
}// end stress tensor computation


