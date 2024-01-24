#include "ref_elem.h"
#include "mesh.h"
#include "state.h"


void update_momentum(DViewCArrayKokkos <double> &node_vel,
                     const size_t stage,
                     const mesh_t &mesh,
                     const CArrayKokkos <double> &A1,
                     const CArrayKokkos <double> &lumped_mass){
    
    FOR_ALL(node_gid, 0, mesh.num_nodes,{
        for (int dim = 0; dim < mesh.num_dims; dim++){
    
            node_vel(stage+1, node_gid, dim ) = node_vel(stage, node_gid, dim) - A1(node_gid, dim)/lumped_mass(node_gid);

        }
    });//end for all
    

}// end update momentum