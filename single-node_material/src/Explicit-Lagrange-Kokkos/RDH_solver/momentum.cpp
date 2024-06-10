#include "ref_elem.h"
#include "mesh.h"
#include "state.h"


void update_momentum(DViewCArrayKokkos <double> &node_vel,
                     const size_t stage,
                     const mesh_t &mesh,
                     const double dt,
                     const CArrayKokkos <double> &L2,
                     const CArrayKokkos <double> &lumped_mass){
    
  
    FOR_ALL(node_gid, 0, mesh.num_nodes,{

        for (int dim = 0; dim < mesh.num_dims; dim++){
 
            node_vel( 1, node_gid, dim ) = node_vel(stage, node_gid, dim) - dt*L2(stage, node_gid, dim)/lumped_mass(node_gid);

        }
        
    });//end for all
    Kokkos::fence();

    // FOR_ALL(elem_gid, 0, mesh.num_elems,{
    //     for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
    //         int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
    //         for (int dim = 0; dim < mesh.num_dims; dim++){
    //             double interp_temp = 0.0;
    //             for (int dof_lid = 0; dof_lid < mesh.num_nodes_in_elem; dof_lid++){
    //                 int dof_gid = mesh.num_nodes_in_elem(elem_gid, dof_lid);
    //                 interp_temp += ref_elem.gauss_leg_basis(gauss_lid, dof_lid)*node_vel(1, dof_gid, dim);
    //             }
    //             node_vel(1, node_gid, dim) = interp_temp;
    //         }
    //     }
    // });
    // Kokkos::fence();
    

}// end update momentum


void get_grad_vel(){

}


void get_sym_grad_vel(){

}


void get_anti_sym_grad_vel(){
    
}