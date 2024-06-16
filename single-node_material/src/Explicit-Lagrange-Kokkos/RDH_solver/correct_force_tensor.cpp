#include "state.h"
#include "mesh.h"
#include <chrono>


void correct_force_tensor(CArrayKokkos <double> &force_tensor,
                         const size_t stage,
                         const mesh_t &mesh,
                         const CArrayKokkos <double> &L2,
                         const CArrayKokkos <double> &M,
                         const CArrayKokkos <double> &m,
                         const CArrayKokkos <double> &F_dot_ones){
        
        CArrayKokkos <double> delta_f(mesh.num_nodes, mesh.num_dims);
        CArrayKokkos <double> M_dot_L2(mesh.num_nodes, mesh.num_dims);
        
        
        FOR_ALL(node_gid, 0, mesh.num_nodes,{
            delta_f(node_gid, 0) = 0.0;
            delta_f(node_gid, 1) = 0.0;
            delta_f(node_gid, 2) = 0.0;

            M_dot_L2(node_gid, 0) = 0.0;
            M_dot_L2(node_gid, 1) = 0.0;
            M_dot_L2(node_gid, 2) = 0.0;
        });// FOR_ALL

        // fill M_dot_L2 
        FOR_ALL(node_gid, 0, mesh.num_nodes,{
            for (int n_gid = 0; n_gid < mesh.num_nodes; n_gid++){
                for(int dim = 0; dim < mesh.num_dims; dim++){
                     M_dot_L2(node_gid, dim) += M(node_gid, n_gid)*L2(stage, n_gid, dim)/m(node_gid);
                }// dim
            }
        });// FOR_ALL

        // fill delta_f 
        FOR_ALL(node_gid, 0, mesh.num_nodes,{
            for(int dim = 0; dim < mesh.num_dims; dim++){
               delta_f(node_gid, dim) = 2.0* M_dot_L2(node_gid, dim) - F_dot_ones(node_gid,dim);
               //printf("delta f : %f \n", delta_f(node_gid, dim));
            }// dim
        });// FOR_ALL

        // fill delta_F
        FOR_ALL(node_gid, 0, mesh.num_nodes,
                zone_gid, 0, mesh.num_zones,
                dim, 0, mesh.num_dims, {
                    
                force_tensor(stage, node_gid, zone_gid, dim) += delta_f(node_gid, dim)/mesh.num_zones;
        });// FOR_ALL

    
}// end correct_force_tensor