#include "matar.h"
#include "mesh.h"
#include "rdh.h"
#include "state.h"
#include "ref_elem.h"
#include <cmath>

#define PI 3.141592653589793

using namespace mtr;


void eval_sie(const DViewCArrayKokkos <double> sie,
              const int elem_gid,
              const int legendre_lid,
              const mesh_t &mesh,
              const fe_ref_elem_t &ref_elem,
              double &val,
              const int stage){
    
    double temp = 0.0;
    for(int i = 0; i < mesh.num_zones_in_elem; i++){
        int zone_gid = mesh.zones_in_elem(elem_gid, i);
        temp += ref_elem.gauss_leg_elem_basis(legendre_lid, i)*sie(stage, zone_gid);
    }
    val = temp;

}// end eval_sie

void eval_vel(const DViewCArrayKokkos <double> vel,
              const int elem_gid,
              const int legendre_lid,
              const mesh_t &mesh,
              const fe_ref_elem_t &ref_elem,
              CArrayKokkos <double> &val,
              const int stage){

    for(int i = 0; i < mesh.num_nodes_in_elem; i++){
        int node_gid = mesh.nodes_in_elem(elem_gid, i);
        for (int dim = 0; dim < mesh.num_dims; dim++){
            val(dim) += ref_elem.gauss_leg_elem_basis(legendre_lid, i)*vel(stage, node_gid, dim);
        }
    }
}// end eval_vel

// compute: \sum_j \nabla \varphi_i(x) \otimes u_j = \grad u (x)
KOKKOS_FUNCTION 
void eval_grad_u(const mesh_t &mesh,
                 const fe_ref_elem_t &ref_elem,
				 const int elem_gid,
                 const int leg_lid,
				 const int leg_gid,
                 const DViewCArrayKokkos <double> &node_vel,
                 const DViewCArrayKokkos <double> &JacInv,
                 double grad_u[3][3],
                 const int stage){
		
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++){
				grad_u[i][j] = 0.0;
			}// j
		}// i

        // FOR_ALL(node_lid, 0, mesh.num_nodes_in_elem, {
        
        for (int dim = 0; dim < mesh.num_dims; dim++){
                
            for (int d = 0; d < mesh.num_dims; d++){

                for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                
                    double JacInvNabla = JacInv( leg_gid, 0, dim)*ref_elem.gauss_leg_grad_basis(leg_lid, node_lid, 0)
                                        + JacInv( leg_gid, 1, dim)*ref_elem.gauss_leg_grad_basis(leg_lid, node_lid, 1)
                                        + JacInv( leg_gid, 2, dim)*ref_elem.gauss_leg_grad_basis(leg_lid, node_lid, 2);
                    // printf("j^{-1}Nabla : %f \n", JacInvNabla);

                    int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
            
                    // printf("node_vel : %f \n", node_vel(stage, node_gid, d));
                    grad_u[dim][d] += JacInvNabla*node_vel(stage, node_gid, d);
                    // printf("grad u : %f \n", grad_u[dim][d]);

                }// end d

            }// end dim

        }//);// end FOR_ALL

}// end eval_grad_u

