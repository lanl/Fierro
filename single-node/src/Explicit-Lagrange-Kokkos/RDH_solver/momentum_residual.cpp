#include "matar.h"
#include "mesh.h"
#include "rdh.h"
#include "state.h"
#include "ref_elem.h"
#include "ref_surf_elem.h"
#include <cmath>

#define PI 3.141592653589793

using namespace mtr;

void get_momentum_rhs(const mesh_t &mesh,
                        const fe_ref_elem_t &ref_elem,
                        const fe_ref_surf_t &ref_surf,
                        const DViewCArrayKokkos <double> &SigmaJacInv,
                        const DViewCArrayKokkos <double> &DetJac,
                        DViewCArrayKokkos <double> &F_u,
                        const int stage,
                        bool &viscosity_cond){
        
        // initialize rhs of momentum equation
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
            for (int i = 0; i<mesh.num_nodes_in_elem; i++){
                for (int j = 0; j<mesh.num_dims; j++){
                    F_u(stage, elem_gid, i, j) = 0.0;
                }
            }
        });// FOR_ALL

        // Do volume integral
        FOR_ALL(elem_gid, 0, mesh.num_elems, {
            for (int node_lid = 0; node_lid<mesh.num_nodes_in_elem; node_lid++){
                int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

                for (int dim = 0; dim<mesh.num_dims; dim++){
                    double temp = 0.0;
                    
                    for (int gauss_lid = 0; gauss_lid < ref_elem.num_gauss_leg_in_elem; gauss_lid++){
                        int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);

                        double SigmaJacInvT_dot_GradPhi = SigmaJacInv(stage, gauss_gid, dim, 0)*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 0)
                                                          + SigmaJacInv(stage, gauss_gid, dim, 1)*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 1)
                                                          + SigmaJacInv(stage, gauss_gid, dim, 2)*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 2);

                        Kokkos::atomic_add(&temp, SigmaJacInvT_dot_GradPhi*ref_elem.gauss_leg_weights(gauss_lid)*DetJac(gauss_gid));
            
                    }

                    // printf("momentum rhs : %f \n", -1.0*temp);
                    F_u(stage, elem_gid, node_lid, dim) = -1.0*temp;

                }

            }

        });// FOR_ALL

        // Do surface integral //

}// end get_momentum_rhs

// For each element E, assemble \PHI^E_i(u) = M^E_V \delta^n u + \int_(t^, t^k) \{ \int_E \nabla \varphi_i \cdot \sigma dx - \int_{\partial E} \varphi_i \sigma \cdot n ds \} dt
void get_momentum_residual(const mesh_t &mesh,
                            const DViewCArrayKokkos <double> &M,
                            const DViewCArrayKokkos <double> &node_vel,
                            const DViewCArrayKokkos <double> &F_u,
                            DViewCArrayKokkos <double> &PHI,
                            const double dt,
                            const int stage,
                            const CArrayKokkos <double> &time_int_weights){

        FOR_ALL(elem_gid, 0, mesh.num_elems, {

            for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                
                // Compute M \cdot \delta^n u //
                double M_dot_u[3];
                M_dot_u[0] = 0.0;
                M_dot_u[1] = 0.0;
                M_dot_u[2] = 0.0;

                for (int i = 0; i < mesh.num_nodes_in_elem; i++){
                    int i_gid = mesh.nodes_in_elem(elem_gid, i);

                    Kokkos::atomic_add(&M_dot_u[0], M(elem_gid, node_lid , i )*(node_vel(stage, i_gid, 0) - node_vel(0, i_gid, 0)) );

                    Kokkos::atomic_add(&M_dot_u[1], M(elem_gid, node_lid , i )*(node_vel(stage, i_gid, 1) - node_vel(0, i_gid, 1)) );

                    Kokkos::atomic_add(&M_dot_u[2], M(elem_gid, node_lid , i )*(node_vel(stage, i_gid, 2) - node_vel(0, i_gid, 2)) );
                }
                // for (int dim = 0; dim < mesh.num_dims; dim++){
                //     printf(" M.du in elem %d, node %d at stage %d and dim %d is : %f \n", elem_gid, node_lid, stage, dim, M_dot_u[dim]);
                // }

                // integrate F_u in time //
                double time_int_F_u[3];
                time_int_F_u[0] = 0.0;
                time_int_F_u[1] = 0.0;
                time_int_F_u[2] = 0.0;

                for (int s = 0; s < stage+1; s++){
                    time_int_F_u[0] += dt*time_int_weights(stage, s)*F_u(s, elem_gid, node_lid, 0);
                    time_int_F_u[1] += dt*time_int_weights(stage, s)*F_u(s, elem_gid, node_lid, 1);
                    time_int_F_u[2] += dt*time_int_weights(stage, s)*F_u(s, elem_gid, node_lid, 2);
                }
                // for (int dim = 0;
                // assemble \PHI^{elem_gid}_{node_lid}(u) //
                PHI(stage, elem_gid, node_lid, 0) = M_dot_u[0] - time_int_F_u[0];
                // printf("PHI in elem %d, node %d and dim 0 is : %f \n", elem_gid, node_lid, PHI(stage, elem_gid, node_lid, 0));

                PHI(stage, elem_gid, node_lid, 1) = M_dot_u[1] - time_int_F_u[1];
                // printf("PHI in elem %d, node %d and dim 1 is : %f \n", elem_gid, node_lid, PHI(stage, elem_gid, node_lid, 1));

                PHI(stage, elem_gid, node_lid, 2) = M_dot_u[2] - time_int_F_u[2];
                // printf("PHI in elem %d, node %d and dim 2 is : %f \n", elem_gid, node_lid, PHI(stage, elem_gid, node_lid, 2));

            }
        });                 
}// end get_momentum_residual