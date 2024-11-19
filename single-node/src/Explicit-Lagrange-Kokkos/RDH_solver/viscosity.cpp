#include "matar.h"
#include "mesh.h"
#include "rdh.h"
#include "state.h"
#include "ref_elem.h"
#include "ref_surf_elem.h"
#include <cmath>
#include "linear_algebra.h"

using namespace mtr;


KOKKOS_INLINE_FUNCTION
void get_viscosity_coefficient(const mesh_t &mesh,
                               const int gauss_gid,
                               double &alpha,
                               const DViewCArrayKokkos <double> &sspd,
                               const DViewCArrayKokkos <double> &den,
                               const DViewCArrayKokkos <double> &J0Inv,
                               const DViewCArrayKokkos <double> &Jac,
                               const DViewCArrayKokkos <double> &h0,
                               const double sym_grad_u[3][3]){
    

    // Compute \| J J0Inv s \|_{\ell^2}
    // compute s (max compression dir)
    double s[3][3];

    double mu[3]; // and min eigenvalue of sym grad u 

	get_spectrum(sym_grad_u, mu, s);

	//for (int i = 0; i < 3; ++i) {
	//    double lambda = mu[i];
	//	  double v[3] = {s[i][0], s[i][1], s[i][2]};
	//    double Av[3] = {0.0, 0.0, 0.0};
	//    // Compute A * v
	//    for (int j = 0; j < 3; ++j) {
	//        Av[j] = sym_grad_u[j][0]*v[0] + sym_grad_u[j][1]*v[1] + sym_grad_u[j][2]*v[2];
	//    }
	//    // Compute lambda * v
	//    double lambda_v[3] = {lambda*v[0], lambda*v[1], lambda*v[2]};
	//    // Compute the difference
	//    double diff = sqrt((Av[0]-lambda_v[0])*(Av[0]-lambda_v[0]) +
	//                       (Av[1]-lambda_v[1])*(Av[1]-lambda_v[1]) +
	//                       (Av[2]-lambda_v[2])*(Av[2]-lambda_v[2]));
	//    if (diff > 1e-6) {
	//        printf("Eigenvector validation failed for eigenvalue %.15f\n", lambda);
	//    } else {
	//        printf("Eigenvector validation passed for eigenvalue %.15f\n", lambda);
	//    }
	//}

	double compression_dir[3];

    compression_dir[0] = 0.0; 
    compression_dir[1] = 0.0; 
    compression_dir[2] = 0.0;

	compression_dir[0] = s[0][0];
	compression_dir[1] = s[0][1];
	compression_dir[2] = s[0][2];

	double min_eig_val = 0.0;
	min_eig_val = mu[0];

    double l2_s;
    compute_l2_norm(compression_dir, l2_s);

    double JJ0Inv[3][3];
    JJ0Inv[0][0] = 0.0; 
    JJ0Inv[0][1] = 0.0; 
    JJ0Inv[0][2] = 0.0;

    JJ0Inv[1][0] = 0.0; 
    JJ0Inv[1][1] = 0.0; 
    JJ0Inv[1][2] = 0.0;

    JJ0Inv[2][0] = 0.0; 
    JJ0Inv[2][1] = 0.0; 
    JJ0Inv[2][2] = 0.0;

    compute_JJ0Inv(Jac, J0Inv, JJ0Inv, gauss_gid, mesh);

    double JJ0Invs[3];
    JJ0Invs[0] = 0; JJ0Invs[1] = 0.0; JJ0Invs[2] = 0.0;
    
    JJ0Invs[0] = JJ0Inv[0][0]*compression_dir[0] + JJ0Inv[0][1]*compression_dir[1] + JJ0Inv[0][2]*compression_dir[2];
    JJ0Invs[1] = JJ0Inv[1][0]*compression_dir[0] + JJ0Inv[1][1]*compression_dir[1] + JJ0Inv[1][2]*compression_dir[2];
    JJ0Invs[2] = JJ0Inv[2][0]*compression_dir[0] + JJ0Inv[2][1]*compression_dir[1] + JJ0Inv[2][2]*compression_dir[2];

    double l2_JJ0Inv_dot_s;
    compute_l2_norm(JJ0Invs, l2_JJ0Inv_dot_s);

    double h = h0(gauss_gid)*l2_JJ0Inv_dot_s/(l2_s + 1.0e-8);

    double coeff = 0.0;
    double eps = 1.0e-12;
    double y = (min_eig_val - eps) / (2.0 * eps);
    if (y < 0.0) 
    { coeff = 0.0; }
    else if (y > 1.0) 
    { coeff = 1.0; }
    else {
        coeff = (3.0 - 2.0 * y) * y * y;
    }

    double div_u = sym_grad_u[0][0] + sym_grad_u[1][1] + sym_grad_u[2][2];
    
	double Fro_norm_grad_u;
	compute_Fro_norm(sym_grad_u, Fro_norm_grad_u);

    double phi_curl = 1.0;

    if (Fro_norm_grad_u > 0.0){
        phi_curl = fabs(div_u/Fro_norm_grad_u);
    }

    alpha = den(gauss_gid)*h*( 0.5*sspd(gauss_gid)*(1.0-coeff)*phi_curl + 2.0*h*fabs( min_eig_val ));
    //if (alpha > 0.0){
	//
	// printf("alpha : %f \n", alpha);
	//
	//}
}


void add_viscosity_momentum(const mesh_t &mesh,
                            const fe_ref_elem_t &ref_elem,
                            const int stage,
                            const DViewCArrayKokkos <double> &node_vel,
                            const DViewCArrayKokkos <double> &sspd,
                            const DViewCArrayKokkos <double> &den,
                            const DViewCArrayKokkos <double> &JacInv,
                            const DViewCArrayKokkos <double> &Jac,
                            const DViewCArrayKokkos <double> &DetJac,
                            const DViewCArrayKokkos <double> &J0Inv,
                            const DViewCArrayKokkos <double> &h0,
                            DViewCArrayKokkos <double> &F){

    FOR_ALL(elem_gid, 0, mesh.num_elems, {
        
        for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
            int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

            for (int gauss_lid = 0; gauss_lid < ref_elem.num_gauss_leg_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);

                double grad_u[3][3];
                double sym_grad_u[3][3];
                
                eval_grad_u(mesh, ref_elem, elem_gid, gauss_lid, gauss_gid, node_vel, JacInv, grad_u, stage);
				symmetrize_matrix(grad_u, sym_grad_u);
				
				double grad_u_dot_Nabla0 = 0.0;
				double grad_u_dot_Nabla1 = 0.0;
				double grad_u_dot_Nabla2 = 0.0;

                grad_u_dot_Nabla0 = sym_grad_u[0][0]*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 0)
                                  + sym_grad_u[0][1]*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 1)
                                  + sym_grad_u[0][2]*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 2);

                grad_u_dot_Nabla1 = sym_grad_u[1][0]*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 0)
                                  + sym_grad_u[1][1]*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 1)
                                  + sym_grad_u[1][2]*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 2);

                grad_u_dot_Nabla2 = sym_grad_u[2][0]*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 0)
                                  + sym_grad_u[2][1]*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 1)
                                  + sym_grad_u[2][2]*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 2);

                double alpha = 0.0;
                
                get_viscosity_coefficient(mesh, gauss_gid, alpha, sspd, den, J0Inv, Jac, h0, sym_grad_u);

                F(stage, elem_gid, node_lid, 0) += alpha
                                                    *grad_u_dot_Nabla0
                                                    *ref_elem.gauss_leg_weights(gauss_lid)
                                                    *DetJac(gauss_gid);

                F(stage, elem_gid, node_lid, 1) += alpha
                                                    *grad_u_dot_Nabla1
                                                    *ref_elem.gauss_leg_weights(gauss_lid)
                                                    *DetJac(gauss_gid);
                
                F(stage, elem_gid, node_lid, 2) += alpha
                                                    *grad_u_dot_Nabla2
                                                    *ref_elem.gauss_leg_weights(gauss_lid)
                                                    *DetJac(gauss_gid);
            }
        }

    });// end FOR_ALL

}


void add_viscosity_energy(const mesh_t &mesh,
                            const fe_ref_elem_t &ref_elem,
                            const int stage,
                            const DViewCArrayKokkos <double> &node_vel,
                            const DViewCArrayKokkos <double> &sspd,
                            const DViewCArrayKokkos <double> &den,
                            const DViewCArrayKokkos <double> &JacInv,
                            const DViewCArrayKokkos <double> &Jac,
                            const DViewCArrayKokkos <double> &DetJac,
                            const DViewCArrayKokkos <double> &J0Inv,
                            const DViewCArrayKokkos <double> &h0,
                            DViewCArrayKokkos <double> &F){

    FOR_ALL(elem_gid, 0, mesh.num_elems, {
        
        for (int zone_lid = 0; zone_lid < mesh.num_zones_in_elem; zone_lid++){

            for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++){
                int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

                for (int gauss_lid = 0; gauss_lid < ref_elem.num_gauss_leg_in_elem; gauss_lid++){
                    int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);

                    double grad_u[3][3];
                    double sym_grad_u[3][3];

                    eval_grad_u(mesh, ref_elem, elem_gid, gauss_lid, gauss_gid, node_vel, JacInv, grad_u, stage);
					symmetrize_matrix(grad_u, sym_grad_u);
					
					double grad_u_dot_Nabla0 = 0.0;
					double grad_u_dot_Nabla1 = 0.0;
					double grad_u_dot_Nabla2 = 0.0;

                    grad_u_dot_Nabla0 = sym_grad_u[0][0]*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 0)
                                      + sym_grad_u[0][1]*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 1)
                                      + sym_grad_u[0][2]*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 2);

                    grad_u_dot_Nabla1 = sym_grad_u[1][0]*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 0)
                                      + sym_grad_u[1][1]*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 1)
                                      + sym_grad_u[1][2]*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 2);

                    grad_u_dot_Nabla2 = sym_grad_u[2][0]*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 0)
                                      + sym_grad_u[2][1]*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 1)
                                      + sym_grad_u[2][2]*ref_elem.gauss_leg_grad_basis(gauss_lid, node_lid, 2);

                    double alpha = 0.0;

                    get_viscosity_coefficient(mesh, gauss_gid, alpha, sspd, den, J0Inv, Jac, h0, sym_grad_u);

                    double grad_u_dot_Nabla_dot_u = grad_u_dot_Nabla0*node_vel(stage, node_gid, 0)
                                                    + grad_u_dot_Nabla1*node_vel(stage, node_gid, 1)
                                                    + grad_u_dot_Nabla2*node_vel(stage, node_gid, 2);

                    F(stage, elem_gid, zone_lid) += alpha
                                                    *grad_u_dot_Nabla_dot_u
                                                    *ref_elem.gauss_leg_elem_basis(gauss_lid, zone_lid)
                                                    *ref_elem.gauss_leg_weights(gauss_lid)
                                                    *DetJac(gauss_gid);

                }
            }
        }

    });// end FOR_ALL

}


