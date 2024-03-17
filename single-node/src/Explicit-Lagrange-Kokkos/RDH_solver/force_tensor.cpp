#include "ref_elem.h"
#include "mesh.h"
#include "state.h"

// This implementation uses a basis expansion of the stress field, thus only the polynomial basis is used to build the "force" tensor.
// See commented out section below for a definition which uses a stress tensor defined at every quadrature point.
// F_{ijlk} = \int \psi_j(x) \psi_l(x) \partial_{x_k} \varphi_i(x) dx , \psi(x) is a thermodynamic basis function and \varphi(x) is a kinematic basis funtion.

// void build_force_tensor( CArrayKokkos <double> &force_tensor,
//                          const size_t stage,
//                          const mesh_t &mesh,
//                          const CArrayKokkos <double> &grad_kine_basis,
//                          const CArrayKokkos <double> &thermo_basis,
//                          const CArrayKokkos <double> &legendre_weights,
//                          const DViewCArrayKokkos <double> &legendre_jacobian_det,
//                          const DViewCArrayKokkos <double> &legendre_jacobian_inverse ){
                    
//     FOR_ALL(elem_gid, 0, mesh.num_elems, {
//         for (int K_dof = 0; K_dof < mesh.num_nodes_in_elem; K_dof++){
//             int global_kinematic_dof = (mesh.num_nodes_in_elem-1)*elem_gid + K_dof;

//             for (int thermodynamic_dof_1 = 0; thermodynamic_dof_1 < mesh.num_zones_in_elem; thermodynamic_dof_1++){
//                 int global_thermodynamic_dof_1 = (mesh.num_zones_in_elem-1)*elem_gid + thermodynamic_dof_1;
                
//                 for (int thermodynamic_dof_2 = 0; thermodynamic_dof_2 < mesh.num_zones_in_elem; thermodynamic_dof_2++){
//                     int global_thermodynamic_dof_2 = (mesh.num_zones_in_elem-1)*elem_gid + thermodynamic_dof_2;

//                     for (int dim = 0; dim < mesh.num_dims; dim++){

//                         for (int leg_lid = 0; leg_lid < mesh.num_leg_gauss_in_elem; leg_lid++){
//                             int leg_gid = mesh.legendre_in_elem(elem_gid, leg_lid);

//                             for (int contraction_dim = 0; contraction_dim < mesh.num_dims; contraction_dim++){
//                                 force_tensor( stage, global_kinematic_dof, global_thermodynamic_dof_1, global_thermodynamic_dof_2, dim) +=  legendre_jacobian_inverse(leg_gid, dim, contraction_dim)
//                                                                                                                 * grad_kine_basis(K_dof, leg_lid, contraction_dim)
//                                                                                                                 * thermo_basis(thermodynamic_dof_1, leg_lid)
//                                                                                                                 * thermo_basis(thermodynamic_dof_2, leg_lid)
//                                                                                                                 * legendre_weights(leg_lid)
//                                                                                                                 * legendre_jacobian_det(leg_gid);
//                             }// end loop over contraction dimension
//                         }// end loop over legendre quadrature nodes
//                     }// end loop over dimension
//                 }// end loop over second set of thermodynamic dofs
//             }// end loop over thermodynamic dofs
//         }// end loop over kinematic dofs
//     });// end FOR_ALL over elements
//     Kokkos::fence();
    
// }// end assemble force tensor routin



// This funciton assembles the global force tensor F_{i,j} = \int (\sigma : \nabla \varphi_{i} ) \psi_{j} dx,
// where \sigma  is the stress tensor, \varphi are kinematic basis functions and \psi are thermodynamic basis functions.

void build_force_tensor( CArrayKokkos <double> &force_tensor,
                         const size_t stage,
                         const mesh_t &mesh,
                         const DViewCArrayKokkos <double> &stress_tensor,
                         const CArrayKokkos <double> &legendre_grad_basis,
                         const CArrayKokkos <double> &bernstein_basis,
                         const CArrayKokkos <double> &legendre_weights,
                         const CArrayKokkos <double> &legendre_jacobian_det,
                         const CArrayKokkos <double> &legendre_jacobian_inverse ){
                    
    FOR_ALL(elem_gid, 0, mesh.num_elems, {
        for (int K_dof = 0; K_dof < mesh.num_nodes_in_elem; K_dof++){
            int global_K_dof = mesh.nodes_in_elem(elem_gid, K_dof);

            for (int T_dof = 0; T_dof < mesh.num_zones_in_elem; T_dof++){
                int global_T_dof = mesh.zones_in_elem(elem_gid, T_dof);

                for (int dim = 0; dim < mesh.num_dims; dim++){

                    for (int leg_lid = 0; leg_lid < mesh.num_leg_gauss_in_elem; leg_lid++){
                        int leg_gid = mesh.legendre_in_elem(elem_gid, leg_lid);

                        for (int i = 0; i < mesh.num_dims; i++){
                            for (int j = 0; j < mesh.num_dims; j++){
                                force_tensor( stage, global_K_dof, global_T_dof, dim) +=  stress_tensor(stage, leg_gid, dim, i)
                                                                                        * legendre_jacobian_inverse(leg_gid, i, j)
                                                                                        * legendre_grad_basis( leg_lid, K_dof, j)
                                                                                        * bernstein_basis(leg_lid, T_dof)
                                                                                        * legendre_weights(leg_lid)
                                                                                        * legendre_jacobian_det(leg_gid);
                            }// end loop over contraction dimension 2
                        }// end loop over contraction dimension 1
                    }// end loop over legendre quadrature nodes
                }// end loop over dimension
            }// end loop over thermodynamic dofs
        }// end loop over kinematic dofs
    });// end FOR_ALL over elements

    Kokkos::fence();
    
}// end function

// }// end assemble force tensor routine
