#include "ref_elem.h"
#include "mesh.h"
#include "state.h"

// This implementation uses a basis expansion of the stress field, thus only the polynomial basis is used to build the "force" tensor.
// See commented out section below for a definition which uses a stress tensor defined at every quadrature point.
// F_{ijlk} = \int \psi_j(x) \psi_l(x) \partial_{x_k} \varphi_i(x) dx , \psi(x) is a thermodynamic basis function and \varphi(x) is a kinematic basis funtion.

void build_force_tensor( CArrayKokkos <double> &force_tensor,
                         const size_t stage,
                         const mesh_t &mesh,
                         const CArrayKokkos <double> &grad_kine_basis,
                         const CArrayKokkos <double> &thermo_basis,
                         const CArrayKokkos <double> &legendre_weights,
                         const DViewCArrayKokkos <double> &legendre_jacobian_det,
                         const DViewCArrayKokkos <double> &legendre_jacobian_inverse ){
                    
    FOR_ALL(elem_gid, 0, mesh.num_elems, {
        for (int kinematic_dof = 0; kinematic_dof < mesh.num_nodes_in_elem; kinematic_dof++){
            int global_kinematic_dof = (mesh.num_nodes_in_elem-1)*elem_gid + kinematic_dof;

            for (int thermodynamic_dof_1 = 0; thermodynamic_dof_1 < mesh.num_zones_in_elem; thermodynamic_dof_1++){
                int global_thermodynamic_dof_1 = (mesh.num_zones_in_elem-1)*elem_gid + thermodynamic_dof_1;
                
                for (int thermodynamic_dof_2 = 0; thermodynamic_dof_2 < mesh.num_zones_in_elem; thermodynamic_dof_2++){
                    int global_thermodynamic_dof_2 = (mesh.num_zones_in_elem-1)*elem_gid + thermodynamic_dof_2;

                    for (int dim = 0; dim < mesh.num_dims; dim++){

                        for (int legendre_node = 0; legendre_node < mesh.num_leg_gauss_in_elem; legendre_node++){
                            int global_legendre_node = mesh.legendre_in_elem(elem_gid, legendre_node);

                            for (int contraction_dim = 0; contraction_dim < mesh.num_dims; contraction_dim++){
                                force_tensor( stage, global_kinematic_dof, global_thermodynamic_dof_1, global_thermodynamic_dof_2, dim) +=  legendre_jacobian_inverse(global_legendre_node, dim, contraction_dim)
                                                                                                                * grad_kine_basis(kinematic_dof, legendre_node, contraction_dim)
                                                                                                                * thermo_basis(thermodynamic_dof_1, legendre_node)
                                                                                                                * thermo_basis(thermodynamic_dof_2, legendre_node)
                                                                                                                * legendre_weights(legendre_node)
                                                                                                                * legendre_jacobian_det(global_legendre_node);
                            }// end loop over contraction dimension
                        }// end loop over legendre quadrature nodes
                    }// end loop over dimension
                }// end loop over second set of thermodynamic dofs
            }// end loop over thermodynamic dofs
        }// end loop over kinematic dofs
    });// end FOR_ALL over elements
    Kokkos::fence();
    
}// end assemble force tensor routin



// This funciton assembles the global force tensor F_{i,j} = \int (\sigma : \nabla \varphi_{i} ) \psi_{j} dx,
// where \sigma  is the stress tensor, \varphi are kinematic basis functions and \psi are thermodynamic basis functions.

// void build_force_tensor( CArrayKokkos <double> &force_tensor,
//                          const size_t stage,
//                          const mesh_t &mesh,
//                          const DViewCArrayKokkos <double> &stress_tensor,
//                          const CArrayKokkos <double> &legendre_grad_basis,
//                          const CArrayKokkos <double> &bernstein_basis,
//                          const CArrayKokkos <double> &legendre_weights,
//                          const DViewCArrayKokkos <double> &legendre_jacobian_det,
//                          const DViewCArrayKokkos <double> &legendre_jacobian_inverse ){
                    
//     FOR_ALL(elem_gid, 0, mesh.num_elems, {
//         for (int kinematic_dof = 0; kinematic_dof < mesh.num_nodes_in_elem; kinematic_dof++){
//             int global_kinematic_dof = (mesh.num_nodes_in_elem-1)*elem_gid + kinematic_dof;

//             for (int thermodynamic_dof = 0; thermodynamic_dof < mesh.num_zones_in_elem; thermodynamic_dof++){
//                 int global_thermodynamic_dof = (mesh.num_zones_in_elem-1)*elem_gid + thermodynamic_dof;

//                 for (int dim = 0; dim < mesh.num_dims; dim++){

//                     for (int legendre_node = 0; legendre_node < mesh.num_leg_gauss_in_elem; legendre_node++){
//                         int global_legendre_node = mesh.legendre_in_elem(elem_gid, legendre_node);

//                         for (int contraction_dim1 = 0; contraction_dim1 < mesh.num_dims; contraction_dim1++){
//                             for (int contraction_dim2 = 0; contraction_dim2 < mesh.num_dims; contraction_dim2++){
//                                 force_tensor( stage, global_kinematic_dof, global_thermodynamic_dof, dim) +=  legendre_jacobian_inverse(global_legendre_node, contraction_dim1, contraction_dim1)
//                                                                                                                 * legendre_grad_basis(kinematic_dof, legendre_node, contraction_dim2)
//                                                                                                                 * stress_tensor(stage, global_legendre_node, dim, contraction_dim1)
//                                                                                                                 * bernstein_basis(thermodynamic_dof, legendre_node)
//                                                                                                                 * legendre_weights(legendre_node)
//                                                                                                                 * legendre_jacobian_det(global_legendre_node);
//                             }// end loop over contraction dimension 2
//                         }// end loop over contraction dimension 1
//                     }// end loop over legendre quadrature nodes
//                 }// end loop over dimension
//             }// end loop over thermodynamic dofs
//         }// end loop over kinematic dofs
//     });// end FOR_ALL over elements
//     Kokkos::fence();
    


// }// end assemble force tensor routine