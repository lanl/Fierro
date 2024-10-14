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
            
            double vel_pred = node_vel(stage, node_gid, dim);
            node_vel( 1, node_gid, dim ) = vel_pred - L2(stage, node_gid, dim)/lumped_mass(node_gid);
            // if (dim == 1 ){
            //     printf("L2 is : %f\n ", L2(stage, node_gid, dim));
            // }

        }
        
    });//end for all
    Kokkos::fence();


}// end update momentum


void get_grad_vel(CArrayKokkos <double> &grad_vel,
                              const DViewCArrayKokkos <double> &vel,
                              const CArrayKokkos <double> &legendre_jacobian_inverse,
                              const mesh_t &mesh,
                              const fe_ref_elem_t &ref_elem,
                              const size_t stage){
    CArrayKokkos <double> J_dot_nabla(mesh.num_elems, mesh.num_leg_gauss_in_elem, mesh.num_nodes_in_elem, mesh.num_dims);

    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for(int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
            for (int dof = 0; dof < mesh.num_nodes_in_elem; dof++){
                for (int k = 0; k < mesh.num_dims; k++){
                    J_dot_nabla(elem_gid, gauss_lid, dof, k) = 0.0;
                }
            }
        }
    });
    Kokkos::fence();


    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
            int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
            for(int i = 0; i < mesh.num_dims; i++){
                
                for (int dof = 0; dof < ref_elem.num_basis; dof++){
                    
                    double temp = 0.0;
                    for (int k = 0; k < mesh.num_dims; k++){
                        temp += legendre_jacobian_inverse(gauss_gid, k, i)
                                                                    *ref_elem.gauss_leg_grad_basis(gauss_lid, dof, k);

                    }// k
                    J_dot_nabla(elem_gid, gauss_lid, dof, i) = temp;
                }
            }
        }
    });
    Kokkos::fence();
    
    // fill grad^s_u
    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
            int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
            for(int i = 0; i < mesh.num_dims; i++){
                for (int j = 0; j < mesh.num_dims; j++){
                    double temp = 0.0;
                    for (int dof = 0; dof < ref_elem.num_basis; dof++){
                        int dof_gid = mesh.nodes_in_elem(elem_gid, dof);
                        
                        temp = J_dot_nabla(elem_gid, gauss_lid, dof, i)*vel(stage, dof_gid, j );

                    }// dof
                    grad_vel(gauss_gid, i, j) = temp;
                }// j
            }// i
        }// gauss_lid
    });// elem_gid
    Kokkos::fence();
}


void get_sym_grad_vel(CArrayKokkos <double> &sym_grad_vel,
                              const DViewCArrayKokkos <double> &vel,
                              const CArrayKokkos <double> &legendre_jacobian_inverse,
                              const mesh_t &mesh,
                              const fe_ref_elem_t &ref_elem,
                              const size_t stage){
        
    CArrayKokkos <double> J_dot_nabla(mesh.num_elems, mesh.num_leg_gauss_in_elem, mesh.num_nodes_in_elem, mesh.num_dims);

    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for(int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
            for (int dof = 0; dof < mesh.num_nodes_in_elem; dof++){
                for (int k = 0; k < mesh.num_dims; k++){
                    J_dot_nabla(elem_gid, gauss_lid, dof, k) = 0.0;
                }
            }
        }
    });
    Kokkos::fence();


    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
            int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
            for(int i = 0; i < mesh.num_dims; i++){
                
                for (int dof = 0; dof < ref_elem.num_basis; dof++){
                    
                    double temp = 0.0;
                    for (int k = 0; k < mesh.num_dims; k++){
                        temp += legendre_jacobian_inverse(gauss_gid, k, i)
                                                                    *ref_elem.gauss_leg_grad_basis(gauss_lid, dof, k);

                    }// k
                    J_dot_nabla(elem_gid, gauss_lid, dof, i) = temp;
                }
            }
        }
    });
    Kokkos::fence();
    
    // fill grad^s_u
    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
            int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
            for(int i = 0; i < mesh.num_dims; i++){
                for (int j = 0; j < mesh.num_dims; j++){
                    double temp = 0.0;
                    for (int dof = 0; dof < ref_elem.num_basis; dof++){
                        int dof_gid = mesh.nodes_in_elem(elem_gid, dof);
                        
                        temp += 0.5*( J_dot_nabla(elem_gid, gauss_lid, dof, i)*vel(stage, dof_gid, j ) );
                        temp += 0.5*( J_dot_nabla(elem_gid, gauss_lid, dof, j)*vel(stage, dof_gid, i ) );

                    }// dof
                    sym_grad_vel(gauss_gid, i, j) = temp;
                }// j
            }// i
        }// gauss_lid
    });// elem_gid
    Kokkos::fence();
}


void get_anti_sym_grad_vel(CArrayKokkos <double> &anti_sym_grad_vel,
                              const DViewCArrayKokkos <double> &vel,
                              const CArrayKokkos <double> &legendre_jacobian_inverse,
                              const mesh_t &mesh,
                              const fe_ref_elem_t &ref_elem,
                              const size_t stage){
    CArrayKokkos <double> J_dot_nabla(mesh.num_elems, mesh.num_leg_gauss_in_elem, mesh.num_nodes_in_elem, mesh.num_dims);

    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for(int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
            for (int dof = 0; dof < mesh.num_nodes_in_elem; dof++){
                for (int k = 0; k < mesh.num_dims; k++){
                    J_dot_nabla(elem_gid, gauss_lid, dof, k) = 0.0;
                }
            }
        }
    });
    Kokkos::fence();


    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
            int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
            for(int i = 0; i < mesh.num_dims; i++){
                
                for (int dof = 0; dof < ref_elem.num_basis; dof++){
                    
                    double temp = 0.0;
                    for (int k = 0; k < mesh.num_dims; k++){
                        temp += legendre_jacobian_inverse(gauss_gid, k, i)
                                                                    *ref_elem.gauss_leg_grad_basis(gauss_lid, dof, k);

                    }// k
                    J_dot_nabla(elem_gid, gauss_lid, dof, i) = temp;
                }
            }
        }
    });
    Kokkos::fence();
    
    // fill grad^s_u
    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
            int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
            for(int i = 0; i < mesh.num_dims; i++){
                for (int j = 0; j < mesh.num_dims; j++){
                    double temp = 0.0;
                    for (int dof = 0; dof < ref_elem.num_basis; dof++){
                        int dof_gid = mesh.nodes_in_elem(elem_gid, dof);
                        
                        temp += 0.5*( J_dot_nabla(elem_gid, gauss_lid, dof, i)*vel(stage, dof_gid, j ) );
                        temp += 0.5*( J_dot_nabla(elem_gid, gauss_lid, dof, j)*vel(stage, dof_gid, i ) );

                    }// dof
                    anti_sym_grad_vel(gauss_gid, i, j) = temp;
                }// j
            }// i
        }// gauss_lid
    });// elem_gid
    Kokkos::fence();
}


void get_div_vel(CArrayKokkos <double> &div_vel,
                              const DViewCArrayKokkos <double> &vel,
                              const CArrayKokkos <double> &legendre_jacobian_inverse,
                              const mesh_t &mesh,
                              const fe_ref_elem_t &ref_elem,
                              const size_t stage){
    CArrayKokkos <double> J_dot_nabla(mesh.num_elems, mesh.num_leg_gauss_in_elem, mesh.num_nodes_in_elem, mesh.num_dims);

    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for(int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
            for (int dof = 0; dof < mesh.num_nodes_in_elem; dof++){
                for (int k = 0; k < mesh.num_dims; k++){
                    J_dot_nabla(elem_gid, gauss_lid, dof, k) = 0.0;
                }
            }
        }
    });
    Kokkos::fence();


    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
            int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
            for(int i = 0; i < mesh.num_dims; i++){
                
                for (int dof = 0; dof < ref_elem.num_basis; dof++){
                    
                    double temp = 0.0;
                    for (int k = 0; k < mesh.num_dims; k++){
                        temp += legendre_jacobian_inverse(gauss_gid, k, i)
                                *ref_elem.gauss_leg_grad_basis(gauss_lid, dof, k);

                    }// k
                    J_dot_nabla(elem_gid, gauss_lid, dof, i) = temp;
                }
            }
        }
    });
    Kokkos::fence();
    
    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
            int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
            double temp = 0.0;
            for(int i = 0; i < mesh.num_dims; i++){
                for (int dof = 0; dof < ref_elem.num_basis; dof++){
                    int dof_gid = mesh.nodes_in_elem(elem_gid, dof);

                    temp += J_dot_nabla(elem_gid, gauss_lid, dof, i)*vel(stage, dof_gid, i);
                }
            }// i
            div_vel(gauss_gid) = temp;
        }// gauss_lid
    });// elem_gid
    Kokkos::fence();
}