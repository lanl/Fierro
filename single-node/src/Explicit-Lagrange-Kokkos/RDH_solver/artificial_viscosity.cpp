#include "ref_elem.h"
#include "mesh.h"
#include "state.h"

void get_artificial_viscosity(CArrayKokkos <double> &sigma_a,
                              const DViewCArrayKokkos <double> &vel,
                              const DViewCArrayKokkos <double> &den,
                              const DViewCArrayKokkos <double> &sspd,
                              const DViewCArrayKokkos <double> &vol,
                              const CArrayKokkos <double> &legendre_jacobian_inverse,
                              const mesh_t &mesh,
                              const fe_ref_elem_t &ref_elem,
                              const size_t stage){

        
        // vargrad_u = 0.5*( \nabla u + (\nabla u)^T ) //
        CArrayKokkos <double> grad_u(mesh.num_elems*mesh.num_leg_gauss_in_elem, mesh.num_dims, mesh.num_dims);

        // Tr(\vargrad_u) //
        CArrayKokkos <double> div_u(mesh.num_elems*mesh.num_leg_gauss_in_elem);

        // mu =  \rho \ell ( c_s + \ell |\Delta u| ) //
        CArrayKokkos <double> mu1(mesh.num_elems*mesh.num_leg_gauss_in_elem);
        CArrayKokkos <double> mu2(mesh.num_elems*mesh.num_leg_gauss_in_elem);
        CArrayKokkos <double> mu(mesh.num_elems*mesh.num_leg_gauss_in_elem);

        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);

                mu(gauss_gid) = 0.0;
                div_u(gauss_gid) = 0.0;

                for(int i = 0; i < mesh.num_dims; i++){
                    for (int j = 0; j < mesh.num_dims; j++){
                        
                        sigma_a(stage, gauss_gid, i, j) = 0.0;
                        grad_u(gauss_gid, i, j) = 0.0;

                    }// j
                }// i

            }// gauss_lid
        });// elem_gid
        Kokkos::fence();

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
                        
                        for (int k = 0; k < mesh.num_dims; k++){
                            J_dot_nabla(elem_gid, gauss_lid, dof, i) += legendre_jacobian_inverse(gauss_gid, k, i)
                                                                        *ref_elem.gauss_leg_grad_basis(gauss_lid, dof, k);

                        }// k
                    }
                }
            }
        });
        Kokkos::fence();
        
        // fill grad_u
        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                for(int i = 0; i < mesh.num_dims; i++){
                    for (int j = 0; j < mesh.num_dims; j++){
                        for (int dof = 0; dof < ref_elem.num_basis; dof++){
                            int dof_gid = mesh.nodes_in_elem(elem_gid, dof);
                            
                            grad_u(gauss_gid, i, j) += 0.5*(J_dot_nabla(elem_gid, gauss_lid, dof, i)*vel(stage, dof_gid, j ) 
                                                        + J_dot_nabla(elem_gid, gauss_lid, dof, j)*vel(stage, dof_gid, i ));
                                                        //J_dot_nabla(elem_gid, gauss_lid, dof, i)*vel(stage, dof_gid, j );//
                        }// dof
                    }// j
                }// i
            }// gauss_lid
        });// elem_gid
        Kokkos::fence();

        // fill div_u
        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                for(int i = 0; i < mesh.num_dims; i++){
                    for (int dof = 0; dof < ref_elem.num_basis; dof++){
                        int dof_gid = mesh.nodes_in_elem(elem_gid, dof);
                        // double J_dot_nabla = 0.0;
                        // for (int  k = 0; k < mesh.num_dims; k++){
                        //     J_dot_nabla += legendre_jacobian_inverse(gauss_gid, k, i)
                        //                     *ref_elem.gauss_leg_grad_basis(gauss_lid, dof, k);
                        // }
                        div_u(gauss_gid) += J_dot_nabla(elem_gid, gauss_lid, dof, i)*vel(stage, dof_gid, i);
                    }
                    
                   // div_u(gauss_gid) += grad_u(gauss_gid, i, i);
                    
                }// i
            }// gauss_lid
        });// elem_gid
        Kokkos::fence();

        // fill mu
        FOR_ALL(elem_gid, 0, mesh.num_elems,{

            double h = pow( vol(elem_gid), 1.0/mesh.num_dims);///mesh.num_leg_gauss_in_elem;
            
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                
                double coeff = 1.0;
                if ( h*div_u(gauss_gid) < 0.0 ){
                    coeff = 0.33;
                }
                mu1(gauss_gid) = den(gauss_gid)*h*( coeff*sspd(gauss_gid) + abs(h*div_u(gauss_gid)) );

                mu2(gauss_gid) = den(gauss_gid)*h*h*pow( abs(h*div_u(gauss_gid)), mesh.Pn-1 );

                mu(gauss_gid) = fmin( mu1(gauss_gid), mu2(gauss_gid) );


            }// gauss_lid
        });// elem_gid
        Kokkos::fence();

        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                for(int i = 0; i < mesh.num_dims; i++){
                    for(int j = 0; j < mesh.num_dims; j++){
                        sigma_a(stage, gauss_gid, i, j) = mu(gauss_gid)*grad_u(gauss_gid, i, j);
                    }
                }
            }
        });
        Kokkos::fence();

        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                for(int i = 0; i < mesh.num_dims; i++){
                    sigma_a(stage, gauss_gid, i, i) += mu(gauss_gid)*div_u(gauss_gid);
                }
            }
        });
        Kokkos::fence();

}// end get_artificial_viscosity

void append_artificial_viscosity(DViewCArrayKokkos <double> &sigma,
                                 const CArrayKokkos <double> &sigma_a,
                                 const mesh_t &mesh,
                                 const size_t stage){

    FOR_ALL(elem_gid, 0, mesh.num_elems,{
        for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
            int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);

            for(int i = 0; i < mesh.num_dims; i++){
                for(int j = 0; j < mesh.num_dims; j++){

                    sigma(stage, gauss_gid, i, j) += sigma_a(stage, gauss_gid, i, j);

                }// j
            }// i

        }// gauss_lid
    });// elem_gid
    Kokkos::fence();

}// end append_artificial_viscosity