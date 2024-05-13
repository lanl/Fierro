#include "ref_elem.h"
#include "mesh.h"
#include "state.h"

void get_artificial_viscosity(CArrayKokkos <double> &sigma_a,
                              const DViewCArrayKokkos <double> &vel,
                              const DViewCArrayKokkos <double> &mat_pt_vel,
                              const DViewCArrayKokkos <double> &den,
                              const DViewCArrayKokkos <double> &sspd,
                              const DViewCArrayKokkos <double> &vol,
                              const CArrayKokkos <double> &legendre_jacobian_inverse,
                              const CArrayKokkos <double> &legendre_jacobian,
                              const CArrayKokkos <double> &legendre_jacobian_inverse_t0,
                              const CArrayKokkos <double> &char_length_t0,
                              const mesh_t &mesh,
                              const fe_ref_elem_t &ref_elem,
                              const size_t stage){

        
        // vargrad_u = 0.5*( \nabla u + (\nabla u)^T ) //
        CArrayKokkos <double> grad_u(mesh.num_elems*mesh.num_leg_gauss_in_elem, mesh.num_dims, mesh.num_dims, "grad_u");
        CArrayKokkos <double> curl_u(mesh.num_elems*mesh.num_leg_gauss_in_elem, mesh.num_dims, "curl_u");
        CArrayKokkos <double> l2_curl_u(mesh.num_elems*mesh.num_leg_gauss_in_elem, "l2_curl_u");
        CArrayKokkos <double> l2_grad_u(mesh.num_elems*mesh.num_leg_gauss_in_elem, "l2_grad_u");

        CArrayKokkos <double> l2_vel_inv(mesh.num_elems*mesh.num_leg_gauss_in_elem, "l2_vel_inv");
        CArrayKokkos <double> l2_JJ0Inv_dot_u(mesh.num_elems*mesh.num_leg_gauss_in_elem, "l2_JJ0Inv_dot_u");
        CArrayKokkos <double> JJ0Inv_dot_u(mesh.num_elems*mesh.num_leg_gauss_in_elem, mesh.num_dims, "JJ0Inv_dot_u");
        CArrayKokkos <double> JJ0Inv(mesh.num_elems*mesh.num_leg_gauss_in_elem, mesh.num_dims, mesh.num_dims, "JJ0Inv");
        

        // Tr(\vargrad_u) //
        CArrayKokkos <double> div_u(mesh.num_elems*mesh.num_leg_gauss_in_elem, "div_u");

        // mu =  \rho \ell ( c_s + \ell |\Delta u| ) //
        CArrayKokkos <double> mu1(mesh.num_elems*mesh.num_leg_gauss_in_elem, "mu_1");
        CArrayKokkos <double> mu2(mesh.num_elems*mesh.num_leg_gauss_in_elem, "mu_2");
        CArrayKokkos <double> mu(mesh.num_elems*mesh.num_leg_gauss_in_elem, "mu");

        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);

                mu(gauss_gid) = 0.0;
                div_u(gauss_gid) = 0.0;
                l2_JJ0Inv_dot_u(gauss_gid) = 0.0;
                l2_vel_inv(gauss_gid) = 0.0;
                l2_curl_u(gauss_gid) = 0.0;
                l2_grad_u(gauss_gid) = 0.0;

                for(int i = 0; i < mesh.num_dims; i++){
                    JJ0Inv_dot_u(gauss_gid, i) = 0.0;
                    for (int j = 0; j < mesh.num_dims; j++){
                        
                        sigma_a(stage, gauss_gid, i, j) = 0.0;
                        grad_u(gauss_gid, i, j) = 0.0;
                        JJ0Inv(gauss_gid, i, j) = 0.0;

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
                            
                            grad_u(gauss_gid, i, j) += 0.5*( J_dot_nabla(elem_gid, gauss_lid, dof, i)*vel(stage, dof_gid, j ) 
                                                        + J_dot_nabla(elem_gid, gauss_lid, dof, j)*vel(stage, dof_gid, i ) );
                                                        //J_dot_nabla(elem_gid, gauss_lid, dof, i)*vel(stage, dof_gid, j );//
                        }// dof
                    }// j
                }// i
            }// gauss_lid
        });// elem_gid
        Kokkos::fence();

        // || nabla^s u ||_{F}
        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                double temp = 0.0;
                for(int i = 0; i < mesh.num_dims; i++){
                    for (int j = 0; j < mesh.num_dims; j++){
                        temp += grad_u(gauss_gid, j, i)*grad_u(gauss_gid, j, i);
                    }// j
                }// i
                l2_grad_u(gauss_gid) = sqrt(temp);
            }// gauss_lid
        });// elem_gid
        Kokkos::fence(); 

        // fill curl u
        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                
                curl_u(gauss_gid, 0) = grad_u(gauss_gid, 2, 1) - grad_u(gauss_gid, 1,2);
                curl_u(gauss_gid, 1) = grad_u(gauss_gid, 0, 2) - grad_u(gauss_gid, 2,0);
                curl_u(gauss_gid, 2) = grad_u(gauss_gid, 1, 0) - grad_u(gauss_gid, 0,1);
                        
            }// gauss_lid
        });// elem_gid
        Kokkos::fence();

        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                for(int i = 0; i < mesh.num_dims; i++){
                  l2_curl_u(gauss_gid) += curl_u(gauss_gid,i)*curl_u(gauss_gid,i);
                }// i
            }// gauss_lid
        });// elem_gid
        Kokkos::fence();

        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                  l2_curl_u(gauss_gid) = sqrt(l2_curl_u(gauss_gid));
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

        // get char_length_t0 perturbation
        // JJ0Inv
        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                for(int i = 0; i < mesh.num_dims; i++){
                    for (int j = 0; j < mesh.num_dims; j++){
                        for (int k = 0; k < mesh.num_dims; k++){
                            JJ0Inv(gauss_gid, i, j) += legendre_jacobian(gauss_gid, i, k)*legendre_jacobian_inverse_t0(gauss_gid, i, k);
                        }
                    }// j
                }// i
            }// gauss_lid
        });// elem_gid
        Kokkos::fence();

        // JJ0Inv_dot_u
        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                for(int i = 0; i < mesh.num_dims; i++){
                    for (int j = 0; j < mesh.num_dims; j++){
                        JJ0Inv_dot_u(gauss_gid, i) += JJ0Inv(gauss_gid, i, j)*mat_pt_vel(gauss_gid, j);
                    }// j
                }// i
            }// gauss_lid
        });// elem_gid
        Kokkos::fence();

        // l2_JJ0Inv_dot_u
        // l2_u
        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                for(int i = 0; i < mesh.num_dims; i++){
                    
                    l2_JJ0Inv_dot_u(gauss_gid) += JJ0Inv_dot_u(gauss_gid, i)*JJ0Inv_dot_u(gauss_gid, i);
                    l2_vel_inv(gauss_gid) += mat_pt_vel(gauss_gid, i)*mat_pt_vel(gauss_gid, i);
                }// i
                // l2_vel_inv(gauss_gid) = 1.0/temp;
            }// gauss_lid
        });// elem_gid
        Kokkos::fence();

        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                for(int i = 0; i < mesh.num_dims; i++){
                    
                    l2_JJ0Inv_dot_u(gauss_gid) = sqrt(l2_JJ0Inv_dot_u(gauss_gid));
                    l2_vel_inv(gauss_gid) = 1.0/sqrt(l2_vel_inv(gauss_gid));
                }// i
                // l2_vel_inv(gauss_gid) = 1.0/temp;
            }// gauss_lid
        });// elem_gid
        Kokkos::fence();

        // fill mu
        FOR_ALL(elem_gid, 0, mesh.num_elems,{

            double h = pow( vol(elem_gid), 1.0/mesh.num_dims );//*l2_JJ0Inv_dot_u(elem_gid)*l2_vel_inv(elem_gid);
            
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);

                //double h = char_length_t0(gauss_gid)*l2_JJ0Inv_dot_u(gauss_gid)*l2_vel_inv(gauss_gid);

                // double delta_u = h*div_u(gauss_gid)*l2_vel_inv(gauss_gid);

                double coeff = 1.0/( 1.0 + exp( 0.50*div_u(gauss_gid) ) );// 0.5*div_u(gauss_gid) ) );//
                //1.0;//0.0;//

                // double coeff2 = 1.0;
                double phi_curl = 1.0;//abs(div_u(gauss_gid)/(l2_grad_u(gauss_gid)+1.0e-06));//fmin(1.0,  fabs( div_u(gauss_gid))/(l2_curl_u(gauss_gid) + 1.0e-15));

                if (l2_grad_u(gauss_gid) > 0.0){
                    phi_curl = abs(div_u(gauss_gid)/l2_grad_u(gauss_gid));
                }

                // if ( div_u(gauss_gid) < 0.0 ){
                //     coeff = 1.0;//1.0/( 1.0 + exp( 100*div_u(gauss_gid) ) );//0.25;//0.33;//0.8;//0.25;//0.3;//1.1;//
                //     //coeff2 = 1.0;//1.0;
                // }

                // mu1(gauss_gid) = den(gauss_gid)*h*( coeff*sspd(gauss_gid) + coeff2*abs(h*div_u(gauss_gid)) );

                // mu2(gauss_gid) = den(gauss_gid)*h*h*pow( abs(h*div_u(gauss_gid)), mesh.Pn-1 );

                // mu(gauss_gid) = fmin( mu1(gauss_gid), mu2(gauss_gid) );

                //mu(gauss_gid) *= phi_curl;

                mu(gauss_gid) = coeff*phi_curl*den(gauss_gid)*h*sspd(gauss_gid);

                mu(gauss_gid) += coeff*2.0*den(gauss_gid)*h*h*abs(div_u(gauss_gid));

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

        // FOR_ALL(elem_gid, 0, mesh.num_elems,{
        //     for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
        //         int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
        //         for(int i = 0; i < mesh.num_dims; i++){
        //             sigma_a(stage, gauss_gid, i, i) += mu(gauss_gid)*div_u(gauss_gid);
        //         }
        //     }
        // });
        // Kokkos::fence();

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