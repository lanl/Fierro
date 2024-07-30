#include "ref_elem.h"
#include "mesh.h"
#include "state.h"

void get_artificial_viscosity(CArrayKokkos <double> &sigma_a,
                              const DViewCArrayKokkos <double> &vel,
                              const DViewCArrayKokkos <double> &mat_pt_vel,
                              const DViewCArrayKokkos <double> &den,
                              const DViewCArrayKokkos <double> &sspd,
                              const DViewCArrayKokkos <double> &vol,
                              DViewCArrayKokkos <double> &mat_pt_h,
                              const CArrayKokkos <double> &legendre_jacobian_inverse,
                              const CArrayKokkos <double> &legendre_jacobian,
                              const CArrayKokkos <double> &legendre_jacobian_inverse_t0,
                              const CArrayKokkos <double> &char_length_t0,
                              const mesh_t &mesh,
                              const fe_ref_elem_t &ref_elem,
                              const size_t stage){

        
        // grad^s u = 0.5*( \nabla u + (\nabla u)^T ) //
        CArrayKokkos <double> grad_u(mesh.num_elems*mesh.num_leg_gauss_in_elem, mesh.num_dims, mesh.num_dims, "grad_u");
        CArrayKokkos <double> curl_u(mesh.num_elems*mesh.num_leg_gauss_in_elem, mesh.num_dims, "curl_u");
        // CArrayKokkos <double> l2_curl_u(mesh.num_elems*mesh.num_leg_gauss_in_elem, "l2_curl_u");
        CArrayKokkos <double> l2_grad_u(mesh.num_elems*mesh.num_leg_gauss_in_elem, "l2_grad_u");
        CArrayKokkos <double> Delta_u(mesh.num_elems*mesh.num_leg_gauss_in_elem, "l2_grad_u");

        CArrayKokkos <double> l2_vel_inv(mesh.num_elems*mesh.num_leg_gauss_in_elem, "l2_vel_inv");
        CArrayKokkos <double> l2_JJ0Inv_dot_u(mesh.num_elems*mesh.num_leg_gauss_in_elem, "l2_JJ0Inv_dot_u");
        CArrayKokkos <double> JJ0Inv_dot_u(mesh.num_elems*mesh.num_leg_gauss_in_elem, mesh.num_dims, "JJ0Inv_dot_u");
        CArrayKokkos <double> JJ0Inv(mesh.num_elems*mesh.num_leg_gauss_in_elem, mesh.num_dims, mesh.num_dims, "JJ0Inv");
        

        // Tr(grad^s u) //
        CArrayKokkos <double> div_u(mesh.num_elems*mesh.num_leg_gauss_in_elem, "div_u");

        CArrayKokkos <double> u_star(mesh.num_elems, mesh.num_dims, "u star");


        CArrayKokkos <double> mu1(mesh.num_elems*mesh.num_leg_gauss_in_elem, "mu_1");
        CArrayKokkos <double> mu2(mesh.num_elems*mesh.num_leg_gauss_in_elem, "mu_2");
        CArrayKokkos <double> mu(mesh.num_elems*mesh.num_leg_gauss_in_elem, "mu");

        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            u_star(elem_gid, 0) = 0.0;
            u_star(elem_gid, 1) = 0.0;
            u_star(elem_gid, 2) = 0.0;

            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);

                mu(gauss_gid) = 0.0;
                mu1(gauss_gid) = 0.0;
                mu2(gauss_gid) = 0.0;
                div_u(gauss_gid) = 0.0;
                Delta_u(gauss_gid) = 0.0;
                l2_JJ0Inv_dot_u(gauss_gid) = 0.0;
                l2_vel_inv(gauss_gid) = 0.0;
                // l2_curl_u(gauss_gid) = 0.0;
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
                        grad_u(gauss_gid, i, j) = temp;
                    }// j
                }// i
            }// gauss_lid
        });// elem_gid
        Kokkos::fence();


        // || div u ||_{F}
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

        // fill u star
        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                for(int i = 0; i < mesh.num_dims; i++){
                    u_star(elem_gid, i) += mat_pt_vel(gauss_gid, i);
                }// i
            }// gauss_lid
            u_star(elem_gid, 0) = u_star(elem_gid, 0)/mesh.num_leg_gauss_in_elem;
            u_star(elem_gid, 1) = u_star(elem_gid, 1)/mesh.num_leg_gauss_in_elem;
            u_star(elem_gid, 2) = u_star(elem_gid, 2)/mesh.num_leg_gauss_in_elem;

        });

        // fill div_u
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
                div_u(gauss_gid) = temp;
            }// gauss_lid
        });// elem_gid
        Kokkos::fence();

        // // get char_length_t0 perturbation
        // JJ0Inv
        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                // printf( " JJ0Inv \n" );
                for(int i = 0; i < mesh.num_dims; i++){
                    for (int j = 0; j < mesh.num_dims; j++){
                        double temp = 0.0;
                        for (int k = 0; k < mesh.num_dims; k++){
                            temp += legendre_jacobian_inverse_t0(gauss_gid, k, j)*legendre_jacobian(gauss_gid, k, i);
                        }
                        JJ0Inv(gauss_gid, i, j) = temp;
                        // printf(" %f ", JJ0Inv(gauss_gid, i, j));
                        //printf(" i = %d, j = %d ", i, j);
                    }// j
                    // printf(" \n ");
                }// i
            }// gauss_lid
        });// elem_gid
        Kokkos::fence();

        // JJ0Inv_dot_u
        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                for(int i = 0; i < mesh.num_dims; i++){
                    double temp = 0.0;
                    for (int j = 0; j < mesh.num_dims; j++){
                        temp += JJ0Inv(gauss_gid, i, j)*(mat_pt_vel(gauss_gid, j)-u_star(elem_gid, j));
                    }// j
                    JJ0Inv_dot_u(gauss_gid, i) = temp;
                    //printf(" JJ0Inv_dot_u : %f \n", mat_pt_vel(gauss_gid, i));
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
                    l2_vel_inv(gauss_gid) += (mat_pt_vel(gauss_gid, i) - u_star(elem_gid, i))*(mat_pt_vel(gauss_gid, i) - u_star(elem_gid, i));

                }// i
            }// gauss_lid
        });// elem_gid
        Kokkos::fence();

        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                    
                l2_vel_inv(gauss_gid) = 1.0/( l2_vel_inv(gauss_gid) + 1.0e-6 );

            }// gauss_lid
        });// elem_gid
        Kokkos::fence();

        // fill Delta_u
        CArrayKokkos <double> Du(mesh.num_elems);
        CArrayKokkos <double> Du_min(1);
        RUN({
            Du_min(0) = 1.0e+12;
        });
        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            Du(elem_gid) = 1.0e+12;
        });// elem_gid

        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                double temp1 =0.0;
                for(int i = 0; i < mesh.num_dims; i++){
                    double temp2 = 0.0;
                    for (int j = 0; j < mesh.num_dims; j++){
                        temp2 += grad_u(gauss_gid, i, j)*(mat_pt_vel(gauss_gid, j)-u_star(elem_gid, j));
                    }// j

                    temp1 += (mat_pt_vel(gauss_gid, i)-u_star(elem_gid, i))*temp2;
                }// i
                Delta_u(gauss_gid) = temp1*l2_vel_inv(gauss_gid);
                Du(elem_gid) = Delta_u(gauss_gid) < Du(elem_gid) ? Delta_u(gauss_gid) : Du(elem_gid);
                Du_min(0) = Delta_u(gauss_gid) < Du_min(0) ? Delta_u(gauss_gid) : Du_min(0);
            }// gauss_lid
        });// elem_gid
        Kokkos::fence();

        // fill mu
        FOR_ALL(elem_gid, 0, mesh.num_elems,{

            double h = pow( vol(elem_gid), 1.0/mesh.num_dims )/mesh.Pn;
            
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                
                double h_pert = 0.0;
                                
                h_pert = h*sqrt( l2_vel_inv(gauss_gid)*l2_JJ0Inv_dot_u(gauss_gid) );
                //printf( "h : %f \n", h_pert);
                
                double coeff = 1.0/( 1.0 + exp( Delta_u(gauss_gid) ));//Du(elem_gid) ) );//Du_min(0) ) );//

                double phi_curl = 1.0;

                if (l2_grad_u(gauss_gid) > 0.0){
                    phi_curl = abs(div_u(gauss_gid)/l2_grad_u(gauss_gid));
                }
                
                mu(gauss_gid) = coeff*phi_curl*den(gauss_gid)*h_pert*sspd(gauss_gid);

                mu(gauss_gid) += coeff*2.0*den(gauss_gid)*h_pert*h_pert*pow( abs( Du(elem_gid) ), 1 );//Du_min(0) ), 1 );//Delta_u(gauss_gid) ),  1);

            }// gauss_lid
        });// elem_gid
        Kokkos::fence();

        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                for(int i = 0; i < mesh.num_dims; i++){
                    for(int j = 0; j < mesh.num_dims; j++){
                        sigma_a(stage, gauss_gid, i, j) = mu(gauss_gid)*grad_u(gauss_gid, i, j);
                        //printf("sigma_a : %f \n", sigma_a(stage, gauss_gid, i, j));
                    }// j
                }// i
            }// gauss_lid
        });
        Kokkos::fence();

       //  FOR_ALL(elem_gid, 0, mesh.num_elems,{
       //      for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
       //          int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
       //          for(int i = 0; i < mesh.num_dims; i++){
       //              sigma_a(stage, gauss_gid, i, i) += mu(gauss_gid)*div_u(gauss_gid);
       //          }// i
       //      }// gauss_lid
       //  });
       //  Kokkos::fence();

        

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


// fill curl u
    // FOR_ALL(elem_gid, 0, mesh.num_elems,{
    //     for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
    //         int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
            
    //         curl_u(gauss_gid, 0) = grad_u(gauss_gid, 2, 1) - grad_u(gauss_gid, 1,2);
    //         curl_u(gauss_gid, 1) = grad_u(gauss_gid, 0, 2) - grad_u(gauss_gid, 2,0);
    //         curl_u(gauss_gid, 2) = grad_u(gauss_gid, 1, 0) - grad_u(gauss_gid, 0,1);
                    
    //     }// gauss_lid
    // });// elem_gid
    // Kokkos::fence();

    // FOR_ALL(elem_gid, 0, mesh.num_elems,{
    //     for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
    //         int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
    //         for(int i = 0; i < mesh.num_dims; i++){
    //           l2_curl_u(gauss_gid) += curl_u(gauss_gid,i)*curl_u(gauss_gid,i);
    //         }// i
    //     }// gauss_lid
    // });// elem_gid
    // Kokkos::fence();

    // FOR_ALL(elem_gid, 0, mesh.num_elems,{
    //     for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
    //         int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
    //           l2_curl_u(gauss_gid) = sqrt(l2_curl_u(gauss_gid));
    //     }// gauss_lid
    // });// elem_gid
    // Kokkos::fence();
