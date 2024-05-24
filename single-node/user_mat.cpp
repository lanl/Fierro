// -----------------------------------------------------------------------------
// This code contains the constitutive relation for a user supplied model
//------------------------------------------------------------------------------
#include "state.h"
#include "mesh.h"


const double fuzz = 1.e-16;


// -----------------------------------------------------------------------------
// This is the user material model function for the equation of state
// An eos function must be supplied or the code will fail to run.
// The pressure and sound speed can be calculated from an analytic eos.
// The pressure can also be calculated using p = -1/3 Trace(Stress)
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
void user_eos_model(const DViewCArrayKokkos <double> &elem_pres,
                    const DViewCArrayKokkos <double> &node_vel,
                    const size_t elem_gid,
                    const size_t mat_id,
                    const DViewCArrayKokkos <double> &elem_state_vars,
                    const DViewCArrayKokkos <double> &elem_sspd,
                    const double den,
                    const double sie){
    
    // -----------------------------------------------------------------------------
    // Required variables are here
    //------------------------------------------------------------------------------
    
    // Mie-Grüneisen eos for metal from Liu et al 2022
    // pres = Ph + gamma*den*(e-Eh)
    // Ph = (den_initial*C0^2*eta)/(1-s*eta)^2
    // Eh = (eta*Ph)/(2*den_init)
    // eta = 1 - den_init/den

    // ---eos
    // statev(0) = gamma
    // statev(1) = c0
    // statev(2) = s0
    // statev(3) = den_ref
    // statev(4) = csmin
    
    // --strength
    // statev(5) = shear modulus
    // statev(6) = yield strength

    
    double gamma = elem_state_vars(elem_gid,0);
    double C0 = elem_state_vars(elem_gid,1);
    double S0 = elem_state_vars(elem_gid,2);
    double den_ref = elem_state_vars(elem_gid,3);
    double csmin = elem_state_vars(elem_gid,4);
    double G = elem_state_vars(elem_gid,5); // shear modulus
    

    
    // calculate the pressure components
    double eta = 1.0 - (den_ref/den);
    double Ph = (den_ref*C0*C0*eta)/((1.0 - S0*eta)*(1.0 - S0*eta));
    double Eh = (eta*Ph)/(2.0*den_ref);

    // full pressure calculation
    // from Liu et al 2022
    elem_pres(elem_gid) = Ph + gamma*den*(sie - Eh);
    
    
    // now calculate the sound speed
    // variables for derivatives and easier readability
    double Vo = 1.0/den_ref;
    double V = 1.0/den;
    double easy = 1.0/(Vo - elem_state_vars(elem_gid,2)*(Vo - V));

    // derivative of Huginot pressure/density with respect to energy
    double dPhdr = (C0*C0)/(den*den) *
                    ( (easy*easy) + 2.0*S0*(Vo - V)*(easy*easy*easy) );
    double dgdr = gamma; //gamma_0 + gamma_1*(elem_state_vars(elem_gid,0)*den_ref)*den*(-1.0/(den*den));
    double dEdr = 0.5 * (dPhdr * (Vo - V) + ((elem_pres(elem_gid) + Ph)/(den*den)));

    // c_EOS is the sound speed squared
    double c_EOS = dPhdr + dgdr*(sie - Eh) + gamma*den*dEdr + (elem_pres(elem_gid)/(den*den))*gamma*den;
    c_EOS = fmax(0.0, c_EOS); //check that sound speed it greater than zero
  
    double B = den*c_EOS*c_EOS; // Bulk Modulus
    
    // final calculation for sound speed for the element
    elem_sspd(elem_gid) = sqrt((B + (4.0/3.0)*G)/den);  // this is adding shear modulus to sound speed
    
    
    if (elem_sspd(elem_gid) < csmin) {
        elem_sspd(elem_gid) = csmin;
    } // end if
    
    return;
    
} // end for user_eos_model





// -----------------------------------------------------------------------------
// This is the user material model function for the stress tensor
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
void user_strength_model(const DViewCArrayKokkos <double> &elem_pres,
                         const DViewCArrayKokkos <double> &elem_stress,
                         const size_t elem_gid,
                         const size_t mat_id,
                         const DViewCArrayKokkos <double> &elem_state_vars,
                         const DViewCArrayKokkos <double> &elem_sspd,
                         const double den,
                         const double sie,
                         const ViewCArrayKokkos <double> &vel_grad,
                         const ViewCArrayKokkos <size_t> &elem_node_gids,
                         const DViewCArrayKokkos <double> &node_coords,
                         const DViewCArrayKokkos <double> &node_vel,
                         const double vol,
                         const double dt,
                         const double rk_alpha){
    

    // statev(0) = var_1
    //   :
    //   :
    //   :
    // statev(N) = var_N

    // --strength
    // statev(5) = shear modulus
    // statev(6) = yield strength
  
    // **WARNING: Assumes single material point per element**

    // shear modulus
    double G_stress = elem_state_vars(elem_gid, 5);

    // yield strength
    double Y_stress = elem_state_vars(elem_gid, 6);
    
    const int num_dims = 3; // for 3D only

    // -----------------------------------------------------------------------------
    // The user must coding goes here
    //------------------------------------------------------------------------------
    // Cauchy deviatoric stress tensor, stress(rk,elem,i,j)

    // create a slice of the stress tensor
    ViewCArrayKokkos <double> stress_n(&elem_stress(0, elem_gid, 0, 0), num_dims, num_dims); // stress at previous timestep
    ViewCArrayKokkos <double> stress(&elem_stress(1, elem_gid, 0, 0), num_dims, num_dims); // stress at this timestep

    double s_dev_n_1D[num_dims * num_dims];
    double s_dev_1D[num_dims * num_dims];

    ViewCArrayKokkos <double> s_dev_n(s_dev_n_1D, num_dims, num_dims); // deviatoric stress at previous timestep
    ViewCArrayKokkos <double> s_dev(s_dev_1D, num_dims, num_dims); // deviatoric stress at this timestep

    for (int i = 0; i < num_dims; i++) {
        for (int j = 0; j < num_dims; j++) {
            s_dev_n(i, j) = stress_n(i, j);
            s_dev(i, j) = stress(i, j);
        }
    }

    // subtract 1/3*trace*I from the stress tensor to get the stress deviators
    //double trace_n = 0.0;
    //double trace = 0.0;
    //for (int i = 0; i < num_dims; i++) {
    //    trace_n += s_dev_n(i, i);
    //    trace += s_dev(i, i);
    //}
    
    
    //for (int i = 0; i < num_dims; i++) {
    //    s_dev_n(i, i) -= trace_n / 3.0;
    //    s_dev(i, i) -= trace / 3.0;
    //}
    
    // Jaumann rate RHS for evolving the stress_dev tensor
    double stress_dev_rhs_1d[num_dims * num_dims];
    ViewCArrayKokkos <double> stress_dev_rhs(stress_dev_rhs_1d, num_dims, num_dims);

    // symmetric part of velocity gradient
    double D_1d[num_dims * num_dims];
    ViewCArrayKokkos <double> D(D_1d, num_dims, num_dims);

    // anti-symmetric part of velocity gradient
    double W_1d[num_dims * num_dims];
    ViewCArrayKokkos <double> W(W_1d, num_dims, num_dims);

    // D = 1/2(vel_grad + vel_grad^T)
    for (int i = 0; i < num_dims; i++) {
      for (int j = 0; j < num_dims; j++) {
        D(i, j) = 0.5 * (vel_grad(i, j) + vel_grad(j, i));
      }
    }
    
    double div = 0;
    for (int i = 0; i < num_dims; i++) {
        div += vel_grad(i, i);
    }
    
    //double diff = fabs(D(0,0)+D(1,1)+D(2,2)-div);
    //if(diff>1e-10){ printf("abs div error = %f \n", diff); }
    //if(elem_gid==1){
    //    printf("div = %f \n", div);
    //}
    
    
    // W = 1/2(vel_grad - vel_grad^T)
    for (int i = 0; i < num_dims; i++) {
      for (int j = 0; j < num_dims; j++) {
          W(i, j) = 0.5 * (vel_grad(i, j) - vel_grad(j, i));
      }
    }

    


    // ----------------------------------------------------------------------
    // Jaumann rate
    //    dsigma'/dt =  2.0G(de'/dt - de_p/dt) + w*sigma' - sigma'*w
    //    de'/dt = D - div/3.0 * I
    //
    // ----------------------------------------------------------------------
    
    // Set up identity matrix
    double eye_size[num_dims * num_dims];
    ViewCArrayKokkos <double> eye(eye_size, num_dims, num_dims);
    // now calculate the identity matrix of s_dev
    for (int i = 0; i < num_dims; i++) {
        for (int j = 0; j < num_dims; j++) {
            
            if (i == j) {
                eye(i,j) = 1;
            }
            else {
                eye(i, j) = 0;
            } // end if
            
        }
    }
    // Now calculate Jaumann for each element
    for (int i = 0; i < num_dims; i++) {
        for (int j = 0; j < num_dims; j++) {
            stress_dev_rhs(i, j) = 2.0 * G_stress * (D(i, j) - ((div / 3.0) * eye(i, j)));
        }
    }
    for (int i = 0; i < num_dims; i++) {
        for (int k = 0; k < num_dims; k++) {
            for (int j = 0; j < num_dims; j++) {
                stress_dev_rhs(i, k) +=
                W(i, j) * s_dev(j, k) - s_dev(i, j) * W(j, k);
            }
        }
    }
    
    // calculate the next stress
    for (int i = 0; i < num_dims; i++) {
        for (int j = 0; j < num_dims; j++) {
            s_dev(i, j) = s_dev_n(i, j) + (rk_alpha * dt * stress_dev_rhs(i, j));
        }
    }
    
    // scale to be on the yield surface, J2 plasticity
    double J2 = 0.0;
    
    for (int i = 0; i < num_dims; i++) {
        for (int j = 0; j < num_dims; j++) {
            J2 += s_dev(i, j) * s_dev(i, j);
        }
    }
    J2 = J2 / 2.0;
    

    // adjust to Mohr's circle
    double factor = Y_stress / sqrt(3.0 * J2 + fuzz);
    
    if (factor < 1.0) {
        for (int i = 0; i < num_dims; i++) {
            for (int j = 0; j < num_dims; j++) {
                s_dev(i, j) *= factor;
            }
        }
    } // end if
    
    
    // set the stress
    for (int i = 0; i < num_dims; i++) {
        for (int j = 0; j < num_dims; j++) {
            stress(i, j) = s_dev(i, j);
        }
    }
    
    return;

} // end of user mat

// -----------------------------------------------------------------------------
// This is the user material model function
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
void user_strength_model_vpsc(const DViewCArrayKokkos <double> &elem_pres,
                              const DViewCArrayKokkos <double> &elem_stress,
                              const size_t elem_gid,
                              const size_t mat_id,
                              const DViewCArrayKokkos <double> &elem_state_vars,
                              const DViewCArrayKokkos <double> &elem_sspd,
                              const double div,
                              const double den,
                              const double sie,
                              const ViewCArrayKokkos <double> &vel_grad,
                              const ViewCArrayKokkos <size_t> &elem_node_gids,
                              const DViewCArrayKokkos <double> &node_coords,
                              const DViewCArrayKokkos <double> &node_vel,
                              const double vol,
                              const double dt,
                              const double rk_alpha){
    

    // statev(0) = var_1
    //   :
    //   :
    //   :
    // statev(N) = var_N

    const int num_dims = 3;
    
    // -----------------------------------------------------------------------------
    // Required variables are here
    //------------------------------------------------------------------------------
    elem_pres(elem_gid) = 1.0e-15;  // pressure
    elem_sspd(elem_gid) = 1.0e-15;  // sound speed

    
    // -----------------------------------------------------------------------------
    // The user must coding goes here
    //------------------------------------------------------------------------------
    
    // For hypo-elastic models
    double D_tensor_values[9];
    double W_tensor_values[9];

    // convert to array syntax with the C-Language access pattern
    ViewCArrayKokkos <double> D_tensor(D_tensor_values, num_dims, num_dims);  // D(i,j)
    ViewCArrayKokkos <double> W_tensor(W_tensor_values, num_dims, num_dims);  // W(i,j)
    
    
    decompose_vel_grad(D_tensor,
                       W_tensor,
                       vel_grad,
                       elem_node_gids,
                       elem_gid,
                       node_coords,
                       node_vel,
                       vol);
    

    // For hypo-elastic models
    double deps_values[9];
    double dW_values[9];
    double drot_values[9];
    
    ViewCMatrixKokkos <double> deps(&deps_values[0], num_dims, num_dims);
    ViewCMatrixKokkos <double> dW(&dW_values[0], num_dims, num_dims);
    ViewCMatrixKokkos <double> drot(&drot_values[0], num_dims, num_dims);
    
    // calculate strain and rotation increments
    for (size_t i = 1; i <= 3; i++) {
        for (size_t j = 1; j <= 3; j++) {
            deps(i, j) = D_tensor(i-1,j-1)*dt; // infinitesimal strain increment
            dW(i, j) = W_tensor(i-1,j-1)*dt;
        }
    } // end for
    
    
    
    return;
    
} // end of ideal_gas

KOKKOS_FUNCTION
void user_strength_model_2D(const DViewCArrayKokkos <double> &elem_pres,
                            const DViewCArrayKokkos <double> &elem_stress,
                            const size_t elem_gid,
                            const size_t mat_id,
                            const DViewCArrayKokkos <double> &elem_state_vars,
                            const DViewCArrayKokkos <double> &elem_sspd,
                            const double den,
                            const double sie,
                            const ViewCArrayKokkos <double> &vel_grad,
                            const ViewCArrayKokkos <size_t> &elem_node_gids,
                            const DViewCArrayKokkos <double> &node_coords,
                            const DViewCArrayKokkos <double> &node_vel,
                            const double vol,
                            const double dt,
                            const double rk_alpha){
    
    /*
     
    const int num_dims = 3; // must be 3 even though its 2D RZ

    // -----------------------------------------------------------------------------
    // The user must coding goes here
    //------------------------------------------------------------------------------
    // Cauchy deviatoric stress tensor, stress(rk,elem,i,j)
    
    
    double div = 0;
    for (int i = 0; i < 2; i++) {
        div += vel_grad(i, i);
    }
    


        // -----------------------------------------------------------------------------
        // The user must coding goes here
        //------------------------------------------------------------------------------
        // Cauchy deviatoric stress tensor, stress(rk,elem,i,j)

        // create a slice of the stress tensor
        ViewCArrayKokkos <double> stress_n(&elem_stress(0, elem_gid, 0, 0), 3, 3); // stress from previous timestep
        ViewCArrayKokkos <double> stress(&elem_stress(1, elem_gid, 0, 0), 3, 3); // stress from this timestep

        double s_dev_n_1D[num_dims*num_dims];
        double s_dev_1D[num_dims*num_dims];

        ViewCArrayKokkos <double> s_dev_n(s_dev_n_1D, num_dims, num_dims); // deviatoric stress from previous timestep
        ViewCArrayKokkos <double> s_dev(s_dev_1D, num_dims, num_dims); // deviatoric stress from this timestep

        for (int i = 0; i < num_dims; i++) {
          for (int j = 0; j < num_dims; j++) {
            s_dev_n(i, j) = stress_n(i, j);
            s_dev(i, j) = stress(i, j);
          }
        }

        // subtract 1/3*trace*I from the stress tensor to get the stress deviators
        double trace_n = 0.0;
        double trace = 0.0;
        for (int i = 0; i < num_dims; i++) {
          trace_n += s_dev_n(i, i);
          trace += s_dev(i, i);
        }
        for (int i = 0; i < num_dims; i++) {
          s_dev_n(i, i) -= trace_n / 3.0;
          s_dev(i, i) -= trace / 3.0;
        }

        // Jaumann rate RHS for evolving the stress_dev tensor
        double stress_dev_rhs_1d[num_dims * num_dims];
        ViewCArrayKokkos <double> stress_dev_rhs(stress_dev_rhs_1d, num_dims, num_dims);

        // symmetric part of velocity gradient
        double D_1d[num_dims * num_dims];
        ViewCArrayKokkos <double> D(D_1d, num_dims, num_dims);

        // anti-symmetric part of velocity gradient
        double W_1d[num_dims * num_dims];
        ViewCArrayKokkos <double> W(W_1d, num_dims, num_dims);

        // D = 1/2(vel_grad + vel_grad^T)
        for (int i = 0; i < num_dims; i++) {
          for (int j = 0; j < num_dims; j++) {
            D(i, j) = 0.5 * (vel_grad(i, j) + vel_grad(j, i));
          }
        }

        // W = 1/2(vel_grad - vel_grad^T)
        for (int i = 0; i < num_dims; i++) {
          for (int j = 0; j < num_dims; j++) {
            W(i, j) = 0.5 * (vel_grad(i, j) - vel_grad(j, i));
          }
        }

        // **WARNING: Assumes single material point per element**

        // shear modulus
    double G_stress = elem_state_vars(elem_gid, 5);

        // yield strength
    double Y_stress = elem_state_vars(elem_gid, 6);

        // ----------------------------------------------------------------------
        // Jaumann rate
        //    dsigma'/dt =  2.0G(de'/dt - de_p/dt) + w*sigma' - sigma'*w
        //    de'/dt = D - div/3.0 * I
        //
        // ----------------------------------------------------------------------
        
        // calculate the RHS for the Jaumann rate for an elastic response in 2D RZ and state
        double term = (s_dev_n(0, 0) - s_dev_n(1, 1))/2.0;

        // setup for the identity matrix
        double eye_size[num_dims * num_dims];
        ViewCArrayKokkos <double> eye(eye_size, num_dims, num_dims);
        // now calculate the identity matrix of s_dev
        for (int i = 0; i < num_dims; i++) {
          for (int j = 0; j < num_dims; j++) {
              
            if (i == j) {
              eye(i,j) = 1;
            }
            else {
                eye(i, j) = 0;
            }
              
          }
        }
        
        // Equations according to Mark Wilkins 1963
        stress_dev_rhs(0, 0) = (2.0 * G_stress * (D(0, 0) - ((div / 3.0)*eye(0, 0)))) + (term*(((D(0, 1)*s_dev(1, 0))+(s_dev(0, 1)*D(1, 0))) + ((D(0, 0)*s_dev(0, 0))+(s_dev(0, 0)*D(0, 0)))))
                                    + s_dev_n(0, 1)*((W(0, 1) * s_dev(1, 0)) - (s_dev(0, 1) * W(1, 0)));
        stress_dev_rhs(0, 1) = (G_stress * D(0, 1)) + ((s_dev_n(0, 1)*((((D(0, 0)*s_dev(0, 1))+(s_dev(0, 0)*D(0, 1))) + ((D(0, 1)*s_dev(1, 1))+(s_dev(0, 1)*D(1, 1))))))
                                  - term*(((W(0, 0) * s_dev(0, 1)) - (s_dev(0, 0) * W(0, 1)))+ ((W(0, 1)*s_dev(1, 1))-(s_dev(0, 1)*W(1, 1)))));
        stress_dev_rhs(0, 2) = 0.0;
        stress_dev_rhs(1, 0) = (G_stress * D(1, 0)) + ((s_dev_n(1, 0)*((((D(1, 0)*s_dev(0, 0))+(s_dev(1, 0)*D(0, 0))) + ((D(1, 1)*s_dev(1, 0))+(s_dev(1, 1)*D(1, 0))))))
                                  - term*(((W(1, 0) * s_dev(0, 0)) - (s_dev(1, 0) * W(0, 0))) + ((W(1, 1)*s_dev(1, 0))-(s_dev(1, 1)*W(1, 0)))));
        stress_dev_rhs(1, 1) = (2.0 * G_stress * (D(1, 1) - ((div / 3.0)*eye(1, 1)))) + ((W(1, 0) * s_dev(0, 1)) - (s_dev(1, 0) * W(0, 1)));
        stress_dev_rhs(1, 2) = 0.0;
        stress_dev_rhs(2, 0) = 0.0;
        stress_dev_rhs(2, 1) = 0.0;
        stress_dev_rhs(2, 2) = (2.0 * G_stress * (D(2, 2) - ((div / 3.0)*eye(2, 2))));

        // calculate the next stress
        for (int i = 0; i < num_dims; i++) {
          for (int j = 0; j < num_dims; j++) {
            s_dev(i, j) = s_dev_n(i, j) + (rk_alpha * dt * stress_dev_rhs(i, j));
          }
        }

        // scale to be on the yield surface, J2 plasticity
        double J2 = 0.0;
        
        for (int i = 0; i < num_dims; i++) {
          for (int j = 0; j < num_dims; j++) {
             J2 += s_dev(i, j) * s_dev(i, j);
           }
         }

        double K = J2 - ((2.0/3.0)*(Y_stress*Y_stress));

        // adjust to Mohr's circle
        double factor = (sqrt(2.0/3.0)*Y_stress) / sqrt(J2 + fuzz);

        // real_t factor = 0.9;  // warning warning testing  WARNING WARNING
        // WARNING testing

        if (K > 0.0) {
          s_dev(0, 0) *= factor;
          s_dev(1, 1) *= factor;
          s_dev(2, 2) *= factor;
          s_dev(0, 1) *= factor;
          s_dev(1, 0) *= factor;
        }
        // set the new stress
        stress(0, 0) = s_dev(0, 0);
        stress(0, 1) = s_dev(0, 1);
        stress(0, 2) = 0.0;
        stress(1, 0) = s_dev(0, 1);
        stress(1, 1) = s_dev(1, 1);
        stress(1, 2) = 0.0;
        stress(2, 0) = 0.0;
        stress(2, 1) = 0.0;
        stress(2, 2) = s_dev(2, 2);
    
    //printf("\n");
    //printf("elem_id = %zu\n", elem_gid);
    //for (int i = 0; i < num_dims; i++) {
    //  for (int j = 0; j < num_dims; j++) {
    //     printf("%f,",s_dev(i, j));
    //   }
    //    printf("\n");
    // }
    //printf("\n");
    
    return;

    */
    
    


    // statev(0) = var_1
    //   :
    //   :
    //   :
    // statev(N) = var_N

    const int num_dims = 3; // must be 3 even though its 2D RZ

    // -----------------------------------------------------------------------------
    // The user must coding goes here
    //------------------------------------------------------------------------------
    // Cauchy deviatoric stress tensor, stress(rk,elem,i,j)

    // create a slice of the stress tensor
    ViewCArrayKokkos <double> stress_n(&elem_stress(0, elem_gid, 0, 0), 3, 3); // stress from previous timestep
    ViewCArrayKokkos <double> stress(&elem_stress(1, elem_gid, 0, 0), 3, 3); // stress from this timestep

    double s_dev_n_1D[3 * 3];
    double s_dev_1D[3 * 3];

    ViewCArrayKokkos <double> s_dev_n(s_dev_n_1D, 3, 3); // deviatoric stress from previous timestep
    ViewCArrayKokkos <double> s_dev(s_dev_1D, 3, 3); // deviatoric stress from this timestep

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        s_dev_n(i, j) = stress_n(i, j);
        s_dev(i, j) = stress(i, j);
      }
    }

    // subtract 1/3*trace*I from the stress tensor to get the stress deviators
    double trace_n = 0.0;
    double trace = 0.0;
    for (int i = 0; i < 3; i++) {
      trace_n += s_dev_n(i, i);
      trace += s_dev(i, i);
    }
    for (int i = 0; i < 3; i++) {
      s_dev_n(i, i) -= trace_n / 3.0;
      s_dev(i, i) -= trace / 3.0;
    }

    // Jaumann rate RHS for evolving the stress_dev tensor
    double stress_dev_rhs_1d[3 * 3];
    ViewCArrayKokkos <double> stress_dev_rhs(stress_dev_rhs_1d, 3, 3);

    // symmetric part of velocity gradient
    double D_1d[3 * 3];
    ViewCArrayKokkos <double> D(D_1d, 3, 3);
    
    // anti-symmetric part of velocity gradient
    double W_1d[num_dims * num_dims];
    ViewCArrayKokkos <double> W(W_1d, num_dims, num_dims);
    
    
    // D = 1/2(vel_grad + vel_grad^T)
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        D(i, j) = 0.5 * (vel_grad(i, j) + vel_grad(j, i));
      }
    }
    
    // W = 1/2(vel_grad - vel_grad^T)
    for (int i = 0; i < num_dims; i++) {
      for (int j = 0; j < num_dims; j++) {
        W(i, j) = 0.5 * (vel_grad(i, j) - vel_grad(j, i));
      }
    }
    
    
    double div = 0;
    for (int i = 0; i < num_dims; i++) {
        div += vel_grad(i, i);
    }

    // vel = [u,v]
    // vel_grad = [du/dx,  du/dy]
    //            [dv/dx,  dv/dy]
    // vel_grad(1,0) = dv/dx
    // vel_grad(0,1) = du/dy
    double curl;
    curl = vel_grad(1,0) - vel_grad(0,1);  // dv/dx - du/dy
    


    // **WARNING: Assumes single material point per element**

    // shear modulus
    double G_stress = elem_state_vars(elem_gid, 5);

    // yield strength
    double Y_stress = elem_state_vars(elem_gid, 6);

    // ----------------------------------------------------------------------
    // Strength used by Maenchen and Sack, and by Wilkins (1964)
    //    dsigma' =  2.0G(de'/dt - de_p/dt)dt - Rotation
    //    de'/dt = D - div/3.0 * I
    // Note: Rotation is only on the first two equations, not hoop stress term
    // ----------------------------------------------------------------------
    
    // calculate the RHS for the Jaumann rate for an elastic response in 2D RZ and state
    double term = (s_dev_n(0, 0) - s_dev_n(1, 1))/2.0;

    // setup for the identity matrix
    double eye_size[9];
    ViewCArrayKokkos <double> eye(eye_size, 3, 3);
    // now calculate the identity matrix of s_dev
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
          
        if (i == j) {
            eye(i,j) = 1;
        }
        else {
            eye(i, j) = 0;
        } // end if
          
      } // end for j
    }  // end for i
    
    
    int approach = 0;
    
    if(approach==1){
        
        // Now calculate Jaumann for each element
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            stress_dev_rhs(i, j) = 2.0 * G_stress * (D(i, j) - ((div / 3.0) * eye(i, j)));
          }
        }
        
        // ----------------------
        // NOTE: the hoop stress doesn't use these Rotation terms
        double sinw = 0.5*curl*rk_alpha*dt;   // the sin of the angle of rotation over dt;
        double cosw = sqrt(1.0 - sinw*sinw);  // sinw^2 + cosw^2 = 1
        
        double Rzz = (s_dev(1, 1) - s_dev(0, 0))*sinw*sinw + 2.0*s_dev(1, 0)*sinw*cosw;
        double Rrr = -Rzz;
        double Rrz = (s_dev(1, 1) - s_dev(0, 0))*sinw*cosw - 2.0*s_dev(1, 0)*sinw*sinw;

        

        // calculate the next stress
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            s_dev(i, j) = s_dev_n(i, j) + (rk_alpha * dt * stress_dev_rhs(i, j));
          } // end for j
        } // end for i
        
        
        // account for the rotations on zz rr and rz values
        s_dev(0, 0) += Rzz;
        s_dev(1, 1) += Rrr;
        s_dev(1, 0) += Rrz;
        s_dev(0, 1) = s_dev(1,0);  // symmetric stress
        

        // scale to be on the yield surface, J2 plasticity
        double J2 = 0.0;

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            J2 += s_dev(i, j) * s_dev(i, j);
          } // end for j
        } // end for i
        J2 = J2 / 2.0;

        // adjust to Mohr's circle
        double factor = Y_stress / sqrt(3.0 * J2 + fuzz);


        if (factor < 1.0) {
          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
              s_dev(i, j) *= factor;
            }
          }
        } // end if

        // set the stress
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            stress(i, j) = s_dev(i, j);
          }
        }
        

        stress(0, 2) = 0.0;
        stress(1, 2) = 0.0;
        stress(2, 0) = 0.0;
        stress(2, 1) = 0.0;
    }
    else {
            
        
        // Equations according to Mark Wilkins 1963
        stress_dev_rhs(0, 0) = (2.0 * G_stress * (D(0, 0) - ((div / 3.0)*eye(0, 0)))) + (term*(((D(0, 1)*s_dev(1, 0))+(s_dev(0, 1)*D(1, 0))) + ((D(0, 0)*s_dev(0, 0))+(s_dev(0, 0)*D(0, 0)))))
                                    + s_dev_n(0, 1)*((W(0, 1) * s_dev(1, 0)) - (s_dev(0, 1) * W(1, 0)));
        stress_dev_rhs(0, 1) = (G_stress * D(0, 1)) + ((s_dev_n(0, 1)*((((D(0, 0)*s_dev(0, 1))+(s_dev(0, 0)*D(0, 1))) + ((D(0, 1)*s_dev(1, 1))+(s_dev(0, 1)*D(1, 1))))))
                                  - term*(((W(0, 0) * s_dev(0, 1)) - (s_dev(0, 0) * W(0, 1)))+ ((W(0, 1)*s_dev(1, 1))-(s_dev(0, 1)*W(1, 1)))));
        stress_dev_rhs(0, 2) = 0.0;
        stress_dev_rhs(1, 0) = (G_stress * D(1, 0)) + ((s_dev_n(1, 0)*((((D(1, 0)*s_dev(0, 0))+(s_dev(1, 0)*D(0, 0))) + ((D(1, 1)*s_dev(1, 0))+(s_dev(1, 1)*D(1, 0))))))
                                  - term*(((W(1, 0) * s_dev(0, 0)) - (s_dev(1, 0) * W(0, 0))) + ((W(1, 1)*s_dev(1, 0))-(s_dev(1, 1)*W(1, 0)))));
        stress_dev_rhs(1, 1) = (2.0 * G_stress * (D(1, 1) - ((div / 3.0)*eye(1, 1)))) + ((W(1, 0) * s_dev(0, 1)) - (s_dev(1, 0) * W(0, 1)));
        stress_dev_rhs(1, 2) = 0.0;
        stress_dev_rhs(2, 0) = 0.0;
        stress_dev_rhs(2, 1) = 0.0;
        stress_dev_rhs(2, 2) = (2.0 * G_stress * (D(2, 2) - ((div / 3.0)*eye(2, 2))));

        // calculate the next stress
        for (int i = 0; i < num_dims; i++) {
          for (int j = 0; j < num_dims; j++) {
            s_dev(i, j) = s_dev_n(i, j) + (rk_alpha * dt * stress_dev_rhs(i, j));
          }
        }

        // scale to be on the yield surface, J2 plasticity
        double J2 = 0.0;
        
        for (int i = 0; i < num_dims; i++) {
          for (int j = 0; j < num_dims; j++) {
             J2 += s_dev(i, j) * s_dev(i, j);
           }
         }

        double K = J2 - ((2.0/3.0)*(Y_stress*Y_stress));

        // adjust to Mohr's circle
        double factor = (sqrt(2.0/3.0)*Y_stress) / sqrt(J2 + fuzz);

        // real_t factor = 0.9;  // warning warning testing  WARNING WARNING
        // WARNING testing

        if (K > 0.0) {
          s_dev(0, 0) *= factor;
          s_dev(1, 1) *= factor;
          s_dev(2, 2) *= factor;
          s_dev(0, 1) *= factor;
          s_dev(1, 0) *= factor;
        }
        // set the new stress
        stress(0, 0) = s_dev(0, 0);
        stress(0, 1) = s_dev(0, 1);
        stress(0, 2) = 0.0;
        stress(1, 0) = s_dev(0, 1);
        stress(1, 1) = s_dev(1, 1);
        stress(1, 2) = 0.0;
        stress(2, 0) = 0.0;
        stress(2, 1) = 0.0;
        stress(2, 2) = s_dev(2, 2);
     
    }

    return;

    
    
} // end of user mat



/*
// -----------------------------------------------------------------------------
// This is the user material model function for the equation of state
// An eos function must be supplied or the code will fail to run.
// The pressure and sound speed can be calculated from an analytic eos.
// The pressure can also be calculated using p = -1/3 Trace(Stress)
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
void foam_eos_model(const DViewCArrayKokkos <double> &elem_pres,
                    const DViewCArrayKokkos <double> &node_vel,
                    const size_t elem_gid,
                    const size_t mat_id,
                    const DViewCArrayKokkos <double> &elem_state_vars,
                    const DViewCArrayKokkos <double> &elem_sspd,
                    const DViewCArrayKokkos <double> &elem_vol,
                    const DViewCArrayKokkos <double> &elem_mass,
                    double den,
                    const double sie,
                    const double dt,
                    const double rk_alpha){
    
    // -----------------------------------------------------------------------------
    // Required variables are here
    //------------------------------------------------------------------------------
    
    // Mie-Grüneisen eos for metal from Liu et al 2022
    // pres = Ph + gamma*den*(e-Eh)
    // Ph = (den_initial*C0^2*eta)/(1-s*eta)^2
    // Eh = (eta*Ph)/(2*den_init)
    // eta = 1 - den_init/den

    // ---eos
    // statev(0) = gamma
    // statev(1) = c0
    // statev(2) = s0
    // statev(3) = den_ref
    // statev(4) = csmin
    
    
    // statev(5) = void_frac
    // statev(6) = sie_mat


    
    // void_frac = void_frac + (V^n+1 - V^n)
    // density = M/V_m
    
    
    // V0_m = M/rho_ref
    
    // den = rho_ref no more void in element
    
    
    
    double gamma = elem_state_vars(elem_gid,0);
    double C0 = elem_state_vars(elem_gid,1);
    double S0 = elem_state_vars(elem_gid,2);
    double den_ref = elem_state_vars(elem_gid,3);
    double csmin = elem_state_vars(elem_gid,4);
    double G = elem_state_vars(elem_gid,5); // shear modulus
    

    
    // calculate the pressure components
    double eta = 1.0 - (den_ref/den);
    double Ph = (den_ref*C0*C0*eta)/((1.0 - S0*eta)*(1.0 - S0*eta));
    double Eh = (eta*Ph)/(2.0*den_ref);

    // full pressure calculation
    // from Liu et al 2022
    elem_pres(elem_gid) = Ph + gamma*den*(sie - Eh);
    
    
    // now calculate the sound speed
    // variables for derivatives and easier readability
    double Vo = 1.0/den_ref;
    double V = 1.0/den;
    double easy = 1.0/(Vo - elem_state_vars(elem_gid,2)*(Vo - V));

    // derivative of Huginot pressure/density with respect to energy
    double dPhdr = (C0*C0)/(den*den) *
                    ( (easy*easy) + 2.0*S0*(Vo - V)*(easy*easy*easy) );
    double dgdr = gamma; //gamma_0 + gamma_1*(elem_state_vars(elem_gid,0)*den_ref)*den*(-1.0/(den*den));
    double dEdr = 0.5 * (dPhdr * (Vo - V) + ((elem_pres(elem_gid) + Ph)/(den*den)));

    // c_EOS is the sound speed squared
    double c_EOS = dPhdr + dgdr*(sie - Eh) + gamma*den*dEdr + (elem_pres(elem_gid)/(den*den))*gamma*den;
    c_EOS = fmax(0.0, c_EOS); //check that sound speed it greater than zero
  
    double B = den*c_EOS*c_EOS; // Bulk Modulus
    
    // final calculation for sound speed for the element
    elem_sspd(elem_gid) = sqrt((B + (4.0/3.0)*G)/den);  // this is adding shear modulus to sound speed
    
    
    if (elem_sspd(elem_gid) < csmin) {
        elem_sspd(elem_gid) = csmin;
    } // end if
    
    return;
    
} // end fcn user_eos_model
*/


/*
void tipton_two_component(const size_t elem_gid,
                          const size_t mat_id,
                          const DViewCArrayKokkos <double> &elem_state_vars,
                          const DViewCArrayKokkos <double> &elem_sspd,
                          const ViewCArrayKokkos <double> &vel_grad,
                          const double vol,
                          const double dt,
                          const double alpha){
    
    
    int num_dims = 3;
    
    //double gamma = elem_state_vars(elem_gid,0);
    //double C0 = elem_state_vars(elem_gid,1);
    //double S0 = elem_state_vars(elem_gid,2);
    //double den_ref = elem_state_vars(elem_gid,3);
    //double csmin = elem_state_vars(elem_gid,4);
    //double G = elem_state_vars(elem_gid,5); // shear modulus
    // Y = elem_state_vars(elem_gid,6)
    
        
    double vol_frac_1 = state_vars(7);  // volume fraction
    double vol_frac_2 = 1 - vol_frac_1;
    
    double den_1 = state_vars(8);
    double den_2 = state_vars(9);
    
    double sspd_1 = state_vars(10);
    double sspd_2 = state_vars(11);
    
    double pres_1 = state_vars(12);
    double pres_2 = state_vars(13);
    
    double sie_1 = state_vars(14);
    double sie_2 = state_vars(15);
    
    auto stress_1 = ViewCArrayKokkos <double> (&state_vars(16),3,3);
    auto stress_2 = ViewCArrayKokkos <double> (&state_vars(24),3,3);
    
    
    // ----
    
    double mass_1 = den_1*vol_frac_1*vol;
    double mass_2 = den_2*vol_frac_2*vol;
    
    double div = 0;
    for (int i = 0; i < num_dims; i++) {
        div += vel_grad(i, i);
    } // end for
    
    double term1 = vol_frac_1/(den_1*sspd_1*sspd_1);
    double term2 = vol_frac_2/(den_2*sspd_2*sspd_2);
    
    double impedance_avg = 1.0/(term1 + term2);
    double pres_avg = (term_1*pres_1 + term_2*pres_2)/(term1 + term2);
    
    double delta_vol_frac_1 = 0;
    double delta_vol_frac_2 = 0;
    
    double term3 = (pres_1 - pres_avg)/(0.25*den_1*sspd_1*sspd_1);
    double term4 = impedance_avg/(den_1*sspd_1*sspd_1) - 1.0;
    
    delta_vol_frac_1 = vol_frac_1*(term3 + term4*div*alpha*dt);
    
    // limit the volume fraction change to ensure realistic vol_fractions
    double min_frac_mag = fmax( fabs(vol_frac_1), fabs(vol_frac_2) );
    if (delta_vol_frac_1<0.0) {
        // do not change the material volume more than exists
        delta_vol_frac_1 = -fmin(fabs(delta_vol_frac_1), min_frac_mag);  // note the -1
    }
    else if (delta_vol_frac_1>0.0) {
        delta_vol_frac_1 = fmin(fabs(delta_vol_frac_1), min_frac_mag);
    }
    
    delta_vol_frac_2 = -delta_vol_frac_1;

    
    // calculate the new volume fractions
    vol_frac_1 += delta_vol_frac_1;
    vol_frac_2 = 1.0 - vol_frac_1;
    
    // error check
    if(vol_frac_1 < 0.0 || vol_frac_1 > 1.0) printf("vol fraction error in component 1 = %f", vol_frac_1);
    if(vol_frac_2 < 0.0 || vol_frac_2 > 1.0) printf("vol fraction error in component 2 = %f", vol_frac_2);
    
    
    // adjust internal energy for the volume change
    //sie_1 -= delta_vol_frac_1*pres_avg*div*alpha*dt;
    //sie_2 -= delta_vol_frac_2*pres_avg*div*alpha*dt;
    
    
    // save the new volume fractions
    state_vars(7) = vol_frac_1;
    
    
    // save the new densities
    state_vars(8) = mass_1/(vol_frac_1*vol);  // den_1
    state_vars(9) = mass_2/(vol_frac_2*vol);  // den_2
    

    
    
    // ensure conservation
    //double delta_elem_sie = elem_sie(elem_gid) - (mass_1*sie_1 + mass_2*sie_2);
    
    
} // end fcn tipton multi-comp model
*/
