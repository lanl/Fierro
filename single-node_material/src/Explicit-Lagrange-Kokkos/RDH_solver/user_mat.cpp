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
                       const DViewCArrayKokkos <double> &elem_stress,
                       const size_t elem_gid,
                       const size_t legendre_gid,
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
    elem_pres(legendre_gid) = Ph + gamma*den*(sie - Eh);
    
    
    // now calculate the sound speed
    // variables for derivatives and easier readability
    double Vo = 1.0/den_ref;
    double V = 1.0/den;
    double easy = 1.0/(Vo - elem_state_vars(elem_gid,2)*(Vo - V));

    // derivative of Huginot pressure/density with respect to energy
    double dPhdr = (C0*C0)/(den*den) *
                    ( (easy*easy) + 2.0*S0*(Vo - V)*(easy*easy*easy) );
    double dgdr = gamma; //gamma_0 + gamma_1*(elem_state_vars(elem_gid,0)*den_ref)*den*(-1.0/(den*den));
    double dEdr = 0.5 * (dPhdr * (Vo - V) + ((elem_pres(legendre_gid) + Ph)/(den*den)));

    // c_EOS is the sound speed squared
    double c_EOS = dPhdr + dgdr*(sie - Eh) + gamma*den*dEdr + (elem_pres(legendre_gid)/(den*den))*gamma*den;
    c_EOS = fmax(0.0, c_EOS); //check that sound speed it greater than zero
  
    double B = den*c_EOS*c_EOS; // Bulk Modulus
    
    // final calculation for sound speed for the element
    elem_sspd(legendre_gid) = sqrt((B + (4.0/3.0)*G)/den);  // this is adding shear modulus to sound speed
    
    
    if (elem_sspd(legendre_gid) < csmin) {
        elem_sspd(legendre_gid) = csmin;
    } // end if
    
    return;
    
} // end for user_eos_model





// -----------------------------------------------------------------------------
// This is the user material model function for the stress tensor
//------------------------------------------------------------------------------


void user_strength_model(CArrayKokkos <double> &deviatoric_stress_rhs,
                         const DViewCArrayKokkos <double> &stress,
                         const CArrayKokkos <double> &mat_pt_state_vars,
                         const CArrayKokkos <double> &sym_grad_vel,
                         const CArrayKokkos <double> &anti_sym_grad_vel,
                         const CArrayKokkos <double> &div_vel,
                         const size_t num_gauss,
                         const size_t stage){
    

    // statev(0) = var_1
    //   :
    //   :
    //   :
    // statev(N) = var_N

    // --strength
    // statev(5) = shear modulus
    // statev(6) = yield strength
  
    // **WARNING: Assumes single material point per element**
    FOR_ALL(gauss_gid, 0, num_gauss,{

        // shear modulus
        double G_stress = 0.286;  //mat_pt_state_vars(gauss_gid, 5);

        // yield strength
        double Y_stress = 0.0026; //mat_pt_state_vars(gauss_gid, 6);

        const int num_dims = 3; // for 3D only

        // -----------------------------------------------------------------------------
        // The user must coding goes here
        //------------------------------------------------------------------------------


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





        // ----------------------------------------------------------------------
        // Jaumann rate
        //    dsigma'/dt =  2.0G(de'/dt - de_p/dt) + w*sigma' - sigma'*w
        //    de'/dt = D - div/3.0 * I
        //
        // ----------------------------------------------------------------------

        // Set up identity matrix
        CArrayKokkos <double> eye(num_dims, num_dims);
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
                deviatoric_stress_rhs(stage, gauss_gid, i, j) = 2.0 * G_stress * (sym_grad_vel(gauss_gid, i, j) - ((div_vel(gauss_gid) / 3.0) * eye(i, j)));
            }
        }
        for (int i = 0; i < num_dims; i++) {
            for (int k = 0; k < num_dims; k++) {
                for (int j = 0; j < num_dims; j++) {
                    deviatoric_stress_rhs(stage, gauss_gid, i, k) += anti_sym_grad_vel(gauss_gid, i, j) * stress(stage, gauss_gid, j, k) - stress(stage, gauss_gid, i, j) * anti_sym_grad_vel(gauss_gid, i, j);
                    //printf("deviatoric stress rhs : %f \n", deviatoric_stress_rhs(stage, gauss_gid, i, j));
                }
            }
        }

    });

    // calculate the next stress
    // for (int i = 0; i < num_dims; i++) {
    //     for (int j = 0; j < num_dims; j++) {
    //         deviatoric_stress(1, gauss_gid, i, j) = deviatoric_stress(0, gauss_gid, i, j) + (dt * stress_dev_rhs(i, j));
    //     }
    // }
    
    // scale to be on the yield surface, J2 plasticity
    // double J2 = 0.0;
    
    // for (int i = 0; i < num_dims; i++) {
    //     for (int j = 0; j < num_dims; j++) {
    //         J2 += s_dev(i, j) * s_dev(i, j);
    //     }
    // }
    // J2 = J2 / 2.0;
    

    // // adjust to Mohr's circle
    // double factor = Y_stress / sqrt(3.0 * J2 + fuzz);
    
    // if (factor < 1.0) {
    //     for (int i = 0; i < num_dims; i++) {
    //         for (int j = 0; j < num_dims; j++) {
    //             s_dev(i, j) *= factor;
    //         }
    //     }
    // } // end if
    
    
    // set the stress
    // for (int i = 0; i < num_dims; i++) {
    //     for (int j = 0; j < num_dims; j++) {
    //         stress(i, j) = s_dev(i, j);
    //     }
    // }
    

} // end of user mat

