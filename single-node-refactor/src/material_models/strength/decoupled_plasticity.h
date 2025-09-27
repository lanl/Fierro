/**********************************************************************************************
� 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and
to permit others to do so.
This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
1.  Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright notice, this list of
conditions and the following disclaimer in the documentation and/or other materials
provided with the distribution.
3.  Neither the name of the copyright holder nor the names of its contributors may be used
to endorse or promote products derived from this software without specific prior
written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********************************************************************************************/
#ifndef DECOUPLED_HYPO_PLASTICITY_H
#define DECOUPLED_HYPO_PLASTICITY_H


/////////////////////////////////////////////////////////////////////////////
///
/// \fn HypoPlasticity
///
/// \brief stress model for hypo plasticity response
///
///  This is a material model function for the deviatoric stress tensor
///
/// \param Material pressure
/// \param Material stress
/// \param Global ID for the material
/// \param Material ID for the element
/// \param Material state variables
/// \param Material Sound speed
/// \param Material density
/// \param Material specific internal energy
/// \param Element velocity gradient
/// \param Element nodes IDs in the element
/// \param Node node coordinates
/// \param Noe velocity 
/// \param Element volume
/// \param Time time step
/// \param Time coefficient in the Runge Kutta time integration step
///
/////////////////////////////////////////////////////////////////////////////
namespace HypoPlasticityModel {

    const double fuzz = 1.e-16;
    
    /**
     * @brief Initialize the strength state variables for the material points.
     *
     * @param MaterialPoints_eos_state_vars   State variables for the equation of state at each material point.
     * @param MaterialPoints_strength_state_vars   State variables for the strength model at each material point.
     * @param eos_global_vars   Global variables for the equation of state.
     * @param strength_global_vars   Global variables for the strength model.
     * @param elem_in_mat_elem   Mapping from material points to mesh elements.
     * @param num_material_points   Number of material points for this material.
     * @param mat_id   Material ID.
     */
    static void init_strength_state_vars(
        const DRaggedRightArrayKokkos <double> &MaterialPoints_eos_state_vars,
        const DRaggedRightArrayKokkos <double> &MaterialPoints_strength_state_vars,
        const RaggedRightArrayKokkos <double> &eos_global_vars,
        const RaggedRightArrayKokkos <double> &strength_global_vars,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const size_t num_material_points,
        const size_t mat_id)
    {


    }  // end of init_strength_state_vars


    /**
     * @brief Calculate the deviatoric stress tensor for the hypo-plasticity response.
     *
     * Evolves the deviatoric stress tensor using the Jaumann rate and applies J2 plasticity yield condition.
     *
     * @param GaussPoints_vel_grad   Velocity gradient at Gauss points.
     * @param node_coords   Coordinates of the nodes.
     * @param node_vel   Velocities of the nodes.
     * @param nodes_in_elem   Node indices in the element.
     * @param MaterialPoints_pres   Pressure at material points.
     * @param MaterialPoints_stress   Stress tensor at material points (output).
     * @param MaterialPoints_stress_n0   Stress tensor at material points from previous time step.
     * @param MaterialPoints_sspd   Sound speed at material points.
     * @param MaterialPoints_eos_state_vars   State variables for the equation of state at material points.
     * @param MaterialPoints_strength_state_vars   State variables for the strength model at material points.
     * @param MaterialPoints_den   Density at the material point.
     * @param MaterialPoints_sie   Specific internal energy at the material point.
     * @param MaterialPoints_shear_modulii   Shear modulus at material points.
     * @param elem_in_mat_elem   Mapping from material points to mesh elements.
     * @param eos_global_vars   Global variables for the equation of state.
     * @param strength_global_vars   Global variables for the strength model.
     * @param vol   Volume of the element.
     * @param dt   Time step size.
     * @param rk_alpha   Runge-Kutta time integration coefficient.
     * @param time   Current simulation time.
     * @param cycle   Current simulation cycle.
     * @param MaterialPoints_lid   Local ID of the material point.
     * @param mat_id   Material ID.
     * @param gauss_gid   Global ID of the Gauss point.
     * @param elem_gid   Global ID of the element.
     */
    KOKKOS_FUNCTION
    static void calc_stress(
        const DCArrayKokkos<double>  &GaussPoints_vel_grad,
        const DistributedDCArray<double> &node_coords,
        const DistributedDCArray<double> &node_vel,
        const DistributedDCArray<size_t>  &nodes_in_elem,
        const DRaggedRightArrayKokkos<double>  &MaterialPoints_pres,
        const DRaggedRightArrayKokkos<double>  &MaterialPoints_stress,
        const DRaggedRightArrayKokkos<double>  &MaterialPoints_stress_n0,
        const DRaggedRightArrayKokkos<double>  &MaterialPoints_sspd,
        const DRaggedRightArrayKokkos <double> &MaterialPoints_eos_state_vars,
        const DRaggedRightArrayKokkos <double> &MaterialPoints_strength_state_vars,
        const double MaterialPoints_den,
        const double MaterialPoints_sie,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_shear_modulii,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const RaggedRightArrayKokkos <double> &eos_global_vars,
        const RaggedRightArrayKokkos <double> &strength_global_vars,
        const double vol,
        const double dt,
        const double rk_alpha,
        const double time,
        const size_t cycle,
        const size_t MaterialPoints_lid,
        const size_t mat_id,
        const size_t gauss_gid,
        const size_t elem_gid)
    {

        // -----------------------------------------------------------------------------
        // Required variables are here
        // ------------------------------------------------------------------------------

        // statev(0) = var_1
        //   :
        //   :
        //   :
        // statev(N) = var_N

        // --strength
        // statev(0) = shear modulus
        // statev(1) = yield strength
    

        // shear modulus
        double G_stress = strength_global_vars(mat_id, 0);

        // yield strength
        double Y_stress = strength_global_vars(mat_id, 1);
        
        const int num_dims = 3; // for 3D only


        // -----------------------------------------------------------------------------
        // The user must coding goes here
        // -----------------------------------------------------------------------------

        // Cauchy deviatoric stress tensor, stress(rk,elem,i,j)

        // create a slice of the stress tensor
        ViewCArrayKokkos <double> stress_n(&MaterialPoints_stress_n0(mat_id, MaterialPoints_lid, 0, 0), num_dims, num_dims); // stress at previous timestep
        ViewCArrayKokkos <double> stress(&MaterialPoints_stress(mat_id, MaterialPoints_lid, 0, 0), num_dims, num_dims); // stress at this timestep

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
                D(i, j) = 0.5 * (GaussPoints_vel_grad(gauss_gid,i, j) + GaussPoints_vel_grad(gauss_gid,j, i));
            }
        }
        
        double div = 0;
        for (int i = 0; i < num_dims; i++) {
            div += GaussPoints_vel_grad(gauss_gid,i, i);
        }
        
        
        // W = 1/2(GaussPoints_vel_grad gauss_gid,- GaussPoints_vel_grad^gauss_gid,T)
        for (int i = 0; i < num_dims; i++) {
            for (int j = 0; j < num_dims; j++) {
                W(i, j) = 0.5 * (GaussPoints_vel_grad(gauss_gid,i, j) - GaussPoints_vel_grad(gauss_gid,j, i));
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

} // end namespace


/////////////////////////////////////////////////////////////////////////////
///
/// \fn HypoElasticPlasticRZ
///
/// \brief stress model for hypo elastic plastic response
///
///  This is a material model function for the deviatoric stress tensor
///
/// \param Material pressure
/// \param Material stress
/// \param Global ID for the material
/// \param Material ID for the element
/// \param Material state variables
/// \param Material Sound speed
/// \param Material density
/// \param Material specific internal energy
/// \param Element velocity gradient
/// \param Element nodes IDs in the element
/// \param Node node coordinates
/// \param Noe velocity 
/// \param Element volume
/// \param Time time step
/// \param Time coefficient in the Runge Kutta time integration step
///
/////////////////////////////////////////////////////////////////////////////
namespace HypoPlasticityRZModel {

    const double fuzz = 1.e-16;

    /**
     * @brief Initialize the strength state variables for the material points (RZ version).
     *
     * @param MaterialPoints_eos_state_vars   State variables for the equation of state at each material point.
     * @param MaterialPoints_strength_state_vars   State variables for the strength model at each material point.
     * @param eos_global_vars   Global variables for the equation of state.
     * @param strength_global_vars   Global variables for the strength model.
     * @param elem_in_mat_elem   Mapping from material points to mesh elements.
     * @param num_material_points   Number of material points.
     * @param mat_id   Material ID.
     */
    static void init_strength_state_vars(
        const DRaggedRightArrayKokkos <double> &MaterialPoints_eos_state_vars,
        const DRaggedRightArrayKokkos <double> &MaterialPoints_strength_state_vars,
        const RaggedRightArrayKokkos <double> &eos_global_vars,
        const RaggedRightArrayKokkos <double> &strength_global_vars,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const size_t num_material_points,
        const size_t mat_id)
    {



    }  // end of init_strength_state_vars


    /**
     * @brief Calculate the deviatoric stress tensor for the hypo-elastic-plastic response in RZ geometry.
     *
     * Evolves the deviatoric stress tensor using the Jaumann rate and applies J2 plasticity yield condition for RZ geometry.
     *
     * @param GaussPoints_vel_grad   Velocity gradient at Gauss points.
     * @param node_coords   Coordinates of the nodes.
     * @param node_vel   Velocities of the nodes.
     * @param nodes_in_elem   Node indices in the element.
     * @param MaterialPoints_pres   Pressure at material points.
     * @param MaterialPoints_stress   Stress tensor at material points (output).
     * @param MaterialPoints_stress_n0   Stress tensor at material points from previous time step.
     * @param MaterialPoints_sspd   Sound speed at material points.
     * @param MaterialPoints_eos_state_vars   State variables for the equation of state at material points.
     * @param MaterialPoints_strength_state_vars   State variables for the strength model at material points.
     * @param MaterialPoints_den   Density at the material point.
     * @param MaterialPoints_sie   Specific internal energy at the material point.
     * @param MaterialPoints_shear_modulii   Shear modulus at material points.
     * @param elem_in_mat_elem   Mapping from material points to mesh elements.
     * @param eos_global_vars   Global variables for the equation of state.
     * @param strength_global_vars   Global variables for the strength model.
     * @param vol   Volume of the element.
     * @param dt   Time step size.
     * @param rk_alpha   Runge-Kutta time integration coefficient.
     * @param time   Current simulation time.
     * @param cycle   Current simulation cycle.
     * @param MaterialPoints_lid   Local ID of the material point.
     * @param mat_id   Material ID.
     * @param gauss_gid   Global ID of the Gauss point.
     * @param elem_gid   Global ID of the element.
     */
    KOKKOS_FUNCTION
    static void calc_stress(
        const DCArrayKokkos<double>  &GaussPoints_vel_grad,
        const DistributedDCArray<double> &node_coords,
        const DistributedDCArray<double> &node_vel,
        const DistributedDCArray<size_t>  &nodes_in_elem,
        const DRaggedRightArrayKokkos<double>  &MaterialPoints_pres,
        const DRaggedRightArrayKokkos<double>  &MaterialPoints_stress,
        const DRaggedRightArrayKokkos<double>  &MaterialPoints_stress_n0,
        const DRaggedRightArrayKokkos<double>  &MaterialPoints_sspd,
        const DRaggedRightArrayKokkos <double> &MaterialPoints_eos_state_vars,
        const DRaggedRightArrayKokkos <double> &MaterialPoints_strength_state_vars,
        const double MaterialPoints_den,
        const double MaterialPoints_sie,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_shear_modulii,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const RaggedRightArrayKokkos <double> &eos_global_vars,
        const RaggedRightArrayKokkos <double> &strength_global_vars,
        const double vol,
        const double dt,
        const double rk_alpha,
        const double time,
        const size_t cycle,
        const size_t MaterialPoints_lid,
        const size_t mat_id,
        const size_t gauss_gid,
        const size_t elem_gid)
    {

        const int num_dims = 3; // must be 3 even though its 2D RZ

        // -----------------------------------------------------------------------------
        // The user must coding goes here
        //------------------------------------------------------------------------------
        // Cauchy deviatoric stress tensor, stress(rk,elem,i,j)

        // create a slice of the stress tensor
        ViewCArrayKokkos <double> stress_n(&MaterialPoints_stress_n0(mat_id, MaterialPoints_lid, 0, 0), 3, 3); // stress from previous timestep
        ViewCArrayKokkos <double> stress(&MaterialPoints_stress(mat_id, MaterialPoints_lid, 0, 0), 3, 3); // stress from this timestep

        double s_dev_n_1D[3 * 3];
        double s_dev_1D[3 * 3];

        ViewCArrayKokkos <double> s_dev_n(&s_dev_n_1D[0], 3, 3); // deviatoric stress from previous timestep
        ViewCArrayKokkos <double> s_dev(&s_dev_1D[0], 3, 3); // deviatoric stress from this timestep

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
        ViewCArrayKokkos <double> stress_dev_rhs(&stress_dev_rhs_1d[0], 3, 3);

        // symmetric part of velocity gradient
        double D_1d[3 * 3];
        ViewCArrayKokkos <double> D(&D_1d[0], 3, 3);
        
        // anti-symmetric part of velocity gradient
        double W_1d[num_dims * num_dims];
        ViewCArrayKokkos <double> W(&W_1d[0], num_dims, num_dims);
        
        
        // D = 1/2(vel_grad + vel_grad^T)
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                D(i, j) = 0.5 * (GaussPoints_vel_grad(gauss_gid,i, j) + GaussPoints_vel_grad(gauss_gid,j, i));
            }
        }
        
        // W = 1/2(vel_grad - vel_grad^T)
        for (int i = 0; i < num_dims; i++) {
            for (int j = 0; j < num_dims; j++) {
                W(i, j) = 0.5 * (GaussPoints_vel_grad(gauss_gid,i, j) - GaussPoints_vel_grad(gauss_gid,j, i));
            }
        }
        
        
        double div = 0;
        for (int i = 0; i < num_dims; i++) {
            div += GaussPoints_vel_grad(gauss_gid,i, i);
        }

        // vel = [u,v]
        // GaussPoints_vel_grad gauss_gid,= [du/dx,  du/dy]
        //            [dv/dx,  dv/dy]
        // GaussPoints_vel_grad(gauss_gid,1,0) = dv/dx
        // GaussPoints_vel_grad(gauss_gid,0,1) = du/dy
        double curl;
        curl = GaussPoints_vel_grad(gauss_gid,1,0) - GaussPoints_vel_grad(gauss_gid,0,1);  // dv/dx - du/dy
        

        // shear modulus
        double G_stress = strength_global_vars(mat_id, 0);

        // yield strength
        double Y_stress = strength_global_vars(mat_id, 1);

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

} // end namespace



#endif // end Header Guard


