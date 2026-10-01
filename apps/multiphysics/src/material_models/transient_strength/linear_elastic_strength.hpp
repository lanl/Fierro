/**********************************************************************************************
© 2020. Triad National Security, LLC. All rights reserved.
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

#ifndef USER_LINEAR_ELASTIC_H
#define USER_LINEAR_ELASTIC_H


/////////////////////////////////////////////////////////////////////////////
///
/// \fn OrthotropicLinearElasticStrengthModel
///
/// \brief user defined strength model
///
///  This is the user material model function for the stress tensor
///
/// \param Element pressure
/// \param Element stress
/// \param Global ID for the element
/// \param Material ID for the element
/// \param Element state variables
/// \param Element Sound speed
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
namespace OrthotropicLinearElasticStrengthModel {

    // ------------------------
    // helper functions 

    // C = A * B  for 3x3 matrices
    template <typename T1>
    KOKKOS_INLINE_FUNCTION
    void matmul3x3(const T1& A,
                   const T1& B,
                   const T1& C) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double sum = 0.0;
                for (int k = 0; k < 3; k++) {
                    sum += A(i,k) * B(k,j);
                }
                C(i,j) = sum;
            }
        }
    } // end function


    // B = A^T for 3x3
    template <typename T1>
    KOKKOS_INLINE_FUNCTION
    void transpose3x3(const T1& A, 
                      const T1& B) {
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++){
                B(i,j) = A(j,i);
        }
    } // end function


    // Generic NxN matrix multiply helper (used for 6x6 here)
    template <typename T1>
    KOKKOS_INLINE_FUNCTION
    void matmulNxN(const T1& A,
                   const T1& B,
                   const T1& C,
                   const int N) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double sum = 0.0;
                for (int k = 0; k < N; k++) {
                    sum += A(i,k) * B(k,j);
                }
                C(i,j) = sum;
            }
        }
    } // end fucntion


    // B = A^T for NxN
    template <typename T1>
    KOKKOS_INLINE_FUNCTION
    void transposeNxN(const T1& A, 
                      const T1& B, 
                      const int N) {
        for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++){
                B(i,j) = A(j,i);
        }
    }   



    /**
    * Symmetric Eigen-Decomposition: sqrt(F^T*F)
    * @brief Jacobi eigenvalue algorithm for a symmetric 3x3 matrix.
    *        Computes eigenvalues (eigval) and eigenvectors (columns of eigvec).
    *        A is overwritten during iteration (pass a copy if you need to keep it).
    */
    KOKKOS_INLINE_FUNCTION
    void jacobiEigenSolve3x3(const ViewCArrayKokkos<double> A,       // 3x3 symmetric, copied
                             const ViewCArrayKokkos<double>& eigval, // size 3
                             const ViewCArrayKokkos<double>& eigvec) // 3x3, columns = eigenvectors
    {
        // Initialize eigvec to identity
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                eigvec(i,j) = (i == j) ? 1.0 : 0.0;

        const int max_sweeps = 50;
        const double tol = 1.0e-14;

        for (int sweep = 0; sweep < max_sweeps; sweep++) {
            // Off-diagonal norm (convergence check)
            double off = fabs(A(0,1)) + fabs(A(0,2)) + fabs(A(1,2));
            if (off < tol) break;

            // Loop over the 3 off-diagonal pairs (p,q)
            const int pairs[3][2] = {{0,1}, {0,2}, {1,2}};

            for (auto& pq : pairs) {
                int p = pq[0], q = pq[1];

                if (fabs(A(p,q)) < 1.0e-18) continue;

                double theta = (A(q,q) - A(p,p)) / (2.0 * A(p,q));
                double t = (theta >= 0.0 ? 1.0 : -1.0) /
                        (fabs(theta) + sqrt(theta*theta + 1.0));
                double c = 1.0 / sqrt(t*t + 1.0);
                double s = t * c;

                double App = A(p,p), Aqq = A(q,q), Apq = A(p,q);

                A(p,p) = App - t * Apq;
                A(q,q) = Aqq + t * Apq;
                A(p,q) = 0.0;
                A(q,p) = 0.0;

                for (int k = 0; k < 3; k++) {
                    if (k != p && k != q) {
                        double Akp = A(k,p), Akq = A(k,q);
                        A(k,p) = c * Akp - s * Akq;
                        A(p,k) = A(k,p);
                        A(k,q) = s * Akp + c * Akq;
                        A(q,k) = A(k,q);
                    }
                }

                // Update eigenvector matrix
                for (int k = 0; k < 3; k++) {
                    double Vkp = eigvec(k,p), Vkq = eigvec(k,q);
                    eigvec(k,p) = c * Vkp - s * Vkq;
                    eigvec(k,q) = s * Vkp + c * Vkq;
                } // end for k

            } // end for pairs
        } // end for sweep

        eigval(0) = A(0,0);
        eigval(1) = A(1,1);
        eigval(2) = A(2,2);
    }  // end function
    

    /**
    * @brief Compute the rotation tensor R from the right polar decomposition
    *        F = R * U, where U = sqrt(F^T F).
    *
    * @param F  Deformation gradient (3x3 MATAR CArray)
    * @param R  Output rotation tensor (3x3 MATAR CArray, must be pre-allocated)
    */
    KOKKOS_INLINE_FUNCTION
    void computeRotationFromF(const ViewCArrayKokkos<double>& F, 
                              const ViewCArrayKokkos<double>& R) {

        double Ft_1D[9];
        ViewCArrayKokkos<double>Ft(&Ft_1D[0],3,3);

        double C_1D[9];
        ViewCArrayKokkos<double>C(&C_1D[0],3,3);

        transpose3x3(F, Ft);
        matmul3x3(Ft, F, C);   // C = F^T F  (right Cauchy-Green tensor)

        double eigval_1D[3];
        ViewCArrayKokkos<double>eigval(&eigval_1D[0],3);

        double eigvec_1D[9];
        ViewCArrayKokkos<double>eigvec(&eigvec_1D[0],3,3); // columns are eigenvectors of C

        jacobiEigenSolve3x3(C, eigval, eigvec);

        // Build U = V * diag(sqrt(lambda)) * V^T
        // Build U_inv = V * diag(1/sqrt(lambda)) * V^T

        double U_1D[9];
        ViewCArrayKokkos<double>U(&U_1D[0],3,3);

        double Uinv_1D[9];
        ViewCArrayKokkos<double>Uinv(&Uinv_1D[0],3,3);

        double V_1D[9];
        ViewCArrayKokkos<double>V(&V_1D[0],3,3);

        double Vt_1D[9];
        ViewCArrayKokkos<double>Vt(&Vt_1D[0],3,3);

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                V(i,j) = eigvec(i,j);

        transpose3x3(V, Vt);


        double sqrtLambda_1D[9];
        ViewCArrayKokkos<double>sqrtLambda(&sqrtLambda_1D[0],3,3);

        double invSqrtLambda_1D[9];
        ViewCArrayKokkos<double>invSqrtLambda(&invSqrtLambda_1D[0],3,3);

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                sqrtLambda(i,j) = 0.0;
                invSqrtLambda(i,j) = 0.0;
            }

        for (int i = 0; i < 3; i++) {
            double lam = fmax(eigval(i), 1.0e-30); // guard against negative/zero
            sqrtLambda(i,i) = sqrt(lam);
            invSqrtLambda(i,i) = 1.0 / sqrt(lam);
        }

        double temporary_1D[9];
        ViewCArrayKokkos<double>temporary(&temporary_1D[0],3,3);

        matmul3x3(V, sqrtLambda, temporary);
        matmul3x3(temporary, Vt, U);

        matmul3x3(V, invSqrtLambda, temporary);
        matmul3x3(temporary, Vt, Uinv);

        // R = F * U^{-1}
        matmul3x3(F, Uinv, R);

    } // end function

    /**
    * @brief Build the 6x6 Bond stress-transformation matrix M from a 3x3 rotation R.
    *        Rotates Voigt-notation stiffness via: C_current = M * C_ref * M^T
    */
    KOKKOS_INLINE_FUNCTION
    void buildBondMatrix(const ViewCArrayKokkos<double>& R, 
                         const ViewCArrayKokkos<double>& M) {

        // Convenience aliases for rotation matrix components
        double R11=R(0,0); 
        double R12=R(0,1); 
        double R13=R(0,2);

        double R21=R(1,0); 
        double R22=R(1,1); 
        double R23=R(1,2);

        double R31=R(2,0); 
        double R32=R(2,1); 
        double R33=R(2,2);

        // Row 1
        M(0,0)=R11*R11; 
        M(0,1)=R12*R12; 
        M(0,2)=R13*R13;
        M(0,3)=2.0*R12*R13; 
        M(0,4)=2.0*R13*R11; 
        M(0,5)=2.0*R11*R12;

        // Row 2
        M(1,0)=R21*R21; 
        M(1,1)=R22*R22; 
        M(1,2)=R23*R23;
        M(1,3)=2.0*R22*R23; 
        M(1,4)=2.0*R23*R21; 
        M(1,5)=2.0*R21*R22;

        // Row 3
        M(2,0)=R31*R31; 
        M(2,1)=R32*R32; 
        M(2,2)=R33*R33;
        M(2,3)=2.0*R32*R33; 
        M(2,4)=2.0*R33*R31; 
        M(2,5)=2.0*R31*R32;

        // Row 4
        M(3,0)=R21*R31; 
        M(3,1)=R22*R32; 
        M(3,2)=R23*R33;
        M(3,3)=R22*R33+R23*R32; 
        M(3,4)=R23*R31+R21*R33; 
        M(3,5)=R21*R32+R22*R31;

        // Row 5
        M(4,0)=R31*R11; 
        M(4,1)=R32*R12; 
        M(4,2)=R33*R13;
        M(4,3)=R32*R13+R33*R12; 
        M(4,4)=R33*R11+R31*R13; 
        M(4,5)=R31*R12+R32*R11;

        // Row 6
        M(5,0)=R11*R21; 
        M(5,1)=R12*R22; 
        M(5,2)=R13*R23;
        M(5,3)=R12*R23+R13*R22; 
        M(5,4)=R13*R21+R11*R23; 
        M(5,5)=R11*R22+R12*R21;
    } // end function



    /**
    * @brief Rotate a 6x6 Voigt-notation stiffness matrix from reference to
    *        current configuration using the deformation gradient F. We seek
    *        C_current = M * C_ref * M^T
    *
    * @param C_ref       Reference stiffness matrix (6x6, Voigt notation)
    * @param F           Deformation gradient (3x3)
    * @param C_current   Output: rotated stiffness matrix (6x6)
    */
    inline void rotateStiffnessFromF(const ViewCArrayKokkos<double>& C_ref,
                                     const ViewCArrayKokkos<double>& F,
                                     const ViewCArrayKokkos<double>& C_current) {

        double R_1D[9];
        ViewCArrayKokkos<double>R(&R_1D[0],3,3);

        computeRotationFromF(F, R);

        double M_1D[36];
        ViewCArrayKokkos<double>M(&M_1D[0],6,6);

        double Mt_1D[36];
        ViewCArrayKokkos<double>Mt(&Mt_1D[0],6,6);

        double temporary_1D[36];
        ViewCArrayKokkos<double>temporary(&temporary_1D[0],6,6);


        buildBondMatrix(R, M);
        transposeNxN(M, Mt, 6);

        matmulNxN(M, C_ref, temporary, 6);          // temp = M * C_ref
        matmulNxN(temporary, Mt, C_current, 6);     // C_current = temp * M^T
    } // end function


    /*
        // routines needed to calculate things like the wave speed in a material

        // Reference stiffness tensor (e.g., cubic crystal, Voigt notation)
        double C_ref_1D[36];
        ViewCArrayKokkos<double> C_ref(&C_ref_1D[0],6,6);
        // ... fill C_ref with your material's reference stiffness values ...

        // Deformation gradient from your solver (already a MATAR RaggeCArray)

        ViewCArrayKokkos<double> F(&Deformation_Grad(...),3,3);
        F(0,0)=1.05; F(0,1)=0.02; F(0,2)=0.00;
        F(1,0)=0.00; F(1,1)=0.98; F(1,2)=0.01;
        F(2,0)=0.01; F(2,1)=0.00; F(2,2)=1.00;

        double C_current_1D[36];
        ViewCArrayKokkos<double> C_current(&C_current_1D[0],6,6);

        rotateStiffnessFromF(C_ref, F, C_current);

        // C_current now holds the stiffness tensor rotated into the
        // current lab-frame orientation, ready for use in the acoustic
        // tensor Q(n) to compute directional sound speed.
    */


    /**
    * @brief Compute Green-Lagrange strain E = 1/2(F^T F - I)
    */
    KOKKOS_INLINE_FUNCTION
    void computeGreenLagrangeStrain(const ViewCArrayKokkos<double>& F,
                                    const ViewCArrayKokkos<double>& E) {
        
        double Ft_1D[9];
        ViewCArrayKokkos<double>Ft(&Ft_1D[0],3,3);

        double C_1D[9];
        ViewCArrayKokkos<double>C(&C_1D[0],3,3);

        transpose3x3(F, Ft);
        matmul3x3(Ft, F, C);   // C = F^T F

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                E(i,j) = 0.5 * (C(i,j) - (i == j ? 1.0 : 0.0));
            }
        }
    } // end function

    /**
    * @brief Convert symmetric 3x3 strain tensor to 6x1 Voigt vector
    *        using engineering strain convention for shear terms.
    */
    KOKKOS_INLINE_FUNCTION
    void strainTensorToVoigt(const ViewCArrayKokkos<double>& E, 
                             const ViewCArrayKokkos<double>& E_voigt) {
        E_voigt(0) = E(0,0);          // E11
        E_voigt(1) = E(1,1);          // E22
        E_voigt(2) = E(2,2);          // E33
        E_voigt(3) = 2.0 * E(1,2);    // 2*E23
        E_voigt(4) = 2.0 * E(0,2);    // 2*E13
        E_voigt(5) = 2.0 * E(0,1);    // 2*E12
    }


    /**
    * @brief Compute 2nd Piola-Kirchhoff stress (Voigt) from reference
    *        stiffness and Green-Lagrange strain (Voigt).
    *           S = C_ref * E
    */
    KOKKOS_INLINE_FUNCTION
    void computeSecondPKStressVoigt(const ViewCArrayKokkos<double>& C_ref,
                                    const ViewCArrayKokkos<double>& E_voigt,
                                    const ViewCArrayKokkos<double>& S_voigt) {

        for (int i = 0; i < 6; i++) {
            double sum = 0.0;
            for (int j = 0; j < 6; j++) {
                sum += C_ref(i,j) * E_voigt(j);
            }
            S_voigt(i) = sum;
        }
    } // end function


    KOKKOS_INLINE_FUNCTION
    void stressVoigtToTensor(const ViewCArrayKokkos<double>& S_voigt, 
                             const ViewCArrayKokkos<double>& S) {
        S(0,0) = S_voigt(0);
        S(1,1) = S_voigt(1);
        S(2,2) = S_voigt(2);
        S(1,2) = S_voigt(3);  S(2,1) = S_voigt(3);
        S(0,2) = S_voigt(4);  S(2,0) = S_voigt(4);
        S(0,1) = S_voigt(5);  S(1,0) = S_voigt(5);
    }


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


    KOKKOS_FUNCTION
    static void calc_stress(
        const DCArrayKokkos<double>  &GaussPoints_vel_grad,
        const MPICArrayKokkos<double> &node_coords,
        const MPICArrayKokkos<double> &node_coords_t0,
        const MPICArrayKokkos <double> &node_vel,
        const DCArrayKokkos<size_t>  &nodes_in_elem,
        const DRaggedRightArrayKokkos<double>  &MaterialPoints_pres,
        const DRaggedRightArrayKokkos<double>  &MaterialPoints_stress,
        const DRaggedRightArrayKokkos<double>  &MaterialPoints_stress_n0,
        const DRaggedRightArrayKokkos<double>  &MaterialPoints_sspd,
        const DRaggedRightArrayKokkos <double> &MaterialPoints_eos_state_vars,
        const DRaggedRightArrayKokkos <double> &MaterialPoints_strength_state_vars,
        const double MaterialPoints_den,
        const double MaterialPoints_sie,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_deformation_grad,
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

        double C[6][6];


        return;
    } // end of user mat



    

    
    static void destroy(
        const DRaggedRightArrayKokkos <double> &MaterialPoints_eos_state_vars,
        const DRaggedRightArrayKokkos <double> &MaterialPoints_strength_state_vars,
        const RaggedRightArrayKokkos <double> &eos_global_vars,
        const RaggedRightArrayKokkos <double> &strength_global_vars,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const size_t num_material_points,
        const size_t mat_ids)
    {

    } // end destory

} // end namespace



#endif // end Header Guard