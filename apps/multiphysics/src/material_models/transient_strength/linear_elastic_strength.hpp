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

#include "cramers_rule.hpp"

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
    // One template parameter per argument: the operands may be views over
    // different tiers (state storage vs. a stack scratch buffer). The
    // accumulator below stays double regardless.
    template <typename T1, typename T2, typename T3>
    KOKKOS_INLINE_FUNCTION
    void matmul3x3(const T1& A,
                   const T2& B,
                   const T3& C) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double sum = real_t(0.0);
                for (int k = 0; k < 3; k++) {
                    sum += A(i,k) * B(k,j);
                }
                C(i,j) = sum;
            }
        }
    } // end function


    // B = A^T for 3x3
    template <typename T1, typename T2>
    KOKKOS_INLINE_FUNCTION
    void transpose3x3(const T1& A, 
                      const T2& B) {
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
                double sum = real_t(0.0);
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
    * Symmetric Eigen-Decomposition: Kokkos::sqrt(F^T*F)
    * @brief Jacobi eigenvalue algorithm for a symmetric 3x3 matrix.
    *        Computes eigenvalues (eigval) and eigenvectors (columns of eigvec).
    *        A is overwritten during iteration (pass a copy if you need to keep it).
    */
    KOKKOS_INLINE_FUNCTION
    void jacobiEigenSolve3x3(const ViewCArrayKokkos<real_t> A,       // 3x3 symmetric, copied
                             const ViewCArrayKokkos<real_t>& eigval, // size 3
                             const ViewCArrayKokkos<real_t>& eigvec) // 3x3, columns = eigenvectors
    {
        // Initialize eigvec to identity
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                eigvec(i,j) = (i == j) ? real_t(1.0) : real_t(0.0);

        const int max_sweeps = 50;
        const real_t tol = real_t(1.0e-14);

        for (int sweep = 0; sweep < max_sweeps; sweep++) {
            // Off-diagonal norm (convergence check)
            real_t off = Kokkos::fabs(A(0,1)) + Kokkos::fabs(A(0,2)) + Kokkos::fabs(A(1,2));
            if (off < tol) break;

            // Loop over the 3 off-diagonal pairs (p,q)
            const int pairs[3][2] = {{0,1}, {0,2}, {1,2}};

            for (auto& pq : pairs) {
                int p = pq[0], q = pq[1];

                if (Kokkos::fabs(A(p,q)) < real_t(1.0e-18)) continue;

                real_t theta = (A(q,q) - A(p,p)) / (real_t(2.0) * A(p,q));
                real_t t = (theta >= real_t(0.0) ? real_t(1.0) : -real_t(1.0)) /
                        (Kokkos::fabs(theta) + Kokkos::sqrt(theta*theta + real_t(1.0)));
                real_t c = real_t(1.0) / Kokkos::sqrt(t*t + real_t(1.0));
                real_t s = t * c;

                real_t App = A(p,p), Aqq = A(q,q), Apq = A(p,q);

                A(p,p) = App - t * Apq;
                A(q,q) = Aqq + t * Apq;
                A(p,q) = real_t(0.0);
                A(q,p) = real_t(0.0);

                for (int k = 0; k < 3; k++) {
                    if (k != p && k != q) {
                        real_t Akp = A(k,p), Akq = A(k,q);
                        A(k,p) = c * Akp - s * Akq;
                        A(p,k) = A(k,p);
                        A(k,q) = s * Akp + c * Akq;
                        A(q,k) = A(k,q);
                    }
                }

                // Update eigenvector matrix
                for (int k = 0; k < 3; k++) {
                    real_t Vkp = eigvec(k,p), Vkq = eigvec(k,q);
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
    *        F = R * U, where U = Kokkos::sqrt(F^T F).
    *
    * @param F  Deformation gradient (3x3 MATAR CArray)
    * @param R  Output rotation tensor (3x3 MATAR CArray, must be pre-allocated)
    */
    KOKKOS_INLINE_FUNCTION
    void computeRotationFromF(const ViewCArrayKokkos<real_t>& F, 
                              const ViewCArrayKokkos<real_t>& R) {

        real_t Ft_1D[9];
        ViewCArrayKokkos<real_t>Ft(&Ft_1D[0],3,3);

        real_t C_1D[9];
        ViewCArrayKokkos<real_t>C(&C_1D[0],3,3);

        transpose3x3(F, Ft);
        matmul3x3(Ft, F, C);   // C = F^T F  (right Cauchy-Green tensor)

        real_t eigval_1D[3];
        ViewCArrayKokkos<real_t>eigval(&eigval_1D[0],3);

        real_t eigvec_1D[9];
        ViewCArrayKokkos<real_t>eigvec(&eigvec_1D[0],3,3); // columns are eigenvectors of C

        jacobiEigenSolve3x3(C, eigval, eigvec);

        // Build U = V * diag(Kokkos::sqrt(lambda)) * V^T
        // Build U_inv = V * diag(1/Kokkos::sqrt(lambda)) * V^T

        real_t U_1D[9];
        ViewCArrayKokkos<real_t>U(&U_1D[0],3,3);

        real_t Uinv_1D[9];
        ViewCArrayKokkos<real_t>Uinv(&Uinv_1D[0],3,3);

        real_t V_1D[9];
        ViewCArrayKokkos<real_t>V(&V_1D[0],3,3);

        real_t Vt_1D[9];
        ViewCArrayKokkos<real_t>Vt(&Vt_1D[0],3,3);

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                V(i,j) = eigvec(i,j);

        transpose3x3(V, Vt);


        real_t sqrtLambda_1D[9];
        ViewCArrayKokkos<real_t>sqrtLambda(&sqrtLambda_1D[0],3,3);

        real_t invSqrtLambda_1D[9];
        ViewCArrayKokkos<real_t>invSqrtLambda(&invSqrtLambda_1D[0],3,3);

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                sqrtLambda(i,j) = real_t(0.0);
                invSqrtLambda(i,j) = real_t(0.0);
            }

        for (int i = 0; i < 3; i++) {
            real_t lam = Kokkos::fmax(eigval(i), real_t(1.0e-30)); // guard against negative/zero
            sqrtLambda(i,i) = Kokkos::sqrt(lam);
            invSqrtLambda(i,i) = real_t(1.0) / Kokkos::sqrt(lam);
        }

        real_t temporary_1D[9];
        ViewCArrayKokkos<real_t>temporary(&temporary_1D[0],3,3);

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
    void buildBondMatrix(const ViewCArrayKokkos<real_t>& R, 
                         const ViewCArrayKokkos<real_t>& M) {

        // Convenience aliases for rotation matrix components
        real_t R11=R(0,0); 
        real_t R12=R(0,1); 
        real_t R13=R(0,2);

        real_t R21=R(1,0); 
        real_t R22=R(1,1); 
        real_t R23=R(1,2);

        real_t R31=R(2,0); 
        real_t R32=R(2,1); 
        real_t R33=R(2,2);

        // Row 1
        M(0,0)=R11*R11; 
        M(0,1)=R12*R12; 
        M(0,2)=R13*R13;
        M(0,3)=real_t(2.0)*R12*R13; 
        M(0,4)=real_t(2.0)*R13*R11; 
        M(0,5)=real_t(2.0)*R11*R12;

        // Row 2
        M(1,0)=R21*R21; 
        M(1,1)=R22*R22; 
        M(1,2)=R23*R23;
        M(1,3)=real_t(2.0)*R22*R23; 
        M(1,4)=real_t(2.0)*R23*R21; 
        M(1,5)=real_t(2.0)*R21*R22;

        // Row 3
        M(2,0)=R31*R31; 
        M(2,1)=R32*R32; 
        M(2,2)=R33*R33;
        M(2,3)=real_t(2.0)*R32*R33; 
        M(2,4)=real_t(2.0)*R33*R31; 
        M(2,5)=real_t(2.0)*R31*R32;

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
    KOKKOS_INLINE_FUNCTION
    void rotateStiffnessFromF(const ViewCArrayKokkos<real_t>& C_ref,
                                     const ViewCArrayKokkos<real_t>& F,
                                     const ViewCArrayKokkos<real_t>& C_current) {

        real_t R_1D[9];
        ViewCArrayKokkos<real_t>R(&R_1D[0],3,3);

        computeRotationFromF(F, R);

        real_t M_1D[36];
        ViewCArrayKokkos<real_t>M(&M_1D[0],6,6);

        real_t Mt_1D[36];
        ViewCArrayKokkos<real_t>Mt(&Mt_1D[0],6,6);

        real_t temporary_1D[36];
        ViewCArrayKokkos<real_t>temporary(&temporary_1D[0],6,6);


        buildBondMatrix(R, M);
        transposeNxN(M, Mt, 6);

        matmulNxN(M, C_ref, temporary, 6);          // temp = M * C_ref
        matmulNxN(temporary, Mt, C_current, 6);     // C_current = temp * M^T
    } // end function


    /*
        // routines needed to calculate things like the wave speed in a material

        // Reference stiffness tensor (e.g., cubic crystal, Voigt notation)
        real_t C_ref_1D[36];
        ViewCArrayKokkos<real_t> C_ref(&C_ref_1D[0],6,6);
        // ... fill C_ref with your material's reference stiffness values ...

        // Deformation gradient from your solver (already a MATAR RaggeCArray)

        ViewCArrayKokkos<real_t> F(&Deformation_Grad(...),3,3);
        F(0,0)=1.05; F(0,1)=0.02; F(0,2)=0.00;
        F(1,0)=0.00; F(1,1)=0.98; F(1,2)=0.01;
        F(2,0)=0.01; F(2,1)=0.00; F(2,2)=1.00;

        real_t C_current_1D[36];
        ViewCArrayKokkos<real_t> C_current(&C_current_1D[0],6,6);

        rotateStiffnessFromF(C_ref, F, C_current);

        // C_current now holds the stiffness tensor rotated into the
        // current lab-frame orientation, ready for use in the acoustic
        // tensor Q(n) to compute directional sound speed.
    */



    /**
    * @brief Compute Green-Lagrange strain E = 1/2(F^T F - I)
    */
    template <typename TF, typename TE>
    KOKKOS_INLINE_FUNCTION
    void computeGreenLagrangeStrain(const TF& F,
                                    const TE& E) {
        
        real_t Ft_1D[9];
        ViewCArrayKokkos<real_t>Ft(&Ft_1D[0],3,3);

        real_t C_1D[9];
        ViewCArrayKokkos<real_t>C(&C_1D[0],3,3);

        transpose3x3(F, Ft);
        matmul3x3(Ft, F, C);   // C = F^T F

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                E(i,j) = real_t(0.5) * (C(i,j) - (i == j ? real_t(1.0) : real_t(0.0)));
            }
        }
    } // end function

    /**
    * @brief Convert symmetric 3x3 strain tensor to 6x1 Voigt vector
    *        using engineering strain convention for shear terms.
    */
    KOKKOS_INLINE_FUNCTION
    void strainTensorToVoigt(const ViewCArrayKokkos<real_t>& E, 
                             const ViewCArrayKokkos<real_t>& E_voigt) {
        E_voigt(0) = E(0,0);          // E11
        E_voigt(1) = E(1,1);          // E22
        E_voigt(2) = E(2,2);          // E33
        E_voigt(3) = real_t(2.0) * E(1,2);    // 2*E23
        E_voigt(4) = real_t(2.0) * E(0,2);    // 2*E13
        E_voigt(5) = real_t(2.0) * E(0,1);    // 2*E12
    }


    /**
    * @brief Compute 2nd Piola-Kirchhoff stress (Voigt) from reference
    *        stiffness and Green-Lagrange strain (Voigt).
    *           S = C_ref * E
    */
    KOKKOS_INLINE_FUNCTION
    void computeSecondPKStressVoigt(const real_t C_ref[6][6],
                                    const ViewCArrayKokkos<real_t>& E_voigt,
                                    const ViewCArrayKokkos<real_t>& S_voigt) {

        for (int i = 0; i < 6; i++) {
            double sum = real_t(0.0);
            for (int j = 0; j < 6; j++) {
                sum += C_ref[i][j] * E_voigt(j);
            }
            S_voigt(i) = sum;
        }
    } // end function


    KOKKOS_INLINE_FUNCTION
    void stressVoigtToTensor(const ViewCArrayKokkos<real_t>& S_voigt, 
                             const ViewCArrayKokkos<real_t>& S) {
        S(0,0) = S_voigt(0);
        S(1,1) = S_voigt(1);
        S(2,2) = S_voigt(2);
        S(1,2) = S_voigt(3);  S(2,1) = S_voigt(3);
        S(0,2) = S_voigt(4);  S(2,0) = S_voigt(4);
        S(0,1) = S_voigt(5);  S(1,0) = S_voigt(5);
    }



    /**
    * @brief Build the 6x6 orthotropic stiffness matrix (Voigt notation)
    *        in the material's reference/principal frame from the 9
    *        standard engineering constants.
    *
    *        Normal-normal 3x3 block is built by inverting the
    *        corresponding compliance sub-block; shear terms are
    *        decoupled and diagonal.
    */
    KOKKOS_INLINE_FUNCTION
    void buildOrthotropicStiffness(const real_t E1, const real_t E2, const real_t E3,
                                   const real_t nu12, const real_t nu13, const real_t nu23,
                                   const real_t G12, const real_t G13, const real_t G23,
                                   real_t C_ref[6][6]) {

        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
                C_ref[i][j] = real_t(0.0);

        // Compliance sub-block (normal-normal), using symmetry nu_ij/E_i = nu_ji/E_j
        const real_t S11 =  real_t(1.0) / E1;
        const real_t S22 =  real_t(1.0) / E2;
        const real_t S33 =  real_t(1.0) / E3;
        const real_t S12 = -nu12 / E1;   // = -nu21/E2
        const real_t S13 = -nu13 / E1;   // = -nu31/E3
        const real_t S23 = -nu23 / E2;   // = -nu32/E3

        // Invert the symmetric 3x3 compliance block analytically
        const real_t detS = S11*(S22*S33 - S23*S23)
                          - S12*(S12*S33 - S23*S13)
                          + S13*(S12*S23 - S22*S13);

        const real_t inv_detS = real_t(1.0) / detS;

        C_ref[0][0] = (S22*S33 - S23*S23) * inv_detS;
        C_ref[1][1] = (S11*S33 - S13*S13) * inv_detS;
        C_ref[2][2] = (S11*S22 - S12*S12) * inv_detS;

        C_ref[0][1] = (S13*S23 - S12*S33) * inv_detS;
        C_ref[1][0] = C_ref[0][1];

        C_ref[0][2] = (S12*S23 - S13*S22) * inv_detS;
        C_ref[2][0] = C_ref[0][2];

        C_ref[1][2] = (S12*S13 - S11*S23) * inv_detS;
        C_ref[2][1] = C_ref[1][2];

        // Shear terms are decoupled and simply the inverse of the shear modulus
        C_ref[3][3] = G23;
        C_ref[4][4] = G13;
        C_ref[5][5] = G12;

    } // end function


    /**
    * @brief Compute K_max = max(C11, C22, C33) directly from the 9
    *        orthotropic engineering constants, without assembling the
    *        full 6x6 stiffness matrix. Only inverts the 3x3 normal block.
    */
    KOKKOS_INLINE_FUNCTION
    double maxBulkModulusFromEngineeringConstants(
        const real_t E1, const real_t E2, const real_t E3,
        const real_t nu12, const real_t nu13, const real_t nu23) {

        const real_t S11 = real_t(1.0)/E1;
        const real_t S22 = real_t(1.0)/E2;
        const real_t S33 = real_t(1.0)/E3;
        const real_t S12 = -nu12/E1;
        const real_t S13 = -nu13/E1;
        const real_t S23 = -nu23/E2;

        const real_t Delta = S11*S22*S33
                           - S11*S23*S23
                           - S22*S13*S13
                           - S33*S12*S12
                           + real_t(2.0)*S12*S13*S23;

        const real_t C11 = (S22*S33 - S23*S23) / (Delta + real_t(1.e-16));
        const real_t C22 = (S11*S33 - S13*S13) / (Delta + real_t(1.e-16));
        const real_t C33 = (S11*S22 - S12*S12) / (Delta + real_t(1.e-16));

        return Kokkos::fmax(C11, Kokkos::fmax(C22, C33));
    } // end function



    static void init_strength_state_vars(
        const DRaggedRightArrayKokkos<real_t> &MaterialPoints_eos_state_vars,
        const DRaggedRightArrayKokkos<real_t> &MaterialPoints_strength_state_vars,
        const RaggedRightArrayKokkos<real_t> &eos_global_vars,
        const RaggedRightArrayKokkos<real_t> &strength_global_vars,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const size_t num_material_points,
        const size_t mat_id)
    {
       

    }  // end of init_strength_state_vars



    KOKKOS_FUNCTION
    static void fill_C_matrix(const RaggedRightArrayKokkos<real_t> &strength_global_vars,
                              real_t C_ref[6][6],
                              const size_t mat_id)
    {
        const real_t E1   = strength_global_vars(mat_id, 0);
        const real_t E2   = strength_global_vars(mat_id, 1);
        const real_t E3   = strength_global_vars(mat_id, 2);
        const real_t nu12 = strength_global_vars(mat_id, 3);
        const real_t nu13 = strength_global_vars(mat_id, 4);
        const real_t nu23 = strength_global_vars(mat_id, 5);
        const real_t G12  = strength_global_vars(mat_id, 6);
        const real_t G13  = strength_global_vars(mat_id, 7);
        const real_t G23  = strength_global_vars(mat_id, 8);

        buildOrthotropicStiffness(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, C_ref);

    } // end of fill_C_matrix


    KOKKOS_FUNCTION
    static void calc_stress(
        const DCArrayKokkos<real_t>  &GaussPoints_vel_grad,
        const MPICArrayKokkos<real_t> &node_coords,
        const MPICArrayKokkos<real_t> &node_coords_t0,
        const MPICArrayKokkos<real_t> &node_vel,
        const DCArrayKokkos<size_t>  &nodes_in_elem,
        const DRaggedRightArrayKokkos<real_t>  &MaterialPoints_pres,
        const DRaggedRightArrayKokkos<real_t>  &MaterialPoints_stress,
        const DRaggedRightArrayKokkos<real_t>  &MaterialPoints_stress_n0,
        const DRaggedRightArrayKokkos<real_t>  &MaterialPoints_sspd,
        const DRaggedRightArrayKokkos<real_t> &MaterialPoints_eos_state_vars,
        const DRaggedRightArrayKokkos<real_t> &MaterialPoints_strength_state_vars,
        const real_t MaterialPoints_den,
        const real_t MaterialPoints_sie,
        const DRaggedRightArrayKokkos<real_t>& MaterialPoints_deformation_grad,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const RaggedRightArrayKokkos<real_t> &eos_global_vars,
        const RaggedRightArrayKokkos<real_t> &strength_global_vars,
        const real_t vol,
        const real_t dt,
        const real_t rk_alpha,
        const real_t time,
        const size_t cycle,
        const size_t MaterialPoints_sid,
        const size_t mat_id,
        const size_t gauss_gid,
        const size_t elem_gid)
    {

        // -----------------------------------------------------------
        // 0) Access deform_grad (input) and stress (output) for this
        //    material point, wrapped as 3x3 views
        // -----------------------------------------------------------
        ViewCArrayKokkos<real_t> F(&MaterialPoints_deformation_grad(mat_id,MaterialPoints_sid,0,0), 3, 3);
        ViewCArrayKokkos<real_t> sigma(&MaterialPoints_stress(mat_id,MaterialPoints_sid,0,0), 3, 3);

        // -----------------------------------------------------------
        // 1) Build the reference (material-frame) orthotropic
        //    stiffness matrix from the engineering constants.
        //    NOTE: if C_ref doesn't change in time, consider caching
        //    this in MaterialPoints_strength_state_vars instead of
        //    rebuilding it every call.
        // -----------------------------------------------------------
        real_t C_ref[6][6];
        fill_C_matrix(strength_global_vars, C_ref, mat_id);


        // -----------------------------------------------------------
        // 2) Green-Lagrange strain:  E = 1/2 (F^T F - I)
        // -----------------------------------------------------------
        real_t strain_1D[9];
        ViewCArrayKokkos <real_t> strain(&strain_1D[0], 3, 3);
        computeGreenLagrangeStrain(F, strain);
        
        // -----------------------------------------------------------
        // 3) Strain tensor -> Voigt -> 2nd PK stress (Voigt) -> tensor
        // -----------------------------------------------------------
        real_t E_voigt_1D[6];
        ViewCArrayKokkos<real_t> E_voigt(&E_voigt_1D[0], 6);
        strainTensorToVoigt(strain, E_voigt);

        real_t S_voigt_1D[6];
        ViewCArrayKokkos<real_t> S_voigt(&S_voigt_1D[0], 6);
        computeSecondPKStressVoigt(C_ref, E_voigt, S_voigt);

        real_t S_1D[9];
        ViewCArrayKokkos<real_t> S(&S_1D[0], 3, 3);
        stressVoigtToTensor(S_voigt, S);

        // -----------------------------------------------------------
        // 4) Push forward 2nd PK stress S to Cauchy stress sigma:
        //        sigma = (1/J) * F * S * F^T
        // -----------------------------------------------------------
        const real_t detF = Kokkos::fmax(det_3x3(F), real_t(1.e-13));

        real_t FS_1D[9];
        ViewCArrayKokkos<real_t> FS(&FS_1D[0], 3, 3);

        real_t Ft_1D[9];
        ViewCArrayKokkos<real_t> Ft(&Ft_1D[0], 3, 3);

        matmul3x3(F, S, FS);
        transpose3x3(F, Ft);
        matmul3x3(FS, Ft, sigma);   // writes directly into MaterialPoints_stress


        for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 3; j++){
            sigma(i,j) /= detF;
        }
        

        return;
    } // end of user mat


    KOKKOS_FUNCTION
    static void calc_sound_speed(const DRaggedRightArrayKokkos<real_t>& MaterialPoints_pres,
                                 const DRaggedRightArrayKokkos<real_t>& MaterialPoints_stress,
                                 const size_t mat_pt_lid,
                                 const size_t mat_id,
                                 const DRaggedRightArrayKokkos<real_t>& MaterialPoints_eos_state_vars,
                                 const DRaggedRightArrayKokkos<real_t>& MaterialPoints_sspd,
                                 const real_t den,
                                 const real_t sie,
                                 const RaggedRightArrayKokkos<real_t> &eos_global_vars,
                                 const RaggedRightArrayKokkos<real_t> &strength_global_vars)
    {

        const real_t E1   = strength_global_vars(mat_id, 0);
        const real_t E2   = strength_global_vars(mat_id, 1);
        const real_t E3   = strength_global_vars(mat_id, 2);
        const real_t nu12 = strength_global_vars(mat_id, 3);
        const real_t nu13 = strength_global_vars(mat_id, 4);
        const real_t nu23 = strength_global_vars(mat_id, 5);
        const real_t G12  = strength_global_vars(mat_id, 6);
        const real_t G13  = strength_global_vars(mat_id, 7);
        const real_t G23  = strength_global_vars(mat_id, 8);
        const real_t csmin  = strength_global_vars(mat_id, 9);

        // sspd = Kokkos::sqrt((B + 4/3G)/rho)
        const real_t B_max = maxBulkModulusFromEngineeringConstants(E1, E2, E3, nu12, nu13, nu23);

        const real_t G_max = Kokkos::fmax(G12, Kokkos::fmax(G13, G23));

        // sound speed
        MaterialPoints_sspd(mat_id, mat_pt_lid) = Kokkos::fmax(csmin, Kokkos::sqrt((B_max + (real_t(4.0)/real_t(3.0))*G_max)/den));

        return;
    } // end func

    

    
    static void destroy(
        const DRaggedRightArrayKokkos<real_t> &MaterialPoints_eos_state_vars,
        const DRaggedRightArrayKokkos<real_t> &MaterialPoints_strength_state_vars,
        const RaggedRightArrayKokkos<real_t> &eos_global_vars,
        const RaggedRightArrayKokkos<real_t> &strength_global_vars,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const size_t num_material_points,
        const size_t mat_ids)
    {

    } // end destory

} // end namespace



#endif // end Header Guard