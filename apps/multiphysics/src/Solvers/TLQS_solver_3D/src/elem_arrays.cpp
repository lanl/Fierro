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

#include "tlqs_solver_3D.hpp"


KOKKOS_FUNCTION
void TLQS3D::get_gradients(
    const double material_matrix[6][6],
    ViewCArrayKokkos <size_t>& nodes_in_elem,
    const DCArrayKokkos <double>& coords_t0,
    const DCArrayKokkos <double>& displacement,
    const CArrayKokkos <double>& displacement_step,
    ViewCArrayKokkos <double>& gauss_point_grad_basis,
    double grad_u[3][3],
    double inv_J[3][3],
    double det_J,
    double PK2_curr_config[6])
{
    // allocate and initialize Jacobian
    double J[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            J[i][j] = 0.0;
        }
    }

    // get Jacobian at the mat point -> J with indices of [ [dx1dxi1 dx2dxi1 dx3dxi1] [dx1dxi2 dx2dxi2 dx3dxi2] [dx1dxi3 dx2dxi3 dx3dxi3] ]
    for (int k = 0; k < gauss_point_grad_basis.dims(0); k++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                J[i][j] += coords_t0(nodes_in_elem(k), j) * gauss_point_grad_basis(k, i);
            }
        }
    }

    // get det(J)
    det_J = J[0][0]*(J[1][1]*J[2][2] - J[2][1]*J[1][2]) - J[0][1]*(J[1][0]*J[2][2] - J[2][0]*J[1][2]) + J[0][2]*(J[1][0]*J[2][1] - J[2][0]*J[1][1]);

    // get inv(J)
    double adjoint[3][3];
    adjoint[0][0] = J[1][1]*J[2][2] - J[2][1]*J[1][2];
    adjoint[0][1] = -(J[0][1]*J[2][2] - J[2][1]*J[0][2]);
    adjoint[0][2] = J[0][1]*J[1][2] - J[1][1]*J[0][2];
    adjoint[1][0] = -(J[1][0]*J[2][2] - J[2][0]*J[1][2]);
    adjoint[1][1] = J[0][0]*J[2][2] - J[2][0]*J[0][2];
    adjoint[1][2] = -(J[0][0]*J[1][2] - J[1][0]*J[0][2]);
    adjoint[2][0] = J[1][0]*J[2][1] - J[2][0]*J[1][1];
    adjoint[2][1] = -(J[0][0]*J[2][1] - J[2][0]*J[0][1]);
    adjoint[2][2] = J[0][0]*J[1][1] - J[1][0]*J[0][1];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            inv_J[i][j] = adjoint[i][j] / det_J;
        }
    }

    // get grad(displacement)
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            grad_u[i][j] = 0.0;
        }
    }

    for (int k = 0; k < gauss_point_grad_basis.dims(0); k++) {
        const size_t node_gid = nodes_in_elem(k);

        const double dpsig_k0 = inv_J[0][0]*gauss_point_grad_basis(k,0) + inv_J[0][1]*gauss_point_grad_basis(k,1) + inv_J[0][2]*gauss_point_grad_basis(k,2);
        const double dpsig_k1 = inv_J[1][0]*gauss_point_grad_basis(k,0) + inv_J[1][1]*gauss_point_grad_basis(k,1) + inv_J[1][2]*gauss_point_grad_basis(k,2);
        const double dpsig_k2 = inv_J[2][0]*gauss_point_grad_basis(k,0) + inv_J[2][1]*gauss_point_grad_basis(k,1) + inv_J[2][2]*gauss_point_grad_basis(k,2);

        for (int j = 0; j < 3; j++) {
            const double u_total = displacement(node_gid, j) + displacement_step(3*node_gid + j);
            grad_u[0][j] += u_total * dpsig_k0;
            grad_u[1][j] += u_total * dpsig_k1;
            grad_u[2][j] += u_total * dpsig_k2;
        }
    }

    // get second PK stress of current configuration
    // ***************************************************
    // WARNING: CURRENTLY ASSUMES ISOTROPIC LINEAR ELASTIC
    // WARNING: NEED TO PUT THIS INTO A SEPARATE FUNCTION
    // WARNING: POINTER CALLED CALC QUASI_STATIC STRESS
    // ***************************************************
    double current_strain[6]; // [Exx Eyy Ezz Eyz Exz Exy]
    current_strain[0] = grad_u[0][0] + 0.5 * (grad_u[0][0]*grad_u[0][0] + grad_u[1][0]*grad_u[1][0] + grad_u[2][0]*grad_u[2][0]);
    current_strain[1] = grad_u[1][1] + 0.5 * (grad_u[0][1]*grad_u[0][1] + grad_u[1][1]*grad_u[1][1] + grad_u[2][1]*grad_u[2][1]);
    current_strain[2] = grad_u[2][2] + 0.5 * (grad_u[0][2]*grad_u[0][2] + grad_u[1][2]*grad_u[1][2] + grad_u[2][2]*grad_u[2][2]);
    current_strain[3] = 0.5 * (grad_u[1][2] + grad_u[2][1] + (grad_u[0][1]*grad_u[0][2] + grad_u[1][1]*grad_u[1][2] + grad_u[2][1]*grad_u[2][2]));
    current_strain[4] = 0.5 * (grad_u[0][2] + grad_u[2][0] + (grad_u[0][0]*grad_u[0][2] + grad_u[1][0]*grad_u[1][2] + grad_u[2][0]*grad_u[2][2]));
    current_strain[5] = 0.5 * (grad_u[0][1] + grad_u[1][0] + (grad_u[0][0]*grad_u[0][1] + grad_u[1][0]*grad_u[1][1] + grad_u[2][0]*grad_u[2][1]));

    // Normal stresses
    PK2_curr_config[0] = material_matrix[0][0] * current_strain[0] + material_matrix[0][1] * current_strain[1] + material_matrix[0][2] * current_strain[2]; // Sxx
    PK2_curr_config[1] = material_matrix[1][0] * current_strain[0] + material_matrix[1][1] * current_strain[1] + material_matrix[1][2] * current_strain[2]; // Syy
    PK2_curr_config[2] = material_matrix[2][0] * current_strain[0] + material_matrix[2][1] * current_strain[1] + material_matrix[2][2] * current_strain[2]; // Szz

    // Shear stresses 
    PK2_curr_config[3] = material_matrix[3][3] * current_strain[3]; // Syz
    PK2_curr_config[4] = material_matrix[4][4] * current_strain[4]; // Sxz
    PK2_curr_config[5] = material_matrix[5][5] * current_strain[5]; // Sxy

} // end get_gradients

KOKKOS_FUNCTION
void TLQS3D::tally_elem_arrays(
    const double material_matrix[6][6],
    const double grad_u[3][3],
    const double inv_J[3][3],
    ViewCArrayKokkos <double>& gauss_point_grad_basis,
    double gauss_point_weight,
    const double PK2_curr_config[6],
    ViewCArrayKokkos <double>& Kel,
    ViewCArrayKokkos <double>& Fel)
{
    const int num_nodes = gauss_point_grad_basis.dims(0);

    // Unpack grad_u for readability
    const double ux = grad_u[0][0];
    const double uy = grad_u[1][0];
    const double uz = grad_u[2][0];
    const double vx = grad_u[0][1];
    const double vy = grad_u[1][1];
    const double vz = grad_u[2][1];
    const double wx = grad_u[0][2];
    const double wy = grad_u[1][2];
    const double wz = grad_u[2][2];

    // Build symmetric 3x3 stress tensor for K2 (extracted from Voigt PK2)
    double S[3][3];
    S[0][0] = PK2_curr_config[0];
    S[0][1] = PK2_curr_config[5];
    S[0][2] = PK2_curr_config[4];
    S[1][0] = PK2_curr_config[5];
    S[1][1] = PK2_curr_config[1];
    S[1][2] = PK2_curr_config[3];
    S[2][0] = PK2_curr_config[4];
    S[2][1] = PK2_curr_config[3];
    S[2][2] = PK2_curr_config[2];

    // temp arrays for forming element matrix
    double B1_a[6][3];              // B1 values for node a
    double CT_matmul_B1_a[6][3];      // C^T * B1_a = C * B1_a (C symmetric), using transpose for better memory access pattern
    double S_mul_glob_grad_a[3];   // only needs to be size 3 because of the sparsity and repeating nature of [S]9x9 and B2
    double B1_b[6][3];              // B1 values for node b

    // looping through each node to avoid dynamic allocations
    // outer loop of node a
    for (int a = 0; a < num_nodes; a++) {

        // dpsig for node a: inv_J * grad_basis(a)
        const double glob_grad_a_0 = inv_J[0][0]*gauss_point_grad_basis(a,0) + inv_J[0][1]*gauss_point_grad_basis(a,1) + inv_J[0][2]*gauss_point_grad_basis(a,2);
        const double glob_grad_a_1 = inv_J[1][0]*gauss_point_grad_basis(a,0) + inv_J[1][1]*gauss_point_grad_basis(a,1) + inv_J[1][2]*gauss_point_grad_basis(a,2);
        const double glob_grad_a_2 = inv_J[2][0]*gauss_point_grad_basis(a,0) + inv_J[2][1]*gauss_point_grad_basis(a,1) + inv_J[2][2]*gauss_point_grad_basis(a,2);

        // columns of B1 for node a
        B1_a[0][0] = glob_grad_a_0*(1+ux);
        B1_a[0][1] = glob_grad_a_0*vx;
        B1_a[0][2] = glob_grad_a_0*wx;

        B1_a[1][0] = glob_grad_a_1*uy;
        B1_a[1][1] = glob_grad_a_1*(1+vy);
        B1_a[1][2] = glob_grad_a_1*wy;

        B1_a[2][0] = glob_grad_a_2*uz;
        B1_a[2][1] = glob_grad_a_2*vz;
        B1_a[2][2] = glob_grad_a_2*(1+wz);

        B1_a[3][0] = glob_grad_a_1*uz + glob_grad_a_2*uy;
        B1_a[3][1] = glob_grad_a_1*vz + glob_grad_a_2*(1+vy);
        B1_a[3][2] = glob_grad_a_1*(1+wz) + glob_grad_a_2*wy;

        B1_a[4][0] = glob_grad_a_2*(1+ux) + glob_grad_a_0*uz;
        B1_a[4][1] = glob_grad_a_0*vz + glob_grad_a_2*vx;
        B1_a[4][2] = glob_grad_a_0*(1+wz) + glob_grad_a_2*wx;

        B1_a[5][0] = glob_grad_a_1*(1+ux) + glob_grad_a_0*uy;
        B1_a[5][1] = glob_grad_a_0*(1+vy) + glob_grad_a_1*vx;
        B1_a[5][2] = glob_grad_a_0*wy + glob_grad_a_1*wx;

        // Precompute C^T * B1_a
        // Used in inner loop as: K1(3a+p, 3b+q) = sum_m CtB1_a[m][p] * B1_b[m][q]
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 3; j++) {
                CT_matmul_B1_a[i][j] = 0.0;
            }
        }
        for (int m = 0; m < 6; m++) {
            for (int k = 0; k < 6; k++) {
                for (int p = 0; p < 3; p++) {
                    CT_matmul_B1_a[m][p] += material_matrix[m][k] * B1_a[k][p];
                }
            }
        }

        // K2 simplification: The geometric stiffness block between nodes a and b is a 
        // diagonal 3x3 matrix where the non-zero terms equal grad(a) · S · grad(b).
        // Precomputing S · grad(a) to save operations in the inner loop.
        // Since S is symmetric, (Sg)⋅h=g⋅S⋅h
        S_mul_glob_grad_a[0] = S[0][0]*glob_grad_a_0 + S[0][1]*glob_grad_a_1 + S[0][2]*glob_grad_a_2;
        S_mul_glob_grad_a[1] = S[1][0]*glob_grad_a_0 + S[1][1]*glob_grad_a_1 + S[1][2]*glob_grad_a_2;
        S_mul_glob_grad_a[2] = S[2][0]*glob_grad_a_0 + S[2][1]*glob_grad_a_1 + S[2][2]*glob_grad_a_2;

        // Fel: Fel(3a+p) -= weight * sum_k B1_a[k][p] * PK2[k]
        // minus because Fel = boundarcy forces - current force in the element aka Fel = F02 - F01
        for (int p = 0; p < 3; p++) {
            double fel_val = 0.0;
            for (int k = 0; k < 6; k++) {
                fel_val += B1_a[k][p] * PK2_curr_config[k];
            }
            Fel(3*a + p) -= gauss_point_weight * fel_val;
        }

        // Inner loop over b for Kel
        for (int b = 0; b < num_nodes; b++) {

            // dpsig for node b
            const double glob_grad_b_0 = inv_J[0][0]*gauss_point_grad_basis(b,0) + inv_J[0][1]*gauss_point_grad_basis(b,1) + inv_J[0][2]*gauss_point_grad_basis(b,2);
            const double glob_grad_b_1 = inv_J[1][0]*gauss_point_grad_basis(b,0) + inv_J[1][1]*gauss_point_grad_basis(b,1) + inv_J[1][2]*gauss_point_grad_basis(b,2);
            const double glob_grad_b_2 = inv_J[2][0]*gauss_point_grad_basis(b,0) + inv_J[2][1]*gauss_point_grad_basis(b,1) + inv_J[2][2]*gauss_point_grad_basis(b,2);

            // B1_b[6][3]: columns of B1 for node b
            double B1_b[6][3];
            B1_b[0][0] = glob_grad_b_0*(1+ux);
            B1_b[0][1] = glob_grad_b_0*vx;
            B1_b[0][2] = glob_grad_b_0*wx;

            B1_b[1][0] = glob_grad_b_1*uy;
            B1_b[1][1] = glob_grad_b_1*(1+vy);
            B1_b[1][2] = glob_grad_b_1*wy;

            B1_b[2][0] = glob_grad_b_2*uz;
            B1_b[2][1] = glob_grad_b_2*vz;
            B1_b[2][2] = glob_grad_b_2*(1+wz);

            B1_b[3][0] = glob_grad_b_1*uz + glob_grad_b_2*uy;
            B1_b[3][1] = glob_grad_b_1*vz + glob_grad_b_2*(1+vy);
            B1_b[3][2] = glob_grad_b_1*(1+wz) + glob_grad_b_2*wy;

            B1_b[4][0] = glob_grad_b_2*(1+ux) + glob_grad_b_0*uz;
            B1_b[4][1] = glob_grad_b_0*vz + glob_grad_b_2*vx;
            B1_b[4][2] = glob_grad_b_0*(1+wz) + glob_grad_b_2*wx;

            B1_b[5][0] = glob_grad_b_1*(1+ux) + glob_grad_b_0*uy;
            B1_b[5][1] = glob_grad_b_0*(1+vy) + glob_grad_b_1*vx;
            B1_b[5][2] = glob_grad_b_0*wy + glob_grad_b_1*wx;

            // K2 scalar: dpsig_a · S3 · dpsig_b (same value along p==q diagonal)
            const double k2_scalar = S_mul_glob_grad_a[0]*glob_grad_b_0 + S_mul_glob_grad_a[1]*glob_grad_b_1 + S_mul_glob_grad_a[2]*glob_grad_b_2;

            // Accumulate 3x3 block into Kel
            for (int p = 0; p < 3; p++) {
                for (int q = 0; q < 3; q++) {
                    double k1_val = 0.0;
                    for (int m = 0; m < 6; m++) {
                        k1_val += CT_matmul_B1_a[m][p] * B1_b[m][q];
                    }
                    const double k2_val = (p == q) ? k2_scalar : 0.0;
                    Kel(3*a+p, 3*b+q) += gauss_point_weight * (k1_val + k2_val);
                }
            }
        }

    }

} // end tally_elem_arrays