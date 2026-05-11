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
void TLQS3D::post_process(
    const double material_matrix[6][6],
    ViewCArrayKokkos <size_t>& nodes_in_elem,
    const DCArrayKokkos <double>& coords_t0,
    const DCArrayKokkos <double>& displacement,
    ViewCArrayKokkos <double>& gauss_point_grad_basis,
    ViewCArrayKokkos <double>& stress,
    ViewCArrayKokkos <double>& strain)
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
    double det_J = J[0][0]*(J[1][1]*J[2][2] - J[2][1]*J[1][2]) - J[0][1]*(J[1][0]*J[2][2] - J[2][0]*J[1][2]) + J[0][2]*(J[1][0]*J[2][1] - J[2][0]*J[1][1]);
    //std::cout << "DET_J: " << det_J << std::endl;
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
    
    double inv_J[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            inv_J[i][j] = adjoint[i][j] / det_J;
        }
    }

    double grad_u[3][3];
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
            const double u_total = displacement(node_gid, j);
            grad_u[0][j] += u_total * dpsig_k0;
            grad_u[1][j] += u_total * dpsig_k1;
            grad_u[2][j] += u_total * dpsig_k2;
        }
    }

    // get second PK stress and Green-Lagrange strain of current configuration
    // ***************************************************
    // WARNING: CURRENTLY ASSUMES ISOTROPIC LINEAR ELASTIC
    // WARNING: NEED TO PUT THIS INTO A SEPARATE FUNCTION
    // WARNING: POINTER CALLED CALC QUASI_STATIC STRESS
    // ***************************************************
    double current_strain[6]; // [Exx Eyy Ezz Eyz Exz Exy]
    strain(0,0) = grad_u[0][0] + 0.5 * (grad_u[0][0]*grad_u[0][0] + grad_u[1][0]*grad_u[1][0] + grad_u[2][0]*grad_u[2][0]);
    strain(1,1) = grad_u[1][1] + 0.5 * (grad_u[0][1]*grad_u[0][1] + grad_u[1][1]*grad_u[1][1] + grad_u[2][1]*grad_u[2][1]);
    strain(2,2) = grad_u[2][2] + 0.5 * (grad_u[0][2]*grad_u[0][2] + grad_u[1][2]*grad_u[1][2] + grad_u[2][2]*grad_u[2][2]);
    strain(1,2) = 0.5 * (grad_u[1][2] + grad_u[2][1] + (grad_u[0][1]*grad_u[0][2] + grad_u[1][1]*grad_u[1][2] + grad_u[2][1]*grad_u[2][2]));
    strain(2,1) = 0.5 * (grad_u[1][2] + grad_u[2][1] + (grad_u[0][1]*grad_u[0][2] + grad_u[1][1]*grad_u[1][2] + grad_u[2][1]*grad_u[2][2]));
    strain(0,2) = 0.5 * (grad_u[0][2] + grad_u[2][0] + (grad_u[0][0]*grad_u[0][2] + grad_u[1][0]*grad_u[1][2] + grad_u[2][0]*grad_u[2][2]));
    strain(2,0) = 0.5 * (grad_u[0][2] + grad_u[2][0] + (grad_u[0][0]*grad_u[0][2] + grad_u[1][0]*grad_u[1][2] + grad_u[2][0]*grad_u[2][2]));
    strain(0,1) = 0.5 * (grad_u[0][1] + grad_u[1][0] + (grad_u[0][0]*grad_u[0][1] + grad_u[1][0]*grad_u[1][1] + grad_u[2][0]*grad_u[2][1]));
    strain(1,0) = 0.5 * (grad_u[0][1] + grad_u[1][0] + (grad_u[0][0]*grad_u[0][1] + grad_u[1][0]*grad_u[1][1] + grad_u[2][0]*grad_u[2][1]));

    // Normal stresses
    stress(0,0) = material_matrix[0][0] * strain(0,0) + material_matrix[0][1] * strain(1,1) + material_matrix[0][2] * strain(2,2); // Sxx
    stress(1,1) = material_matrix[1][0] * strain(0,0) + material_matrix[1][1] * strain(1,1) + material_matrix[1][2] * strain(2,2); // Syy
    stress(2,2) = material_matrix[2][0] * strain(0,0) + material_matrix[2][1] * strain(1,1) + material_matrix[2][2] * strain(2,2); // Szz

    // Shear stresses 
    stress(1,2) = material_matrix[3][3] * strain(1,2); // Syz
    stress(2,1) = material_matrix[3][3] * strain(1,2); // Syz

    stress(0,2) = material_matrix[4][4] * strain(0,2); // Sxz
    stress(2,0) = material_matrix[4][4] * strain(0,2); // Sxz

    stress(0,1) = material_matrix[5][5] * strain(0,1); // Sxy
    stress(1,0) = material_matrix[5][5] * strain(0,1); // Sxy

} // end get_gradients