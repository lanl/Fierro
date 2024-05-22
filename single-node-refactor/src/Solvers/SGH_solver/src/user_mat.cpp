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

#include "sgh_solver.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn user_eos_model
///
/// \brief <insert brief description>
///
/// REMOVE TO MATERIAL MODEL LIBRARY
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void user_eos_model(const DCArrayKokkos<double>& elem_pres,
    const DCArrayKokkos<double>& elem_stress,
    const size_t elem_gid,
    const size_t mat_id,
    const DCArrayKokkos<double>& elem_state_vars,
    const DCArrayKokkos<double>& elem_sspd,
    const double den,
    const double sie)
{
    const int num_dims = 3;

    // -----------------------------------------------------------------------------
    // Required variables are here
    // ------------------------------------------------------------------------------
    elem_pres(elem_gid) = 0.0;  // pressure
    elem_sspd(elem_gid) = 1.0e-8;  // sound speed

    // pressure = 1/3tr(stress)
    for (int i = 0; i < num_dims; i++) {
        elem_pres(elem_gid) -= elem_stress(i, i);
    }
    elem_pres(elem_gid) *= 1.0 / 3.0;

    return;
} // end for user_eos_model

/////////////////////////////////////////////////////////////////////////////
///
/// \fn user_strength_model
///
/// REMOVE TO MATERIAL MODEL LIBRARY
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void user_strength_model(const DCArrayKokkos<double>& elem_pres,
    const DCArrayKokkos<double>& elem_stress,
    const size_t elem_gid,
    const size_t mat_id,
    const DCArrayKokkos<double>& elem_state_vars,
    const DCArrayKokkos<double>& elem_sspd,
    const double den,
    const double sie,
    const ViewCArrayKokkos<double>& vel_grad,
    const ViewCArrayKokkos<size_t>& elem_node_gids,
    const DCArrayKokkos<double>&    node_coords,
    const DCArrayKokkos<double>&    node_vel,
    const double vol,
    const double dt,
    const double rk_alpha)
{
    // statev(0) = var_1
    //   :
    //   :
    //   :
    // statev(N) = var_N

    // int num_dims = 3;

    // -----------------------------------------------------------------------------
    // The user must coding goes here
    // ------------------------------------------------------------------------------

    return;
} // end of user mat

/////////////////////////////////////////////////////////////////////////////
///
/// \fn user_strength_model_vpsc
///
/// REMOVE TO MATERIAL MODEL LIBRARY
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void user_strength_model_vpsc(const DCArrayKokkos<double>& elem_pres,
    const DCArrayKokkos<double>& elem_stress,
    const size_t elem_gid,
    const size_t mat_id,
    const DCArrayKokkos<double>& elem_state_vars,
    const DCArrayKokkos<double>& elem_sspd,
    const double den,
    const double sie,
    const ViewCArrayKokkos<double>& vel_grad,
    const ViewCArrayKokkos<size_t>& elem_node_gids,
    const DCArrayKokkos<double>&    node_coords,
    const DCArrayKokkos<double>&    node_vel,
    const double vol,
    const double dt,
    const double rk_alpha)
{
    // statev(0) = var_1
    //   :
    //   :
    //   :
    // statev(N) = var_N

    const int num_dims = 3;

    // -----------------------------------------------------------------------------
    // Required variables are here
    // ------------------------------------------------------------------------------
    elem_pres(elem_gid) = 1.0e-15;  // pressure
    elem_sspd(elem_gid) = 1.0e-15;  // sound speed

    // -----------------------------------------------------------------------------
    // The user must coding goes here
    // ------------------------------------------------------------------------------

    // For hypo-elastic models
    double D_tensor_values[9];
    double W_tensor_values[9];

    // convert to array syntax with the C-Language access pattern
    ViewCArrayKokkos<double> D_tensor(D_tensor_values, num_dims, num_dims);  // D(i,j)
    ViewCArrayKokkos<double> W_tensor(W_tensor_values, num_dims, num_dims);  // W(i,j)

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

    ViewCMatrixKokkos<double> deps(&deps_values[0], num_dims, num_dims);
    ViewCMatrixKokkos<double> dW(&dW_values[0], num_dims, num_dims);
    ViewCMatrixKokkos<double> drot(&drot_values[0], num_dims, num_dims);

    // calculate strain and rotation increments
    for (size_t i = 1; i <= 3; i++) {
        for (size_t j = 1; j <= 3; j++) {
            deps(i, j) = D_tensor(i - 1, j - 1) * dt; // infinitesimal strain increment
            dW(i, j)   = W_tensor(i - 1, j - 1) * dt;
        }
    } // end for

    return;
} // end of user_strength_model_vpsc
