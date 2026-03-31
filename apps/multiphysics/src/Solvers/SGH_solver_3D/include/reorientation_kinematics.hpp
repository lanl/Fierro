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

#ifndef REORIENTATION_KINEMATICS_H
#define REORIENTATION_KINEMATICS_H

#include "matar.h"
#include "ELEMENTS.h"
#include "mesh_io.hpp"
#include "state.hpp"
#include <cmath>

namespace ReorientationKinematics
{

/////////////////////////////////////////////////////////////////////////////
///
/// \fn rotation_matrix_y
///
/// \brief Rotation matrix about y-axis (x2)
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void rotation_matrix_y(double angle, double R[3][3]) {
    double c = cos(angle);
    double s = sin(angle);
    R[0][0] = c;   R[0][1] = 0.0; R[0][2] = s;
    R[1][0] = 0.0; R[1][1] = 1.0; R[1][2] = 0.0;
    R[2][0] = -s;  R[2][1] = 0.0; R[2][2] = c;
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn rotation_matrix_z
///
/// \brief Rotation matrix about z-axis (x3)
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void rotation_matrix_z(double angle, double R[3][3]) {
    double c = cos(angle);
    double s = sin(angle);
    R[0][0] = c;   R[0][1] = -s;  R[0][2] = 0.0;
    R[1][0] = s;   R[1][1] = c;   R[1][2] = 0.0;
    R[2][0] = 0.0; R[2][1] = 0.0; R[2][2] = 1.0;
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn mat_mult_3x3
///
/// \brief Multiply two 3x3 matrices: C = A * B; used for computing rotations
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
void mat_mult_3x3(const double A[3][3], const double B[3][3], double C[3][3]) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            C[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn prescribe_reorientation_kinematics
///
/// \brief Prescribes rigid-body rotation + cohesive zone opening
///        Overrides solver's position/velocity update for validation testing
///
/////////////////////////////////////////////////////////////////////////////
inline void prescribe_reorientation_kinematics(
    const swage::Mesh& mesh,
    State_t& State,
    const CArrayKokkos<double>& initial_coords,
    const CArrayKokkos<int>&    cz_b_side_flag,
    double time_value,
    double dt_stage,
    double omega_y,
    double omega_z,
    double cz_opening_rate)
{
    const double angle_y_t   = omega_y * time_value;
    const double angle_z_t   = omega_z * time_value;
    const double angle_y_tdt = omega_y * (time_value + dt_stage);
    const double angle_z_tdt = omega_z * (time_value + dt_stage);

    const double s_t   = cz_opening_rate * time_value;
    const double s_tdt = cz_opening_rate * (time_value + dt_stage);

    double Ry_t[3][3], Rz_t[3][3], R_t[3][3];
    double Ry_tdt[3][3], Rz_tdt[3][3], R_tdt[3][3];

    rotation_matrix_y(angle_y_t,   Ry_t);
    rotation_matrix_z(angle_z_t,   Rz_t);
    mat_mult_3x3(Rz_t, Ry_t, R_t);

    rotation_matrix_y(angle_y_tdt, Ry_tdt);
    rotation_matrix_z(angle_z_tdt, Rz_tdt);
    mat_mult_3x3(Rz_tdt, Ry_tdt, R_tdt);

    // n(t) = R(t)*[1,0,0] => column 0
    const double n_t[3]   = { R_t[0][0],   R_t[1][0],   R_t[2][0]   };
    const double n_tdt[3] = { R_tdt[0][0], R_tdt[1][0], R_tdt[2][0] };

    FOR_ALL(node_gid, 0, mesh.num_nodes, {
        const double Xx = initial_coords(node_gid,0);
        const double Xy = initial_coords(node_gid,1);
        const double Xz = initial_coords(node_gid,2);

        // x(t) = R(t)*X0
        double xt   = R_t[0][0]*Xx + R_t[0][1]*Xy + R_t[0][2]*Xz;
        double yt   = R_t[1][0]*Xx + R_t[1][1]*Xy + R_t[1][2]*Xz;
        double zt   = R_t[2][0]*Xx + R_t[2][1]*Xy + R_t[2][2]*Xz;

        // x(t+dt_stage) = R(t+dt_stage)*X0
        double xtdt = R_tdt[0][0]*Xx + R_tdt[0][1]*Xy + R_tdt[0][2]*Xz;
        double ytdt = R_tdt[1][0]*Xx + R_tdt[1][1]*Xy + R_tdt[1][2]*Xz;
        double ztdt = R_tdt[2][0]*Xx + R_tdt[2][1]*Xy + R_tdt[2][2]*Xz;

        // Add opening displacement for B-side nodes
        if (cz_b_side_flag(node_gid)) {
            xt   += s_t   * n_t[0];   yt   += s_t   * n_t[1];   zt   += s_t   * n_t[2];
            xtdt += s_tdt * n_tdt[0]; ytdt += s_tdt * n_tdt[1]; ztdt += s_tdt * n_tdt[2];
        }

        State.node.coords(node_gid,0) = xt;
        State.node.coords(node_gid,1) = yt;
        State.node.coords(node_gid,2) = zt;

        State.node.vel(node_gid,0) = (xtdt - xt) / dt_stage;
        State.node.vel(node_gid,1) = (ytdt - yt) / dt_stage;
        State.node.vel(node_gid,2) = (ztdt - zt) / dt_stage;
    });

    Kokkos::fence();
}

} // end namespace ReorientationKinematics

#endif // REORIENTATION_KINEMATICS_H