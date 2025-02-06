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

#ifndef BOUNDARY_STRESS_TIME_VARYING_H
#define BOUNDARY_STRESS_TIME_VARYING_H

#include "boundary_conditions.h"

struct BoundaryConditionEnums_t;

namespace TimeVaryingStressBC
{
// add an enum for boundary statevars and global vars
enum BCVars
{
    sig_00 = 0,
    sig_11 = 1,
    sig_22 = 2,
    sig_12 = 3,
    sig_02 = 4,
    sig_01 = 5,
    var_decay = 6,
    time_start = 7,
    time_end = 8
};

/////////////////////////////////////////////////////////////////////////////
///
/// \fn stress
///
/// \brief This is a function to force the boundary to move with a specified
///        pressure that varies in time.  The function form is an exponential
///        decay, given by:
///         if (t_end > time > t_start) then
///              sigma(t) = sigma0 exp(-varDecay*(time - time_start) )
///
/// \param Mesh object
/// \param Boundary condition enums to select options
/// \param Boundary condition global variables array
/// \param Boundary condition state variables array
/// \param corner stress
/// \param Time of the simulation
/// \param Boundary global index for the surface node
/// \param Boundary set local id
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
static void stress(const Mesh_t& mesh,
    const DCArrayKokkos<BoundaryConditionEnums_t>& BoundaryConditionEnums,
    const RaggedRightArrayKokkos<double>& stress_bc_global_vars,
    const DCArrayKokkos<double>& bc_state_vars,
    const ViewCArrayKokkos <double>& corner_surf_force,
    const ViewCArrayKokkos <double>& corner_surf_normal,
    const double time_value,
    const size_t rk_stage,
    const size_t bdy_node_gid,
    const size_t bdy_set)
{

    double sigma_1D[9];
    ViewCArrayKokkos<double> sigma(&sigma_1D[0], 3,3);  // its 3D even in 2D
   
    // Cauchy stress is symmetric
    sigma(0,0) = stress_bc_global_vars(bdy_set,BCVars::sig_00);
    sigma(1,1) = stress_bc_global_vars(bdy_set,BCVars::sig_11);
    sigma(2,2) = stress_bc_global_vars(bdy_set,BCVars::sig_22);

    sigma(1,2) = stress_bc_global_vars(bdy_set,BCVars::sig_12);
    sigma(2,1) = stress_bc_global_vars(bdy_set,BCVars::sig_12);

    sigma(0,2) = stress_bc_global_vars(bdy_set,BCVars::sig_02);
    sigma(2,0) = stress_bc_global_vars(bdy_set,BCVars::sig_02);

    sigma(0,1) = stress_bc_global_vars(bdy_set,BCVars::sig_01);
    sigma(1,0) = stress_bc_global_vars(bdy_set,BCVars::sig_01);

    const double var_decay = stress_bc_global_vars(bdy_set, BCVars::var_decay);
    const double time_start = stress_bc_global_vars(bdy_set, BCVars::time_start);
    const double time_end   = stress_bc_global_vars(bdy_set, BCVars::time_end);

    // directions are as follows:
    // x_plane  = 0,
    // y_plane  = 1,
    // z_plane  = 2,

    // Set pressure in that direction to specified value
    // if t_end > time > t_start
    // simga(t) = simga0 exp(-varDecay*(time - time_start) )
    if (time_value >= time_start && time_value <= time_end) {
        // the time difference
        const double time_delta = time_value - time_start;

        // the desired mag
        const double mag = exp(-var_decay * time_delta);

        
        // sigma(0,0)*normal(0) + sigma(0,1)*normal(1) + sigma(0,2)*normal(2)
        // sigma(1,0)*normal(0) + sigma(1,1)*normal(1) + sigma(1,2)*normal(2)
        // sigma(2,0)*normal(0) + sigma(2,1)*normal(1) + sigma(2,2)*normal(2)
        for(size_t i=0; i<mesh.num_dims; i++){
            for(size_t j=0; j<mesh.num_dims; j++){
                corner_surf_force(i) += mag*sigma(i,j)*corner_surf_normal(j);
            } // end j
        } // end i

    } // end if within the time frame

    return;
} // end func
} // end namespace

#endif // end Header Guard