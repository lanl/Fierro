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

#ifndef BOUNDARY_DISP_PISTON_H
#define BOUNDARY_DISP_PISTON_H

#include "boundary_conditions.hpp"

struct BoundaryConditionEnums_t;

namespace PistonDisplacementBC
{
    // add an enum for boundary statevars and global vars
    enum BCVars
    {
        x_comp = 0,
        y_comp = 1,
        z_comp = 2,
        hydro_bc_disp_0 = 3,
        hydro_bc_disp_1 = 4,
        hydro_bc_disp_2 = 5,
        hydro_bc_disp_t_start = 6,
        hydro_bc_disp_t_end = 7
    };

/////////////////////////////////////////////////////////////////////////////
///
/// \fn displacement
///
/// \brief This is a function to set the displacement in one direction to a
///        specified displacement.  The other components can freely slide on
///        the piston
///
/// \param Mesh object
/// \param Boundary condition enums to select options
/// \param Boundary condition global variables array
/// \param Boundary condition state variables array
/// \param Node displacement
/// \param Time of the simulation
/// \param Boundary global index for the surface node
/// \param Boundary set local id
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
static void displacement(const swage::Mesh& mesh,
    const DCArrayKokkos<BoundaryConditionEnums_t>& BoundaryConditionEnums,
    const RaggedRightArrayKokkos<double>& disp_bc_global_vars,
    const DCArrayKokkos<double>& bc_state_vars,
    const DCArrayKokkos<double>& node_disp,
    const double time_value,
    const size_t rk_stage,
    const size_t bdy_node_gid,
    const size_t bdy_set)
{


    // directions are:
    // x_plane  = 0,
    // y_plane  = 1,
    // z_plane  = 2,

    const double x_comp = disp_bc_global_vars(bdy_set, BCVars::x_comp);
    const double y_comp = disp_bc_global_vars(bdy_set, BCVars::y_comp);
    const double z_comp = disp_bc_global_vars(bdy_set, BCVars::z_comp);
    const double hydro_bc_disp_0 = disp_bc_global_vars(bdy_set, BCVars::hydro_bc_disp_0);
    const double hydro_bc_disp_1 = disp_bc_global_vars(bdy_set, BCVars::hydro_bc_disp_1);
    const double hydro_bc_disp_2 = disp_bc_global_vars(bdy_set, BCVars::hydro_bc_disp_2);
    const double hydro_bc_disp_t_start = disp_bc_global_vars(bdy_set, BCVars::hydro_bc_disp_t_start);
    const double hydro_bc_disp_t_end   = disp_bc_global_vars(bdy_set, BCVars::hydro_bc_disp_t_end);

    // Set displacement to that direction to specified value
    // if t_end > time > t_start
    // v(t) = v0 exp(-v1*(time - time_start) )
    if (time_value >= hydro_bc_disp_t_start && time_value <= hydro_bc_disp_t_end) {
        // the time difference
        const double time_delta = time_value - hydro_bc_disp_t_start;

        // the desired displacement
        const double disp = hydro_bc_disp_0 + hydro_bc_disp_1*time_delta + 0.5*hydro_bc_disp_2*time_delta*time_delta;

        // magnitude of the normal in the specified direction
        double mag = 0.0;
        for (size_t dim = 0; dim<mesh.num_dims; dim++){
            mag += disp_bc_global_vars(bdy_set,dim)*disp_bc_global_vars(bdy_set,dim);
        } // will make sure it's a unit vector
        mag = sqrt(mag);

        for (size_t dim = 0; dim<mesh.num_dims; dim++){
            // remove the displacement in the specified direction
            //  direction normal = disp_bc_global_vars(bdy_set,dim)/mag
            node_disp(bdy_node_gid, dim) -= node_disp(bdy_node_gid, dim) * (disp_bc_global_vars(bdy_set,dim)/mag);

            // add the desired displacement in the specified direction
            node_disp(bdy_node_gid, dim) += disp * (disp_bc_global_vars(bdy_set,dim)/mag);
        }
        
    } // end if on time

    return;
}     // end func
} // end namespace

#endif // end Header Guard