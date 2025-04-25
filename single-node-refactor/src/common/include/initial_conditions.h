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
#ifndef FIERRO_IC_H
#define FIERRO_IC_H

#include <map>

namespace init_conds
{
    // initial conditions for a vector field
    enum init_vector_conds
    {
        noICsVec = 0,  // do nothing

        stationary = 1,   // set components in the vector to zero

        // uniform
        cartesian = 2,       // cart vector
        radialVec = 3,          // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
        sphericalVec = 4,       // spherical

        // linear variation
        radialLinearVec = 5,         // linear variation from 0,0,0
        sphericalLinearVec = 6,      // linear variation from 0,0,0

        // vortical initial conditions
        tgVortexVec = 7

        // user defined here
    };

    // initial conditions for a scalar field
    enum init_scalar_conds
    {
        noICsScalar = 0,  // do nothing

        // uniform
        uniform = 1,  // same value everywhere
        radialScalar = 2,     // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
        sphericalScalar = 3,  // spherical

        // linear variations
        xlinearScalar = 4, // linear variations from 0,0,0
        ylinearScalar = 5,
        zlinearScalar = 6,
        // xy_linear = 5,
        // xz_linear = 6,
        // yz_linear = 7,
        // xyz_linear = 8,
        // rad_linear = 9,      // linear variation from 0,0,0
        // sph_linear = 10,     // linear variation from 0,0,0

        // vortical initial conditions
        tgVortexScalar = 11
    };

} // end of initial conditions namespace




static std::map<std::string, init_conds::init_vector_conds> vector_ics_type_map
{
    { "static", init_conds::stationary },
    { "cartesian", init_conds::cartesian },
    { "radial", init_conds::radialVec },
    { "spherical", init_conds::sphericalVec },
    { "radial_linear", init_conds::radialLinearVec },
    { "spherical_linear", init_conds::sphericalLinearVec },
    { "tg_vortex", init_conds::tgVortexVec }
};

static std::map<std::string, init_conds::init_scalar_conds> scalar_ics_type_map
{
    { "uniform", init_conds::uniform },
    { "radial", init_conds::radialScalar },
    { "spherical", init_conds::sphericalScalar },
    { "x_linear", init_conds::xlinearScalar },
    { "y_linear", init_conds::ylinearScalar },
    { "z_linear", init_conds::zlinearScalar },
    { "tg_vortex", init_conds::tgVortexScalar }
};

#endif // end Header Guard