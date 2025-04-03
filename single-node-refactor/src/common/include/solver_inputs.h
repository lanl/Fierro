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

#ifndef FIERRO_SOLVER_INPUT_OPTIONS_H
#define FIERRO_SOLVER_INPUT_OPTIONS_H
#include <stdio.h>
#include "matar.h"

namespace solver_input
{
// solver method
enum method
{
    NONE = 0,
    SGH3D = 1,
    SGHRZ = 2,
    SGTM3D = 3,
};
} // end of namespace

static std::map<std::string, solver_input::method> solver_map
{
    { "dynx_FE", solver_input::SGH3D },
    { "dynx_FE_rz", solver_input::SGHRZ },
    { "thrmex_FE", solver_input::SGTM3D }
};
// quasi-static mechanics FE (qz-FE)
// quasi-static thermal-mechanical FE  (qz-thmec-FE)
// quasi-static mechanical GF (qz-GF)
// quasi-static mechanical large-strain GF 

/////////////////////////////////////////////////////////////////////////////
///
/// \structsolver_input_t
///
/// \brief Struct for holding metadata on which solvers are used.
///
/////////////////////////////////////////////////////////////////////////////
struct solver_input_t
{
    solver_input::method method = solver_input::NONE;
}; // solver_input_t

// ----------------------------------
// valid inputs for solver options
// ----------------------------------
static std::vector<std::string> str_solver_inps
{
    "method",
    "id"
};

// ----------------------------------
// required inputs for solver options
// ----------------------------------
static std::vector<std::string> solver_required_inps
{
    "method",
    "id"
};

#endif // end Header Guard