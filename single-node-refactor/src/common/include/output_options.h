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

#ifndef FIERRO_OUTPUT_OPTIONS_H
#define FIERRO_OUTPUT_OPTIONS_H
#include <stdio.h>
#include "matar.h"

namespace output_options
{
// output file options
enum format
{
    vtk = 0,
    ensight = 1,
    state = 2
};

// timer output level
enum timer_output_level
{
    thorough = 0,
};
} // end of namespace

static std::map<std::string, output_options::format> output_format_map
{
    { "vtk", output_options::vtk },
    { "ensight", output_options::ensight },
    { "state", output_options::state }
};

static std::map<std::string, output_options::timer_output_level> timer_output_level_map
{
    { "thorough", output_options::thorough }
};

/////////////////////////////////////////////////////////////////////////////
///
/// \struct output_options_t
///
/// \brief Output related options for a Fierro simulation
///
/////////////////////////////////////////////////////////////////////////////
struct output_options_t
{
    output_options::format format;  ///< Format for the output files
    output_options::timer_output_level timer_level; ///< How often to output mesh WARNING: CURRENTLY UNUSED

    double graphics_time_step   = 1.0;  ///< How often to write a graphics dump in time
    int graphics_iteration_step = 2000000;  ///< How often to write a graphics dump by iteration count
}; // output_options_t

// ----------------------------------
// valid inputs for output options
// ----------------------------------
static std::vector<std::string> str_output_options_inps
{
    "timer_output_level",
    "output_file_format",
    "graphics_time_step",
    "graphics_iteration_step"
};

// ----------------------------------
// required inputs for output options
// ----------------------------------
static std::vector<std::string> output_options_required_inps
{
};

#endif // end Header Guard