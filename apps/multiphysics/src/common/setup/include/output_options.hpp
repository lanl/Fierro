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

#ifndef FIERRO_OUTPUT_OPTIONS_H
#define FIERRO_OUTPUT_OPTIONS_H
#include "matar.h"
#include "state.hpp"

namespace output_options
{
// output file options
enum format
{
    viz = 0,
    state = 1,
    viz_and_state = 2,
    ensight = 3
};

// timer output level
enum timer_output_level
{
    thorough = 0,
};
} // end of namespace

extern std::map<std::string, output_options::format> output_format_map;

extern std::map<std::string, output_options::timer_output_level> timer_output_level_map;


// nodal state variables writen to file
extern std::map<std::string, node_state> node_outputs_map;

// gauss point state variables writen to file
extern std::map<std::string, gauss_pt_state> gauss_pt_outputs_map;

// material point state variables writen to file
extern std::map<std::string, material_pt_state> mat_pt_outputs_map;

// element average state variables writen to file
extern std::map<std::string, material_pt_state> elem_outputs_map;


/////////////////////////////////////////////////////////////////////////////
///
/// \struct OutputOptions_t
///
/// \brief Output related options for a Fierro simulation
///
/////////////////////////////////////////////////////////////////////////////
struct OutputOptions_t
{
    output_options::format format;  ///< Format for the output files
    output_options::timer_output_level timer_level; ///< How often to output mesh WARNING: CURRENTLY UNUSED

    real_t graphics_time_step   = real_t(1.0);  ///< How often to write a graphics dump in time
    int graphics_iteration_step = 2000000;  ///< How often to write a graphics dump by iteration count

    std::vector<material_pt_state> output_elem_state;
    std::vector<node_state> output_node_state;
    std::vector<gauss_pt_state> output_gauss_pt_state;
    std::vector<material_pt_state> output_mat_pt_state;

}; // OutputOptions_t

// ----------------------------------
// valid inputs for output options
// ----------------------------------
extern std::vector<std::string> str_output_options_inps;

// ----------------------------------
// required inputs for output options
// ----------------------------------
extern std::vector<std::string> output_options_required_inps;

#endif // end Header Guard