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
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>
#include <vector>
#include <variant>
#include <algorithm>
#include <map>

#include "string_utils.hpp"

#include "matar.h"
#include "parse_tools.hpp"
#include "parse_dynamic_inputs.hpp"

// simulation parameters contains:
//   mesh_input
//   output_options
//   dynamic_options
//   solver_inputs
//   region_setups
#include "simulation_parameters.hpp"



// ==============================================================================
//   Function Definitions
// ==============================================================================



// =================================================================================
//    Parse Dynamic Options regions
// =================================================================================
void parse_dynamic_options(Yaml::Node& root, dynamic_options_t& dynamic_options)
{
    Yaml::Node& yaml = root["dynamic_options"];

    // get the material variables names set by the user
    std::vector<std::string> user_dynamic_inps;

    // extract words from the input file and validate they are correct
    validate_inputs(yaml, user_dynamic_inps, str_dyn_opts_inps, dyn_opts_required_inps);

    // loop over the words in the material input definition
    for (auto& a_word : user_dynamic_inps) {

        // Start time
        if (a_word.compare("time_initial") == 0) {
            double time_initial = yaml[a_word].As<double>();
            dynamic_options.time_initial = time_initial;
        } // start time
        // End time
        else if (a_word.compare("time_final") == 0) {
            double time_final = yaml[a_word].As<double>();
            dynamic_options.time_final = time_final;
        } // end time
        // Minimum time step
        else if (a_word.compare("dt_min") == 0) {
            double dt_min = yaml[a_word].As<double>();
            dynamic_options.dt_min = dt_min;
        }
        // Maximum time step
        else if (a_word.compare("dt_max") == 0) {
            double dt_max = yaml[a_word].As<double>();
            dynamic_options.dt_max = dt_max;
        }
        // Initial time step
        else if (a_word.compare("dt_start") == 0) {
            double dt_start = yaml[a_word].As<double>();
            dynamic_options.dt_start = dt_start;
        }
        // CFL valid timestep
        else if (a_word.compare("dt_cfl") == 0) {
            double dt_cfl = yaml[a_word].As<double>();
            dynamic_options.dt_cfl = dt_cfl;
        }
        // End cycle count
        else if (a_word.compare("cycle_stop") == 0) {
            int cycle_stop = yaml[a_word].As<int>();
            dynamic_options.cycle_stop = cycle_stop;
        }
        // Machine precision small value
        else if (a_word.compare("fuzz") == 0) {
            double fuzz = yaml[a_word].As<double>();
            dynamic_options.fuzz = fuzz;
        }
        // Very small value
        else if (a_word.compare("tiny") == 0) {
            double tiny = yaml[a_word].As<double>();
            dynamic_options.tiny = tiny;
        }
        // Single precision value
        else if (a_word.compare("small") == 0) {
            double small = yaml[a_word].As<double>();
            dynamic_options.small = small;
        }
        //  Number of RK stages
        else if (a_word.compare("rk_num_stages") == 0) {
            int rk_num_stages = yaml[a_word].As<int>();
            dynamic_options.rk_num_stages = rk_num_stages;
        }
        else {
            std::cout << "ERROR: invalid input: " << a_word << std::endl;
            std::cout << "Valid options are: " << std::endl;
            for (const auto& element : str_dyn_opts_inps) {
                std::cout << element << std::endl;
            }
            throw std::runtime_error("**** User Dynamics Not Understood ****");
        }
    } // end for words in dynamic options
} // end of function to parse region
