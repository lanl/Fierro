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

#include "string_utils.h"

#include "matar.h"
#include "parse_tools.hpp"
#include "parse_solver_inputs.hpp"

// simulation parameters contains:
//   mesh_input
//   output_options
//   dynamic_options
//   solver_inputs
//   region_setups
#include "simulation_parameters.h"



// ==============================================================================
//   Function Definitions
// ==============================================================================


// =================================================================================
//    Parse Solver options
// =================================================================================
void parse_solver_input(Yaml::Node& root, std::vector<solver_input_t>& solver_inputs)
{
    Yaml::Node& solver_yaml = root["solver_options"];

    size_t num_solvers = solver_yaml.Size();

    solver_inputs = std::vector<solver_input_t>(num_solvers);

    // a check on solverl_id not being specified more than once or not at all
    CArray <bool> check_solver_ids(num_solvers);
    check_solver_ids.set_values(false);

    // loop over the solvers specified in the YAML file
    for (int s_id = 0; s_id < num_solvers; s_id++) {


        // read the variables names
        Yaml::Node& inps_yaml = root["solver_options"][s_id]["solver"];

        // get the solver variables names set by the user
        std::vector<std::string> user_inputs;

        // extract words from the input file and validate they are correct
        validate_inputs(inps_yaml, user_inputs, str_solver_inps, solver_required_inps);


        // loop over the words in the solver input definition and find the solver id
        int solver_id = -1;
        for (auto& a_word : user_inputs) {

            Yaml::Node& solver_inps_yaml = root["solver_options"][s_id]["solver"][a_word];

            if (a_word.compare("id") == 0) {
                solver_id = root["solver_options"][s_id]["solver"]["id"].As<int>();

                if (solver_id<0 || solver_id>=num_solvers){
                    std::cout << "ERROR: invalid solver id specified in the solver definition " << std::endl;
            
                    throw std::runtime_error("**** Solver id is out of bounds ****");
                } // end check on solver_id range

                if (check_solver_ids(solver_id) == true){
                    std::cout << "ERROR: solver id = " << solver_id << " was already specified "<< std::endl;
                    throw std::runtime_error("**** Multiple solvers used the same solver_id ****");
                }
                else {
                    check_solver_ids(solver_id) = true;
                } // end check on solver_id

            } // end id
        } // end loop over all solver inputs for this solver

        if (solver_id<0){
            std::cout << "ERROR: solver id must be specified in the solver definition " << std::endl;
            
            throw std::runtime_error("**** Solver id is missing ****");
        } // end check on solver_id specified


        // loop over the words in the input
        for (auto& a_word : user_inputs) {

            // get solver method
            if (a_word.compare("method") == 0) {
                // input order is s_id, but we save things in the array using solver_id
                std::string method = root["solver_options"][s_id]["solver"][a_word].As<std::string>();

                auto map = solver_map;

                // set the method
                if (map.find(method) != map.end()) {
                    solver_inputs[solver_id].method = map[method];  // save it to solver_id value, input order may differ
                }
                else{
                    std::cout << "ERROR: invalid method option input in YAML file: " << method << std::endl;
                    std::cout << "Valid options are: " << std::endl;

                    for (const auto& pair : map) {
                        std::cout << "\t" << pair.first << std::endl;
                    }
                    throw std::runtime_error("**** Solver Input Method Not Understood ****");
                } // end if
            } // method
            else if (a_word.compare("id") == 0) {
                // do nothing, we already got the id
            }
            else if (a_word.compare("time_end") == 0){
                double t_end = root["solver_options"][s_id]["solver"]["time_end"].As<double>();

                if (t_end<0){
                    std::cout << "ERROR: invalid ending time specified in the solver definition " << std::endl;
            
                    throw std::runtime_error("**** Solver time is negative ****");
                } // end check on time range

                solver_inputs[solver_id].time_end = t_end;
            }
            // add solver_vars parsing here
            else {
                std::cout << "ERROR: invalid input: " << a_word << std::endl;

                std::cout << "Valid options are: " << std::endl;

                for (const auto& element : str_solver_inps) {
                    std::cout << element << std::endl;
                }

                throw std::runtime_error("**** Solver Inputs Not Understood ****");
            }
        } // end loop over solver options
    } // end loop over solvers
} // end of parse solver options
