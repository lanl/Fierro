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

#include "parse_yaml.hpp"

#include "parse_tools.hpp"
#include "parse_output_options.hpp"

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
//    Parse Output options
// =================================================================================
void parse_output_options(Yaml::Node& root, 
                          output_options_t& output_options)
{
    Yaml::Node& out_opts = root["output_options"];

    // get the mesh variables names set by the user
    std::vector<std::string> user_inputs;

    // extract words from the input file and validate they are correct
    validate_inputs(out_opts, user_inputs, str_output_options_inps, output_options_required_inps);

    // loop over the output options
    for (auto& a_word : user_inputs) {

        // get output format
        if (a_word.compare("output_file_format") == 0) {
            std::string format = root["output_options"][a_word].As<std::string>();

            auto map = output_format_map;

            // set the output format
            if (map.find(format) != map.end()) {
                output_options.format = map[format];
            }
            else{
                std::cout << "ERROR: invalid output option input in YAML file: " << format << std::endl;
                std::cout << "Valid options are: " << std::endl;

                for (const auto& pair : map) {
                    std::cout << "\t" << pair.first << std::endl;
                }
                throw std::runtime_error("**** Output Format Not Understood ****");
            } // end if
        } // output format
        // get timer_output_level
        else if (a_word.compare("timer_output_level") == 0) {
            std::string timer_level = root["output_options"][a_word].As<std::string>();

            auto map = timer_output_level_map;

            // set the timer_output_level
            if (map.find(timer_level) != map.end()) {
                output_options.timer_level = map[timer_level];
            }
            else{
                std::cout << "ERROR: invalid timer output option input in YAML file: " << timer_level << std::endl;
                std::cout << "Valid options are: " << std::endl;

                for (const auto& pair : map) {
                    std::cout << "\t" << pair.first << std::endl;
                }
                throw std::runtime_error("**** Time Output Level Syntax Not Understood ****");
            } // end if
        } // timer_level
        // Graphics time step
        else if (a_word.compare("graphics_time_step") == 0) {
            real_t graphics_time_step = root["output_options"][a_word].As<real_t>();

            output_options.graphics_time_step = graphics_time_step;
        } // graphics_time_step
        // Graphics iteration step
        else if (a_word.compare("graphics_iteration_step") == 0) {
            int graphics_iteration_step = root["output_options"][a_word].As<int>();

            output_options.graphics_iteration_step = graphics_iteration_step;
        } // graphics_iteration_step
        else if (a_word.compare("elem_field_outputs") == 0) {
                Yaml::Node & output_vars_yaml = root["output_options"][a_word];

                size_t num_output_vars = output_vars_yaml.Size();

                // loop over the output vars
                for (int var_id = 0; var_id < num_output_vars; var_id++){
                    std::string var_name = root["output_options"][a_word][var_id].As<std::string>();
                    
                    // save the enum fields to the outputs
                    if (elem_outputs_map.find(var_name) != elem_outputs_map.end()) {
                        auto field_var = elem_outputs_map[var_name]; // get the enum for this field variabled
                        
                        output_options.output_elem_state.push_back(field_var);
                    }
                    else{
                        throw std::runtime_error("**** Element Field Ouput Variable Name Not Understood ****");
                    } // end if

                } // end for over variables
        } // end of elem fields outputs
        else if (a_word.compare("node_field_outputs") == 0) {
                Yaml::Node & output_vars_yaml = root["output_options"][a_word];

                size_t num_output_vars = output_vars_yaml.Size();

                // loop over the output vars
                for (int var_id = 0; var_id < num_output_vars; var_id++){
                    std::string var_name = root["output_options"][a_word][var_id].As<std::string>();
                    
                    // save the enum fields to the outputs
                    if (node_outputs_map.find(var_name) != node_outputs_map.end()) {
                        auto field_var = node_outputs_map[var_name]; // get the enum for this field variabled
                        
                        output_options.output_node_state.push_back(field_var);
                    }
                    else{
                        throw std::runtime_error("**** Node Field Ouput Variable Name Not Understood ****");
                    } // end if

                } // end for over variables
        } // end of nodal fields outputs
        else if (a_word.compare("gauss_pt_field_outputs") == 0) {
                Yaml::Node & output_vars_yaml = root["output_options"][a_word];

                size_t num_output_vars = output_vars_yaml.Size();

                // loop over the output vars
                for (int var_id = 0; var_id < num_output_vars; var_id++){
                    std::string var_name = root["output_options"][a_word][var_id].As<std::string>();
                    
                    // save the enum fields to the outputs
                    if (gauss_pt_outputs_map.find(var_name) != gauss_pt_outputs_map.end()) {
                        auto field_var = gauss_pt_outputs_map[var_name]; // get the enum for this field variabled
                        
                        output_options.output_gauss_pt_state.push_back(field_var);
                    }
                    else{
                        throw std::runtime_error("**** Gauss Pnt Field Ouput Variable Name Not Understood ****");
                    } // end if

                } // end for over variables
        } // end of gauss_pt fields outputs
        else if (a_word.compare("mat_pt_field_outputs") == 0) {
                Yaml::Node & output_vars_yaml = root["output_options"][a_word];

                size_t num_output_vars = output_vars_yaml.Size();

                // loop over the output vars
                for (int var_id = 0; var_id < num_output_vars; var_id++){
                    std::string var_name = root["output_options"][a_word][var_id].As<std::string>();
                    
                    // save the enum fields to the outputs
                    if (mat_pt_outputs_map.find(var_name) != mat_pt_outputs_map.end()) {
                        auto field_var = mat_pt_outputs_map[var_name]; // get the enum for this field variabled
                        
                        output_options.output_mat_pt_state.push_back(field_var);
                    }
                    else{
                        throw std::runtime_error("**** Mat Pnt Field Ouput Variable Name Not Understood ****");
                    } // end if

                } // end for over variables
        } // end of elem fields outputs
        else {
            std::cout << "ERROR: invalid input: " << a_word << std::endl;

            std::cout << "Valid options are: " << std::endl;
            for (const auto& element : str_output_options_inps) {
                std::cout << element << std::endl;
            }
            throw std::runtime_error("**** Graphics Not Understood ****");
        }
    } // end user_inputs
} // end of parse mesh options



