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

#include "Yaml.hpp"
#include "matar.h"

// simulation parameters contains:
//   mesh_input
//   output_options
//   dynamic_options
//   solver_inputs
//   region_setups
#include "simulation_parameters.h"

// boundary conditions
#include "boundary_conditions.h"


// =================================================================================
//    A function to print out all possible inputs when --help flag is given
// =================================================================================
void print_inputs()
{

    //std::vector<std::string>& str_dyn_opts_inps, 
    // mesh_input.h defines str_mesh_inps
    //std::vector<std::string>& str_output_options_inps,
    //std::vector<std::string>& str_bc_inps,
    //std::vector<std::string>& str_bc_surface_inps
    
    std::cout << "\n";
    std::cout << " Fierro: Forging the future of engineering analysis \n\n";

     std::cout << " USE \n";
    std::cout << "    To run Fierro, supply an input file writen in yaml  \n";
    std::cout << "       ./Fierro input.yaml \n\n";

    std::cout << " YAML INPUT OPTIONS \n";

// ---
    std::cout << "    dynamic_options \n";
    for (auto field : str_dyn_opts_inps){
        std::cout << "       "<< field << "\n";
    }
// ---
    std::cout << "    mesh_options \n";
    for (auto field : str_mesh_inps){
        std::cout << "       "<< field << "\n";
    }
// ---
    std::cout << "    output_options \n";
    for (auto field : str_output_options_inps){
        if(field.compare("elem_field_outputs") == 0){
            std::cout << "       elem_field_outputs \n";
            for (const auto& pair : elem_outputs_map) {
                std::cout << "         - "<< pair.first << "\n";
            }
        }
        else {
            std::cout << "       "<< field << "\n";
        }
    }
// ---
    std::cout << "    solver_options \n";
    for (auto field : str_solver_inps){
        std::cout << "       "<< field << "\n";
    }
// ---
    std::cout << "    boundary_conditions \n";
    for (auto field : str_bc_inps){
        if(field.compare("surface") == 0){
            std::cout << "       surface \n";
            for (auto subfield : str_bc_surface_inps){
                std::cout << "          "<< subfield << "\n";
            }
        }
        else {
            std::cout << "       "<< field << "\n";
        }
    } 
    
    std::cout << "  " << std::endl;

    return;


} // end print_inputs


// =================================================================================
//    Extract words from the input file and validate they are correct
// =================================================================================
void validate_inputs(
    Yaml::Node& yaml, 
    std::vector<std::string>& user_inputs, 
    std::vector<std::string>& str_valid_inputs,
    std::vector<std::string>& str_required_inputs)
{
    for (auto item = yaml.Begin(); item != yaml.End(); item++) {
        std::string var_name = (*item).first;


        user_inputs.push_back(var_name);

        // validate input: user_inputs match words in the str_valid_inputs
        if (std::find(str_valid_inputs.begin(), str_valid_inputs.end(), var_name) == str_valid_inputs.end()) {
            std::cout << "ERROR: invalid input: " << var_name << std::endl;
        } // end if variable exists
    } // end for item in this yaml input

    // Add checks for required inputs here
    bool valid = false;
    // Use std::all_of to check if all elements of str_required_inputs are found in user_inputs
    valid = std::all_of(str_required_inputs.begin(), str_required_inputs.end(),[&user_inputs](const std::string& str) {
                return std::find(user_inputs.begin(), user_inputs.end(), str)!= user_inputs.end();
            });

    if (valid == false){
        std::cout << "ERROR: Missing required YAML inputs "<< std::endl;
        std::cout << "Required inputs are:" << std::endl;
        for (const auto& inp : str_required_inputs) {
            std::cout << inp << std::endl;
        }
        throw std::runtime_error("**** Missing required inputs ****");
    }

} // end validate inputs




// =================================================================================
//    Print a yaml file to 6 levels
// =================================================================================
/*
void print_yaml(Yaml::Node root)
{
    if (!VERBOSE) {
        return;
    }

    Yaml::Node& layer0_items = root;

    if (layer0_items.Size() != 0) {
        std::cout << "\n";

        if (VERBOSE) {
            std::cout << "Layer 0 member size = " << layer0_items.Size() << "\n";
        }

        for (auto layer0_item = layer0_items.Begin(); layer0_item != layer0_items.End(); layer0_item++) {
            // print the outlayer variable
            std::cout << (*layer0_item).first << "\n";

            Yaml::Node& layer1_items = (*layer0_item).second;

            // layer 1
            if (layer1_items.Size() != 0) {
                size_t count_layer_1 = 0;

                if (VERBOSE) {
                    std::cout << "\t num_items in this layer = " << layer1_items.Size() << "\n";
                }

                for (auto layer1_item = layer1_items.Begin(); layer1_item != layer1_items.End(); layer1_item++) {
                    std::string text_here1 = (*layer1_item).first;
                    if (text_here1.size() > 0) {
                        std::cout << "\t " << text_here1 << "\n";
                    }
                    else{
                        std::cout << "\t " << count_layer_1 << "\n";
                    }

                    // layer 2
                    Yaml::Node& layer2_items = (*layer1_item).second;

                    if (layer2_items.Size() != 0) {
                        size_t count_layer_2 = 0;

                        if (VERBOSE) {
                            std::cout << "\t\t num_items in this layer = " << layer2_items.Size() << "\n";
                        }

                        for (auto layer2_item = layer2_items.Begin(); layer2_item != layer2_items.End(); layer2_item++) {
                            std::string text_here2 = (*layer2_item).first;
                            if (text_here2.size() > 0) {
                                std::cout << "\t\t " << text_here2 << std::endl;
                            }
                            else{
                                std::cout << "\t\t " << count_layer_2 << "\n";
                            }

                            // layer 3
                            Yaml::Node& layer3_items = (*layer2_item).second;

                            if (layer3_items.Size() != 0) {
                                size_t count_layer_3 = 0;

                                if (VERBOSE) {
                                    std::cout << "\t\t\t num_items in this layer = " << layer3_items.Size() << "\n";
                                }

                                for (auto layer3_item = layer3_items.Begin(); layer3_item != layer3_items.End(); layer3_item++) {
                                    std::string text_here3 = (*layer3_item).first;
                                    if (text_here3.size() > 0) {
                                        std::cout << "\t\t\t " << text_here3 << std::endl;
                                    }
                                    else{
                                        std::cout << "\t\t\t " << count_layer_3 << "\n";
                                    }

                                    // layer 4
                                    Yaml::Node& layer4_items = (*layer3_item).second;

                                    if (layer4_items.Size() != 0) {
                                        size_t count_layer_4 = 0;

                                        if (VERBOSE) {
                                            std::cout << "\t\t\t\t num_items in layer 4 = " << layer4_items.Size() << "\n";
                                        }

                                        for (auto layer4_item = layer4_items.Begin(); layer4_item != layer4_items.End(); layer4_item++) {
                                            std::string text_here4 = (*layer4_item).first;
                                            if (text_here4.size() > 0) {
                                                std::cout << "\t\t\t\t " << text_here4 << std::endl;
                                            }
                                            else{
                                                std::cout << "\t\t\t\t " << count_layer_4 << "\n";
                                            }

                                            // layer 5
                                            Yaml::Node& layer5_items = (*layer4_item).second;

                                            if (layer5_items.Size() != 0) {
                                                size_t count_layer_5 = 0;
                                                if (VERBOSE) {
                                                    std::cout << "\t\t\t\t\t num_items in layer 5 = " << layer5_items.Size() << "\n";
                                                }

                                                for (auto layer5_item = layer5_items.Begin(); layer5_item != layer5_items.End(); layer5_item++) {
                                                    std::string text_here5 = (*layer5_item).first;
                                                    if (text_here5.size() > 0) {
                                                        std::cout << "\t\t\t\t\t " << text_here5 << std::endl;
                                                    }
                                                    else{
                                                        std::cout << "\t\t\t\t\t " << count_layer_5 << "\n";
                                                    }

                                                    // layer 6
                                                    Yaml::Node& layer6_items = (*layer5_item).second;

                                                    if (layer6_items.Size() != 0) {
                                                        size_t count_layer_6 = 0;
                                                        if (VERBOSE) {
                                                            std::cout << "\t\t\t\t\t\t num_items in layer 6 = " << layer6_items.Size() << "\n";
                                                        }

                                                        for (auto layer6_item = layer6_items.Begin(); layer6_item != layer6_items.End(); layer6_item++) {
                                                            std::string text_here6 = (*layer6_item).first;
                                                            if (text_here6.size() > 0) {
                                                                std::cout << "\t\t\t\t\t\t layer 6 = " << text_here6 << std::endl;
                                                            }
                                                            else{
                                                                std::cout << "\t\t\t\t\t\t " << count_layer_6 << "\n";
                                                            }

                                                            count_layer_6++;
                                                        } // end loop over layer 6
                                                    } // end of layer6

                                                    count_layer_5++;
                                                } // end loop over layer 5
                                            } // end if layer5 exists

                                            count_layer_4++;
                                        } // end loop over layer4
                                    } // end if layer4 exists

                                    count_layer_3++;
                                } // end loop over layer 3 items
                            } // end if layer 3 exists

                            count_layer_2++;
                        } // end if loop over layer 2
                    } // end if layer 2 exists

                    count_layer_1++;
                } // end loop over layer 1
            } // end if layer 1 exists
        } // end loop over layer 0
    } // end if layer0 exists
} // end print yaml function
*/
