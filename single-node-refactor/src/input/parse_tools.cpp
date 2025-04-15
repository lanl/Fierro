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


// check to see if a word is in the string, the sentence is a string
bool contains_word(const std::string& text, const std::string& word) {
  return text.find(word) != std::string::npos;
}

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
    std::cout << "    Fierro: Forging the future of engineering analysis                       \n\n";
    
    std::cout << " DESCRIPTION                                                                 \n";
    std::cout << "    A multi-solver multi-physics multi-material code that runs in parallel   \n";
    std::cout << "    across CPU and GPU computer architectures.   \n\n";

    std::cout << " SOLVERS                                                                     \n";
    std::cout << "    dynx_FE --    A conservative transient finite element solver for high    \n";
    std::cout << "                  strain rate, shock-driven problems in 3D Cartesian         \n";
    std::cout << "                  coordinates. This solver is for unstructured meshes with   \n";
    std::cout << "                  linear hexahedral elements.                                \n";
    std::cout << "    dynx_FE_rz -- A conservative transient finite element solver for high    \n";
    std::cout << "                  strain rate, shock-driven problems in 2D axisymmetric      \n";
    std::cout << "                  coordinates. The axis of symmetry (i.e., axis of rotation) \n";
    std::cout << "                  is the horizontal-axis. This solver preservers symmetry on \n";
    std::cout << "                  on polar meshes with uniform angular spacing.  This solver \n"; 
    std::cout << "                  is for unstructured meshes with linear quadralateral       \n";
    std::cout << "                  elements.                                                  \n";
    std::cout << "    thrmex_FE  -- A conservative transient finite element solver for thermal \n";
    std::cout << "                  mechanical problems in 3D Cartesian coordinates. This      \n";
    std::cout << "                  solver is for unstructured meshes with linear hexahedral   \n";
    std::cout << "                  elements.                                                  \n\n";

    std::cout << " USE                                                                         \n";
    std::cout << "    To run Fierro, supply an input file written in yaml                      \n";
    std::cout << "       ./Fierro input.yaml                                                   \n\n";

    std::cout << " YAML INPUT                                                                  \n";
    std::cout << "    Example input files are provided in,                                     \n";
    std::cout << "        Fierro/single-node-refactor/example_inputs                           \n\n";

    std::cout << "...press Enter to continue or 'q' to quit...                                 \n"; 
    char letter = std::cin.get();  // Waits for user to press Enter
    if (letter == 'q' || letter == 'Q') {
        return;
    } // end if

    std::cout << "    All input key words are listed below here. In some cases, the options, a \n";
    std::cout << "    number or blank fields are provided to show the input structure (e.g. a  \n";
    std::cout << "    vector input). Words like <word> are possible string inputs to use.      \n";
    std::cout << "    Double and interger values are also denoted. The user will need to input \n";
    std::cout << "    a double or integer value.  Inputs left blank are doubles or integers.   \n\n";

// ---
    std::cout << "    dynamic_options: \n";
    for (auto field : str_dyn_opts_inps){
        std::cout << "       "<< field << ":\n";
    }
// ---
    std::cout << "    mesh_options: \n";
    for (auto field : str_mesh_inps){
        if(field.compare("origin") == 0){
            std::cout << "       "<< field << ": [double, double, double]\n";
        }
        else if(field.compare("source") == 0){
            std::cout << "       source:";
            for (const auto& pair : mesh_input_source_map) {
                std::cout << " <" << pair.first << ">";
            }
            std::cout << "\n";
        }
        else if(field.compare("type") == 0){
            std::cout << "       type:";
            for (const auto& pair : mesh_input_type_map) {
                std::cout << " <" << pair.first << ">";
            }
            std::cout << "\n";
        }
        else if(field.compare("length") == 0){
            std::cout << "       "<< field << ": [double, double, double]\n";
        }
        else if(field.compare("num_elems") == 0){
            std::cout << "       "<< field << ": [int, int, int]\n";
        }
        else {
            std::cout << "       "<< field << ":\n";
        }
    }
// ---
    std::cout << "    output_options: \n";
    for (auto field : str_output_options_inps){
        if(field.compare("elem_field_outputs") == 0){
            std::cout << "       elem_field_outputs: \n";
            for (const auto& pair : elem_outputs_map) {
                std::cout << "         - "<< pair.first << "\n";
            }
        }
        else if(field.compare("timer_output_level:") == 0){
            std::cout << "       timer_output_level: \n";
            for (const auto& pair :  timer_output_level_map) {
                std::cout << "         - "<< pair.first << "\n";
            }
        }
        else if(field.compare("output_file_format:") == 0){
            std::cout << "       output_file_format: \n";
            for (const auto& pair :  output_format_map) {
                std::cout << "         - "<< pair.first << "\n";
            }
        }
        else if(field.compare("node_field_outputs") == 0){
            std::cout << "       node_field_outputs: \n";
            for (const auto& pair : node_outputs_map) {
                std::cout << "         - "<< pair.first << "\n";
            }
        }
        else if(field.compare("mat_pt_field_outputs") == 0){
            std::cout << "       mat_pt_field_outputs: \n";
            for (const auto& pair : mat_pt_outputs_map) {
                std::cout << "         - "<< pair.first << "\n";
            }
        }
        else if(field.compare("gauss_pt_field_outputs") == 0){
            std::cout << "       gauss_pt_field_outputs: \n";
            for (const auto& pair : gauss_pt_outputs_map) {
                std::cout << "         - "<< pair.first << "\n";
            }
        }
        else {
            std::cout << "       "<< field << ":\n";
        }
    }
// ---
    std::cout << "    solver_options: \n";
    std::cout << "       #...as many solvers as you need... \n";
    std::cout << "       - solver: \n";
    for (auto field : str_solver_inps){
        if(field.compare("method") == 0){
            std::cout << "          method:";
            for (const auto& pair : solver_map) {
                std::cout << " <" << pair.first << ">";
            }
            std::cout << "\n";
        }
        else {
            std::cout << "          "<< field << ":\n";
        }
    }
// ---
    std::cout << "    boundary_conditions: \n";
    std::cout << "       #...as many conditions as you need... \n";
    std::cout << "       - boundary_condition: \n";
    for (auto field : str_bc_inps){
        if(field.compare("surface") == 0){
            std::cout << "          surface: \n";

            for (auto subfield : str_bc_surface_inps){

                if(subfield.compare("type") == 0){

                    std::cout << "             type:";
                    for (const auto& pair : bc_surface_map) {
                        std::cout << " <" << pair.first << ">";
                    }
                    std::cout << "\n";

                } // end if type
                else if(subfield.compare("origin") == 0){
                    std::cout << "             "<< subfield << ": [double, double, double]\n";
                }
                else {
                    std::cout << "             "<< subfield << ":\n";
                } // end if type

            } // end for

        }
        else if(field.compare("velocity_model") == 0){
            std::cout << "          velocity_model:";
            size_t count = 0;
            for (const auto& pair : bc_velocity_model_map) {
                if(count==4){
                    std::cout << "\n"; // new line
                    std::cout << "                         "; // tab in
                    count=0;
                }
                std::cout << " <" << pair.first << ">";
                count++;
            }
            std::cout << "\n";
        }
        else if(field.compare("stress_model") == 0){
            std::cout << "          stress_model:";
            for (const auto& pair : bc_stress_model_map) {
                std::cout << " <" << pair.first << ">";
            }
            std::cout << "\n";
        }
        else if(field.compare("temperature_model") == 0){
            std::cout << "          temperature_model:";
            for (const auto& pair : bc_temperature_model_map) {
                std::cout << " <" << pair.first << ">";
            }
            std::cout << "\n";
        }
        else if(contains_word(field, "global_vars")){
            std::cout << "          " << field <<":\n";
            std::cout << "             - double/int \n";
            std::cout << "             - double/int \n";
            std::cout << "             #...as many values as you need... \n";
            std::cout << "             - double/int \n";
        }
        else {
            std::cout << "          "<< field << ":\n";
        }
    } // end bcs
 // ---
    std::cout << "    materials: \n";
    std::cout << "       #...as many materials as you need... \n";
    std::cout << "       - material: \n";
    for (auto field : str_material_inps){
        
        if(field.compare("eos_model") == 0){
            std::cout << "          eos_model:";
            size_t count=0;
            for (const auto& pair : eos_models_map) {
                if(count==2){
                    std::cout << "\n"; // new line
                    std::cout << "                    "; // tab in
                    count=0;
                }
                std::cout << " <" << pair.first << ">";
                count++;
            }
            std::cout << "\n";
        }
        else if(field.compare("eos_model_type") == 0){
            std::cout << "          eos_model_type:";
            for (const auto& pair : eos_type_map) {
                std::cout << " <" << pair.first << ">";
            }
            std::cout << "\n";
        }
        else if(field.compare("strength_model") == 0){
            std::cout << "          strength_model:";
            size_t count=0;
            for (const auto& pair : strength_models_map) {
                if(count==2){
                    std::cout << "\n"; // new line
                    std::cout << "                         "; // tab in
                    count=0;
                }
                std::cout << " <" << pair.first << ">";
                count++;
            }
            std::cout << "\n";
        }
        else if(field.compare("strength_model_type") == 0){
            std::cout << "          strength_model_type:";
            for (const auto& pair : strength_type_map) {
                std::cout << " <" << pair.first << ">";
            }
            std::cout << "\n";
        }
        else if(field.compare("dissipation_model") == 0){
            std::cout << "          dissipation_model:";
            size_t count=0;
            for (const auto& pair : dissipation_model_map) {
                if(count==2){
                    std::cout << "\n"; // new line
                    std::cout << "                            "; // tab in
                    count=0;
                }
                std::cout << " <" << pair.first << ">";
                count++;
            }
            std::cout << "\n";
        }
        else if(field.compare("erosion_model") == 0){
            std::cout << "          erosion_model:";
            for (const auto& pair : erosion_model_map) {
                std::cout << " <" << pair.first << ">";
            }
            std::cout << "\n";
        }
        else if(contains_word(field, "global_vars")){
            std::cout << "          " << field <<":\n";
            std::cout << "             - double/int \n";
            std::cout << "             - double/int \n";
            std::cout << "             #...as many values as you need... \n";
            std::cout << "             - double/int \n";
        }
        else {
            std::cout << "          "<< field << ":\n";
        }
    } // end mats
// ---
    std::cout << "    regions: \n";
    std::cout << "       #...as many regions as you need... \n";
    std::cout << "       - region: \n";
    for (auto field : str_region_inps){
        if(field.compare("volume") == 0){
            std::cout << "          volume: \n";

            for (auto subfield : str_region_volume_inps){

                if(subfield.compare("type") == 0){

                    std::cout << "             type:";
                    for (const auto& pair : region_type_map) {
                        std::cout << " <" << pair.first << ">";
                    }
                    std::cout << "\n";

                } // end if type
                else if(subfield.compare("origin") == 0){
                    std::cout << "             "<< subfield << ": [double, double, double]\n";
                }
                else {
                    std::cout << "             "<< subfield << ":\n";
                } // end if type

            } // end for
        } // end if volume
        else if(field.compare("volume_fraction") == 0){
            std::cout << "          volume_fraction: \n";

            for (auto subfield : str_region_volfrac_inps){

                if(subfield.compare("type") == 0){

                    std::cout << "             type:";
                    for (const auto& pair : scalar_ics_type_map) {
                        std::cout << " <" << pair.first << ">";
                    }
                    std::cout << "\n";

                } // end if type
                else if(subfield.compare("origin") == 0){
                    std::cout << "             "<< subfield << ": [double, double, double]\n";
                }
                else {
                    std::cout << "             "<< subfield << ":\n";
                } // end if type

            } // end for
        } // end if volfrac
        else if(field.compare("density") == 0){
            std::cout << "          density: \n";

            for (auto subfield : str_region_den_inps){

                if(subfield.compare("type") == 0){

                    std::cout << "             type:";
                    for (const auto& pair : scalar_ics_type_map) {
                        std::cout << " <" << pair.first << ">";
                    }
                    std::cout << "\n";

                } // end if type
                else if(subfield.compare("origin") == 0){
                    std::cout << "             "<< subfield << ": [double, double, double]\n";
                }
                else {
                    std::cout << "             "<< subfield << ":\n";
                } // end if type

            } // end for
        } // end if density
        else if(field.compare("specific_internal_energy") == 0){
            std::cout << "          specific_internal_energy: \n";

            for (auto subfield : str_region_sie_inps){

                if(subfield.compare("type") == 0){

                    std::cout << "             type:";
                    for (const auto& pair : scalar_ics_type_map) {
                        std::cout << " <" << pair.first << ">";
                    }
                    std::cout << "\n";

                } // end if type
                else if(subfield.compare("origin") == 0){
                    std::cout << "             "<< subfield << ": [double, double, double]\n";
                }
                else {
                    std::cout << "             "<< subfield << ":\n";
                } // end if type

            } // end for
        } // end if sie
        else if(field.compare("internal_energy") == 0){
            std::cout << "          internal_energy: \n";

            for (auto subfield : str_region_ie_inps){

                if(subfield.compare("type") == 0){

                    std::cout << "             type:";
                    for (const auto& pair : scalar_ics_type_map) {
                        std::cout << " <" << pair.first << ">";
                    }
                    std::cout << "\n";

                } // end if type
                else if(subfield.compare("origin") == 0){
                    std::cout << "             "<< subfield << ": [double, double, double]\n";
                }
                else {
                    std::cout << "             "<< subfield << ":\n";
                } // end if type

            } // end for
        } // end if ie
        else if(field.compare("specific_heat") == 0){
            std::cout << "          specific_heat: \n";

            for (auto subfield : str_region_specific_heat_inps){

                if(subfield.compare("type") == 0){

                    std::cout << "             type:";
                    for (const auto& pair : scalar_ics_type_map) {
                        std::cout << " <" << pair.first << ">";
                    }
                    std::cout << "\n";

                } // end if type
                else if(subfield.compare("origin") == 0){
                    std::cout << "             "<< subfield << ": [double, double, double]\n";
                }
                else {
                    std::cout << "             "<< subfield << ":\n";
                } // end if type

            } // end for
        } // end if specific heat
        else if(field.compare("thermal_conductivity") == 0){
            std::cout << "          thermal_conductivity: \n";

            for (auto subfield : str_region_thermal_conductivity_inps){

                if(subfield.compare("type") == 0){

                    std::cout << "             type:";
                    for (const auto& pair : scalar_ics_type_map) {
                        std::cout << " <" << pair.first << ">";
                    }
                    std::cout << "\n";

                } // end if type
                else if(subfield.compare("origin") == 0){
                    std::cout << "             "<< subfield << ": [double, double, double]\n";
                }
                else {
                    std::cout << "             "<< subfield << ":\n";
                } // end if type

            } // end for
        } // end if thermal_conductivity
        else if(field.compare("temperature") == 0){
            std::cout << "          temperature: \n";

            for (auto subfield : str_region_temperature_inps){

                if(subfield.compare("type") == 0){

                    std::cout << "             type:";
                    for (const auto& pair : scalar_ics_type_map) {
                        std::cout << " <" << pair.first << ">";
                    }
                    std::cout << "\n";

                } // end if type
                else if(subfield.compare("origin") == 0){
                    std::cout << "             "<< subfield << ": [double, double, double]\n";
                }
                else {
                    std::cout << "             "<< subfield << ":\n";
                } // end if type

            } // end for
        } // end if temperature
        else if(field.compare("velocity") == 0){
            std::cout << "          velocity: \n";
            

            for (auto subfield : str_region_vel_inps){

                if(subfield.compare("type") == 0){
                    std::cout << "             type:";

                    size_t count=0;
                    for (const auto& pair : vector_ics_type_map) {
                        if(count==4){
                            std::cout << "\n"; // new line
                            std::cout << "                  "; // tab in
                            count=0;
                        }
                        std::cout << " <" << pair.first << ">";
                        count++;
                    }
                    std::cout << "\n";

                } // end if type
                else if(subfield.compare("origin") == 0){
                    std::cout << "             "<< subfield << ": [double, double, double]\n";
                }
                else {
                    std::cout << "             "<< subfield << ":\n";
                } // end if type

            } // end for
        } // end if sie
        else{
            std::cout << "          "<< field << ":\n";
        } // end if

    } // end for


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
