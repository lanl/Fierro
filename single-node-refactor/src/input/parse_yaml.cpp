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
#include "parse_yaml.h"
#include "simulation_parameters.h"
#include "boundary_conditions.h"

// velocity bc files
#include "constant_velocity_bc.h"
#include "no_velocity_bc.h"
#include "piston_velocity_bc.h"
#include "reflected_velocity_bc.h"
#include "time_varying_velocity_bc.h"
#include "user_defined_velocity_bc.h"
#include "zero_velocity_bc.h"

// temperature bc files
#include "constant_temp_bc.h"


// eos files
#include "gamma_law_eos.h"
#include "no_eos.h"
#include "user_defined_eos.h"
#include "void_eos.h"

// strength
#include "no_strength.h"
#include "user_defined_strength.h"

// erosion files
#include "basic_erosion.h"
#include "no_erosion.h"

// fracture files
#include "user_defined_fracture.h"


#define PI 3.141592653589793

bool VERBOSE = false;

// ==============================================================================
//   Function Definitions
// ==============================================================================

// =================================================================================
//    Print a yaml file to 6 levels
// =================================================================================
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

// =================================================================================
//    Parse YAML file
// =================================================================================
void parse_yaml(Yaml::Node& root, SimulationParameters_t& SimulationParamaters, Material_t& Materials, BoundaryCondition_t& Boundary)
{
    // if (VERBOSE) {
    //     printf("\n");
    //     std::cout << "Printing YAML Input file:" << std::endl;
    // }
    // print the input file
    print_yaml(root);

    if (VERBOSE) {
        printf("\n");
        std::cout << "Parsing YAML meshing options:" << std::endl;
    }
    parse_mesh_input(root, SimulationParamaters.mesh_input);

    if (VERBOSE) {
        printf("\n");
        std::cout << "Parsing YAML dynamic options:" << std::endl;
    }
    parse_dynamic_options(root, SimulationParamaters.dynamic_options);

    if (VERBOSE) {
        printf("\n");
        std::cout << "Parsing YAML output options:" << std::endl;
    }
    parse_output_options(root, SimulationParamaters.output_options);

    if (VERBOSE) {
        printf("\n");
        std::cout << "Parsing YAML solver options:" << std::endl;
    }
    parse_solver_input(root, SimulationParamaters.solver_inputs);

    if (VERBOSE) {
        printf("\n");
        std::cout << "Parsing YAML boundary condition options:" << std::endl;
    }
    parse_bcs(root, Boundary);

    if (VERBOSE) {
        printf("\n");
        std::cout << "Parsing YAML regions:" << std::endl;
    }
    // parse the region yaml text into a vector of region_fills
    parse_regions(root, 
                  SimulationParamaters.region_fills,
                  SimulationParamaters.region_fills_host);


    if (VERBOSE) {
        printf("\n");
        std::cout << "Parsing YAML materials:" << std::endl;
    }
    // parse the material yaml text into a vector of materials
    parse_materials(root, Materials);
}

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

        // print the variable
        if (VERBOSE) {
            std::cout << "This is var name = " << var_name << "\n";
        }

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
//    Parse Solver options
// =================================================================================
void parse_solver_input(Yaml::Node& root, std::vector<solver_input_t>& solver_input)
{
    Yaml::Node& solver_yaml = root["solver_options"];

    size_t num_solvers = solver_yaml.Size();

    if (VERBOSE) {
        std::cout << "Num solvers = " << num_solvers << std::endl;
    }

    solver_input = std::vector<solver_input_t>(num_solvers);

    // loop over the solvers specified in the YAML file
    for (int solver_id = 0; solver_id < num_solvers; solver_id++) {
        // read the variables names
        Yaml::Node& inps_yaml = root["solver_options"][solver_id]["solver"];

        // get the solver variables names set by the user
        std::vector<std::string> user_inputs;

        // extract words from the input file and validate they are correct
        validate_inputs(inps_yaml, user_inputs, str_solver_inps, solver_required_inps);

        // loop over the words in the input
        for (auto& a_word : user_inputs) {
            if (VERBOSE) {
                std::cout << a_word << std::endl;
            }

            // get solver method
            if (a_word.compare("method") == 0) {
                std::string method = root["solver_options"][solver_id]["solver"][a_word].As<std::string>();

                auto map = solver_map;

                // set the method
                if (map.find(method) != map.end()) {
                    solver_input[solver_id].method = map[method];
                    if (VERBOSE) {
                        std::cout << "\tmethod = " << method << std::endl;
                    }
                }
                else{
                    std::cout << "ERROR: invalid method option input in YAML file: " << method << std::endl;
                    std::cout << "Valid options are: " << std::endl;

                    for (const auto& pair : map) {
                        std::cout << "\t" << pair.first << std::endl;
                    }
                } // end if
            } // method
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
        if (VERBOSE) {
            std::cout << a_word << std::endl;
        }

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
        //  Number of RK bins
        else if (a_word.compare("rk_num_bins") == 0) {
            int rk_num_bins = yaml[a_word].As<int>();
            dynamic_options.rk_num_bins = rk_num_bins;
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

// =================================================================================
//    Parse Mesh options
// =================================================================================
void parse_mesh_input(Yaml::Node& root, mesh_input_t& mesh_input)
{
    Yaml::Node& mesh_yaml = root["mesh_options"];

    // get the mesh variables names set by the user
    std::vector<std::string> user_mesh_inputs;


    // extract words from the input file and validate they are correct
    validate_inputs(mesh_yaml, user_mesh_inputs, str_mesh_inps, mesh_required_inps);

    // loop over the words in the material input definition
    for (auto& a_word : user_mesh_inputs) {
        if (VERBOSE) {
            std::cout << a_word << std::endl;
        }

        Yaml::Node& material_inps_yaml = root["mesh_options"][a_word];

        // get mesh source [generate or from file]
        if (a_word.compare("source") == 0) {
            std::string source = root["mesh_options"][a_word].As<std::string>();

            auto map = mesh_input_source_map;

            // set the mesh source
            if (map.find(source) != map.end()) {
                mesh_input.source = map[source];
                if (VERBOSE) {
                    std::cout << "\tsource = " << source << std::endl;
                }
            }
            else{
                std::cout << "ERROR: invalid mesh option input in YAML file: " << source << std::endl;
                std::cout << "Valid options are: " << std::endl;

                for (const auto& pair : map) {
                    std::cout << "\t" << pair.first << std::endl;
                }
                throw std::runtime_error("**** Region Not Understood ****");
            } // end if
        } // source
        // Number of dimensions for the mesh
        else if (a_word.compare("num_dims") == 0) {
            int num_dim = root["mesh_options"][a_word].As<int>();
            if (VERBOSE) {
                std::cout << "\tNum dimensions = " << num_dim << std::endl;
            }

            mesh_input.num_dims = num_dim;
        } // Number of dimensions
        // get mesh type for generation
        else if (a_word.compare("type") == 0) {
            std::string type = root["mesh_options"][a_word].As<std::string>();

            auto map = mesh_input_type_map;

            // set the mesh type
            if (map.find(type) != map.end()) {
                mesh_input.type = map[type];
                if (VERBOSE) {
                    std::cout << "\ttype = " << type << std::endl;
                }
            }
            else{
                std::cout << "ERROR: invalid mesh option input in YAML file: " << type << std::endl;
                std::cout << "Valid options are: " << std::endl;

                for (const auto& pair : map) {
                    std::cout << "\t" << pair.first << std::endl;
                }
                throw std::runtime_error("**** Region Not Understood ****");
            } // end if
        } // type
        // Get mesh file path
        else if (a_word.compare("file_path") == 0) {
            std::string path = root["mesh_options"][a_word].As<std::string>();
            if (VERBOSE) {
                std::cout << "\tfile_path = " << path << std::endl;
            }

            // absolute path to file or local to the director where exe is run
            mesh_input.file_path = path;  

        } // end file path
        // Origin for the mesh
        else if (a_word.compare("origin") == 0) {
            std::string origin = root["mesh_options"][a_word].As<std::string>();
            if (VERBOSE) {
                std::cout << "\torigin = " << origin << std::endl;
            }

            // get the origin numbers, values are words
            std::vector<std::string> numbers = exact_array_values(origin, ",");

            double x1 = std::stod(numbers[0]);
            double y1 = std::stod(numbers[1]);
            double z1;

            if(numbers.size()==3){ 
                // 3D
                z1 = std::stod(numbers[2]);
            }
            else {
                // 2D
                z1 = 0.0;
            } // end if

            if (VERBOSE) {
                std::cout << "\tx1 = " << x1 << std::endl;
                std::cout << "\ty1 = " << y1 << std::endl;
                std::cout << "\tz1 = " << z1 << std::endl;
            }

            // storing the origin values as
            mesh_input.origin[0] = x1;
            mesh_input.origin[1] = y1;
            mesh_input.origin[2] = z1;
        }
        // Extents of the mesh
        else if (a_word.compare("length") == 0) {
            std::string length = root["mesh_options"][a_word].As<std::string>();
            if (VERBOSE) {
                std::cout << "\tlength = " << length << std::endl;
            }

            std::vector<std::string> numbers = exact_array_values(length, ",");

            double l1 = std::stod(numbers[0]);
            double l2 = std::stod(numbers[1]);
            double l3;

            if(numbers.size()==3){ 
                // 3D
                l3 = std::stod(numbers[2]);
            }
            else {
                // 2D
                l3 = 0.0;
            } // end if


            if (VERBOSE) {
                std::cout << "\tl1 = " << l1 << std::endl;
                std::cout << "\tl2 = " << l2 << std::endl;
                std::cout << "\tl3 = " << l3 << std::endl;
            }

            // storing the length values
            mesh_input.length[0] = l1;
            mesh_input.length[1] = l2;
            mesh_input.length[2] = l3;
        }
        // Number of elements per direction
        else if (a_word.compare("num_elems") == 0) {
            std::string elem_count = root["mesh_options"][a_word].As<std::string>();
            if (VERBOSE) {
                std::cout << "\tnum_elems = " << elem_count << std::endl;
            }

            // get the elem_count numbers, values are words
            std::vector<std::string> numbers = exact_array_values(elem_count, ",");

            int n1 = std::stod(numbers[0]);
            int n2 = std::stod(numbers[1]);
            int n3;

            if(numbers.size()==3){ 
                // 3D
                n3 = std::stod(numbers[2]);
            }
            else {
                // 2D
                n3 = 0;
            } // end if


            if (VERBOSE) {
                std::cout << "\tn1 = " << n1 << std::endl;
                std::cout << "\tn2 = " << n2 << std::endl;
                std::cout << "\tn3 = " << n3 << std::endl;
            }

            // storing the number of elements
            mesh_input.num_elems[0] = n1;
            mesh_input.num_elems[1] = n2;
            mesh_input.num_elems[2] = n3;
        }
        // Polynomial order for the mesh
        else if (a_word.compare("polynomial_order") == 0) {
            int p_order = root["mesh_options"][a_word].As<int>();
            if (VERBOSE) {
                std::cout << "\tPoly order = " << p_order << std::endl;
            }

            mesh_input.p_order = p_order;
        } // polynomial order
        // inner radius for 2D RZ meshes
        else if (a_word.compare("inner_radius") == 0) {
            double inner_radius = root["mesh_options"][a_word].As<double>();
            if (VERBOSE) {
                std::cout << "\tInner Radius = " << inner_radius << std::endl;
            }

            mesh_input.inner_radius = inner_radius;
        } // inner radius for 2D RZ meshes
        // outer radius for 2D RZ meshes
        else if (a_word.compare("outer_radius") == 0) {
            double outer_radius = root["mesh_options"][a_word].As<double>();
            if (VERBOSE) {
                std::cout << "\tOuter Radius = " << outer_radius << std::endl;
            }

            mesh_input.outer_radius = outer_radius;
        } // outer radius for 2D RZ meshes
        // starting angle for 2D RZ meshes
        else if (a_word.compare("starting_angle") == 0) {
            double starting_angle = root["mesh_options"][a_word].As<double>();
            if (VERBOSE) {
                std::cout << "\tStarting angle = " << starting_angle << std::endl;
            }

            mesh_input.starting_angle = starting_angle;
        } // starting angle for 2D RZ meshes
        // ending angle for 2D RZ meshes
        else if (a_word.compare("ending_angle") == 0) {
            double ending_angle = root["mesh_options"][a_word].As<double>();
            if (VERBOSE) {
                std::cout << "\tEnding angle = " << ending_angle << std::endl;
            }

            mesh_input.ending_angle = ending_angle;
        } // ending angle for 2D RZ meshes
        // Number of radial elements for 2D RZ meshes
        else if (a_word.compare("num_radial_elems") == 0) {
            int num_radial_elems = root["mesh_options"][a_word].As<int>();
            if (VERBOSE) {
                std::cout << "\tNumber of radial elements = " << num_radial_elems << std::endl;
            }

            mesh_input.num_radial_elems = num_radial_elems;
        } // Number of radial elements for 2D RZ meshes
        // Number of angular elements for 2D RZ meshes
        else if (a_word.compare("num_angular_elems") == 0) {
            int num_angular_elems = root["mesh_options"][a_word].As<int>();
            if (VERBOSE) {
                std::cout << "\tNumber of angular elements = " << num_angular_elems << std::endl;
            }

            mesh_input.num_angular_elems = num_angular_elems;
        } // Number of angular elements for 2D RZ meshes
        else {
            std::cout << "ERROR: invalid input: " << a_word << std::endl;
            std::cout << "Valid options are: " << std::endl;

            for (const auto& element : str_mesh_inps) {
                std::cout << element << std::endl;
            }
            throw std::runtime_error("**** Mesh Not Understood ****");
        }



        // -----------------------------------------------
        // check for consistency in input settings

        if (mesh_input.source == mesh_input::file && mesh_input.file_path.empty()) {
            std::cout << "ERROR: When the mesh source is a file, a file_path must be set to point to the mesh file" << std::endl;
            std::cout << "A mesh can either be generated or read in from a file, but not both" << std::endl;
        }

        if (mesh_input.source == mesh_input::generate) {

            // check to see if a file path was set
            if(mesh_input.file_path.size()>0){
                std::cout << "ERROR: When the mesh source is set to generate, a mesh file cannot be passed in" << std::endl;
                exit(0);
            }

        }

        // -----------------------------------------------


    } // end user_mesh_inputs
} // end of parse mesh options

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
        if (VERBOSE) {
            std::cout << "Word : " << a_word << std::endl;
        }

        // get output format
        if (a_word.compare("output_file_format") == 0) {
            std::string format = root["output_options"][a_word].As<std::string>();

            auto map = output_format_map;

            // set the output format
            if (map.find(format) != map.end()) {
                output_options.format = map[format];
                if (VERBOSE) {
                    std::cout << "\tformat = " << format << std::endl;
                }
            }
            else{
                std::cout << "ERROR: invalid output option input in YAML file: " << format << std::endl;
                std::cout << "Valid options are: " << std::endl;

                for (const auto& pair : map) {
                    std::cout << "\t" << pair.first << std::endl;
                }
            } // end if
        } // output format
        // get timer_output_level
        else if (a_word.compare("timer_output_level") == 0) {
            std::string timer_level = root["output_options"][a_word].As<std::string>();

            auto map = timer_output_level_map;

            // set the timer_output_level
            if (map.find(timer_level) != map.end()) {
                output_options.timer_level = map[timer_level];
                if (VERBOSE) {
                    std::cout << "\ttimer_level = " << timer_level << std::endl;
                }
            }
            else{
                std::cout << "ERROR: invalid output option input in YAML file: " << timer_level << std::endl;
                std::cout << "Valid options are: " << std::endl;

                for (const auto& pair : map) {
                    std::cout << "\t" << pair.first << std::endl;
                }
            } // end if
        } // timer_level
        // Graphics time step
        else if (a_word.compare("graphics_time_step") == 0) {
            real_t graphics_time_step = root["output_options"][a_word].As<real_t>();
            if (VERBOSE) {
                std::cout << "\tgraphics_time_step = " << graphics_time_step << std::endl;
            }

            output_options.graphics_time_step = graphics_time_step;
        } // graphics_time_step
        // Graphics iteration step
        else if (a_word.compare("graphics_iteration_step") == 0) {
            int graphics_iteration_step = root["output_options"][a_word].As<int>();
            if (VERBOSE) {
                std::cout << "\tgraphics_iteration_step = " << graphics_iteration_step << std::endl;
            }

            output_options.graphics_iteration_step = graphics_iteration_step;
        } // graphics_iteration_step
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

// =================================================================================
//    Parse Fill regions
// =================================================================================
void parse_regions(Yaml::Node& root, 
                   CArrayKokkos<RegionFill_t>& region_fills, 
                   CArray<RegionFill_host_t>&  region_fills_host)

{
    Yaml::Node& region_yaml = root["regions"];

    size_t num_regions = region_yaml.Size();

    region_fills = CArrayKokkos<RegionFill_t>(num_regions , "sim_param.region_fills");
    region_fills_host = CArray<RegionFill_host_t>(num_regions); 

    // loop over the fill regions specified
    for (int reg_id = 0; reg_id < num_regions; reg_id++) {
        // read the variables names
        Yaml::Node& inps_yaml = root["regions"][reg_id]["fill_volume"];

        // get the material variables names set by the user
        std::vector<std::string> user_str_region_inps;

        // extract words from the input file and validate they are correct
        validate_inputs(inps_yaml, user_str_region_inps, str_region_inps, region_required_inps);

        // loop over the words in the material input definition
        for (auto& a_word : user_str_region_inps) {
            if (VERBOSE) {
                std::cout << a_word << std::endl;
            }

            Yaml::Node& material_inps_yaml = root["regions"][reg_id]["fill_volume"][a_word];

            // set the values
            if (a_word.compare("material_id") == 0) {
                int id = root["regions"][reg_id]["fill_volume"][a_word].As<int>();

                RUN({
                    region_fills(reg_id).material_id = id;
                });
            } // mat_id
            else if (a_word.compare("den") == 0) {
                double den = root["regions"][reg_id]["fill_volume"]["den"].As<double>();

                // check for a valid density else save it
                if (den < 0.0) {
                    std::cout << "ERROR: density is negative: " << den << std::endl;
                }
                else {
                    RUN({
                        region_fills(reg_id).den = den;  // NOTE: GPUs will require a RUN({})
                    });
                }
            } // den
            else if (a_word.compare("sie") == 0) {
                // specific internal energy

                double sie = root["regions"][reg_id]["fill_volume"]["sie"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tsie = " << sie << std::endl;
                }

                RUN({
                    region_fills(reg_id).sie = sie;
                });
            } // sie
            else if (a_word.compare("ie") == 0) {
                // extensive internal energy

                double ie = root["regions"][reg_id]["fill_volume"]["ie"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tie = " << ie << std::endl;
                }

                RUN({
                    region_fills(reg_id).ie = ie;
                });
            } // ie
            else if (a_word.compare("speed") == 0) {
                double speed = root["regions"][reg_id]["fill_volume"]["speed"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tspeed = " << speed << std::endl;
                }
                RUN({
                    region_fills(reg_id).speed = speed;
                });
            } // speed
            else if (a_word.compare("temperature") == 0) {
                double temperature = root["regions"][reg_id]["fill_volume"]["temperature"].As<double>();
                if (VERBOSE) {
                    std::cout << "\ttemperature = " << temperature << std::endl;
                }
                RUN({
                    region_fills(reg_id).temperature = temperature;
                });
            } // temperature
            else if (a_word.compare("specific_heat") == 0) {
                double specific_heat = root["regions"][reg_id]["fill_volume"]["specific_heat"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tspecific_heat = " << specific_heat << std::endl;
                }
                RUN({
                    region_fills(reg_id).specific_heat = specific_heat;
                });
            } // specific_heat
            else if (a_word.compare("thermal_conductivity") == 0) {
                double thermal_conductivity = root["regions"][reg_id]["fill_volume"]["thermal_conductivity"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tthermal_conductivity = " << thermal_conductivity << std::endl;
                }
                RUN({
                    region_fills(reg_id).thermal_conductivity = thermal_conductivity;
                });
            } // thermal_conductivity
            else if (a_word.compare("u") == 0) {
                // x-component of velocity
                double u = root["regions"][reg_id]["fill_volume"]["u"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tu = " << u << std::endl;
                }
                RUN({
                    region_fills(reg_id).u = u;
                });
            } // u
            else if (a_word.compare("v") == 0) {
                // y-component of velocity
                double v = root["regions"][reg_id]["fill_volume"]["v"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tv = " << v << std::endl;
                }

                RUN({
                    region_fills(reg_id).v = v;
                });
            } // v
            else if (a_word.compare("w") == 0) {
                // z-component of velocity

                double w = root["regions"][reg_id]["fill_volume"]["w"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tw = " << w << std::endl;
                }

                RUN({
                    region_fills(reg_id).w = w;
                });
            } // w
            else if (a_word.compare("radius1") == 0) {
                // inner radius of sphere/cylinder

                double radius1 = root["regions"][reg_id]["fill_volume"]["radius1"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tradius1 = " << radius1 << std::endl;
                }

                RUN({
                    region_fills(reg_id).radius1 = radius1;
                });
            } // radius1
            else if (a_word.compare("radius2") == 0) {
                // outer radius of sphere/cylinder

                double radius2 = root["regions"][reg_id]["fill_volume"]["radius2"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tradius2 = " << radius2 << std::endl;
                }

                RUN({
                    region_fills(reg_id).radius2 = radius2;
                });
            } // radius2
            else if (a_word.compare("x1") == 0) {
                // inner plane

                double x1 = root["regions"][reg_id]["fill_volume"]["x1"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tx1 = " << x1 << std::endl;
                }

                RUN({
                    region_fills(reg_id).x1 = x1;
                });
            } // x1
            else if (a_word.compare("x2") == 0) {
                // outer plane

                double x2 = root["regions"][reg_id]["fill_volume"]["x2"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tx2 = " << x2 << std::endl;
                }

                RUN({
                    region_fills(reg_id).x2 = x2;
                });
            } // x2
            else if (a_word.compare("y1") == 0) {
                // inner plane

                double y1 = root["regions"][reg_id]["fill_volume"]["y1"].As<double>();
                if (VERBOSE) {
                    std::cout << "\ty1 = " << y1 << std::endl;
                }

                RUN({
                    region_fills(reg_id).y1 = y1;
                });
            } // y1
            else if (a_word.compare("y2") == 0) {
                // outer plane

                double y2 = root["regions"][reg_id]["fill_volume"]["y2"].As<double>();
                if (VERBOSE) {
                    std::cout << "\ty2 = " << y2 << std::endl;
                }

                RUN({
                    region_fills(reg_id).y2 = y2;
                });
            } // y2
            else if (a_word.compare("z1") == 0) {
                // inner plane

                double z1 = root["regions"][reg_id]["fill_volume"]["z1"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tz1 = " << z1 << std::endl;
                }

                RUN({
                    region_fills(reg_id).z1 = z1;
                });
            } // z1
            else if (a_word.compare("z2") == 0) {
                // outer plane

                double z2 = root["regions"][reg_id]["fill_volume"]["z2"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tz2 = " << z2 << std::endl;
                }

                RUN({
                    region_fills(reg_id).z2 = z2;
                });
            } // z2
            else if (a_word.compare("scale_x") == 0) {
                // outer plane

                double scale_x = root["regions"][reg_id]["fill_volume"]["scale_x"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tscale_x = " << scale_x << std::endl;
                }

                region_fills_host(reg_id).scale_x = scale_x;

            } // scale_x
            else if (a_word.compare("scale_y") == 0) {
                // outer plane

                double scale_y = root["regions"][reg_id]["fill_volume"]["scale_y"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tscale_y = " << scale_y << std::endl;
                }

                region_fills_host(reg_id).scale_y = scale_y;

            } // scale_y
            else if (a_word.compare("scale_z") == 0) {
                // outer plane

                double scale_z = root["regions"][reg_id]["fill_volume"]["scale_z"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tscale_z = " << scale_z << std::endl;
                }

                region_fills_host(reg_id).scale_z = scale_z;

            } // scale_z
            else if (a_word.compare("velocity") == 0) {

                // velocity fill region type
                std::string type = root["regions"][reg_id]["fill_volume"]["velocity"].As<std::string>();

                if (VERBOSE) {
                    std::cout << "\tvelocity = " << type << std::endl;
                }
                // set the volume tag type
                if (velocity_type_map.find(type) != velocity_type_map.end()) {
                 
                    // velocity_type_map[type] returns enum value, e.g., init_conds::velocity 
                    switch(velocity_type_map[type]){

                        case init_conds::cartesian:
                            std::cout << "Setting velocity initial conditions type to cartesian " << std::endl;
                            RUN({
                                region_fills(reg_id).velocity = init_conds::cartesian;
                            });
                            break;

                         case init_conds::radial:
                            std::cout << "Setting velocity initial conditions type to radial " << std::endl;
                            RUN({
                                region_fills(reg_id).velocity = init_conds::radial;
                            });
                            break;

                         case init_conds::spherical:
                            std::cout << "Setting velocity initial conditions type to spherical " << std::endl;
                            RUN({
                                region_fills(reg_id).velocity = init_conds::spherical;
                            });
                            break;

                         case init_conds::radial_linear:
                            std::cout << "Setting velocity initial conditions type to radial_linear " << std::endl;
                            RUN({
                                region_fills(reg_id).velocity = init_conds::radial_linear;
                            });
                            break;

                         case init_conds::spherical_linear:
                            std::cout << "Setting velocity initial conditions type to spherical_linear " << std::endl;
                            RUN({
                                region_fills(reg_id).velocity = init_conds::spherical_linear;
                            });
                            break;

                         case init_conds::tg_vortex:
                            std::cout << "Setting velocity initial conditions type to tg_vortex " << std::endl;
                            RUN({
                                region_fills(reg_id).velocity = init_conds::tg_vortex;
                            });
                            break;

                         case init_conds::no_ic_vel:
                            std::cout << "Setting velocity initial conditions type to no velocity" << std::endl;
                            RUN({ 
                                region_fills(reg_id).velocity = init_conds::no_ic_vel;
                            });
                            break;

                        default:

                            RUN({ 
                                region_fills(reg_id).velocity = init_conds::no_ic_vel;
                            });

                            std::cout << "ERROR: No valid velocity intial conditions type input " << std::endl;
                            std::cout << "Valid IC types are: " << std::endl;
                            
                            for (const auto& pair : velocity_type_map) {
                                std::cout << pair.second << std::endl;
                            }

                            throw std::runtime_error("**** Velocity Initial Conditions Type Not Understood ****");
                            break;
                    } // end switch

                    if (VERBOSE) {
                        std::cout << "\tvolume_fill = " << type << std::endl;
                    } // end if

                }
                else{
                    std::cout << "ERROR: invalid input: " << type << std::endl;
                    throw std::runtime_error("**** Velocity IC Not Understood ****");
                } // end if
            } // end velocity fill type

            else if (a_word.compare("temperature_distribution") == 0) {

                // temperature_distribution fill region type
                std::string type = root["regions"][reg_id]["fill_volume"]["temperature_distribution"].As<std::string>();

                if (VERBOSE) {
                    std::cout << "\ttemperature_distribution = " << type << std::endl;
                }
                // set the volume tag type NOTE: rename to remove reference to velocity, change to distribution
                if (velocity_type_map.find(type) != velocity_type_map.end()) {
                 
                    // velocity_type_map[type] returns enum value, e.g., init_conds::velocity 
                    switch(velocity_type_map[type]){

                        case init_conds::cartesian:
                            std::cout << "Setting temperature_distribution initial conditions type to cartesian " << std::endl;
                            RUN({
                                region_fills(reg_id).temp_distribution = init_conds::cartesian;
                            });
                            break;

                         case init_conds::radial:
                            std::cout << "Setting temperature_distribution initial conditions type to radial " << std::endl;
                            RUN({
                                region_fills(reg_id).temp_distribution = init_conds::radial;
                            });
                            break;

                         case init_conds::spherical:
                            std::cout << "Setting temperature_distribution initial conditions type to spherical " << std::endl;
                            RUN({
                                region_fills(reg_id).temp_distribution = init_conds::spherical;
                            });
                            break;

                         case init_conds::radial_linear:
                            std::cout << "Setting temperature_distribution initial conditions type to radial_linear " << std::endl;
                            RUN({
                                region_fills(reg_id).temp_distribution = init_conds::radial_linear;
                            });
                            break;

                         case init_conds::spherical_linear:
                            std::cout << "Setting temperature_distribution initial conditions type to spherical_linear " << std::endl;
                            RUN({
                                region_fills(reg_id).temp_distribution = init_conds::spherical_linear;
                            });
                            break;

                         case init_conds::tg_vortex:
                            std::cout << "Setting temperature_distribution initial conditions type to tg_vortex " << std::endl;
                            RUN({
                                region_fills(reg_id).temp_distribution = init_conds::tg_vortex;
                            });
                            break;

                         case init_conds::no_ic_vel:
                            std::cout << "Setting temperature_distribution initial conditions type to no temperature" << std::endl;
                            RUN({ 
                                region_fills(reg_id).temp_distribution = init_conds::no_ic_vel;
                            });
                            break;

                        default:

                            RUN({ 
                                region_fills(reg_id).temp_distribution = init_conds::no_ic_vel;
                            });

                            std::cout << "ERROR: No valid temperature_distribution intial conditions type input " << std::endl;
                            std::cout << "Valid IC types are: " << std::endl;
                            
                            for (const auto& pair : velocity_type_map) {
                                std::cout << pair.second << std::endl;
                            }

                            throw std::runtime_error("**** Temperature Initial Conditions Type Not Understood ****");
                            break;
                    } // end switch

                    if (VERBOSE) {
                        std::cout << "\tvolume_fill = " << type << std::endl;
                    } // end if

                }
                else{
                    std::cout << "ERROR: invalid input: " << type << std::endl;
                    throw std::runtime_error("**** Temperature IC Not Understood ****");
                } // end if
            } // end velocity fill type

            //
            else if (a_word.compare("type") == 0) {

                // region volume fill type
                std::string type = root["regions"][reg_id]["fill_volume"]["type"].As<std::string>();

                if (VERBOSE) {
                    std::cout << "\ttype = " << type << std::endl;
                }

                // set the velocity tag type
                if (region_type_map.find(type) != region_type_map.end()) {
                 
                    // region_type_map[type] returns enum value, e.g., init_conds::velocity 
                    switch(region_type_map[type]){

                        case region::global:
                            std::cout << "Setting volume fill type to global " << std::endl;
                            RUN({
                                region_fills(reg_id).volume = region::global;
                            });
                            break;

                        case region::box:
                            std::cout << "Setting volume fill type to box " << std::endl;
                            RUN({
                                region_fills(reg_id).volume = region::box;
                            });
                            break;

                        case region::cylinder:
                            std::cout << "Setting volume fill type to cylinder " << std::endl;
                            RUN({
                                region_fills(reg_id).volume = region::cylinder;
                            });
                            break;

                        case region::sphere:
                            std::cout << "Setting volume fill type to sphere " << std::endl;
                            RUN({
                                region_fills(reg_id).volume = region::sphere;
                            });
                            break;

                        case region::readVoxelFile:
                            std::cout << "Setting volume fill type to readVoxelFile " << std::endl;
                            RUN({
                                region_fills(reg_id).volume = region::readVoxelFile;
                            });
                            break;

                        case region::no_volume:
                            std::cout << "Setting volume fill type to none " << std::endl;
                            RUN({
                                region_fills(reg_id).volume = region::no_volume;
                            });
                            break;
                        default:

                            RUN({ 
                                region_fills(reg_id).volume = region::no_volume;
                            });

                            std::cout << "ERROR: No valid region volume fill type input " << std::endl;
                            std::cout << "Valid IC volume fill types are: " << std::endl;
                            
                            for (const auto& pair : region_type_map) {
                                std::cout << pair.second << std::endl;
                            }

                            throw std::runtime_error("**** Region Volume Fill Type Not Understood ****");
                            break;
                    } // end switch
                }
                 
                else{
                    std::cout << "ERROR: invalid input: " << type << std::endl;
                    throw std::runtime_error("**** Volume Fill Not Understood ****");
                } // end if

            } // end volume fill type
            // Get mesh file path
            else if (a_word.compare("file_path") == 0) {
                // region volume fill type
                std::string path = root["regions"][reg_id]["fill_volume"]["file_path"].As<std::string>();

                if (VERBOSE) {
                    std::cout << "\tfile_path = " << path << std::endl;
                }

                // absolute path to file or local to the director where exe is run
                region_fills_host(reg_id).file_path = path;   // saving the absolute file path
   

            } // end file path
            //
            else if (a_word.compare("origin") == 0) {
                std::string origin = root["regions"][reg_id]["fill_volume"]["origin"].As<std::string>();
                if (VERBOSE) {
                    std::cout << "\torigin = " << origin << std::endl;
                }

                // get the origin numbers, values are words
                std::vector<std::string> numbers = exact_array_values(origin, ",");

                double x1 = std::stod(numbers[0]);
                double y1 = std::stod(numbers[1]);
                double z1;

                if(numbers.size()==3){ 
                    // 3D
                    z1 = std::stod(numbers[2]);
                }
                else {
                    // 2D
                    z1 = 0.0;
                } //

                if (VERBOSE) {
                    std::cout << "\tx1 = " << x1 << std::endl;
                    std::cout << "\ty1 = " << y1 << std::endl;
                    std::cout << "\tz1 = " << z1 << std::endl;
                }

                // storing the origin values as (x1,y1,z1)
                RUN({
                    region_fills(reg_id).origin[0] = x1;
                    region_fills(reg_id).origin[1] = y1;
                    region_fills(reg_id).origin[2] = z1;
                });
            } // origin
            else {
                std::cout << "ERROR: invalid input: " << a_word << std::endl;
                std::cout << "Valid options are: " << std::endl;
                for (const auto& element : str_region_inps) {
                    std::cout << element << std::endl;
                }
                throw std::runtime_error("**** Region Not Understood ****");
            }
        } // end for words in material

        // -----------------------------------------------
        // check for consistency in input settings

        // check to see if a file path is empty
        if(region_fills_host(reg_id).file_path.empty()){

            RUN({
                // if the following is true, stop simulation; must add all mesh read options
                if (region_fills(reg_id).volume == region::readVoxelFile) {
                    Kokkos::abort("\n********************************************************************************************\n"
                                    "ERROR: \n"
                                    "When using a file to initialize a region, a file_path must be set to point to the mesh file\n"
                                    "********************************************************************************************\n");
                }
            });
        } // end if check

        // check to see if a file path was set
        if(region_fills_host(reg_id).file_path.size()>0){
            RUN({
                if (region_fills(reg_id).volume != region::readVoxelFile){  
                    // this means it is a geometric definition of the region
                    Kokkos::abort("\n********************************************************************************************\n"
                                    "ERROR: \n"
                                    "When a geometric entity defines the region, a mesh file cannot be passed to set the region\n"
                                    "********************************************************************************************\n");
                }
            });

        }  // end if        

        // -----------------------------------------------


    } // end loop over regions
} // end of function to parse region

// =================================================================================
//    Parse Material Definitions
// =================================================================================
void parse_materials(Yaml::Node& root, Material_t& Materials)
{
    Yaml::Node& material_yaml = root["materials"];

    size_t num_materials = material_yaml.Size();
    std::cout << "Number of materials =  "<< num_materials << std::endl;

    Materials.num_mats = num_materials;

    // --- allocate memory for arrays inside material struct ---

    // setup
    Materials.MaterialSetup = DCArrayKokkos<MaterialSetup_t>(num_materials, "material_setup");

    // function pointers to material models
    Materials.MaterialFunctions = CArrayKokkos<MaterialFunctions_t>(num_materials, "material_functions");

    // enums
    Materials.MaterialEnums = DCArrayKokkos<MaterialEnums_t>(num_materials, "material_enums");

    // these are temp arrays to store global variables given in the yaml input file for each material, 100 vars is the max allowable
    DCArrayKokkos<double> tempGlobalEOSVars(num_materials, 100, "temp_array_eos_vars");
    DCArrayKokkos<double> tempGlobalStrengthVars(num_materials, 100, "temp_array_strength_vars");

    Materials.num_eos_global_vars      =  CArrayKokkos <size_t> (num_materials, "num_eos_global_vars");
    Materials.num_strength_global_vars =  CArrayKokkos <size_t> (num_materials, "num_strength_global_vars");


    // initialize the num of global vars to 0 for all models
    FOR_ALL(mat_id, 0, num_materials, {
        Materials.num_eos_global_vars(mat_id) = 0;
        Materials.num_strength_global_vars(mat_id) = 0;
    });

    // loop over the materials specified
    for (int mat_id = 0; mat_id < num_materials; mat_id++) {
        // read the variables names
        Yaml::Node& inps_yaml = root["materials"][mat_id]["material"];

        size_t num_vars_set = inps_yaml.Size();

        std::cout << "Number of vars set =  "<< num_vars_set << std::endl;

        // get the material variables names set by the user
        std::vector<std::string> user_str_material_inps;

        // extract words from the input file and validate they are correct
        validate_inputs(inps_yaml, user_str_material_inps, str_material_inps, material_hydrodynamics_required_inps);

        // loop over the words in the material input definition
        for (auto& a_word : user_str_material_inps) {
            if (VERBOSE) {
                std::cout << a_word << std::endl;
            }

            Yaml::Node& material_inps_yaml = root["materials"][mat_id]["material"][a_word];

            // set the values in the input for this word

            // set the values
            if (a_word.compare("q1") == 0) {
                double q1 = root["materials"][mat_id]["material"]["q1"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tq1 = " << q1 << std::endl;
                }
                RUN({
                    Materials.MaterialFunctions(mat_id).q1 = q1;
                });
            } // q1
            else if (a_word.compare("q1ex") == 0) {
                double q1ex = root["materials"][mat_id]["material"]["q1ex"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tq1ex = " << q1ex << std::endl;
                }
                RUN({
                    Materials.MaterialFunctions(mat_id).q1ex = q1ex;
                });
            } // q1ex
            else if (a_word.compare("q2") == 0) {
                // outer plane

                double q2 = root["materials"][mat_id]["material"]["q2"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tq2 = " << q2 << std::endl;
                }

                RUN({
                    Materials.MaterialFunctions(mat_id).q2 = q2;
                });
            } // q2
            else if (a_word.compare("q2ex") == 0) {
                // outer plane

                double q2ex = root["materials"][mat_id]["material"]["q2ex"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tq2ex = " << q2ex << std::endl;
                }
                RUN({
                    Materials.MaterialFunctions(mat_id).q2ex = q2ex;
                });
            } // q1ex
            else if (a_word.compare("id") == 0) {
                int m_id = root["materials"][mat_id]["material"]["id"].As<int>();
                if (VERBOSE) {
                    std::cout << "\tid = " << m_id << std::endl;
                }
                RUN({
                    Materials.MaterialFunctions(mat_id).id = m_id;
                });
            } // id
            //extract eos model
            else if (a_word.compare("eos_model_type") == 0) {
                std::string type = root["materials"][mat_id]["material"]["eos_model_type"].As<std::string>();

                // set the eos type
                if (eos_type_map.find(type) != eos_type_map.end()) {

                    // eos_type_map[type] returns enum value, e.g., model::decoupled
                    switch(eos_type_map[type]){
                        case model::decoupledEOSType:
                            std::cout << "Setting EOS type to decoupled " << std::endl;
                            RUN({
                                Materials.MaterialEnums(mat_id).EOSType = model::decoupledEOSType;
                            });

                            Materials.MaterialEnums.host(mat_id).EOSType = model::decoupledEOSType;
                            break;

                        case model::coupledEOSType:
                            std::cout << "Setting EOS type to coupled " << std::endl;
                            RUN({
                                Materials.MaterialEnums(mat_id).EOSType = model::coupledEOSType;
                            });
                            Materials.MaterialEnums.host(mat_id).EOSType = model::coupledEOSType;
                            break;

                        default:
                            RUN({ 
                                Materials.MaterialEnums(mat_id).EOSType = model::noEOSType;
                            });
                            Materials.MaterialEnums.host(mat_id).EOSType = model::noEOSType;
                            std::cout << "ERROR: No valid EOS type input " << std::endl;
                            std::cout << "Valid EOS types are: " << std::endl;
                            
                            for (const auto& pair : eos_type_map) {
                                std::cout << pair.second << std::endl;
                            }

                            throw std::runtime_error("**** EOS Type Not Understood ****");
                            break;
                    } // end switch

                    if (VERBOSE) {
                        std::cout << "\teos = " << type << std::endl;
                    }
                } 
                else{
                    std::cout << "ERROR: invalid eos type input: " << type << std::endl;
                } // end if
            }


            else if (a_word.compare("eos_model") == 0) {
                std::string eos = root["materials"][mat_id]["material"]["eos_model"].As<std::string>();

                // set the EOS
                if (eos_models_map.find(eos) != eos_models_map.end()) {
                    
                    switch(eos_models_map[eos]){

                        case model::noEOS:
                            RUN({
                                Materials.MaterialFunctions(mat_id).calc_pressure    = &NoEOSModel::calc_pressure;
                                Materials.MaterialFunctions(mat_id).calc_sound_speed = &NoEOSModel::calc_sound_speed;
                            });
                            if (VERBOSE) {
                                std::cout << "\teos_model = " << eos << std::endl;
                            }
                            break;

                        case model::gammaLawGasEOS:
                            RUN({
                                Materials.MaterialFunctions(mat_id).calc_pressure    = &GammaLawGasEOSModel::calc_pressure;
                                Materials.MaterialFunctions(mat_id).calc_sound_speed = &GammaLawGasEOSModel::calc_sound_speed;
                            });
                            if (VERBOSE) {
                                std::cout << "\teos_model = " << eos << std::endl;
                            }
                            break;

                        case model::voidEOS:
                            RUN({
                                Materials.MaterialFunctions(mat_id).calc_pressure    = &VoidEOSModel::calc_pressure;
                                Materials.MaterialFunctions(mat_id).calc_sound_speed = &VoidEOSModel::calc_sound_speed;
                            });
                            if (VERBOSE) {
                                std::cout << "\teos_model = " << eos << std::endl;
                            }
                            break;

                        case model::userDefinedEOS:
                            RUN({
                                Materials.MaterialFunctions(mat_id).calc_pressure    = &UserDefinedEOSModel::calc_pressure;
                                Materials.MaterialFunctions(mat_id).calc_sound_speed = &UserDefinedEOSModel::calc_sound_speed;
                            });
                            if (VERBOSE) {
                                std::cout << "\teos_model = " << eos << std::endl;
                            }
                            break;
                        default:
                            std::cout << "ERROR: invalid input: " << eos << std::endl;
                            throw std::runtime_error("**** EOS Not Understood ****");
                            break;
                    } // end switch on EOS type
                }
                else{
                    std::cout << "ERROR: invalid EOS input: " << eos << std::endl;
                    throw std::runtime_error("**** EOS Not Understood ****");
                } // end if
            } // EOS model

            // Type of strength model
            else if (a_word.compare("strength_model_type") == 0) {
                std::string strength_model_type = root["materials"][mat_id]["material"]["strength_model_type"].As<std::string>();

                // set the EOS
                if (strength_type_map.find(strength_model_type) != strength_type_map.end()) {
                    
                    switch(strength_type_map[strength_model_type]){

                        case model::noStrengthType:
                            RUN({
                                Materials.MaterialEnums(mat_id).StrengthType = model::noStrengthType;
                            });
                            Materials.MaterialEnums.host(mat_id).StrengthType = model::noStrengthType;
                            if (VERBOSE) {
                                std::cout << "\tstrength_model_type_type = " << strength_model_type << std::endl;
                            }
                            break;

                        case model::incrementBased:
                            RUN({
                                Materials.MaterialEnums(mat_id).StrengthType = model::incrementBased;
                            });
                            Materials.MaterialEnums.host(mat_id).StrengthType = model::incrementBased;
                            
                            if (VERBOSE) {
                                std::cout << "\tstrength_model_type = " << strength_model_type << std::endl;
                            }
                            break;
                        case model::stateBased:
                            RUN({
                                Materials.MaterialEnums(mat_id).StrengthType = model::stateBased;
                            });
                            Materials.MaterialEnums.host(mat_id).StrengthType = model::stateBased;
                            
                            std::cout << "ERROR: state_based models not yet defined: " << std::endl;
                            throw std::runtime_error("**** ERROR: state_based models not yet defined ****");
                            if (VERBOSE) {
                                std::cout << "\tstrength_model_type = " << strength_model_type << std::endl;
                            }
                            break;
                        default:
                            std::cout << "ERROR: invalid strength type input: " << strength_model_type << std::endl;
                            throw std::runtime_error("**** Strength Model Type Not Understood ****");
                            break;
                    } // end switch on EOS type
                }
                else{
                    std::cout << "ERROR: Invalid strength model type input: " << strength_model_type << std::endl;
                    throw std::runtime_error("**** Strength model type not understood ****");
                } // end if
            } // Strength model type
            
            // Set specific strength model
            else if (a_word.compare("strength_model") == 0) {
                std::string strength_model = root["materials"][mat_id]["material"]["strength_model"].As<std::string>();

                // set the EOS
                if (strength_models_map.find(strength_model) != strength_models_map.end()) {
                    
                    switch(strength_models_map[strength_model]){

                        case model::noStrengthModel:
                            RUN({
                                Materials.MaterialFunctions(mat_id).calc_stress = &NoStrengthModel::calc_stress;
                            });
                            if (VERBOSE) {
                                std::cout << "\tstrength_model = " << strength_model << std::endl;
                            }
                            break;

                        case model::userDefinedStrength:
                            
                            RUN({
                                Materials.MaterialFunctions(mat_id).calc_stress = &UserDefinedStrengthModel::calc_stress;
                            });
                            if (VERBOSE) {
                                std::cout << "\tstrength_model = " << strength_model << std::endl;
                            }
                            break;
                        default:
                            std::cout << "ERROR: invalid strength input: " << strength_model << std::endl;
                            throw std::runtime_error("**** Strength model Not Understood ****");
                            break;
                    } // end switch on EOS type
                }
                else{
                    std::cout << "ERROR: invalid Strength model input: " << strength_model << std::endl;
                    throw std::runtime_error("**** Strength model Not Understood ****");
                } // end if
            } // Strength model
            
            //extract erosion model
            else if (a_word.compare("erosion_model") == 0) {
                std::string erosion_model = root["materials"][mat_id]["material"]["erosion_model"].As<std::string>();

                // set the erosion model
                if (erosion_model_map.find(erosion_model) != erosion_model_map.end()) {

                    // erosion_model_map[erosion_model] returns enum value, e.g., model::erosion
                    switch(erosion_model_map[erosion_model]){
                        case model::basicErosion:
                            Materials.MaterialEnums.host(mat_id).ErosionModels = model::basicErosion;
                            RUN({
                                Materials.MaterialEnums(mat_id).ErosionModels = model::basicErosion;
                                Materials.MaterialFunctions(mat_id).erode = &BasicErosionModel::erode;
                            });
                            break;
                        case model::noErosion:
                            Materials.MaterialEnums.host(mat_id).ErosionModels = model::noErosion;
                            RUN({
                                Materials.MaterialEnums(mat_id).ErosionModels = model::noErosion;
                                Materials.MaterialFunctions(mat_id).erode = &NoErosionModel::erode;
                            });
                            break;
                        default:
                            std::cout << "ERROR: invalid erosion input: " << erosion_model << std::endl;
                            throw std::runtime_error("**** Erosion model Not Understood ****");
                            break;
                    } // end switch

                    if (VERBOSE) {
                        std::cout << "\terosion = " << erosion_model << std::endl;
                    }

                } 
                else{
                    std::cout << "ERROR: invalid erosion type input: " << erosion_model<< std::endl;
                    throw std::runtime_error("**** Erosion model Not Understood ****");
                    break;
                } // end if

            } // erosion model variables
            else if (a_word.compare("erode_tension_val") == 0) {
                double erode_tension_val = root["materials"][mat_id]["material"]["erode_tension_val"].As<double>();
                if (VERBOSE) {
                    std::cout << "\terode_tension_val = " << erode_tension_val << std::endl;
                }
                RUN({
                    Materials.MaterialFunctions(mat_id).erode_tension_val = erode_tension_val;
                });
            } // erode_tension_val
            else if (a_word.compare("erode_density_val") == 0) {
                double erode_density_val = root["materials"][mat_id]["material"]["erode_density_val"].As<double>();
                if (VERBOSE) {
                    std::cout << "\terode_density_val = " << erode_density_val << std::endl;
                }
                RUN({
                    Materials.MaterialFunctions(mat_id).erode_density_val = erode_density_val;
                });
            } // erode_density_val
            
            // exact the eos_global_vars
            else if (a_word.compare("eos_global_vars") == 0) {
                Yaml::Node & mat_global_vars_yaml = root["materials"][mat_id]["material"][a_word];

                size_t num_global_vars = mat_global_vars_yaml.Size();
                
                RUN({ 
                    Materials.num_eos_global_vars(mat_id) = num_global_vars;
                });

                if(num_global_vars>100){
                    throw std::runtime_error("**** Per material, the code only supports up to 100 global vars in the input file ****");
                } // end check on num_global_vars
               
                if (VERBOSE) {
                    std::cout << "num global eos vars = " << num_global_vars << std::endl;
                }

                // store the global eos model parameters
                for (int global_var_id = 0; global_var_id < num_global_vars; global_var_id++) {
                    double eos_var = root["materials"][mat_id]["material"]["eos_global_vars"][global_var_id].As<double>();
                    

                    RUN({
                        tempGlobalEOSVars(mat_id, global_var_id) = eos_var;
                    });

                    if (VERBOSE) {
                        std::cout << "\t var = " << eos_var << std::endl;
                    }
                } // end loop over global vars
            } // "eos_global_vars"
            
            // exact the strength_global_vars
            else if (a_word.compare("strength_global_vars") == 0) {
                Yaml::Node & mat_global_vars_yaml = root["materials"][mat_id]["material"][a_word];

                size_t num_global_vars = mat_global_vars_yaml.Size();
                
                RUN({ 
                    Materials.num_strength_global_vars(mat_id) = num_global_vars;
                });


                if (VERBOSE) {
                    std::cout << "num global strength vars = " << num_global_vars << std::endl;
                }

                // store the global eos model parameters
                for (int global_var_id = 0; global_var_id < num_global_vars; global_var_id++) {
                    double strength_var = root["materials"][mat_id]["material"]["strength_global_vars"][global_var_id].As<double>();
                    
                    RUN({
                        tempGlobalStrengthVars(mat_id,global_var_id) = strength_var;
                    });

                    if (VERBOSE) {
                        std::cout << "\t var = " << strength_var << std::endl;
                    }
                } // end loop over global vars
            } // "eos_global_vars"
            
            else {
                std::cout << "ERROR: invalid input: " << a_word << std::endl;
                std::cout << "Valid options are: " << std::endl;
            
                for (const auto& element : str_material_inps) {
                    std::cout << element << std::endl;
                }
                throw std::runtime_error("**** Material Input Not Understood ****");
            }
        } // end for words in material
    } // end loop over materials

    // allocate ragged rigght memory to hold the model global variables
    Materials.eos_global_vars = RaggedRightArrayKokkos <double> (Materials.num_eos_global_vars, "Materials.eos_global_vars");
    Materials.strength_global_vars = RaggedRightArrayKokkos <double> (Materials.num_strength_global_vars, "Materials.strength_global_vars");



    // save the global variables
    FOR_ALL(mat_id, 0, num_materials, {
        
        for (size_t var_lid=0; var_lid<Materials.num_eos_global_vars(mat_id); var_lid++){
            Materials.eos_global_vars(mat_id, var_lid) = tempGlobalEOSVars(mat_id, var_lid);
        } // end for eos var_lid

        for (size_t var_lid=0; var_lid<Materials.num_strength_global_vars(mat_id); var_lid++){
            Materials.strength_global_vars(mat_id, var_lid) = tempGlobalStrengthVars(mat_id, var_lid);
        } // end for strength var_lid

    }); // end for loop over materials

} // end of function to parse material information



// =================================================================================
//    Parse Boundary Conditions
// =================================================================================
void parse_bcs(Yaml::Node& root, BoundaryCondition_t& BoundaryConditions)
{
    Yaml::Node& bc_yaml = root["boundary_conditions"];

    size_t num_bcs = bc_yaml.Size();

    BoundaryConditions.num_bcs = num_bcs;

    BoundaryConditions.BoundaryConditionSetup = CArrayKokkos <BoundaryConditionSetup_t>(num_bcs, "bc_setup_vars");

    // device functions
    BoundaryConditions.BoundaryConditionFunctions = CArrayKokkos <BoundaryConditionFunctions_t>(num_bcs, "bc_fcns");

    // enums to select options with boundary conditions
    BoundaryConditions.BoundaryConditionEnums  = DCArrayKokkos<BoundaryConditionEnums_t> (num_bcs,"bc_enums");  

    // the state for boundary conditions
    BoundaryConditions.bc_global_vars = DCArrayKokkos<double>(num_bcs, 5, "bc_global_values"); // increase 4 for more params
    BoundaryConditions.bc_state_vars  = DCArrayKokkos<double>(num_bcs, 4, "bc_state_values");  // WARNING a place holder


    // loop over the BC specified
    for (int bc_id = 0; bc_id < num_bcs; bc_id++) {
        // read the variables names
        Yaml::Node& inps_yaml = bc_yaml[bc_id]["boundary_condition"];

        // get the material variables names set by the user
        std::vector<std::string> user_str_bc_inps;

        // extract words from the input file and validate they are correct
        validate_inputs(inps_yaml, user_str_bc_inps, str_bc_inps, bc_required_inps);

        // loop over the words in the material input definition
        for (auto& a_word : user_str_bc_inps) {
            if (VERBOSE) {
                std::cout << a_word << std::endl;
            }

            Yaml::Node& inps_yaml = bc_yaml[bc_id]["boundary_condition"][a_word];

            // get solver for this boundary condition
            if (a_word.compare("solver") == 0) {
                std::string solver = bc_yaml[bc_id]["boundary_condition"][a_word].As<std::string>();

                auto map = solver_map;

                // set the solver
                if (map.find(solver) != map.end()) {
                    solver_input::method bc_solver = map[solver];

                    RUN({
                        BoundaryConditions.BoundaryConditionEnums(bc_id).solver = bc_solver;
                    });

                    if (VERBOSE) {
                        std::cout << "\tsolver = " << solver << std::endl;
                    }
                }
                else{
                    std::cout << "ERROR: invalid boundary condition option input in YAML file: " << solver << std::endl;
                    std::cout << "Valid options are: " << std::endl;

                    for (const auto& pair : map) {
                        std::cout << "\t" << pair.first << std::endl;
                    }
                } // end if
            } // solver
            // get boundary condition type
            else if (a_word.compare("type") == 0) {
                std::string type = bc_yaml[bc_id]["boundary_condition"][a_word].As<std::string>();

                auto map = bc_type_map;

                // set the type
                if (map.find(type) != map.end()) {
                    auto bc_type = map[type];

                    // bc_type_map[type] returns enum value, e.g., boundary_conditions::velocity_constant
                    switch(map[type]){

                        case boundary_conditions::constantVelocityBC :
                            std::cout << "Setting velocity bc " << std::endl;
                            
                            RUN({
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCHydroType = boundary_conditions::constantVelocityBC ;
                                BoundaryConditions.BoundaryConditionFunctions(bc_id).velocity = &ConstantVelocityBC::velocity;
                            });
                            break;

                        case boundary_conditions::timeVaringVelocityBC:
                            std::cout << "Setting velocity bc " << std::endl;
                            
                            RUN({
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCHydroType = boundary_conditions::timeVaringVelocityBC;
                                BoundaryConditions.BoundaryConditionFunctions(bc_id).velocity = &TimeVaryingVelocityBC::velocity;
                            });
                            break;
                        
                        case boundary_conditions::reflectedVelocityBC:
                            std::cout << "Setting velocity bc " << std::endl;
                            
                            RUN({
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCHydroType = boundary_conditions::reflectedVelocityBC;
                                BoundaryConditions.BoundaryConditionFunctions(bc_id).velocity = &ReflectedVelocityBC::velocity;
                            });
                            break;

                        case boundary_conditions::zeroVelocityBC:
                            std::cout << "Setting velocity bc " << std::endl;
                            
                            RUN({
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCHydroType = boundary_conditions::zeroVelocityBC;
                                BoundaryConditions.BoundaryConditionFunctions(bc_id).velocity = &ZeroVelocityBC::velocity;
                            });
                            break;
                        case boundary_conditions::userDefinedVelocityBC:
                            std::cout << "Setting velocity bc " << std::endl;
                            
                            RUN({
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCHydroType = boundary_conditions::userDefinedVelocityBC;
                                BoundaryConditions.BoundaryConditionFunctions(bc_id).velocity = &UserDefinedVelocityBC::velocity;
                            });
                            break;
                        case boundary_conditions::pistonVelocityBC:
                            std::cout << "Setting velocity bc " << std::endl;
                            
                            RUN({
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCHydroType = boundary_conditions::pistonVelocityBC;
                                BoundaryConditions.BoundaryConditionFunctions(bc_id).velocity = &UserDefinedVelocityBC::velocity;
                            });
                            break;   
                        case boundary_conditions::temperature:
                            std::cout << "Setting temperature bc " << std::endl;
                            
                            RUN({
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCHydroType = boundary_conditions::temperature;
                                BoundaryConditions.BoundaryConditionFunctions(bc_id).temperature = &ConstantTempBC::temperature;
                            });
                            break;                        
                        default:
                            
                            std::cout << "Setting velocity bc " << std::endl;
                            RUN({
                                BoundaryConditions.BoundaryConditionEnums(bc_id).BCHydroType = boundary_conditions::noVelocityBC;
                                BoundaryConditions.BoundaryConditionFunctions(bc_id).velocity = &NoVelocityBC::velocity;
                            });
                            // no velocity specified is default
                            break;
                        
                    } // end switch

                    if (VERBOSE) {
                        std::cout << "\ttype = " << type << std::endl;
                    }
                }
                else{
                    std::cout << "ERROR: invalid boundary condition option input in YAML file: " << type << std::endl;
                    std::cout << "Valid options are: " << std::endl;

                    for (const auto& pair : map) {
                        std::cout << "\t" << pair.first << std::endl;
                    }
                } // end if
            } // type
            // get boundary condition direction
            else if (a_word.compare("location") == 0) {
                std::string location = bc_yaml[bc_id]["boundary_condition"][a_word].As<std::string>();

                auto map = bc_location_map;

                // set the direction
                if (map.find(location) != map.end()) {
                    auto bc_location = map[location];
                    RUN({
                        BoundaryConditions.BoundaryConditionEnums(bc_id).Location = bc_location;
                    });
                    if (VERBOSE) {
                        std::cout << "\tlocation = " << location << std::endl;
                    }
                }
                else{
                    std::cout << "ERROR: invalid boundary condition option input in YAML file: " << location << std::endl;
                    std::cout << "Valid options are: " << std::endl;

                    for (const auto& pair : map) {
                        std::cout << "\t" << pair.first << std::endl;
                    }
                } // end if
            } // direction
            // get boundary condition geometry
            // get boundary condition direction
            else if (a_word.compare("direction") == 0) {
                std::string direction = bc_yaml[bc_id]["boundary_condition"][a_word].As<std::string>();

                auto map = bc_direction_map;

                // set the direction
                if (map.find(direction) != map.end()) {
                    auto bc_direction = map[direction];
                    RUN({
                        BoundaryConditions.BoundaryConditionEnums(bc_id).Direction = bc_direction;
                    });
                    if (VERBOSE) {
                        std::cout << "\tdirection = " << direction << std::endl;
                    }
                }
                else{
                    std::cout << "ERROR: invalid boundary condition option input in YAML file: " << direction << std::endl;
                    std::cout << "Valid options are: " << std::endl;

                    for (const auto& pair : map) {
                        std::cout << "\t" << pair.first << std::endl;
                    }
                } // end if
            } // direction
            // get boundary condition geometry
            else if (a_word.compare("geometry") == 0) {
                std::string geometry = bc_yaml[bc_id]["boundary_condition"][a_word].As<std::string>();

                auto map = bc_geometry_map;

                // set the geometry
                if (map.find(geometry) != map.end()) {
                    auto bc_geometry = map[geometry];
                    RUN({
                        BoundaryConditions.BoundaryConditionSetup(bc_id).geometry = bc_geometry;
                    });

                    if (VERBOSE) {
                        std::cout << "\tgeometry = " << geometry << std::endl;
                    }
                }
                else{
                    std::cout << "ERROR: invalid boundary condition option input in YAML file: " << geometry << std::endl;
                    std::cout << "Valid options are: " << std::endl;

                    for (const auto& pair : map) {
                        std::cout << "\t" << pair.first << std::endl;
                    }
                } // end if
            } // geometry
            // set the value
            else if (a_word.compare("value") == 0) {
                double value = bc_yaml[bc_id]["boundary_condition"][a_word].As<double>();
                RUN({
                    BoundaryConditions.BoundaryConditionSetup(bc_id).value = value;
                });
            } // value
            // set the u
            else if (a_word.compare("u") == 0) {
                double u = bc_yaml[bc_id]["boundary_condition"][a_word].As<double>();
                RUN({
                    BoundaryConditions.bc_global_vars(bc_id,0) = u;
                });
            } // u
            // set the v
            else if (a_word.compare("v") == 0) {
                double v = bc_yaml[bc_id]["boundary_condition"][a_word].As<double>();
                RUN({
                    BoundaryConditions.bc_global_vars(bc_id,1) = v;
                });
            } // v
            // set the w
            else if (a_word.compare("w") == 0) {
                double w = bc_yaml[bc_id]["boundary_condition"][a_word].As<double>();
                RUN({
                    BoundaryConditions.bc_global_vars(bc_id,2) = w;
                });
            } // w
            else if (a_word.compare("temperature") == 0) {
                double temperature = bc_yaml[bc_id]["boundary_condition"][a_word].As<double>();
                RUN({
                    BoundaryConditions.bc_global_vars(bc_id,4) = temperature;
                });
            } // temperature
            // set the bc_vel_0
            else if (a_word.compare("hydro_bc_vel_0") == 0) {
                double hydro_bc_vel_0 = bc_yaml[bc_id]["boundary_condition"][a_word].As<double>();
                RUN({
                    BoundaryConditions.bc_global_vars(bc_id,0) = hydro_bc_vel_0;
                });
            } // 
            // set the bc_vel_1
            else if (a_word.compare("hydro_bc_vel_1") == 0) {
                double hydro_bc_vel_1 = bc_yaml[bc_id]["boundary_condition"][a_word].As<double>();
                RUN({
                    BoundaryConditions.bc_global_vars(bc_id,1) = hydro_bc_vel_1;
                });
            } // 
            // set the bc_vel_start
            else if (a_word.compare("hydro_bc_vel_t_start") == 0) {
                double hydro_bc_vel_t_start = bc_yaml[bc_id]["boundary_condition"][a_word].As<double>();
                RUN({
                    BoundaryConditions.bc_global_vars(bc_id,2) = hydro_bc_vel_t_start;
                });
            } // 
            // set the bc_vel_end
            else if (a_word.compare("hydro_bc_vel_t_end") == 0) {
                double hydro_bc_vel_t_end = bc_yaml[bc_id]["boundary_condition"][a_word].As<double>();
                RUN({
                    BoundaryConditions.bc_global_vars(bc_id,3) = hydro_bc_vel_t_end;
                });
            } // 
            else if (a_word.compare("origin") == 0) {
                std::string origin = bc_yaml[bc_id]["boundary_condition"][a_word].As<std::string>();
                if (VERBOSE) {
                    std::cout << "\torigin = " << origin << std::endl;
                }

                // get the origin numbers, values are words
                std::vector<std::string> numbers = exact_array_values(origin, ",");

                double x1 = std::stod(numbers[0]);
                double y1 = std::stod(numbers[1]);
                double z1;

                if(numbers.size()==3){ 
                    // 3D
                    z1 = std::stod(numbers[2]);
                }
                else {
                    // 2D
                    z1 = 0.0;
                } //

                if (VERBOSE) {
                    std::cout << "\tx1 = " << x1 << std::endl;
                    std::cout << "\ty1 = " << y1 << std::endl;
                    std::cout << "\tz1 = " << z1 << std::endl;
                }
                // storing the origin values as (x1,y1,z1)

                RUN({
                    BoundaryConditions.BoundaryConditionSetup(bc_id).origin[0] = x1;
                    BoundaryConditions.BoundaryConditionSetup(bc_id).origin[1] = y1;
                    BoundaryConditions.BoundaryConditionSetup(bc_id).origin[2] = z1;
                });

            } // origin
            else {
                std::cout << "ERROR: invalid input: " << a_word << std::endl;
                std::cout << "Valid options are: " << std::endl;
                for (const auto& element : str_bc_inps) {
                    std::cout << element << std::endl;
                }
                throw std::runtime_error("**** Boundary Conditions Not Understood ****");
            }
        } // end for words in boundary conditions

    } // end loop over BCs specified

    // copy the enum values to the host 
    BoundaryConditions.BoundaryConditionEnums.update_host();

} // end of function to parse region