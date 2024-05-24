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

#include "matar.h"
#include "parse_yaml.h"
#include "simulation_parameters.h"

#define PI 3.141592653589793

bool VERBOSE = false;

// ==============================================================================
//   Function Definitions
// ==============================================================================
// modified code from stack over./flow for string delimiter parsing
std::vector<std::string> exact_array_values(std::string s, std::string delimiter)
{
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    // remove first and last char in the string, which are [ and ]
    s.erase(s.begin());
    s.erase(s.end() - 1);

    // now parse the values in the array into a vector
    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token     = s.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back(token);
    }

    res.push_back(s.substr(pos_start));
    return res;
} // end of extract_array_values

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
void parse_yaml(Yaml::Node& root, simulation_parameters_t& sim_param)
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
    parse_mesh_input(root, sim_param.mesh_input);

    if (VERBOSE) {
        printf("\n");
        std::cout << "Parsing YAML dynamic options:" << std::endl;
    }
    parse_dynamic_options(root, sim_param.dynamic_options);

    if (VERBOSE) {
        printf("\n");
        std::cout << "Parsing YAML output options:" << std::endl;
    }
    parse_output_options(root, sim_param.output_options);

    if (VERBOSE) {
        printf("\n");
        std::cout << "Parsing YAML solver options:" << std::endl;
    }
    parse_solver_input(root, sim_param.solver_inputs);

    if (VERBOSE) {
        printf("\n");
        std::cout << "Parsing YAML boundary condition options:" << std::endl;
    }
    parse_bcs(root, sim_param.boundary_conditions);

    if (VERBOSE) {
        printf("\n");
        std::cout << "Parsing YAML regions:" << std::endl;
    }
    // parse the region yaml text into a vector of region_fills
    parse_regions(root, sim_param.region_fills);

    if (VERBOSE) {
        printf("\n");
        std::cout << "Parsing YAML materials:" << std::endl;
    }
    // parse the material yaml text into a vector of materials
    parse_materials(root, sim_param.materials);
}

// =================================================================================
//    Extract words from the input file and validate they are correct
// =================================================================================
void validate_inputs(Yaml::Node& yaml, std::vector<std::string>& user_inputs, std::vector<std::string>& str_valid_inputs)
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
        validate_inputs(inps_yaml, user_inputs, str_solver_inps);

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
    validate_inputs(yaml, user_dynamic_inps, str_dyn_opts_inps);

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

    mesh_input.origin = DCArrayKokkos<double> (3, "mesh_input.origin");
    mesh_input.length = DCArrayKokkos<double> (3, "mesh_input.length");
    mesh_input.num_elems = DCArrayKokkos<int> (3, "mesh_input.num_elems");


    // extract words from the input file and validate they are correct
    validate_inputs(mesh_yaml, user_mesh_inputs, str_mesh_inps);

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

                if (mesh_input.source == mesh_input::generate && !mesh_input.file_path.empty()) {
                    std::cout << "ERROR: When the mesh source is set to generate, a mesh file cannot be passed in" << std::endl;
                    exit(0);
                }
            }
            else{
                std::cout << "ERROR: invalid mesh option input in YAML file: " << source << std::endl;
                std::cout << "Valid options are: " << std::endl;

                for (const auto& pair : map) {
                    std::cout << "\t" << pair.first << std::endl;
                }
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
            } // end if
        } // type
        // Get mesh file path
        else if (a_word.compare("file_path") == 0) {
            std::string path = root["mesh_options"][a_word].As<std::string>();
            if (VERBOSE) {
                std::cout << "\tfile_path = " << path << std::endl;
            }

            mesh_input.file_path = path;

            if (mesh_input.source == mesh_input::file && mesh_input.file_path.empty()) {
                std::cout << "ERROR: When the mesh source is a file, a file_path must be set to point to the mesh file" << std::endl;
                std::cout << "A mesh can either be generated or read in from a file, but not both" << std::endl;
            }

            if (mesh_input.source == mesh_input::generate) {
                std::cout << "ERROR: When the mesh source is set to generate, a mesh file cannot be passed in" << std::endl;
                exit(0);
            }
        } // file path
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
            double z1 = std::stod(numbers[2]);

            if (VERBOSE) {
                std::cout << "\tx1 = " << x1 << std::endl;
                std::cout << "\ty1 = " << y1 << std::endl;
                std::cout << "\tz1 = " << z1 << std::endl;
            }

            // storing the origin values as
            RUN({
                mesh_input.origin(0) = x1;
                mesh_input.origin(1) = y1;
                mesh_input.origin(2) = z1;
            });
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
            double l3 = std::stod(numbers[2]);

            if (VERBOSE) {
                std::cout << "\tl1 = " << l1 << std::endl;
                std::cout << "\tl2 = " << l2 << std::endl;
                std::cout << "\tl3 = " << l3 << std::endl;
            }

            // storing the length values
            RUN({
                mesh_input.length(0) = l1;
                mesh_input.length(1) = l2;
                mesh_input.length(2) = l3;
            });
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
            int n3 = std::stod(numbers[2]);

            if (VERBOSE) {
                std::cout << "\tn1 = " << n1 << std::endl;
                std::cout << "\tn2 = " << n2 << std::endl;
                std::cout << "\tn3 = " << n3 << std::endl;
            }

            // storing the number of elements
            RUN({
                mesh_input.num_elems(0) = n1;
                mesh_input.num_elems(1) = n2;
                mesh_input.num_elems(2) = n3;
            });
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
        }
    } // end user_mesh_inputs
} // end of parse mesh options

// =================================================================================
//    Parse Output options
// =================================================================================
void parse_output_options(Yaml::Node& root, output_options_t& output_options)
{
    Yaml::Node& out_opts = root["output_options"];

    // get the mesh variables names set by the user
    std::vector<std::string> user_inputs;

    // extract words from the input file and validate they are correct
    validate_inputs(out_opts, user_inputs, str_output_options_inps);

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
        }
    } // end user_inputs
} // end of parse mesh options

// =================================================================================
//    Parse Fill regions
// =================================================================================
void parse_regions(Yaml::Node& root, DCArrayKokkos<reg_fill_t>& region_fills)
{
    Yaml::Node& region_yaml = root["regions"];

    size_t num_regions = region_yaml.Size();

    region_fills = DCArrayKokkos<reg_fill_t>(num_regions , "sim_param.region_fills");

    for(int i=0; i< num_regions; i++){
        region_fills.host(i).origin = DCArrayKokkos<double> (3, "region_fills.origin");
        region_fills.host(i).origin.update_device();
    }
    region_fills.update_device();

    // loop over the fill regions specified
    for (int reg_id = 0; reg_id < num_regions; reg_id++) {
        // read the variables names
        Yaml::Node& inps_yaml = root["regions"][reg_id]["fill_volume"];

        // get the material variables names set by the user
        std::vector<std::string> user_str_region_inps;

        // extract words from the input file and validate they are correct
        validate_inputs(inps_yaml, user_str_region_inps, str_region_inps);

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
                    std::cout << "\tie = " << v << std::endl;
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
            else if (a_word.compare("type") == 0) {
                std::string type = root["regions"][reg_id]["fill_volume"]["type"].As<std::string>();
                if (VERBOSE) {
                    std::cout << "\ttype = " << type << std::endl;
                }

                // set the volume tag type
                if (region_type_map.find(type) != region_type_map.end()) {
                    auto vol = region_type_map[type];

                    RUN({
                        region_fills(reg_id).volume = vol;
                    });
                    if (VERBOSE) {
                        std::cout << "\tvolume_fill = " << type << std::endl;
                    }
                    if (VERBOSE) {
                        std::cout << vol << std::endl;
                    }
                }
                else{
                    std::cout << "ERROR: invalid input: " << type << std::endl;
                } // end if
            } // end volume fill type
            else if (a_word.compare("velocity") == 0) {
                std::string type = root["regions"][reg_id]["fill_volume"]["velocity"].As<std::string>();
                if (VERBOSE) {
                    std::cout << "\tvelocity = " << type << std::endl;
                }

                // set the volume tag type
                if (velocity_type_map.find(type) != velocity_type_map.end()) {
                    auto vel = velocity_type_map[type];

                    RUN({
                        region_fills(reg_id).velocity = vel;
                    });

                    if (VERBOSE) {
                        std::cout << "\tvelocity_fill = " << type << std::endl;
                        std::cout << vel << std::endl;
                    }
                }
                else{
                    std::cout << "ERROR: invalid input: " << type << std::endl;
                } // end if
            } // end velocity
            else if (a_word.compare("origin") == 0) {
                std::string origin = root["regions"][reg_id]["fill_volume"]["origin"].As<std::string>();
                if (VERBOSE) {
                    std::cout << "\torigin = " << origin << std::endl;
                }

                // get the origin numbers, values are words
                std::vector<std::string> numbers = exact_array_values(origin, ",");

                double x1 = std::stod(numbers[0]);
                double y1 = std::stod(numbers[1]);
                double z1 = std::stod(numbers[2]);

                if (VERBOSE) {
                    std::cout << "\tx1 = " << x1 << std::endl;
                    std::cout << "\ty1 = " << y1 << std::endl;
                    std::cout << "\tz1 = " << z1 << std::endl;
                }

                // storing the origin values as (x1,y1,z1)
                RUN({
                    region_fills(reg_id).origin(0) = x1;
                    region_fills(reg_id).origin(1) = y1;
                    region_fills(reg_id).origin(2) = z1;
                });
            } // origin
            else {
                std::cout << "ERROR: invalid input: " << a_word << std::endl;
                std::cout << "Valid options are: " << std::endl;
                for (const auto& element : str_region_inps) {
                    std::cout << element << std::endl;
                }
            }
        } // end for words in material
    } // end loop over regions
} // end of function to parse region

// =================================================================================
//    Parse Material Definitions
// =================================================================================
void parse_materials(Yaml::Node& root, DCArrayKokkos<material_t>& materials)
{
    Yaml::Node& material_yaml = root["materials"];

    size_t num_materials = material_yaml.Size();

    materials = DCArrayKokkos<material_t>(num_materials, "sim_param.materials");

    // allocate room for each material to store eos_global_vars

    // loop over the materials specified
    for (int mat_id = 0; mat_id < num_materials; mat_id++) {
        // read the variables names
        Yaml::Node& inps_yaml = root["materials"][mat_id]["material"];

        size_t num_vars_set = inps_yaml.Size();

        // get the material variables names set by the user
        std::vector<std::string> user_str_material_inps;

        // extract words from the input file and validate they are correct
        validate_inputs(inps_yaml, user_str_material_inps, str_material_inps);

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
                    materials(mat_id).q1 = q1;
                });
            } // q1
            else if (a_word.compare("q1ex") == 0) {
                double q1ex = root["materials"][mat_id]["material"]["q1ex"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tq1ex = " << q1ex << std::endl;
                }
                RUN({
                    materials(mat_id).q1ex = q1ex;
                });
            } // q1ex
            else if (a_word.compare("q2") == 0) {
                // outer plane

                double q2 = root["materials"][mat_id]["material"]["q2"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tq2 = " << q2 << std::endl;
                }

                RUN({
                    materials(mat_id).q2 = q2;
                });
            } // q2
            else if (a_word.compare("q2ex") == 0) {
                // outer plane

                double q2ex = root["materials"][mat_id]["material"]["q2ex"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tq2ex = " << q2ex << std::endl;
                }
                RUN({
                    materials(mat_id).q2ex = q2ex;
                });
            } // q1ex
            else if (a_word.compare("id") == 0) {
                int m_id = root["materials"][mat_id]["material"]["id"].As<int>();
                if (VERBOSE) {
                    std::cout << "\tid = " << m_id << std::endl;
                }
                RUN({
                    materials(mat_id).id = m_id;
                });
            } // id
            else if (a_word.compare("elastic_modulus") == 0) {
                double elastic_modulus = root["materials"][mat_id]["material"]["elastic_modulus"].As<double>();
                if (VERBOSE) {
                    std::cout << "\telastic_modulus = " << elastic_modulus << std::endl;
                }
                RUN({
                    materials(mat_id).elastic_modulus = elastic_modulus;
                });
            } // elastic_modulus
            else if (a_word.compare("poisson_ratio") == 0) {
                double poisson_ratio = root["materials"][mat_id]["material"]["poisson_ratio"].As<double>();
                if (VERBOSE) {
                    std::cout << "\tpoisson_ratio = " << poisson_ratio << std::endl;
                }
                RUN({
                    materials(mat_id).poisson_ratio = poisson_ratio;
                });
            } // poisson_ratio
            else if (a_word.compare("eos_model") == 0) {
                std::string eos = root["materials"][mat_id]["material"]["eos_model"].As<std::string>();

                // set the EOS
                if (eos_map.find(eos) != eos_map.end()) {
                    
                    switch(eos_map[eos]){

                        case model::no_eos:
                            RUN({
                                materials(mat_id).eos_type = model::no_eos;
                                materials(mat_id).eos_model = no_eos;
                            });
                            if (VERBOSE) {
                                std::cout << "\teos_model = " << eos << std::endl;
                            }
                            break;

                        case model::ideal_gas:
                            RUN({
                                materials(mat_id).eos_type = model::ideal_gas;
                                materials(mat_id).eos_model = ideal_gas;
                            });
                            if (VERBOSE) {
                                std::cout << "\teos_model = " << eos << std::endl;
                            }
                            break;
                        case model::user_defined_eos:
                            RUN({
                                materials(mat_id).eos_type = model::user_defined_eos;
                                materials(mat_id).eos_model = user_eos_model;
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

                        case model::no_strength_type:
                            RUN({
                                materials(mat_id).strength_type = model::no_strength_type;
                            });
                            if (VERBOSE) {
                                std::cout << "\tstrength_model_type_type = " << strength_model_type << std::endl;
                            }
                            break;

                        case model::increment_based:
                            RUN({
                                materials(mat_id).strength_type = model::increment_based;
                            });
                            
                            if (VERBOSE) {
                                std::cout << "\tstrength_model_type = " << strength_model_type << std::endl;
                            }
                            break;
                        case model::state_based:
                            RUN({
                                materials(mat_id).strength_type = model::state_based;
                            });
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

                        case model::no_strength_model:
                            RUN({
                                materials(mat_id).strength_model = no_strength;
                            });
                            if (VERBOSE) {
                                std::cout << "\tstrength_model = " << strength_model << std::endl;
                            }
                            break;

                        case model::user_defined_strength:
                            RUN({
                                materials(mat_id).strength_model = user_strength_model;
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
            // exact the eos_global_vars
            else if (a_word.compare("eos_global_vars") == 0) {
                Yaml::Node & mat_global_vars_yaml = root["materials"][mat_id]["material"][a_word];

                size_t num_global_vars = mat_global_vars_yaml.Size();

                std::cout << "*** parsing num global eos vars = " << num_global_vars << std::endl;
                

                materials.host(mat_id).eos_global_vars = DCArrayKokkos<double>(num_global_vars, "material.eos_global_vars");
                materials.host(mat_id).eos_global_vars.update_device();
                RUN({
                    materials(mat_id).num_eos_global_vars = num_global_vars;
                });
                materials.update_device();

                if (VERBOSE) {
                    std::cout << "num global eos vars = " << num_global_vars << std::endl;
                }

                for (int global_var_id = 0; global_var_id < num_global_vars; global_var_id++) {
                    double eos_var = root["materials"][mat_id]["material"]["eos_global_vars"][global_var_id].As<double>();
                    

                    RUN({
                        materials(mat_id).eos_global_vars(global_var_id) = eos_var;
                    });

                    if (VERBOSE) {
                        std::cout << "\t var = " << eos_var << std::endl;
                    }
                } // end loop over global vars
            } // "eos_global_vars"
            else {
                std::cout << "ERROR: invalid input: " << a_word << std::endl;
                std::cout << "Valid options are: " << std::endl;
                for (const auto& element : str_material_inps) {
                    std::cout << element << std::endl;
                }
            }
        } // end for words in material
    } // end loop over materials
} // end of function to parse material information

// =================================================================================
//    Parse Boundary Conditions
// =================================================================================
void parse_bcs(Yaml::Node& root, DCArrayKokkos<boundary_condition_t>& boundary_conditions)
{
    Yaml::Node& bc_yaml = root["boundary_conditions"];

    size_t num_bcs = bc_yaml.Size();

    boundary_conditions = DCArrayKokkos<boundary_condition_t>(num_bcs, "sim_param.boundary_conditions");

    for(int i=0; i< num_bcs; i++){
        boundary_conditions.host(i).origin = DCArrayKokkos<double> (3, "boundary_conditions.origin");
        boundary_conditions.host(i).origin.update_device();
    }
    boundary_conditions.update_device();

    // loop over the fill regions specified
    for (int bc_id = 0; bc_id < num_bcs; bc_id++) {
        // read the variables names
        Yaml::Node& inps_yaml = bc_yaml[bc_id]["boundary_condition"];

        // get the material variables names set by the user
        std::vector<std::string> user_str_bc_inps;

        // extract words from the input file and validate they are correct
        validate_inputs(inps_yaml, user_str_bc_inps, str_bc_inps);

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
                        boundary_conditions(bc_id).solver = bc_solver;
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
                    RUN({
                        boundary_conditions(bc_id).type = bc_type;
                    });

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
            else if (a_word.compare("direction") == 0) {
                std::string direction = bc_yaml[bc_id]["boundary_condition"][a_word].As<std::string>();

                auto map = bc_direction_map;

                // set the direction
                if (map.find(direction) != map.end()) {
                    auto bc_direction = map[direction];
                    RUN({
                        boundary_conditions(bc_id).direction = bc_direction;
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
                        boundary_conditions(bc_id).geometry = bc_geometry;
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
                    boundary_conditions(bc_id).value = value;
                });
            } // value
            // set the u
            else if (a_word.compare("u") == 0) {
                double u = bc_yaml[bc_id]["boundary_condition"][a_word].As<double>();
                RUN({
                    boundary_conditions(bc_id).u = u;
                });
            } // u
            // set the v
            else if (a_word.compare("v") == 0) {
                double v = bc_yaml[bc_id]["boundary_condition"][a_word].As<double>();
                RUN({
                    boundary_conditions(bc_id).v = v;
                });
            } // v
            // set the w
            else if (a_word.compare("w") == 0) {
                double w = bc_yaml[bc_id]["boundary_condition"][a_word].As<double>();
                RUN({
                    boundary_conditions(bc_id).w = w;
                });
            } // w
            else if (a_word.compare("origin") == 0) {
                std::string origin = bc_yaml[bc_id]["boundary_condition"][a_word].As<std::string>();
                if (VERBOSE) {
                    std::cout << "\torigin = " << origin << std::endl;
                }

                // get the origin numbers, values are words
                std::vector<std::string> numbers = exact_array_values(origin, ",");

                double x1 = std::stod(numbers[0]);
                double y1 = std::stod(numbers[1]);
                double z1 = std::stod(numbers[2]);
                if (VERBOSE) {
                    std::cout << "\tx1 = " << x1 << std::endl;
                    std::cout << "\ty1 = " << y1 << std::endl;
                    std::cout << "\tz1 = " << z1 << std::endl;
                }
                // storing the origin values as (x1,y1,z1)

                RUN({
                    boundary_conditions(bc_id).origin(0) = x1;
                    boundary_conditions(bc_id).origin(1) = y1;
                    boundary_conditions(bc_id).origin(2) = z1;
                });
            } // origin
            else {
                std::cout << "ERROR: invalid input: " << a_word << std::endl;
                std::cout << "Valid options are: " << std::endl;
                for (const auto& element : str_bc_inps) {
                    std::cout << element << std::endl;
                }
            }
        } // end for words in material
    } // end loop over regions
} // end of function to parse region