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
#include "parse_mesh_inputs.hpp"

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
//    Parse Mesh options
// =================================================================================
void parse_mesh_inputs(Yaml::Node& root, mesh_input_t& mesh_input)
{
    Yaml::Node& mesh_yaml = root["mesh_options"];

    // get the mesh variables names set by the user
    std::vector<std::string> user_mesh_inputs;


    // extract words from the input file and validate they are correct
    validate_inputs(mesh_yaml, user_mesh_inputs, str_mesh_inps, mesh_required_inps);

    // loop over the words in the material input definition
    for (auto& a_word : user_mesh_inputs) {


        Yaml::Node& material_inps_yaml = root["mesh_options"][a_word];

        // get mesh source [generate or from file]
        if (a_word.compare("source") == 0) {
            std::string source = root["mesh_options"][a_word].As<std::string>();

            auto map = mesh_input_source_map;

            // set the mesh source
            if (map.find(source) != map.end()) {
                mesh_input.source = map[source];
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

            mesh_input.num_dims = num_dim;
        } // Number of dimensions
        // get mesh type for generation
        else if (a_word.compare("type") == 0) {
            std::string type = root["mesh_options"][a_word].As<std::string>();

            auto map = mesh_input_type_map;

            // set the mesh type
            if (map.find(type) != map.end()) {
                mesh_input.type = map[type];
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

            // absolute path to file or local to the director where exe is run
            mesh_input.file_path = path;  

        } // end file path
        // Origin for the mesh
        else if (a_word.compare("origin") == 0) {
            std::string origin = root["mesh_options"][a_word].As<std::string>();

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

            // storing the origin values as
            mesh_input.origin[0] = x1;
            mesh_input.origin[1] = y1;
            mesh_input.origin[2] = z1;
        }
        // Extents of the mesh
        else if (a_word.compare("length") == 0) {
            std::string length = root["mesh_options"][a_word].As<std::string>();

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

            // storing the length values
            mesh_input.length[0] = l1;
            mesh_input.length[1] = l2;
            mesh_input.length[2] = l3;
        }
        // Number of elements per direction
        else if (a_word.compare("num_elems") == 0) {
            std::string elem_count = root["mesh_options"][a_word].As<std::string>();

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

            // storing the number of elements
            mesh_input.num_elems[0] = n1;
            mesh_input.num_elems[1] = n2;
            mesh_input.num_elems[2] = n3;
        }
        // Polynomial order for the mesh
        else if (a_word.compare("polynomial_order") == 0) {
            int p_order = root["mesh_options"][a_word].As<int>();

            mesh_input.p_order = p_order;
        } // polynomial order
        // inner radius for 2D RZ meshes
        else if (a_word.compare("inner_radius") == 0) {
            double inner_radius = root["mesh_options"][a_word].As<double>();

            mesh_input.inner_radius = inner_radius;
        } // inner radius for 2D RZ meshes
        // outer radius for 2D RZ meshes
        else if (a_word.compare("outer_radius") == 0) {
            double outer_radius = root["mesh_options"][a_word].As<double>();

            mesh_input.outer_radius = outer_radius;
        } // outer radius for 2D RZ meshes
        // starting angle for 2D RZ meshes
        else if (a_word.compare("starting_angle") == 0) {
            double starting_angle = root["mesh_options"][a_word].As<double>();

            mesh_input.starting_angle = starting_angle;
        } // starting angle for 2D RZ meshes
        // ending angle for 2D RZ meshes
        else if (a_word.compare("ending_angle") == 0) {
            double ending_angle = root["mesh_options"][a_word].As<double>();

            mesh_input.ending_angle = ending_angle;
        } // ending angle for 2D RZ meshes
        // Number of radial elements for 2D RZ meshes
        else if (a_word.compare("num_radial_elems") == 0) {
            int num_radial_elems = root["mesh_options"][a_word].As<int>();

            mesh_input.num_radial_elems = num_radial_elems;
        } // Number of radial elements for 2D RZ meshes
        // Number of angular elements for 2D RZ meshes
        else if (a_word.compare("num_angular_elems") == 0) {
            int num_angular_elems = root["mesh_options"][a_word].As<int>();

            mesh_input.num_angular_elems = num_angular_elems;
        } // Number of angular elements for 2D RZ meshes
        else if (a_word.compare("scale_x") == 0) {

            double scale_x = root["mesh_options"][a_word].As<double>();

            mesh_input.scale_x = scale_x;

        } // end scale_x
        else if (a_word.compare("scale_y") == 0) {

            double scale_y = root["mesh_options"][a_word].As<double>();

            mesh_input.scale_y = scale_y;

        } // end scale_y
        else if (a_word.compare("scale_z") == 0) {

            double scale_z = root["mesh_options"][a_word].As<double>();

            mesh_input.scale_z = scale_z;

        } // end scale_z
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
