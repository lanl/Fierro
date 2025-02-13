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

#ifndef FIERRO_PARSE_YAML_H
#define FIERRO_PARSE_YAML_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <sys/stat.h>

#include "matar.h"

#include "Yaml.hpp"

struct SimulationParameters_t;
struct solver_input_t;
struct mesh_input_t;
struct RegionFill_t;
struct RegionFill_host_t;
struct output_options_t;

struct Material_t;
struct MaterialSetup_t;
struct MaterialFunctions_t;
struct MaterialEnums_t;

struct BoundaryCondition_t;
struct BoundaryConditionSetup_t;
struct BoundaryConditionEnums_t;
struct BoundaryConditionFunctions_t;

struct dynamic_options_t;

using namespace mtr;

// checks to see if a path exists
static bool DoesPathExist(const std::string& s)
{
    struct stat buffer;
    return (stat(s.c_str(), &buffer) == 0);
}

// prints the contents of a parsed yaml file
void print_yaml(Yaml::Node root);

// Read and validate user inputs
void validate_inputs(
    Yaml::Node& yaml, 
    std::vector<std::string>& user_inputs, 
    std::vector<std::string>& str_valid_inputs,
    std::vector<std::string>& str_required_inputs);

// utility function for parsing YAML file
void parse_yaml(Yaml::Node& root, SimulationParameters_t& SimulationParamaters, Material_t& Materials, BoundaryCondition_t& Boundary);

// Parse the solver related data
void parse_solver_input(Yaml::Node& root, std::vector<solver_input_t>& solver_input);

// Parse dynamic time related options
void parse_dynamic_options(Yaml::Node& root, dynamic_options_t& dynamic_options);

// Parse the mesh related data
void parse_mesh_input(Yaml::Node& root, mesh_input_t& mesh_input);

// Parse output options
void parse_output_options(Yaml::Node& root, output_options_t& output_options);

// parse the region text
void parse_regions(Yaml::Node& root, 
                   CArrayKokkos<RegionFill_t>& region_fills,
                   CArray<RegionFill_host_t>& region_fills_host);

// parse the region text
void parse_materials(Yaml::Node& root, Material_t& Materials, const size_t num_dims);

// parse the boundary condition text
void parse_bcs(Yaml::Node& root, BoundaryCondition_t& BoundaryConditions, const size_t num_solvers);

#endif // end Header Guard