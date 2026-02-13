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

#include "parse_yaml.hpp"

#include "parse_tools.hpp"
#include "parse_regions.hpp"
#include "parse_mesh_inputs.hpp"
#include "parse_solver_inputs.hpp"
#include "parse_material_inputs.hpp"
#include "parse_bdy_conds_inputs.hpp"
#include "parse_dynamic_inputs.hpp"
#include "parse_output_options.hpp"

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
//    Parse YAML file
// =================================================================================
void parse_yaml(Yaml::Node& root, SimulationParameters_t& SimulationParamaters, Material_t& Materials, BoundaryCondition_t& Boundary)
{

    parse_mesh_inputs(root, SimulationParamaters.mesh_input);

    parse_dynamic_options(root, SimulationParamaters.dynamic_options);

    parse_output_options(root, SimulationParamaters.output_options);

    parse_solver_input(root, SimulationParamaters.solver_inputs);

    // parse the region yaml text into a vector of boundary conditions
    size_t num_solvers = SimulationParamaters.solver_inputs.size();
    parse_bcs(root, Boundary, num_solvers);


    // parse the region yaml text into a vector of region_fills
    parse_regions(root, 
                  SimulationParamaters.region_setups.reg_fills_in_solver,
                  SimulationParamaters.region_setups.num_reg_fills_in_solver,
                  SimulationParamaters.region_setups.region_fills,
                  SimulationParamaters.region_setups.region_fills_host,
                  SimulationParamaters.region_setups.fill_gauss_states,
                  SimulationParamaters.region_setups.fill_node_states,
                  num_solvers);

    // parse the material yaml text into a vector of materials
    parse_materials(root, Materials, SimulationParamaters.mesh_input.num_dims);
    parse_multimaterial_options(root, Materials);
}




