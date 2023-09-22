/**********************************************************************************************
 © 2020. Triad National Security, LLC. All rights reserved.
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
#pragma once
#ifndef SIMULATION_PARAMETERS_DYNAMIC_OPTIMIZATION_H
#define SIMULATION_PARAMETERS_DYNAMIC_OPTIMIZATION_H

#include "yaml-serializable.h"
#include "Simulation_Parameters.h"

SERIALIZABLE_ENUM(FUNCTION_TYPE,
  OBJECTIVE, 
  MULTI_OBJECTIVE_TERM, 
  EQUALITY_CONSTRAINT, 
  INEQUALITY_CONSTRAINT, 
  VECTOR_EQUALITY_CONSTRAINT, 
  VECTOR_INEQUALITY_CONSTRAINT
)

SERIALIZABLE_ENUM(TO_MODULE_TYPE,
  Kinetic_Energy_Minimize,
  Multi_Objective,
  Heat_Capacity_Potential_Minimize,
  Strain_Energy_Minimize,
  Mass_Constraint,
  Moment_of_Inertia_Constraint,
  Heat_Capacity_Potential_Constraint
)

SERIALIZABLE_ENUM(OPTIMIZATION_PROCESS, none, topology_optimization, shape_optimization)
SERIALIZABLE_ENUM(OPTIMIZATION_OBJECTIVE, minimize_kinetic_energy)
SERIALIZABLE_ENUM(CONSTRAINT_TYPE, mass)
SERIALIZABLE_ENUM(RELATION, equality)
struct Optimization_Constraint {
  std::optional<double> value;
  CONSTRAINT_TYPE type;
  RELATION relation;
};
IMPL_YAML_SERIALIZABLE_FOR(Optimization_Constraint, value, type, relation)

struct Optimization_Options : Yaml::ValidatedYaml, Yaml::DerivedFields {
  OPTIMIZATION_PROCESS optimization_process = OPTIMIZATION_PROCESS::none;
  OPTIMIZATION_OBJECTIVE optimization_objective;
  std::vector<Optimization_Constraint> constraints;
  bool method_of_moving_asymptotes;
  double simp_penalty_power;
  double density_epsilon;

  void validate() { }
};
IMPL_YAML_SERIALIZABLE_FOR(Optimization_Options, 
  optimization_process, optimization_objective, 
  constraints, method_of_moving_asymptotes,
  simp_penalty_power, density_epsilon
)

struct Simulation_Parameters_Dynamic_Optimization : public Simulation_Parameters {
  // --- Mesh regions and material fills ---
  int NB  = 6; // number of boundary patch sets to tag
  int NBD = 2; //number of density boundary conditions

  std::optional<Optimization_Options> optimization_options;
  
  //When on, all element nodes connected to a boundary condition patch will have their density constrained
  bool thick_condition_boundary = true;

  //file output parameters
  int optimization_output_freq;

  //Topology Optimization parameters
  double penalty_power = 3.0;

  bool nodal_density_flag = true;
  
  // Non-serialized fields
  bool helmholtz_filter  = false;
  double density_epsilon = 0.0001;
  //list of TO functions needed by problem
  std::vector<TO_MODULE_TYPE> TO_Module_List;
  std::vector<FUNCTION_TYPE> TO_Function_Type;
  std::vector<int> TO_Module_My_FEA_Module;
  std::vector<std::vector<int>> FEA_Module_My_TO_Modules;
  std::vector<std::vector<double>> Function_Arguments;

  //Topology Optimization flags
  bool topology_optimization_on = false;
  bool shape_optimization_on    = false;

  void derive_from_optimization_options() {
    if (!optimization_options.has_value())
      return;
    auto options = optimization_options.value();

    shape_optimization_on = options.optimization_process == OPTIMIZATION_PROCESS::shape_optimization;
    topology_optimization_on = options.optimization_process == OPTIMIZATION_PROCESS::topology_optimization;

    TO_Module_List.resize(options.constraints.size());
    TO_Function_Type.resize(options.constraints.size());
    Function_Arguments.resize(options.constraints.size());

    for (size_t i = 0; i < options.constraints.size(); i++) {
      auto constraint = options.constraints[i];

      // TODO: This whole thing is pretty messed up.
      // Both FEA Modules and TO_Modules really need to be structs
      // of their own. ATM we are assuming a lot about the input.
      // If the constraints aren't set correctly, we will end up with a lot
      // of duplicate/ill defined TO module specifications.
      if (constraint.type == CONSTRAINT_TYPE::mass)
        TO_Module_List[i] = TO_MODULE_TYPE::Mass_Constraint;
      if (constraint.relation == RELATION::equality)
        TO_Function_Type[i] = FUNCTION_TYPE::EQUALITY_CONSTRAINT;
      if (constraint.value.has_value())
        Function_Arguments[i] = { constraint.value.value() };
    }

    switch (options.optimization_objective) {
      case OPTIMIZATION_OBJECTIVE::minimize_kinetic_energy:
        add_TO_module(TO_MODULE_TYPE::Kinetic_Energy_Minimize, FUNCTION_TYPE::OBJECTIVE, {});
        break;
      default:
        throw Yaml::ConfigurationException("Unsupported optimization objective " 
          + to_string(options.optimization_objective)
        );
    }
  }

  void map_TO_to_FEA() {
    // Now we allocate the vectors for each of the currently identified modules.
    TO_Module_My_FEA_Module.resize(TO_Module_List.size());
    FEA_Module_My_TO_Modules.resize(FEA_Modules_List.size());
    
    // Finally we can set up the maps.
    // 
    // TODO: This should really use two `std::map<int, std::set<int>>`s
    // instead of this vector stuff.
    for (size_t to_index = 0; to_index < TO_Module_List.size(); to_index++) {
      auto to_module = TO_Module_List[to_index];
      size_t fea_index = find_TO_module_dependency(to_module);

      TO_Module_My_FEA_Module[to_index] = fea_index;
      FEA_Module_My_TO_Modules[fea_index].push_back(to_index);
    }
  }

  void derive() {
    derive_from_optimization_options();
    map_TO_to_FEA();
  }

  void validate() {
    // Check that the FEA module dependencies are satisfied.
    for (auto to_module : TO_Module_List)
      find_TO_module_dependency(to_module); 
  }

  void add_TO_module(TO_MODULE_TYPE type, FUNCTION_TYPE function_type, std::vector<double> arguments) {
    if (std::find(TO_Module_List.begin(), TO_Module_List.end(), type) != TO_Module_List.end())
      return; // Already have it.
    
    TO_Module_List.push_back(type);
    TO_Function_Type.push_back(function_type);
    Function_Arguments.push_back(arguments);
  }

  /**
   * Find the TO module dependency in the FEA_Module_List list.
   * Returns the index of the dependency satisfier. 
   * 
   * Throws a Yaml::ConfigurationException if one does not exist.
  */
  size_t find_TO_module_dependency(TO_MODULE_TYPE type) {
    switch (type) {
      case TO_MODULE_TYPE::Kinetic_Energy_Minimize:
        validate_one_of_modules_are_specified({FEA_MODULE_TYPE::SGH, FEA_MODULE_TYPE::Dynamic_Elasticity});
        return find_one_of_module({FEA_MODULE_TYPE::SGH, FEA_MODULE_TYPE::Dynamic_Elasticity});
      case TO_MODULE_TYPE::Heat_Capacity_Potential_Minimize:
        validate_module_is_specified(FEA_MODULE_TYPE::Heat_Conduction);
        return find_module(FEA_MODULE_TYPE::Heat_Conduction);
      case TO_MODULE_TYPE::Heat_Capacity_Potential_Constraint:
        validate_module_is_specified(FEA_MODULE_TYPE::Heat_Conduction);
        return find_module(FEA_MODULE_TYPE::Heat_Conduction);
      case TO_MODULE_TYPE::Mass_Constraint:
        validate_module_is_specified(FEA_MODULE_TYPE::Inertial);
        return find_module(FEA_MODULE_TYPE::Inertial);
      case TO_MODULE_TYPE::Moment_of_Inertia_Constraint:
        validate_module_is_specified(FEA_MODULE_TYPE::Inertial);
        return find_module(FEA_MODULE_TYPE::Inertial);
      default:
        throw Yaml::ConfigurationException(
          "Unsupported optimization module type " + to_string(type)
        );
    }
  }
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Simulation_Parameters_Dynamic_Optimization, Simulation_Parameters,
  optimization_options, nodal_density_flag, thick_condition_boundary,
  optimization_output_freq
)

#endif // end HEADER_H
