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

#ifndef SIMULATION_PARAMETERS_TOPOLOGY_OPTIMIZATION_H
#define SIMULATION_PARAMETERS_TOPOLOGY_OPTIMIZATION_H

#include "utilities.h"
#include "Simulation_Parameters.h"
#include "yaml-serializable.h"

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
  Heat_Capacity_Potential_Constraint,
  MULTI_OBJECTIVE_TERM,
  Thermo_Elastic_Strain_Energy_Minimize,
  Strain_Energy_Constraint
)

SERIALIZABLE_ENUM(OPTIMIZATION_PROCESS, none, topology_optimization, shape_optimization)
SERIALIZABLE_ENUM(OPTIMIZATION_OBJECTIVE,
  minimize_compliance, minimize_kinetic_energy,
  minimize_thermal_resistance, multi_objective
)
SERIALIZABLE_ENUM(CONSTRAINT_TYPE, mass, moment_of_inertia)
SERIALIZABLE_ENUM(CONSTRAINT_COMPONENT, xx, yy, zz, xy, xz, yz)
inline int component_to_int(CONSTRAINT_COMPONENT c) {
  switch(c) {
    case CONSTRAINT_COMPONENT::xx:
      return 0;
    case CONSTRAINT_COMPONENT::yy:
      return 1;
    case CONSTRAINT_COMPONENT::zz:
      return 2;
    case CONSTRAINT_COMPONENT::xy:
      return 3;
    case CONSTRAINT_COMPONENT::xz:
      return 4;
    case CONSTRAINT_COMPONENT::yz:
      return 5;
    default:
      throw std::runtime_error("Unsupported component " + to_string(c));
  }
}

SERIALIZABLE_ENUM(RELATION, equality)
SERIALIZABLE_ENUM(DENSITY_FILTER, none, hemlholtz_filter)
SERIALIZABLE_ENUM(MULTI_OBJECTIVE_STRUCTURE, linear)

struct Optimization_Constraint : Yaml::ValidatedYaml {
  double value;
  CONSTRAINT_TYPE type;
  RELATION relation;
  std::optional<CONSTRAINT_COMPONENT> component;
  bool inertia_center_x = false;
  bool inertia_center_y = false;
  bool inertia_center_z = false;

  void validate() {
    if (type == CONSTRAINT_TYPE::moment_of_inertia) {
      if (!component.has_value())
        throw Yaml::ConfigurationException("`component` field required for constraint type " + to_string(type));
    }
  }
};
YAML_ADD_REQUIRED_FIELDS_FOR(Optimization_Constraint, value, type)
IMPL_YAML_SERIALIZABLE_FOR(Optimization_Constraint, 
  value, type, relation, component,
  inertia_center_x, inertia_center_y, inertia_center_z
)


struct MultiObjectiveModule {
  OPTIMIZATION_OBJECTIVE type;
  double weight_coefficient;
};
YAML_ADD_REQUIRED_FIELDS_FOR(MultiObjectiveModule, type, weight_coefficient)
IMPL_YAML_SERIALIZABLE_FOR(MultiObjectiveModule, type, weight_coefficient)

struct Optimization_Options : Yaml::ValidatedYaml, Yaml::DerivedFields {
  OPTIMIZATION_PROCESS optimization_process = OPTIMIZATION_PROCESS::none;
  OPTIMIZATION_OBJECTIVE optimization_objective;
  std::vector<Optimization_Constraint> constraints;
  bool method_of_moving_asymptotes;
  double simp_penalty_power;
  double density_epsilon;
  DENSITY_FILTER density_filter = DENSITY_FILTER::none;

  MULTI_OBJECTIVE_STRUCTURE multi_objective_structure = MULTI_OBJECTIVE_STRUCTURE::linear;
  std::vector<MultiObjectiveModule> multi_objective_modules;
};
IMPL_YAML_SERIALIZABLE_FOR(Optimization_Options, 
  optimization_process, optimization_objective, 
  constraints, method_of_moving_asymptotes,
  simp_penalty_power, density_epsilon, density_filter,
  multi_objective_structure, multi_objective_modules
)

struct Simulation_Parameters_Topology_Optimization : public Simulation_Parameters {
  // --- Mesh regions and material fills ---
  int NB  = 6; // number of boundary patch sets to tag
  int NBD = 2; //number of density boundary conditions

  Optimization_Options optimization_options;
  
  //When on, all element nodes connected to a boundary condition patch will have their density constrained
  bool thick_condition_boundary = true;

  //file output parameters
  int optimization_output_freq = 2000;

  // Non-serialized fields
  //list of TO functions needed by problem
  std::vector<TO_MODULE_TYPE> TO_Module_List {};
  std::vector<FUNCTION_TYPE> TO_Function_Type {};
  std::vector<int> TO_Module_My_FEA_Module {};
  std::vector<std::vector<int>> FEA_Module_My_TO_Modules {};
  std::vector<std::vector<double>> Function_Arguments {};

  std::vector<int> Multi_Objective_Modules;
  std::vector<double> Multi_Objective_Weights;

  //Topology Optimization flags
  bool topology_optimization_on = false;
  bool shape_optimization_on    = false;
  bool nodal_density_flag       = true;
  bool helmholtz_filter         = false;

  void derive_objective_module() {
    switch (optimization_options.optimization_objective) {
      case OPTIMIZATION_OBJECTIVE::minimize_compliance:
        add_TO_module(TO_MODULE_TYPE::Strain_Energy_Minimize, FUNCTION_TYPE::OBJECTIVE, {});
        break;
      case OPTIMIZATION_OBJECTIVE::minimize_thermal_resistance:
        add_TO_module(TO_MODULE_TYPE::Heat_Capacity_Potential_Minimize, FUNCTION_TYPE::OBJECTIVE, {});
        break;
      case OPTIMIZATION_OBJECTIVE::multi_objective:
        add_TO_module(TO_MODULE_TYPE::Multi_Objective, FUNCTION_TYPE::OBJECTIVE, {});
        derive_multi_objectives();
        break;
      default:
        throw Yaml::ConfigurationException("Unsupported optimization objective " + to_string(optimization_options.optimization_objective));
    }
  }
  void derive_multi_objectives() {
    if (optimization_options.optimization_objective != OPTIMIZATION_OBJECTIVE::multi_objective)
      return;
    
    for (auto mod : optimization_options.multi_objective_modules) {
      TO_MODULE_TYPE to_type;
      switch (mod.type) {
        case OPTIMIZATION_OBJECTIVE::minimize_compliance:
          to_type = TO_MODULE_TYPE::Strain_Energy_Minimize;
          break;
        case OPTIMIZATION_OBJECTIVE::minimize_thermal_resistance:
          to_type = TO_MODULE_TYPE::Heat_Capacity_Potential_Minimize;
          break;
        default:
          throw Yaml::ConfigurationException("Unsupported sub-objective " + to_string(mod.type));
      }
      Multi_Objective_Modules.push_back(TO_Module_List.size());
      Multi_Objective_Weights.push_back(mod.weight_coefficient);
      add_TO_module(to_type, FUNCTION_TYPE::MULTI_OBJECTIVE_TERM, {});
    }
  }
  void derive_constraint_modules() {
    for (auto constraint : optimization_options.constraints) {
      FUNCTION_TYPE f_type;
      switch (constraint.relation) {
        case RELATION::equality:
          f_type = FUNCTION_TYPE::EQUALITY_CONSTRAINT;
          break;
        default:
          throw Yaml::ConfigurationException("Unsupported relation " + to_string(constraint.relation));
      }

      switch (constraint.type) {
        case CONSTRAINT_TYPE::mass:
          add_TO_module(
            TO_MODULE_TYPE::Mass_Constraint, 
            f_type, 
            {constraint.value}
          );
          break;

        case CONSTRAINT_TYPE::moment_of_inertia:
          add_TO_module(
            TO_MODULE_TYPE::Moment_of_Inertia_Constraint, 
            f_type, 
            {constraint.value, (double)component_to_int(constraint.component.value())}
          );
          break;

        default:
          throw Yaml::ConfigurationException("Unsupported constraint type " + to_string(constraint.type));
      }
    }
  }
  void derive() {
    if (optimization_options.density_filter == DENSITY_FILTER::hemlholtz_filter)
      helmholtz_filter = true;

    shape_optimization_on = optimization_options.optimization_process == OPTIMIZATION_PROCESS::shape_optimization;
    topology_optimization_on = optimization_options.optimization_process == OPTIMIZATION_PROCESS::topology_optimization;

    derive_objective_module();
    derive_constraint_modules();
  
    // Now we allocate the vectors for each of the currently identified modules.
    TO_Module_My_FEA_Module.resize(TO_Module_List.size());
    FEA_Module_My_TO_Modules.resize(FEA_Modules_List.size());
    
    // Finally we can set up the maps.
    // 
    // TODO: This should really use two `std::map<int, std::set<int>>`s
    // instead of this vector stuff.
    for (size_t to_index = 0; to_index < TO_Module_List.size(); to_index++) {
      auto to_module = TO_Module_List[to_index];
      auto fea_module = get_TO_module_dependency(to_module);
      if (!fea_module.has_value())
        continue;

      size_t fea_index = find_module(fea_module.value());

      TO_Module_My_FEA_Module[to_index] = fea_index;
      FEA_Module_My_TO_Modules[fea_index].push_back(to_index);
    }
  }

  void validate() {
    validate_unique_vector(TO_Module_List, "Duplicate Topology Optimization Module listed: ");
    validate_unique_vector(FEA_Modules_List, "Duplicate FEA Module listed: ");
  }

  void add_TO_module(TO_MODULE_TYPE type, FUNCTION_TYPE function_type, std::vector<double> arguments) {
    TO_Module_List.push_back(type);
    TO_Function_Type.push_back(function_type);
    Function_Arguments.push_back(arguments);
    
    auto fea_module = get_TO_module_dependency(type);
    if (fea_module.has_value())
      ensure_module(fea_module.value());
  }
  std::optional<FEA_MODULE_TYPE> get_TO_module_dependency(TO_MODULE_TYPE type) {
    switch (type) {
      case TO_MODULE_TYPE::Heat_Capacity_Potential_Minimize:
        return FEA_MODULE_TYPE::Heat_Conduction;
      case TO_MODULE_TYPE::Heat_Capacity_Potential_Constraint:
        return FEA_MODULE_TYPE::Heat_Conduction;
      case TO_MODULE_TYPE::Mass_Constraint:
        return FEA_MODULE_TYPE::Inertial;
      case TO_MODULE_TYPE::Moment_of_Inertia_Constraint:
        return FEA_MODULE_TYPE::Inertial;
      case TO_MODULE_TYPE::Strain_Energy_Minimize:
        return FEA_MODULE_TYPE::Elasticity;
      case TO_MODULE_TYPE::Thermo_Elastic_Strain_Energy_Minimize:
        return FEA_MODULE_TYPE::Heat_Conduction;
      case TO_MODULE_TYPE::Strain_Energy_Constraint:
        return FEA_MODULE_TYPE::Elasticity;
      default:
        return {};
    }
  }
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Simulation_Parameters_Topology_Optimization, Simulation_Parameters,
  optimization_options, nodal_density_flag, thick_condition_boundary,
  optimization_output_freq
)


// class Simulation_Parameters_Topology_Optimization : public Simulation_Parameters
// {
//  public:
//   Simulation_Parameters_Topology_Optimization(Implicit_Solver *solver_pointer);
//   virtual ~Simulation_Parameters_Topology_Optimization();
//   virtual void input();
//   virtual void FEA_module_setup();
//   virtual void apply_settings();
//   //==============================================================================
//   //   Mesh Variables
//   //==============================================================================

//   // --- Mesh regions and material fills ---
//   int NB; // number of boundary patch sets to tag
//   int NBD; //number of density boundary conditions

//   //Topology Optimization flags
//   bool topology_optimization_on, shape_optimization_on, nodal_density_flag, helmholtz_filter;
//   std::string multi_objective_structure;
  
//   //When on, all element nodes connected to a boundary condition patch will have their density constrained
//   bool thick_condition_boundary;

//   //method of moving asymptotes enabled for the optimization algorithm
//   bool mma_on;

//   //file output parameters
//   int optimization_output_freq;

//   //Topology Optimization parameters
//   real_t penalty_power, density_epsilon;

//   //pointer to Solver object (just used to consolidate error handling for now)
//   Implicit_Solver *solver_pointer_;

//   //volumes to hold density constant
  
//   //types of TO functions
//   enum function_type {OBJECTIVE, MULTI_OBJECTIVE_TERM, EQUALITY_CONSTRAINT, INEQUALITY_CONSTRAINT, VECTOR_EQUALITY_CONSTRAINT, VECTOR_INEQUALITY_CONSTRAINT};

//   //list of TO functions needed by problem
//   std::vector<std::string> TO_Module_List;
//   std::vector<function_type> TO_Function_Type;
//   std::vector<int> TO_Module_My_FEA_Module;
//   std::vector<int> Multi_Objective_Modules;
//   std::vector<real_t> Multi_Objective_Weights;
//   std::vector<std::vector<real_t>> Function_Arguments;
//   int nTO_modules, nmulti_objective_modules;
// };

#endif // end HEADER_H
