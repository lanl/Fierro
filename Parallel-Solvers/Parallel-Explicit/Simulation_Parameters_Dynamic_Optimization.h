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

#ifndef SIMULATION_PARAMETERS_DYNAMIC_OPTIMIZATION_H
#define SIMULATION_PARAMETERS_DYNAMIC_OPTIMIZATION_H

//#include "utilities.h"
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
};
IMPL_YAML_SERIALIZABLE_FOR(Optimization_Options, optimization_process, optimization_objective, constraints)

struct Simulation_Parameters_Dynamic_Optimization : public Simulation_Parameters {
  // --- Mesh regions and material fills ---
  int NB  = 6; // number of boundary patch sets to tag
  int NBD = 2; //number of density boundary conditions

  Optimization_Options optimization_options;
  
  //When on, all element nodes connected to a boundary condition patch will have their density constrained
  bool thick_condition_boundary = true;

  //file output parameters
  int optimization_output_freq;

  //Topology Optimization parameters
  real_t penalty_power = 3.0;

  
  //Non-serialized fields
  //list of TO functions needed by problem
  std::vector<TO_MODULE_TYPE> TO_Module_List {
    TO_MODULE_TYPE::Kinetic_Energy_Minimize
  };
  std::vector<FUNCTION_TYPE> TO_Function_Type {
    FUNCTION_TYPE::OBJECTIVE
  };
  std::vector<int> TO_Module_My_FEA_Module;
  std::vector<std::vector<int>> FEA_Module_My_TO_Modules;
  std::vector<std::vector<real_t>> Function_Arguments;

  // TODO: Implement a real structure here.
  // Then we can support this stuff.
  // std::vector<int> Multi_Objective_Modules;
  // std::vector<real_t> Multi_Objective_Weights;
  int nTO_modules, nmulti_objective_modules;

  //Topology Optimization flags
  bool topology_optimization_on = false;
  bool shape_optimization_on = false;
  bool nodal_density_flag;;

  void derive() {
    shape_optimization_on = optimization_options.optimization_process == OPTIMIZATION_PROCESS::shape_optimization;
    topology_optimization_on = optimization_options.optimization_process == OPTIMIZATION_PROCESS::topology_optimization;

    if (optimization_options.constraints.size() > 0) {
      TO_Module_List.push_back(TO_MODULE_TYPE::Kinetic_Energy_Minimize);
      TO_Function_Type.push_back(FUNCTION_TYPE::OBJECTIVE);
      Function_Arguments.push_back({});
    }
    TO_Module_List.resize(optimization_options.constraints.size());
    TO_Function_Type.resize(optimization_options.constraints.size());
    Function_Arguments.resize(optimization_options.constraints.size());
    for (size_t i = 0; i < optimization_options.constraints.size(); i++) {
      auto constraint = optimization_options.constraints[i];

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
  
    // Take a pass first to ensure that all necessary modules are loaded.
    for (auto to_module : TO_Module_List)
      ensure_module(get_TO_module_dependency(to_module)); 

    // Now we allocate the vectors for each of the currently identified modules.
    TO_Module_My_FEA_Module.resize(TO_Module_List.size());
    FEA_Module_My_TO_Modules.resize(FEA_Modules_List.size());
    
    // Finally we can set up the maps.
    // 
    // TODO: This should really use two `std::map<int, std::set<int>>`s
    // instead of this vector stuff.
    for (size_t to_index = 0; to_index < TO_Module_List.size(); to_index++) {
      auto to_module = TO_Module_List[to_index];
      size_t fea_index = find_module(get_TO_module_dependency(to_module));

      TO_Module_My_FEA_Module[to_index] = fea_index;
      FEA_Module_My_TO_Modules[fea_index].push_back(to_index);
    }
  }
  void validate() { }

  FEA_MODULE_TYPE get_TO_module_dependency(TO_MODULE_TYPE type) {
    switch (type) {
      case TO_MODULE_TYPE::Kinetic_Energy_Minimize:
        return FEA_MODULE_TYPE::SGH;
      case TO_MODULE_TYPE::Heat_Capacity_Potential_Minimize:
        return FEA_MODULE_TYPE::Heat_Conduction;
      case TO_MODULE_TYPE::Heat_Capacity_Potential_Constraint:
        return FEA_MODULE_TYPE::Heat_Conduction;
      case TO_MODULE_TYPE::Mass_Constraint:
        return FEA_MODULE_TYPE::Inertial;
      case TO_MODULE_TYPE::Moment_of_Inertia_Constraint:
        return FEA_MODULE_TYPE::Inertial;
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

// //forward declare
// class Solver;

// class Simulation_Parameters_Dynamic_Optimization : public Simulation_Parameters
// {
//  public:
//   Simulation_Parameters_Dynamic_Optimization(Solver *solver_pointer);
//   virtual ~Simulation_Parameters_Dynamic_Optimization();
//   virtual void input();
//   virtual void apply_settings();
//   virtual void FEA_module_setup();
//   //==============================================================================
//   //   Mesh Variables
//   //==============================================================================

//   // --- Mesh regions and material fills ---
//   int NB; // number of boundary patch sets to tag
//   int NBD; //number of density boundary conditions

//   //Topology Optimization flags
//   bool topology_optimization_on, shape_optimization_on, nodal_density_flag;
  
//   //When on, all element nodes connected to a boundary condition patch will have their density constrained
//   bool thick_condition_boundary;

//   //file output parameters
//   int optimization_output_freq;

//   //Topology Optimization parameters
//   real_t penalty_power;

//   //pointer to Solver object (just used to consolidate error handling for now)
//   Solver *solver_pointer_;

//   //volumes to hold density constant
  
//   //types of TO functions
//   enum function_type {OBJECTIVE, MULTI_OBJECTIVE_TERM, EQUALITY_CONSTRAINT, INEQUALITY_CONSTRAINT, VECTOR_EQUALITY_CONSTRAINT, VECTOR_INEQUALITY_CONSTRAINT};

//   //list of TO functions needed by problem
//   std::vector<std::string> TO_Module_List;
//   std::vector<function_type> TO_Function_Type;
//   std::vector<int> TO_Module_My_FEA_Module;
//   std::vector<std::vector<int>> FEA_Module_My_TO_Modules;
//   std::vector<int> Multi_Objective_Modules;
//   std::vector<real_t> Multi_Objective_Weights;
//   std::vector<std::vector<real_t>> Function_Arguments;
//   int nTO_modules, nmulti_objective_modules;
// };

#endif // end HEADER_H
