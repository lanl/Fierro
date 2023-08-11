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
#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H

#include "yaml-serializable.h"
#include "configuration-validation.h"
#include <set>
#include <vector>
#include <string>
#include <stdio.h>
#include <optional>
#include <algorithm>
#include <filesystem>

template<typename T>
inline void validate_unique_vector(const std::vector<T>& vec, std::string err_msg) {
  std::set<T> seen;
  for (auto v : vec) {
    if (seen.find(v) != seen.end()) {
      throw Yaml::ConfigurationException(
        err_msg + to_string(v)
      );
      seen.insert(v);
    }
  }
}

SERIALIZABLE_ENUM(SOLVER_TYPE, SGH, Implicit)
SERIALIZABLE_ENUM(MESH_FORMAT,
    ensight,
    tecplot,
    vtk,
    ansys_dat
)

SERIALIZABLE_ENUM(OUTPUT_FORMAT, vtk, vtu)
SERIALIZABLE_ENUM(TIMER_VERBOSITY, standard, thorough)

SERIALIZABLE_ENUM(FEA_MODULE_TYPE,
  Elasticity,
  Heat_Conduction,
  SGH,
  Inertial,
  Thermo_Elasticity,
  Eulerian
)

SERIALIZABLE_ENUM(ELEMENT_TYPE, 
  quad4, quad8, quad12,
  hex8, hex20, hex32
)

struct Input_Options : Yaml::ValidatedYaml, Yaml::DerivedFields {
  std::string mesh_file_name;
  MESH_FORMAT mesh_file_format;

  ELEMENT_TYPE element_type = ELEMENT_TYPE::hex8;
  bool zero_index_base = false;

  // Non-serialized fields
  int words_per_line;
  int elem_words_per_line;

  /**
   * Determine a couple of file parsing parameters from the specified filetype.
  */
  void derive() {
    if (mesh_file_format == MESH_FORMAT::ansys_dat) {
      words_per_line = 4;
      elem_words_per_line = 11;
    } else {
      switch (mesh_file_format) {
        case MESH_FORMAT::ensight:
          words_per_line = 1;
          break;
        case MESH_FORMAT::vtk:
        case MESH_FORMAT::tecplot:
          words_per_line = 3;
          break;
        default:
          break;
      }

      switch (element_type) {
        case ELEMENT_TYPE::hex8:
          elem_words_per_line = 8;
          break;
        case ELEMENT_TYPE::quad4:
          elem_words_per_line = 4;
          break;
        // TODO: Implement handling for other element types
        default:
          throw Yaml::ConfigurationException("Unsupported element type `" + to_string(element_type) + "`.");
          break;
      }
    }

    mesh_file_name = std::filesystem::absolute(mesh_file_name).string();
  }
  
  /**
   * Ensures that the provided filepath is valid.
  */
  void validate() {
    Yaml::validate_filepath(mesh_file_name);
  }
};
IMPL_YAML_SERIALIZABLE_FOR(Input_Options, mesh_file_name, mesh_file_format, element_type, zero_index_base)

struct Output_Options {
  int graphics_step_frequency;
  double graphics_step;
  OUTPUT_FORMAT output_file_format;
  size_t max_num_user_output_vars=0;
};
IMPL_YAML_SERIALIZABLE_FOR(Output_Options, graphics_step_frequency, graphics_step, output_file_format, max_num_user_output_vars)


SERIALIZABLE_ENUM(BOUNDARY_TAG, 
    x_plane,   // tag an x-plane
    y_plane,   // tag an y-plane
    z_plane,   // tag an z-plane
    cylinder,  // tag an cylindrical surface
    sphere,    // tag a spherical surface
    readFile   // read from a file
)

SERIALIZABLE_ENUM(BOUNDARY_FEA_CONDITION, fixed_displacement, fixed_temperature)
SERIALIZABLE_ENUM(LOADING_CONDITION_TYPE, surface_traction, surface_heat_flux)
SERIALIZABLE_ENUM(LOADING_SPECIFICATION, normal, coordinated)

struct Loading_Condition : Yaml::ValidatedYaml {
  std::string id;
  BOUNDARY_TAG surface;
  LOADING_CONDITION_TYPE condition_type;
  std::optional<double> plane_position {};
  std::optional<double> flux_value {};
  std::optional<double> component_x {};
  std::optional<double> component_y {};
  std::optional<double> component_z {};
  std::optional<LOADING_SPECIFICATION> specification {};

  void validate_surface_heat_flux() {
    std::string type_name = to_string(LOADING_CONDITION_TYPE::surface_heat_flux);
    if (component_x.has_value() || component_y.has_value() || component_z.has_value())
      throw Yaml::ConfigurationException("Do not specify xyz components for " + type_name);

    if (!flux_value.has_value())
      throw Yaml::ConfigurationException("`flux_value` required for " + type_name);
    if (!specification.has_value())
      throw Yaml::ConfigurationException("`specification` required for " + type_name);
  }

  void validate_surface_traction() {
    std::string type_name = to_string(LOADING_CONDITION_TYPE::surface_traction);
    
    if (flux_value.has_value())
      throw Yaml::ConfigurationException("Do not specify `flux_value` for " + type_name);
    if (specification.has_value())
      throw Yaml::ConfigurationException("Do not provide `specification` for " + type_name);
    
    if (!component_x.has_value() || !component_y.has_value() || !component_z.has_value())
      throw Yaml::ConfigurationException("`component_[x,y,z]` values required for " + type_name);
  }

  void validate() {
    switch (condition_type) {
      case LOADING_CONDITION_TYPE::surface_heat_flux:
        validate_surface_heat_flux();
        break;
      case LOADING_CONDITION_TYPE::surface_traction:
        validate_surface_traction();
        break;
      default:
        throw Yaml::ConfigurationException(
          "Unhandled loading condition type: " + to_string(condition_type)
        );
    }
  }
};
YAML_ADD_REQUIRED_FIELDS_FOR(Loading_Condition, condition_type, surface)
IMPL_YAML_SERIALIZABLE_FOR(Loading_Condition, 
  id, surface, plane_position, condition_type,
  flux_value, component_x, component_y, component_z,
  specification
)

struct FEA_Boundary_Condition {
    std::string id;
    BOUNDARY_TAG surface;
    double value;
    BOUNDARY_FEA_CONDITION condition_type;
    // TODO: these should probably just be value
    std::optional<double> temperature_value;
    std::optional<double> displacement_value;
    std::optional<double> plane_position;
};
YAML_ADD_REQUIRED_FIELDS_FOR(FEA_Boundary_Condition, surface, condition_type)
IMPL_YAML_SERIALIZABLE_FOR(FEA_Boundary_Condition, 
  id, surface, value, condition_type, 
  temperature_value, plane_position, displacement_value
)

struct FEA_Module_Config {
  FEA_MODULE_TYPE type;
  std::vector<FEA_Boundary_Condition> boundary_conditions;
  std::vector<Loading_Condition> loading_conditions;
};
IMPL_YAML_SERIALIZABLE_FOR(FEA_Module_Config, type, boundary_conditions, loading_conditions)

struct Simulation_Parameters : Yaml::ValidatedYaml, Yaml::DerivedFields {
  SOLVER_TYPE solver_type;
  TIMER_VERBOSITY timer_output_level;
  int num_dims = 3;
  bool restart_file = false;
  Input_Options input_options;
  Output_Options output_options;
  bool report_runtime = true;

  // fea_modules holds configuration for a subset of
  // the modules listed in FEA_Modules_List
  std::vector<FEA_Module_Config> fea_modules;

  // Non-serialized fields
  int p_order = 0;
  double unit_scaling = 1.0;
  std::vector<FEA_MODULE_TYPE> FEA_Modules_List;
  std::set<FEA_MODULE_TYPE> fea_module_must_read;
  std::vector<bool> enable_inertia_center {false, false, false};
  std::vector<double> moment_of_inertia_center {0.0, 0.0, 0.0};


  void derive() {
    for (auto& spec : fea_modules) {
      FEA_Modules_List.push_back(spec.type);
    }
  }

  void validate_element_type() {
    auto et = input_options.element_type;
    bool invalid_et = (et == ELEMENT_TYPE::quad4 || et == ELEMENT_TYPE::quad8 || et == ELEMENT_TYPE::quad12) && num_dims == 3;
    invalid_et = invalid_et ||
      (et == ELEMENT_TYPE::hex8 || et == ELEMENT_TYPE::hex20 || et == ELEMENT_TYPE::hex32) && num_dims == 2;
    
    if (invalid_et) 
      throw Yaml::ConfigurationException(
        "Invalid element type " + to_string(et) + " for number of dimensions " + std::to_string(num_dims)
      );
  }
  void validate() {
    validate_element_type();
    validate_unique_vector(FEA_Modules_List, "Duplicate FEA Module Found: ");
  }

  /**
   * Given a module type, return the optionally present
   * configuration associated with it.
  */
  std::optional<FEA_Module_Config> get_module_config(FEA_MODULE_TYPE type) {
    for (auto spec : fea_modules) 
      if (spec.type == type) return spec;
    return {};
  }

  /**
   * If a module with the provided type is not present,
   * add it to the list of modules without any configuration.
  */
  size_t ensure_module(FEA_MODULE_TYPE type) {
    size_t i = find_module(type);
    if (i == FEA_Modules_List.size())
      FEA_Modules_List.push_back(type);
    return i;
  }

  /**
   * Ensure that the module is provided.
   * 
   * If module configuration of this type is not found,
   * add it to the module configurations.
   * 
   * If there is one present already, don't do anything.
  */
  size_t ensure_module(FEA_Module_Config default_spec) {
    for (size_t i = 0; i < fea_modules.size(); i++) {
      if (fea_modules[i].type == default_spec.type) 
        return i;
    }

    fea_modules.push_back(default_spec);
    return ensure_module(default_spec.type);
  }

  /**
   * Checks to see if a module of this type is already loaded,
   * with or without additional configuration present.
  */
  bool has_module(FEA_MODULE_TYPE type) {
    return find_module(type) != FEA_Modules_List.size();
  }

  /**
   * Returns the index of the module in the list of modules.
   * Returns FEA_Modules_List.size() if it isn't present.
  */
  size_t find_module(FEA_MODULE_TYPE type) {
    size_t i = 0;
    for (; i < FEA_Modules_List.size(); i++) {
      if (FEA_Modules_List[i] == type)
        break;
    }
    return i;
  }
};
YAML_ADD_REQUIRED_FIELDS_FOR(Simulation_Parameters,
  solver_type, num_dims, input_options
)
IMPL_YAML_SERIALIZABLE_FOR(Simulation_Parameters, 
  solver_type, restart_file, input_options, timer_output_level,
  num_dims, output_options, fea_modules,
  report_runtime
)
#endif
