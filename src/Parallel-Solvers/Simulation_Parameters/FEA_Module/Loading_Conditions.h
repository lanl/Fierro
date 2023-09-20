#pragma once
#include "yaml-serializable.h"
#include "Simulation_Parameters/FEA_Module/Boundary_Conditions.h"

SERIALIZABLE_ENUM(LOADING_CONDITION_TYPE, surface_traction, surface_heat_flux, body_force)
SERIALIZABLE_ENUM(LOADING_SPECIFICATION, normal, coordinated)

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

struct Loading_Condition : Yaml::ValidatedYaml {
  BOUNDARY_TYPE surface;
  LOADING_CONDITION_TYPE condition_type;
  std::optional<double> plane_position {};
  std::optional<double> flux_value  {};
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
  surface, plane_position, condition_type,
  flux_value, component_x, component_y, component_z,
  specification
)