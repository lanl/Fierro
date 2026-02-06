#pragma once

#include "yaml-serializable.h"
#include <string>
#include <cmath>
#include "Simulation_Parameters/Geometry.h"

struct Optimization_Bound_Constraint_Region : Yaml::DerivedFields, Yaml::ValidatedYaml {
    size_t id;
    Volume volume;
    VELOCITY_TYPE velocity = VELOCITY_TYPE::cartesian;
    std::optional<double> lower_density_bound;
    std::optional<double> upper_density_bound;
    double set_lower_density_bound;
    double set_upper_density_bound;

    void validate() {
      if (lower_density_bound.has_value()){
        if(lower_density_bound.value()<0)
        throw Yaml::ConfigurationException("Specify values for lower relative density greater than or equal to 0 but less than or equal to 1");
      }

      if (upper_density_bound.has_value()){
        if(upper_density_bound.value()<0)
        throw Yaml::ConfigurationException("Specify values for upper relative density greater than 0 but less than or equal to 1");
      }
    }

    void derive() {
      set_lower_density_bound = lower_density_bound.value_or(0);
      set_upper_density_bound = upper_density_bound.value_or(1);
    }
};
IMPL_YAML_SERIALIZABLE_FOR(Optimization_Bound_Constraint_Region, 
    id, volume, lower_density_bound, upper_density_bound
)