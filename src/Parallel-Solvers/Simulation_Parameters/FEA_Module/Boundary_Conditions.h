#pragma once
#include "yaml-serializable.h"

SERIALIZABLE_ENUM(BOUNDARY_TYPE, 
    x_plane,   // tag an x-plane
    y_plane,   // tag an y-plane
    z_plane,   // tag an z-plane
    cylinder,  // tag an cylindrical surface
    sphere,    // tag a spherical surface
    readFile   // read from a file
)

SERIALIZABLE_ENUM(CONDITION_TYPE, 
    fixed_displacement, 
    fixed_temperature, 
    fixed_velocity, 
    reflected
)

struct Boundary_Condition {
    BOUNDARY_TYPE surface;
    CONDITION_TYPE condition_type;
    double value;
    // TODO: these should probably just be value
    std::optional<double> temperature_value;
    std::optional<double> displacement_value;
    std::optional<double> plane_position;
};
YAML_ADD_REQUIRED_FIELDS_FOR(Boundary_Condition, surface, condition_type)
IMPL_YAML_SERIALIZABLE_FOR(Boundary_Condition,
  surface, value, condition_type, 
  temperature_value, plane_position, displacement_value
)