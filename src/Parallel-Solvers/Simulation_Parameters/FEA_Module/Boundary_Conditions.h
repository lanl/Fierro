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

SERIALIZABLE_ENUM(BOUNDARY_CONDITION_TYPE, 
    displacement, 
    temperature, 
    velocity,
    reflected,
    fixed_position,
    pressure,
    acceleration,
    contact
)

struct Surface {
    BOUNDARY_TYPE type;
    double plane_position;
    
    KOKKOS_FUNCTION
    size_t planar_surface_index() {
      switch (type) {
        case BOUNDARY_TYPE::x_plane:
          return 0;
        case BOUNDARY_TYPE::y_plane:
          return 1;
        case BOUNDARY_TYPE::z_plane:
          return 2;
        default:
          // Not sure what to do about this error case, since we can't
          // throw on the device.
          return -1;
      }
    }
};
YAML_ADD_REQUIRED_FIELDS_FOR(Surface, type)
IMPL_YAML_SERIALIZABLE_FOR(Surface, type, plane_position)


/**
 * Device struct representation of Boundary_Condition.
*/
struct boundary_t {
    BOUNDARY_CONDITION_TYPE type = BOUNDARY_CONDITION_TYPE::fixed_position;
    Surface surface;
    double value = 0;
    double u = 0;
    double v = 0;
    double w = 0;
};

// TODO: Consider using type discrimination to get different validation 
// for different types of boundaries.
struct Boundary_Condition : boundary_t { };
YAML_ADD_REQUIRED_FIELDS_FOR(Boundary_Condition, surface, type)
IMPL_YAML_SERIALIZABLE_FOR(Boundary_Condition, surface, value, type, u, v, w)
