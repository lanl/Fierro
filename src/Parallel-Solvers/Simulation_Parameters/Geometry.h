#pragma once

#include "yaml-serializable.h"
#include "matar.h"
#include <assert.h>

using namespace mtr;

SERIALIZABLE_ENUM(VOLUME_TYPE,
    global,   // tag every cell in the mesh
    box,      // tag all cells inside a box
    cylinder, // tag all cells inside a cylinder
    sphere    // tag all cells inside a sphere
)

struct Volume : Yaml::ValidatedYaml {
    VOLUME_TYPE type = VOLUME_TYPE::sphere;
    double radius1 = 0;
    double radius2 = 0;

    double x1, x2, y1, y2, z1, z2;

    KOKKOS_FUNCTION
    bool contains(const double* elem_coords) { 
      double radius;

      switch(type) {
      case VOLUME_TYPE::global:
        return true;

      case VOLUME_TYPE::box:
        return ( elem_coords[0] >= x1 && elem_coords[0] <= x2
              && elem_coords[1] >= y1 && elem_coords[1] <= y2
              && elem_coords[2] >= z1 && elem_coords[2] <= z2 );

      case VOLUME_TYPE::cylinder:
        radius = sqrt( elem_coords[0]*elem_coords[0] +
                        elem_coords[1]*elem_coords[1] ); 
        return ( radius >= radius1
              && radius <= radius2 );

      case VOLUME_TYPE::sphere:
        radius = sqrt( elem_coords[0]*elem_coords[0] +
                        elem_coords[1]*elem_coords[1] +
                        elem_coords[2]*elem_coords[2] );
        return ( radius >= radius1
              && radius <= radius2 );

      default:
        return false;
      }
    }

    KOKKOS_FUNCTION
    double get_volume() {
      switch(type) {
      case VOLUME_TYPE::global:
        assert(0); // Cannot evaluate the volume fo an unbounded region.

      case VOLUME_TYPE::box:
        return (x2 - x1) * (y2 - y1) * (z2 - z1);

      case VOLUME_TYPE::cylinder:
        // 2D cylinder
        return M_PI * (std::pow(radius2, 2) - std::pow(radius1, 2));

      case VOLUME_TYPE::sphere:
        return (4.0 / 3.0) * M_PI * (std::pow(radius2, 3) - std::pow(radius1, 3));
      
      default:
        assert(0); // Unsupported volume type.
        return 0; // Compiler complains that execution can fall out of void func.
      }
    }
    
    void validate() {
      if (type == VOLUME_TYPE::cylinder || type == VOLUME_TYPE::sphere) {
        if (radius2 < 0) {
          throw Yaml::ConfigurationException(
            "The inner radius cannot be less than 0 for volume type: "
              + to_string(type)
          );
        }
        if (radius2 <= radius1) {
          throw Yaml::ConfigurationException(
            "The outer radius must be greater than the inner radius for volume type: "
              + to_string(type)
          );
        }
      }
    }
};
IMPL_YAML_SERIALIZABLE_FOR(Volume, type, radius1, radius2, x1, x2, y1, y2, z1, z2)