#pragma once
#include "yaml-serializable.h"
#include <string>
#include <cmath>

SERIALIZABLE_ENUM(VOLUME_TYPE,
    global,   // tag every cell in the mesh
    box,      // tag all cells inside a box
    cylinder, // tag all cells inside a cylinder
    sphere    // tag all cells inside a sphere
)

SERIALIZABLE_ENUM(VELOCITY_TYPE,
    // uniform
    cartesian,          // cart velocity
    radial,             // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
    spherical,          // spherical
    // linear variation
    radial_linear,      // linear variation from 0,0,0
    spherical_linear,   // linear variation from 0,0,0
    // vortical initial conditions
    tg_vortex
)

struct mat_fill_t {
    std::string material_id;
    VOLUME_TYPE volume = VOLUME_TYPE::sphere;
    VELOCITY_TYPE velocity = VELOCITY_TYPE::cartesian;
    double radius1, radius2;
    double u,v,w;
    double speed;
    double sie;
    double den;
    
    // TODO: These aren't set from the YAML at all
    // They are only set in the test problems.
    double x1, x2, y1, y2, z1, z2;
    
    // KOKKOS_FUNCTION
    // bool contains(const double* elem_coords) { 
    //   double radius;

    //   switch(volume)
    //   {
    //     case VOLUME_TAG::global:
    //       return true;

    //     case VOLUME_TAG::box:
    //       return ( elem_coords[0] >= x1 && elem_coords[0] <= x2
    //             && elem_coords[1] >= y1 && elem_coords[1] <= y2
    //             && elem_coords[2] >= z1 && elem_coords[2] <= z2 );

    //     case VOLUME_TAG::cylinder:
    //       radius = sqrt( elem_coords[0]*elem_coords[0] +
    //                      elem_coords[1]*elem_coords[1] ); 
    //       return ( radius >= radius1
    //             && radius <= radius2 );

    //     case VOLUME_TAG::sphere:
    //       radius = sqrt( elem_coords[0]*elem_coords[0] +
    //                      elem_coords[1]*elem_coords[1] +
    //                      elem_coords[2]*elem_coords[2] );
    //       return ( radius >= radius1
    //             && radius <= radius2 );
        
    //     default:
    //       return false;
    //   }
    // }
};

struct Region : Yaml::DerivedFields, Yaml::ValidatedYaml, mat_fill_t {
    std::string id;
    std::optional<double> sie;
    std::optional<double> u;
    std::optional<double> v;
    std::optional<double> w;
    std::optional<double> radius1;
    std::optional<double> radius2;

    void validate() {
      if ( (u.has_value() || v.has_value() || w.has_value())
          && velocity != VELOCITY_TYPE::cartesian) {
        std::cerr << "Warning: velocity values (u, v, w) ignored for velocity type: " << to_string(velocity) << "." << std::endl;
      }
      if (volume == VOLUME_TAG::cylinder || volume == VOLUME_TAG::sphere) {
        if (!(radius1.has_value() && radius2.has_value())) {
          throw Yaml::ConfigurationException(
            "Both the inner and outer radii (radius1, radius2) must be specified for volume type "
              + to_string(volume)
          );
        }
      }
    }

    void derive() {
      if (volume == VOLUME_TAG::sphere) {
        if (radius2.has_value()) {
          if (sie.has_value())
            std::cerr << "Warning: sie value being overwritten" << std::endl;
          sie = (963.652344 * std::pow((1.2 / 30.0), 3)) / std::pow(radius2.value(), 3);
        }
      }

      mat_fill_t::radius1 = radius1.value_or(0);
      mat_fill_t::radius2 = radius2.value_or(0);
      mat_fill_t::sie = sie.value_or(0);
      mat_fill_t::u = u.value_or(0);
      mat_fill_t::v = v.value_or(0);
      mat_fill_t::w = w.value_or(0);
    }
};
IMPL_YAML_SERIALIZABLE_FOR(Region, 
    id, volume, material_id,
    den, sie, velocity, u, v, w, 
    radius1, radius2, speed
)