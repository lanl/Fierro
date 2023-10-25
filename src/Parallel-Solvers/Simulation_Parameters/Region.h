#pragma once
#include "yaml-serializable.h"
#include <string>
#include <cmath>
#include "Simulation_Parameters/Geometry.h"

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
    size_t material_id;
    Volume volume;
    VELOCITY_TYPE velocity = VELOCITY_TYPE::cartesian;
    double radius1, radius2;
    double u,v,w;
    double speed;
    double sie;
    double den;
};

struct Region : Yaml::DerivedFields, Yaml::ValidatedYaml, mat_fill_t {
    size_t id;
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
      if (volume.type == VOLUME_TYPE::cylinder || volume.type == VOLUME_TYPE::sphere) {
        if (!(radius1.has_value() && radius2.has_value())) {
          throw Yaml::ConfigurationException(
            "Both the inner and outer radii (radius1, radius2) must be specified for volume type "
              + to_string(volume.type)
          );
        }
      }
    }

    void derive() {
      if (volume.type == VOLUME_TYPE::sphere) {
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