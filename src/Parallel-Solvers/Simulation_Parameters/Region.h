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
    double u,v,w;
    double speed;
    double sie;
    double den;
    bool extensive_energy_setting = false;
};

struct Region : Yaml::DerivedFields, Yaml::ValidatedYaml, mat_fill_t {
    size_t id;
    std::optional<double> sie;
    std::optional<double> ie;
    std::optional<double> u;
    std::optional<double> v;
    std::optional<double> w;

    void validate() {
      if ( (u.has_value() || v.has_value() || w.has_value())
          && velocity != VELOCITY_TYPE::cartesian) {
        std::cerr << "Warning: velocity values (u, v, w) ignored for velocity type: " << to_string(velocity) << "." << std::endl;
      }
    }

    void derive() {
      if (sie.has_value() == ie.has_value()){
        throw Yaml::ConfigurationException("Specify values for exactly one of: energy (ie) or specific energy (sie).");
      }
      
      mat_fill_t::u = u.value_or(0);
      mat_fill_t::v = v.value_or(0);
      mat_fill_t::w = w.value_or(0);

      if (ie.has_value()){
        mat_fill_t::extensive_energy_setting = true;
        sie = ie.value()  / (den * volume.get_volume());
      }

      mat_fill_t::sie = sie.value_or(0);
      // Remove ie now that we have set sie.
      // It doesn't really make sense to have both, 
      // and having both means we can't deserialize -> reserialize.
      ie = {};
    }
};
IMPL_YAML_SERIALIZABLE_FOR(Region, 
    id, volume, material_id,
    den, sie, velocity, u, v, w, ie,
    speed
)