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
#ifndef SIMULATION_PARAMETERS_BASE_H
#define SIMULATION_PARAMETERS_BASE_H

#include "matar.h"
#include <cmath>
#include <stdexcept>

#include "mesh.h"
#include "state.h"
#include "yaml-serializable.h"
#include "Simulation_Parameters.h"
#include "material_t.h"

using namespace mtr;


SERIALIZABLE_ENUM(VOLUME_TAG,
    global,   // tag every cell in the mesh
    box,      // tag all cells inside a box
    cylinder, // tag all cells inside a cylinder
    sphere    // tag all cells inside a sphere
)

SERIALIZABLE_ENUM(TIME_OUTPUT_LEVEL,
    none,   // output no time sequence information of forward solve
    low,      // output just final timestep of forward solve
    high, // output time sequence for forward solve
    extreme    // output time sequence for all phases of solver
)

SERIALIZABLE_ENUM(BOUNDARY_HYDRO_CONDITION,
    fixed,        // zero velocity
    reflected,    // reflected or wall condition
    velocity,     // constant velocity
    pressure,     // constant pressure
    acceleration, // constant acceleration
    contact       // contact surface
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

struct Time_Variables : Yaml::DerivedFields {
    double time_final = 1.0;
    double time_initial = 0.0;
    double dt_min     = 1e-8;
    double dt_max     = 1e-2;
    double dt_start   = 1e-5;
    double dt_cfl     = 0.4;
    int cycle_stop    = 2000000;
    double fuzz       = 1e-16; // machine precision
    double tiny       = 1e-12; // very very small (between real_t and single)
    double small      = 1e-8;  // single precision
    TIME_OUTPUT_LEVEL output_time_sequence_level = TIME_OUTPUT_LEVEL::high;

    // Non-serialized Fields
    double dt;
    void derive() {
      dt = dt_start;
    }
};
IMPL_YAML_SERIALIZABLE_FOR(Time_Variables, output_time_sequence_level,
  time_initial, time_final, dt_min, dt_max, dt_start, dt_cfl,
  cycle_stop, fuzz, tiny, small
)

struct mat_fill_t {
    size_t mat_id;
    VOLUME_TAG volume = VOLUME_TAG::sphere;
    VELOCITY_TYPE velocity = VELOCITY_TYPE::cartesian;
    double radius1, radius2;
    double u,v,w;
    double speed;
    double sie;
    double den;
    
    // TODO: These aren't set from the YAML at all
    // They are only set in the test problems.
    double x1, x2, y1, y2, z1, z2;
    
    KOKKOS_FUNCTION
    bool contains(const double* elem_coords) { 
      double radius;

      switch(volume)
      {
        case VOLUME_TAG::global:
          return true;

        case VOLUME_TAG::box:
          return ( elem_coords[0] >= x1 && elem_coords[0] <= x2
                && elem_coords[1] >= y1 && elem_coords[1] <= y2
                && elem_coords[2] >= z1 && elem_coords[2] <= z2 );

        case VOLUME_TAG::cylinder:
          radius = sqrt( elem_coords[0]*elem_coords[0] +
                         elem_coords[1]*elem_coords[1] ); 
          return ( radius >= radius1
                && radius <= radius2 );

        case VOLUME_TAG::sphere:
          radius = sqrt( elem_coords[0]*elem_coords[0] +
                         elem_coords[1]*elem_coords[1] +
                         elem_coords[2]*elem_coords[2] );
          return ( radius >= radius1
                && radius <= radius2 );
        
        default:
          return false;
      }
    }
};
struct MaterialFill : Yaml::DerivedFields, Yaml::ValidatedYaml, mat_fill_t {
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
IMPL_YAML_SERIALIZABLE_FOR(
  MaterialFill, volume, mat_id, 
  den, sie, velocity, u, v, w, 
  radius1, radius2, speed
)

struct boundary_t {
    BOUNDARY_HYDRO_CONDITION condition_type = BOUNDARY_HYDRO_CONDITION::fixed;
    BOUNDARY_TAG surface = BOUNDARY_TAG::sphere;
    double value;
    double u, v, w;
    
    KOKKOS_FUNCTION
    size_t planar_surface_index() {
      switch (surface) {
        case BOUNDARY_TAG::x_plane:
          return 0;
        case BOUNDARY_TAG::y_plane:
          return 1;
        case BOUNDARY_TAG::z_plane:
          return 2;
        default:
          // Not sure what to do about this error case, since we can't
          // throw on the device.
          return -1;
      }
    }
};

struct Boundary : Yaml::ValidatedYaml, Yaml::DerivedFields, boundary_t {
    std::string id;
    std::optional<double> u;
    std::optional<double> v;
    std::optional<double> w;

    void derive() {
      boundary_t::u = u.value_or(0);
      boundary_t::v = v.value_or(0);
      boundary_t::w = w.value_or(0);
    }

    void validate() {
      if ( (u.has_value() || v.has_value() || w.has_value())
          && condition_type != BOUNDARY_HYDRO_CONDITION::velocity) {
        std::cerr << "Warning: velocity values (u, v, w) ignored for boundary condition type " << to_string(condition_type) << "." << std::endl;
      }

      if (condition_type == BOUNDARY_HYDRO_CONDITION::reflected) {
        switch (surface) {
          case BOUNDARY_TAG::x_plane:
          case BOUNDARY_TAG::y_plane:
          case BOUNDARY_TAG::z_plane:
            break;
          default:
            throw Yaml::ConfigurationException(
              "Invalid surface type `" + to_string(surface) + 
              "` with boundary condition type of `" + to_string(condition_type) + "`."
            );
        }
      }
    }
};
IMPL_YAML_SERIALIZABLE_FOR(Boundary, id, surface, value, condition_type, u, v, w)

struct Graphics_Options : Yaml::DerivedFields {
  bool output_velocity_flag = true;
  bool output_strain_flag   = true;
  bool output_stress_flag   = false;

  bool displaced_mesh_flag  = true;
  bool strain_max_flag      = false;

  size_t graphics_cyc_ival  = 1000000;
  double graphics_dt_ival   = 0.25;

  // Non-serialized
  size_t graphics_id = 0;
  double graphics_time;  // the times for writing graphics dump
  CArray <double> graphics_times;

  void derive() {
    graphics_time     = graphics_dt_ival;
    graphics_times    = CArray<double>(2000);
    graphics_times(0) = 0.0;
  }
};
IMPL_YAML_SERIALIZABLE_FOR(Graphics_Options, 
  output_velocity_flag, output_stress_flag, output_strain_flag,
  strain_max_flag, displaced_mesh_flag, graphics_cyc_ival, graphics_dt_ival
)

#endif
