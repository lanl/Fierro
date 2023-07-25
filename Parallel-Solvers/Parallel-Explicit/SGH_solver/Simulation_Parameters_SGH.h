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

#ifndef SIMULATION_PARAMETERS_SGH_H
#define SIMULATION_PARAMETERS_SGH_H

#include "matar.h"
#include <cmath>
#include <stdexcept>
#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>

#include "mesh.h"
#include "state.h"
#include "yaml-serializable.h"
#include "Simulation_Parameters.h"
#include "user_material_functions.h"

using namespace mtr;

// TODO: This should be in some header or something.
//eos forward declaration
typedef void eos_function_type(
  const DViewCArrayKokkos <double> &elem_pres,
  const DViewCArrayKokkos <double> &elem_stress,
  const size_t elem_gid,
  const size_t mat_id,
  const DViewCArrayKokkos <double> &elem_state_vars,
  const DCArrayKokkos <double> &global_vars,
  const DViewCArrayKokkos <double> &elem_sspd,
  const double den,
  const double sie);
KOKKOS_FUNCTION eos_function_type ideal_gas;

typedef void strength_function_type(
  const DViewCArrayKokkos <double> &elem_pres,
  const DViewCArrayKokkos <double> &elem_stress,
  const size_t elem_gid,
  const size_t mat_id,
  const DViewCArrayKokkos <double> &elem_state_vars,
  const DCArrayKokkos <double> &global_vars,
  const DViewCArrayKokkos <double> &elem_sspd,
  const double den,
  const double sie,
  const ViewCArrayKokkos <double> &vel_grad,
  const ViewCArrayKokkos <size_t>  &elem_node_gids,
  const DViewCArrayKokkos <double> &node_coords,
  const DViewCArrayKokkos <double> &node_vel,
  const double vol,
  const double dt,
  const double alpha,
  const size_t cycle);
KOKKOS_FUNCTION strength_function_type user_strength_model;


SERIALIZABLE_ENUM(VOLUME_TAG,
    global,   // tag every cell in the mesh
    box,      // tag all cells inside a box
    cylinder, // tag all cells inside a cylinder
    sphere    // tag all cells inside a sphere
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

SERIALIZABLE_ENUM(EOS_MODEL, ideal_gas, user_eos_model)
SERIALIZABLE_ENUM(STRENGTH_MODEL, user_strength_model)
SERIALIZABLE_ENUM(STRENGTH_TYPE, hypo, hyper)
SERIALIZABLE_ENUM(STRENGTH_SETUP, input, user_input)
SERIALIZABLE_ENUM(RUN_LOCATION, host, device)

struct Time_Variables {
    double time_final = 1.0;
    double dt_min     = 1e-8;
    double dt_max     = 1e-2;
    double dt_start   = 1e-5;
    double dt         = dt_start;
    double dt_cfl     = 0.4;
    int cycle_stop    = 2000000;
    double fuzz       = 1e-16; // machine precision
    double tiny       = 1e-12; // very very small (between real_t and single)
    double small      = 1e-8;  // single precision
};
IMPL_YAML_SERIALIZABLE_FOR(Time_Variables, 
  time_final, dt_min, dt_max, dt_start, dt_cfl, dt,
  cycle_stop, fuzz, tiny, small
)

struct material_t : Yaml::DerivedFields {
  EOS_MODEL eos_model_type           = EOS_MODEL::ideal_gas;
  STRENGTH_TYPE strength_type        = STRENGTH_TYPE::hypo;
  STRENGTH_SETUP strength_setup      = STRENGTH_SETUP::user_input;
  STRENGTH_MODEL strength_model_type = STRENGTH_MODEL::user_strength_model;
  RUN_LOCATION strength_run_location = RUN_LOCATION::host;

  double q1;
  double q2;
  double q1ex;
  double q2ex;
  std::vector<double> state_vars;
  std::vector<double> global_vars;

  // Non-serialized fields
  size_t num_state_vars = 0;
  size_t num_global_vars = 0;
      
  eos_function_type* eos_model = NULL;
  strength_function_type* strength_model = NULL;

  void derive_function_pointers() {
    switch (eos_model_type) {
      case EOS_MODEL::ideal_gas:
        eos_model = ideal_gas;
        break;
      default:
        throw Yaml::ConfigurationException("Unsupported EOS Model Type " + to_string(eos_model_type));
        break;
    }

    switch (strength_model_type) {
      case STRENGTH_MODEL::user_strength_model:
        strength_model = user_strength_model;
        break;
      default:
        strength_model = NULL;
        break;
    }
  }

  void derive() {
    num_state_vars = state_vars.size();
    num_global_vars = global_vars.size();
    derive_function_pointers();
  }
};
IMPL_YAML_SERIALIZABLE_FOR(material_t, 
  eos_model_type, strength_model_type, strength_type, strength_setup, 
  strength_run_location,
  q1, q2, q1ex, q2ex, 
  state_vars, global_vars
)

struct mat_fill_t : Yaml::DerivedFields, Yaml::ValidatedYaml {
    VOLUME_TAG volume = VOLUME_TAG::sphere;
    size_t mat_id;
    double den;
    std::optional<double> sie;
    VELOCITY_TYPE velocity = VELOCITY_TYPE::cartesian;
    std::optional<double> u;
    std::optional<double> v;
    std::optional<double> w;
    double speed;

    std::optional<double> radius1;
    std::optional<double> radius2;

    // Non-serialized
    
    // TODO: These aren't set from the YAML at all
    // They are only set in the test problems.
    double x1;
    double x2;
    double y1;
    double y2;
    double z1;
    double z2;

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
    }

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
          return ( radius >= radius1.value()
                && radius <= radius2.value() );

        case VOLUME_TAG::sphere:
          radius = sqrt( elem_coords[0]*elem_coords[0] +
                         elem_coords[1]*elem_coords[1] +
                         elem_coords[2]*elem_coords[2] );
          return ( radius >= radius1.value()
                && radius <= radius2.value() );
        
        default:
          throw std::runtime_error("Unsupported volume type: " + to_string(volume));
      }
    }
};
IMPL_YAML_SERIALIZABLE_FOR(
  mat_fill_t, volume, mat_id, 
  den, sie, velocity, u, v, w, 
  radius1, radius2, speed
)

struct boundary_t : Yaml::ValidatedYaml {
    std::string id;
    BOUNDARY_TAG surface = BOUNDARY_TAG::sphere;
    double value;
    BOUNDARY_HYDRO_CONDITION condition_type = BOUNDARY_HYDRO_CONDITION::fixed;
    std::optional<double> u;
    std::optional<double> v;
    std::optional<double> w;

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

    size_t planar_surface_index() {
      switch (surface) {
        case BOUNDARY_TAG::x_plane:
          return 0;
        case BOUNDARY_TAG::y_plane:
          return 1;
        case BOUNDARY_TAG::z_plane:
          return 2;
        default:
          throw std::runtime_error("Attempted to get surface index with invalid boundary type");
      }
    }
};
IMPL_YAML_SERIALIZABLE_FOR(boundary_t, id, surface, value, condition_type, u, v, w)

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

struct Simulation_Parameters_SGH : Simulation_Parameters {
  Time_Variables time_variables;
  std::vector<material_t> material_options;
  std::vector<mat_fill_t> region_options;
  std::vector<boundary_t> boundary_conditions;
  Graphics_Options graphics_options;

  bool gravity_flag   = false;
  bool report_runtime = true;

  size_t rk_num_stages = 2;
  int NB   = 6; // number of boundaries
  int NBSF = 4; //number of surface density force conditions
  int NBV  = 2; //number of surface sets used to specify a fixed displacement on nodes belonging to respective surfaces

  //Non-serialized fields
  size_t max_num_state_vars;
  size_t max_num_global_vars;
  size_t rk_num_bins;
  double time_value = 0.0;
  DCArrayKokkos<double> state_vars;
  DCArrayKokkos<double> global_vars;

  DCArrayKokkos<mat_fill_t> mat_fill;
  DCArrayKokkos<material_t> material;
  DCArrayKokkos<boundary_t> boundary; 
  std::vector<double> gravity_vector {9.81, 0., 0.};

  void init_material_variable_arrays(size_t nstate_vars, size_t nglobal_vars) {
    state_vars = DCArrayKokkos <double> (material_options.size(), nstate_vars);
    global_vars = DCArrayKokkos <double> (material_options.size(), nglobal_vars);

    for (size_t i = 0; i < material_options.size(); i++) {
      auto mat = material_options[i];

      for (size_t j = 0; j < mat.state_vars.size(); j++)
        state_vars.host(i, j) = mat.state_vars[j];

      for (size_t j = 0; j < mat.global_vars.size(); j++)
        global_vars.host(i, j) = mat.global_vars[j];
    }
  }

  template<typename T> void from_vector(DCArrayKokkos<T>& array, const std::vector<T>& vec) {
    array = DCArrayKokkos<T>(vec.size());
    for (size_t i = 0; i < vec.size(); i++)
      array.host(i) = vec[i];
  }
  void derive_kokkos_arrays() {
    max_num_state_vars = 0;
    for (auto mo : material_options) 
      max_num_state_vars = std::max(max_num_state_vars, mo.state_vars.size());

    max_num_global_vars = 0;
    for (auto mo : material_options) 
      max_num_global_vars = std::max(max_num_global_vars, mo.global_vars.size());

    init_material_variable_arrays(max_num_state_vars, max_num_global_vars);

    from_vector(mat_fill, region_options);
    from_vector(material, material_options);
    from_vector(boundary, boundary_conditions);

    // Send to device.
    mat_fill.update_device();
    boundary.update_device();
    material.update_device();
    state_vars.update_device();
    global_vars.update_device();

    // Re-derive function pointers after move to device.
    for (size_t i = 0; i < material.size(); i++) {
      RUN_CLASS({ material(i).derive_function_pointers(); });
    }
  }

  void derive() {
    derive_kokkos_arrays();

    rk_num_bins = rk_num_stages;

    ensure_module(FEA_MODULE_TYPE::SGH);
  }
  void validate() { }
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Simulation_Parameters_SGH, Simulation_Parameters, 
  time_variables, material_options, region_options, 
  boundary_conditions, gravity_flag, report_runtime, rk_num_stages,
  NB, NBSF, NBV,
  graphics_options
)

// class Simulation_Parameters_SGH : public Simulation_Parameters
// {
//  public:
//   Simulation_Parameters_SGH();
//   virtual ~Simulation_Parameters_SGH();
//   virtual void input();
//   virtual void apply_settings();
//   virtual void FEA_module_setup();
    
//     // applying initial conditions
//   enum setup
//   {
//       none = 0,
//       Sedov3D = 1,
//       SedovRZ = 2,
        
//       Noh3D = 3,
//       NohRZ = 4,
        
//       SodZ = 5,
//       Sod3DX = 6,
//       Sod3DY = 7,
//       Sod3DZ = 8,
        
//       TriplePoint = 9,
//       TaylorAnvil = 10,
//   };

//   void select_problem(setup problem_selector);
    
//   // end of initial conditions enum

//   setup test_problem;

//   //==============================================================================
//   //   Mesh Variables
//   //==============================================================================

//   // --- Mesh regions and material fills ---
//   int NB; // number of boundary patch sets to tag
//   int NBSF; //number of surface force density boundary conditions
//   int NBV; //number of velocity boundary conditions


//   // --- Graphics output variables ---
//   bool output_velocity_flag, output_stress_flag, output_strain_flag, strain_max_flag, displaced_mesh_flag;

//   DCArrayKokkos <material_t> material;
//   DCArrayKokkos <double> state_vars; // array to hold init model variables
//   DCArrayKokkos <double> global_vars; // array to hold global varibles for user material model

//   DCArrayKokkos <mat_fill_t> mat_fill;
//   DCArrayKokkos <boundary_t> boundary;

//   // --- num vars ----
//   size_t num_dims;

//   size_t num_materials;
//   size_t max_num_state_vars;
//   size_t max_num_global_vars;

//   size_t num_fills;
//   size_t num_bcs;

//   // --- Graphics output variables ---
//   size_t graphics_id;
//   size_t graphics_cyc_ival;

//   CArray <double> graphics_times;
//   double graphics_dt_ival;
//   double graphics_time;  // the times for writing graphics dump


//   // --- Time and cycling variables ---
//   double time_value;
//   double time_final;
//   double dt;
//   double dt_max;
//   double dt_min;
//   double dt_cfl;
//   double dt_start;

//   size_t rk_num_stages;
//   size_t rk_num_bins;

//   size_t cycle;
//   size_t cycle_stop;

//   // --- Precision variables ---
//   double fuzz;  // machine precision
//   double tiny;  // very very small (between real_t and single)
//   double small;   // single precision

//   // -- Integration rule
//   int num_gauss_points;

//   //Body force parameters
//   bool gravity_flag;
//   real_t gravity_vector[3];

//   //SGH Solver Flags

//   //Possible options to parse; setting marked arbitrary if the options are not limited

//   options_multimap sgh_possible_options
//   {
//     { "solver_type", "SGH"}
//   }; // end std::map

//   nested_options_multimap sgh_possible_options_nested2
//   {
//     { "solver_options",
//         {
//             { "test_problem", "Sedov3D"},
//             { "mesh_file_name", "arbitrary"},
//             { "num_dims", "arbitrary"},
//             { "mesh_file_format", "ensight"},
//             { "mesh_file_format", "tecplot"},
//             { "mesh_file_format", "vtk"}
//         }
//     }
//   }; // end std::map

//   doubly_nested_options_multimap sgh_possible_options_nested3
//   {

//     //mesh options
//     { "mesh",
//         {
//             { "create",
//                 {
//                     {"type" , "3D-Box"}
//                 }
//             },
//             { "format",
//                 {
//                     {"type" , "geo"}
//                 }
//             },
//             { "parameters",
//                 {
//                     {"length_x" , "1.0"},
//                     {"length_y" , "1.0"},
//                     {"length_z" , "1.0"},
//                     {"num_x_elems", "10"},
//                     {"num_y_elems", "10"},
//                     {"num_z_elems", "10"},
//                     {"inner_radius" , "1.0"},
//                     {"outer_radius" , "1.0"},
//                     {"starting_angle" , "0.0"},
//                     {"ending_angle" , "180.0"},
//                     {"num_radial_elems" , "10"},
//                     {"num_angular_elems" , "10"},
//                     {"origin", "[0,0,0]"},
//                     {"order", "1"}
//                 }
//             }
//         }
//     }
//   }; // end std::map
  
// };

#endif // end HEADER_H
