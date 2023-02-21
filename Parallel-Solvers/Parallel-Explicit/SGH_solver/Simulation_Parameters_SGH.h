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

#include "utilities.h"
#include "state.h"
#include "matar.h"
#include "mesh.h"
#include "Simulation_Parameters.h"
using namespace utils;

class Simulation_Parameters_SGH : public Simulation_Parameters
{
 public:
  Simulation_Parameters_SGH();
  virtual ~Simulation_Parameters_SGH();
  virtual void input();
  virtual void apply_settings();
  virtual void FEA_module_setup();
  std::string yaml_input(std::string filename);
    
    // applying initial conditions
  enum setup
  {
      none = 0,
      Sedov3D = 1,
      SedovRZ = 2,
        
      Noh3D = 3,
      NohRZ = 4,
        
      SodZ = 5,
      Sod3DX = 6,
      Sod3DY = 7,
      Sod3DZ = 8,
        
      TriplePoint = 9,
      TaylorAnvil = 10,
  };

  
  void select_problem(setup problem_selector);
    
  // end of initial conditions enum

  setup test_problem;

  //==============================================================================
  //   Mesh Variables
  //==============================================================================

  // --- Mesh regions and material fills ---
  int NB; // number of boundary patch sets to tag
  int NBSF; //number of surface force density boundary conditions
  int NBV; //number of velocity boundary conditions


  // --- Graphics output variables ---
  bool output_velocity_flag, output_stress_flag, output_strain_flag, strain_max_flag, displaced_mesh_flag;

  CArrayKokkos <material_t> material;
  CArrayKokkos <double> state_vars; // array to hold init model variables

  CArrayKokkos <mat_fill_t> mat_fill;
  CArrayKokkos <boundary_t> boundary;

  // --- num vars ----
  size_t num_dims;

  size_t num_materials;
  size_t max_num_state_vars;

  size_t num_fills;
  size_t num_bcs;

  // --- Graphics output variables ---
  size_t graphics_id;
  size_t graphics_cyc_ival;

  CArray <double> graphics_times;
  double graphics_dt_ival;
  double graphics_time;  // the times for writing graphics dump


  // --- Time and cycling variables ---
  double time_value;
  double time_final;
  double dt;
  double dt_max;
  double dt_min;
  double dt_cfl;
  double dt_start;

  size_t rk_num_stages;
  size_t rk_num_bins;

  size_t cycle;
  size_t cycle_stop;

  // --- Precision variables ---
  double fuzz;  // machine precision
  double tiny;  // very very small (between real_t and single)
  double small;   // single precision

  // -- Integration rule
  int num_gauss_points;

  //Body force parameters
  bool gravity_flag;
  real_t gravity_vector[3];

  //SGH Solver Flags

  //Possible options to parse; setting marked arbitrary if the options are not limited

  options_multimap sgh_possible_options
  {
    { "solver_type", "SGH"}
  }; // end std::map

  nested_options_multimap sgh_possible_options_nested2
  {
    { "solver_options",
        {
            { "test_problem", "Sedov3D"},
            { "mesh_file_name", "arbitrary"},
            { "num_dims", "arbitrary"},
            { "mesh_file_format", "ensight"},
            { "mesh_file_format", "tecplot"},
            { "mesh_file_format", "vtk"}
        }
    }
  }; // end std::map

  doubly_nested_options_multimap sgh_possible_options_nested3
  {

    //mesh options
    { "mesh",
        {
            { "create",
                {
                    {"type" , "3D-Box"}
                }
            },
            { "format",
                {
                    {"type" , "geo"}
                }
            },
            { "parameters",
                {
                    {"length_x" , "1.0"},
                    {"length_y" , "1.0"},
                    {"length_z" , "1.0"},
                    {"num_x_elems", "10"},
                    {"num_y_elems", "10"},
                    {"num_z_elems", "10"},
                    {"inner_radius" , "1.0"},
                    {"outer_radius" , "1.0"},
                    {"starting_angle" , "0.0"},
                    {"ending_angle" , "180.0"},
                    {"num_radial_elems" , "10"},
                    {"num_angular_elems" , "10"},
                    {"origin", "[0,0,0]"},
                    {"order", "1"}
                }
            }
        }
    }
  }; // end std::map
  
};

//eos forward declaration
KOKKOS_FUNCTION
void ideal_gas(const DViewCArrayKokkos <double> &elem_pres,
               const DViewCArrayKokkos <double> &elem_stress,
               const size_t elem_gid,
               const size_t mat_id,
               const DViewCArrayKokkos <double> &elem_state_vars,
               const DViewCArrayKokkos <double> &elem_sspd,
               const double den,
               const double sie);

#endif // end HEADER_H
