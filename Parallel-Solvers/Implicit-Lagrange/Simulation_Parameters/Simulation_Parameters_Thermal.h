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

#ifndef SIMULATION_PARAMETERS_THERMAL_H
#define SIMULATION_PARAMETERS_THERMAL_H

#include "utilities.h"
#include "yaml-serializable.h"
#include "Simulation_Parameters.h"
using namespace utils;

struct Simulation_Parameters_Thermal: public Simulation_Parameters {
  // --- Boundary Conditions ---
  int NB   = 6; // number of boundary patch sets to tag
  int NBSF = 4; //number of surface heat flux boundary conditions
  int NBT  = 2; //number of temperature boundary conditions

  // --- Graphics output variables ---
  bool output_temperature_flag     = true; 
  bool output_heat_flux_flag       = true;
  bool output_temperature_gradient_flag = false;
  bool flux_max_flag = false;

  // --- Constitutive Parameters ---
  double Thermal_Conductivity = 10;

  // --- Integration Scheme
  int num_gauss_points = 2;

  //debug and performance reporting flags
  bool multigrid_timers = false;
  bool direct_solver_flag = false;

  //Body flux parameters
  bool thermal_flag = false;
  //bool electric_flag;
  double specific_internal_energy_rate = 1;

  //Linear Solver Flags
  bool equilibrate_matrix_flag = false;

  //Topology Optimization parameters
  //double maximum_strain, maximum_strain_energy;
};
IMPL_YAML_SERIALIZABLE_WITH_BASE(Simulation_Parameters_Thermal, Simulation_Parameters)

// class Simulation_Parameters_Thermal: public Simulation_Parameters
// {
//  public:
//   Simulation_Parameters_Thermal();
//   virtual ~Simulation_Parameters_Thermal();
//   virtual void input();
//   virtual void apply_settings() {}
//   //==============================================================================
//   //   Thermal FEA problem parameters
//   //==============================================================================

//   // --- Boundary Conditions ---
//   int NB; // number of boundary patch sets to tag
//   int NBSF; //number of surface heat flux boundary conditions
//   int NBT; //number of temperature boundary conditions

//   // --- Graphics output variables ---
//   bool output_temperature_flag, output_temperature_gradient_flag, output_heat_flux_flag, flux_max_flag;

//   // --- Constitutive Parameters ---
//   real_t Thermal_Conductivity;

//   // --- Integration Scheme
//   int num_gauss_points;

//   //debug and performance reporting flags
//   bool multigrid_timers, direct_solver_flag;

//   //Body flux parameters
//   bool thermal_flag, electric_flag;
//   real_t specific_internal_energy_rate;

//   //Linear Solver Flags
//   bool equilibrate_matrix_flag;

//   //Topology Optimization parameters
//   real_t maximum_strain, maximum_strain_energy;
// };

#endif // end HEADER_H
