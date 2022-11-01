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

#include "utilities.h"
#include "Simulation_Parameters_Thermal.h"

using namespace utils;

Simulation_Parameters_Thermal::Simulation_Parameters_Thermal() : Simulation_Parameters(){

  //initialize data and flags to defaults
  output_temperature_flag = false;
  output_temperature_gradient_flag = false;
  output_heat_flux_flag = false;
  report_runtime_flag = false;
  unit_scaling = 1;
  flux_max_flag = false;
  direct_solver_flag = false;
  thermal_flag = false;
  multigrid_timers = false;
  equilibrate_matrix_flag = false;
  // ---- boundary conditions ---- //
  NB = 0;
  NBSF = 0;
  NBT = 0;
}

Simulation_Parameters_Thermal::~Simulation_Parameters_Thermal(){
}

void Simulation_Parameters_Thermal::input(){
  
  Simulation_Parameters::input();
  //multigrid_timers = true;
  equilibrate_matrix_flag = false;

  //simulation spatial dimension
  num_dim = 3;
  unit_scaling = 1;

  //polynomial interpolation order
  p_order = 0;

  // ---- graphics information ---- //
  output_temperature_flag = true;
  output_heat_flux_flag = true;
  
  //Isotropic Conductivity parameters to move into a child class later
  Thermal_Conductivity = 10;

  //Gauss-Legendre parameters
  num_gauss_points = 2;

  //debug and performance report flags
  report_runtime_flag = true;

  // ---- boundary conditions ---- //
  NB = 6; // number of boundary conditions for this module
  NBSF = 4; //number of surface heat flux conditions
  NBT = 2; //number of surface sets used to specify a fixed temperature on nodes belonging to respective surfaces

  //apply body forces
  thermal_flag = false;
  specific_internal_energy_rate = 1;

}