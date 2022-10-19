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
#include "Simulation_Parameters_Elasticity.h"

using namespace utils;

Simulation_Parameters_Elasticity::Simulation_Parameters_Elasticity() : Simulation_Parameters(){

  //initialize data and flags to defaults
  output_displacement_flag = false;
  output_strain_flag = false;
  output_stress_flag = false;
  displaced_mesh_flag = false;
  report_runtime_flag = false;
  unit_scaling = 1;
  strain_max_flag = false;
  direct_solver_flag = false;
  gravity_flag = false;
  multigrid_timers = false;
  equilibrate_matrix_flag = false;
  // ---- boundary conditions ---- //
  NB = 0; 
  NBSF = 0; 
  NBD = 0;
}

Simulation_Parameters_Elasticity::~Simulation_Parameters_Elasticity(){
}

void Simulation_Parameters_Elasticity::input(){
  
  Simulation_Parameters::input();
  //output settings
  output_displacement_flag = true;
  //requires displacement flag to be true
  displaced_mesh_flag = true;
  
  output_strain_flag = true;
  output_stress_flag = false;

  //multigrid_timers = true;
  equilibrate_matrix_flag = false;

  //simulation spatial dimension
  num_dim = 3;
  unit_scaling = 1;

  //polynomial interpolation order
  p_order = 0;
  
  //Static isotropic parameters to move into a child class later
  Elastic_Modulus = 200000000000;
  Poisson_Ratio = 0.3;
  
  //Gauss-Legendre integration order
  num_gauss_points = 2;

  //debug and performance report flags
  report_runtime_flag = true;

  // ---- boundary conditions ---- //
  NB = 6; // number of boundaries
  NBSF = 4; //number of surface density force conditions
  NBD = 2; //number of surface sets used to specify a fixed displacement on nodes belonging to respective surfaces

  //apply body forces
  gravity_flag = false;
  gravity_vector[0] = 9.81;
  gravity_vector[1] = 0;
  gravity_vector[2] = 0;

}