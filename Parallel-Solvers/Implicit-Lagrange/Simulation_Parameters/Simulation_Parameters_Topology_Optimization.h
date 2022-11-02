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

#ifndef SIMULATION_PARAMETERS_TOPOLOGY_OPTIMIZATION_H
#define SIMULATION_PARAMETERS_TOPOLOGY_OPTIMIZATION_H

#include "utilities.h"
#include "Simulation_Parameters.h"
using namespace utils;

//forward declare
class Implicit_Solver;

class Simulation_Parameters_Topology_Optimization : public Simulation_Parameters
{
 public:
  Simulation_Parameters_Topology_Optimization(Implicit_Solver *solver_pointer);
  virtual ~Simulation_Parameters_Topology_Optimization();
  virtual void input();
  virtual void FEA_module_setup();
  //==============================================================================
  //   Mesh Variables
  //==============================================================================

  // --- Mesh regions and material fills ---
  int NB; // number of boundary patch sets to tag
  int NBD; //number of density boundary conditions
  
  // --- Dimensional and mesh constants ---
  int num_dim;
  int p_order;

  //debug and performance reporting flags
  bool report_runtime_flag;

  //Topology Optimization flags
  bool nodal_density_flag;
  
  //When on, all element nodes connected to a boundary condition patch will have their density constrained
  bool thick_condition_boundary;

  //file output parameters
  int optimization_output_freq;

  //Topology Optimization parameters
  real_t penalty_power;

  //pointer to Solver object (just used to consolidate error handling for now)
  Implicit_Solver *solver_pointer_;

  //volumes to hold density constant
  
  //types of TO functions
  enum function_type {OBJECTIVE, MULTI_OBJECTIVE_TERM, EQUALITY_CONSTRAINT, INEQUALITY_CONSTRAINT, VECTOR_EQUALITY_CONSTRAINT, VECTOR_INEQUALITY_CONSTRAINT};

  //list of TO functions needed by problem
  std::vector<std::string> TO_Module_List;
  std::vector<function_type> TO_Function_Type;
  std::vector<int> TO_Module_My_FEA_Module;
  std::vector<int> Multi_Objective_Modules;
  std::vector<real_t> Multi_Objective_Weights;
  std::vector<std::vector<real_t>> Function_Arguments;
  int nTO_modules, nmulti_objective_modules;
};

#endif // end HEADER_H
