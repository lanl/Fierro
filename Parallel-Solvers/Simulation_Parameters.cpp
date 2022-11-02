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
#include "Simulation_Parameters.h"
#include "elements.h"
#include "swage.h"

using namespace utils;

Simulation_Parameters::Simulation_Parameters(){

  //initialize data and flags to defaults
  report_runtime_flag = 0;
  unit_scaling = 1;
  restart_file = false;
  tecplot_input = ansys_dat_input = false;
  nfea_modules = 0;
  element_type = "Hex8";
}

Simulation_Parameters::~Simulation_Parameters(){
}

void Simulation_Parameters::input(){
  //file input flags
  tecplot_input = false;
  restart_file = false;
  ansys_dat_input = false;

  //simulation spatial dimension
  num_dim = 3;
  unit_scaling = 1;
  
  //file readin parameters
  words_per_line = 1;
  tecplot_words_per_line = 3;
  ansys_dat_node_words_per_line = 4;
  elem_words_per_line = 8;
  ansys_dat_elem_words_per_line = 11;
  element_type = "Hex8";

  //debug and performance report flags
  report_runtime_flag = 1;

}

void Simulation_Parameters::FEA_module_setup(){
  
  //initial buffer size for FEA module list storage
  int buffer_size = 10 + nfea_modules;
  FEA_Module_List.resize(buffer_size);
  fea_module_must_read.resize(buffer_size);
  int start_module = nfea_modules;

  //decides which FEA modules to setup based on user decided implicit solves
  FEA_Module_List[nfea_modules] = "Elasticity";
  nfea_modules++;
  FEA_Module_List[nfea_modules] = "Heat_Conduction";
  nfea_modules++;
  //example for later
  if(nfea_modules==buffer_size){
    buffer_size += 10;
    FEA_Module_List.resize(buffer_size);
    fea_module_must_read.resize(buffer_size);
  }

  //initialize
  for(int imodule = start_module; imodule < nfea_modules; imodule++){
    fea_module_must_read[imodule] = false;
  }
}
