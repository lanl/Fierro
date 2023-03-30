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

#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H

#include "utilities.h"
#include "Yaml.hpp"
#include "mpi.h"
#include <stdio.h>
#include <vector>
#include <string>
using namespace utils;

class Simulation_Parameters
{
 public:
  Simulation_Parameters();
  virtual ~Simulation_Parameters();
  virtual void input(); //typically sets default problem parameters
  virtual void apply_settings();
  virtual std::string yaml_input(std::string filename); //reads in user defined parameters
  virtual void FEA_module_setup();
  virtual void yaml_FEA_module_setup();

  //==============================================================================
  //   Mesh Variables
  //==============================================================================

  // --- Mesh regions and material fills ---
  int NR; // number of Regions
  int NC; // number of contours
  int NF; // number of fill
  int NB; // number of boundary patch sets to tag

  // --- Dimensional and mesh constants ---
  int num_dim;
  int p_order;
  
  //file input parameters 
  int words_per_line, elem_words_per_line, tecplot_words_per_line, vtk_words_per_line, ansys_dat_node_words_per_line, ansys_dat_elem_words_per_line;
  char *format_specification;  //per line file format when reading dofs
  real_t unit_scaling;
  bool restart_file;
  bool tecplot_input, ansys_dat_input, vtk_input, zero_index_base;
  std::string element_type;
  std::string solver_type, mesh_file_name, mesh_file_format;

  //debug and performance reporting flags
  int report_runtime_flag;
 
  //necessary FEA modules
  std::vector<std::string> FEA_Module_List;
  std::vector<bool> fea_module_must_read;
  int nfea_modules;

  //====================================================================================================
  //  possible values in dictionary of options (our yaml parser supports three nested levels of options)
  //====================================================================================================
  typedef std::multimap<std::string,std::string> options_multimap;

  typedef std::multimap<std::string,std::multimap<std::string,std::string>> nested_options_multimap;

  typedef std::multimap<std::string,std::multimap<std::string,std::multimap<std::string,std::string>>> doubly_nested_options_multimap;

  options_multimap possible_options;

  nested_options_multimap possible_options_nested2;

  doubly_nested_options_multimap possible_options_nested3;

  std::map<std::string,std::string> set_options;

  std::map<std::string,std::map<std::string,std::string>> set_options_nested2;

  std::map<std::string,std::map<std::string,std::map<std::string,std::string>>> set_options_nested3;

  typedef std::pair<std::string,std::string> option_setting_pair;

  typedef options_multimap::iterator multimap_iterator;

  typedef nested_options_multimap::iterator multimap_iterator_nested2;
  
  typedef doubly_nested_options_multimap::iterator multimap_iterator_nested3;

  //MPI data
  int myrank; //index of this mpi rank in the world communicator
  int nranks; //number of mpi ranks in the world communicator
  MPI_Comm world; //stores the default communicator object (MPI_COMM_WORLD)

};

#endif // end HEADER_H
