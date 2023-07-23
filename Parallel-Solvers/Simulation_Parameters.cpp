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

// Simulation_Parameters::Simulation_Parameters(){

//   //initialize data and flags to defaults
//   report_runtime_flag = 0;
//   unit_scaling = 1;
//   restart_file = false;
//   tecplot_input = ansys_dat_input = false;
//   zero_index_base = false;
//   nfea_modules = 0;
//   element_type = "Hex8";
//   filtered_density = false;
//   enable_inertia_center.push_back(false);
//   enable_inertia_center.push_back(false);
//   enable_inertia_center.push_back(false);
//   moment_of_inertia_center = std::vector<double>(3);
//   output_file_format = "vtk";
//   timer_output_level = "standard";

//   //MPI info
//   world = MPI_COMM_WORLD; //used for convenience to represent all the ranks in the job
//   MPI_Comm_rank(world,&myrank);
//   MPI_Comm_size(world,&nranks);
// }

// Simulation_Parameters::~Simulation_Parameters(){
// }

// void Simulation_Parameters::input(){
//   //file input flags
//   tecplot_input = false;
//   restart_file = false;
//   ansys_dat_input = false;
//   vtk_input = false;
//   zero_index_base = false;

//   //simulation spatial dimension
//   num_dim = 3; //simulation spatial dimension
//   p_order = 0; //polynomial interpolation order
//   unit_scaling = 1;
  
//   //file readin parameters
//   words_per_line = 1;
//   vtk_words_per_line = tecplot_words_per_line = 3;
//   ansys_dat_node_words_per_line = 4;
//   ansys_dat_elem_words_per_line = 11;

//   //default element types
//   if(num_dim==3)
//     element_type = "Hex8";
//   else if(num_dim==2)
//     element_type = "Quad4";

//   if(element_type == "Hex8") 
//     elem_words_per_line = 8;
//   else if(element_type == "Quad4") 
//     elem_words_per_line = 4;

//   //debug and performance report flags
//   report_runtime_flag = 1;

// }

// //==============================================================================
// //    Communicate user settings from YAML file and apply to class members
// //==============================================================================

// void Simulation_Parameters::apply_settings(){
//   std::string current_option;

//   current_option = "solver_type"; //string for the parameter to find
//   if(set_options.find(current_option)!=set_options.end()){
//       //set parameter here
//       solver_type = set_options[current_option];
//       set_options.erase(current_option);
//   }

//   if(myrank==0)
//     std::cout << "Solver Type is " << solver_type << std::endl;

//   current_option = "solver_options:mesh_file_name"; //string for the parameter to find
//   if(set_options.find(current_option)!=set_options.end()){
//       //set parameter here
//       mesh_file_name = set_options[current_option];
//       set_options.erase(current_option);
//   }
  
//   current_option = "solver_options:mesh_file_format"; //string for the parameter to find
//   if(set_options.find(current_option)!=set_options.end()){
//       //set parameter here
//       mesh_file_format = set_options[current_option];
//       set_options.erase(current_option);
//   }

//   current_option = "solver_options:output_file_format"; //string for the parameter to find
//   if(set_options.find(current_option)!=set_options.end()){
//       //set parameter here
//       output_file_format = set_options[current_option];
//       set_options.erase(current_option);
//   }

//   current_option = "solver_options:timer_output_level"; //string for the parameter to find
//   if(set_options.find(current_option)!=set_options.end()){
//       //set parameter here
//       timer_output_level = set_options[current_option];
//       set_options.erase(current_option);
//   }

//   current_option = "solver_options:num_dims"; //string for the parameter to find
//   if(set_options.find(current_option)!=set_options.end()){
//       //set parameter here
//       num_dim = std::stoi(set_options[current_option]);
//       set_options.erase(current_option);
//   }

//   current_option = "output_options:graphics_step_frequency"; //string for the parameter to find
//   if(set_options.find(current_option)!=set_options.end()){
//       //set parameter here
//       file_output_frequency = std::stoi(set_options[current_option]);
//       set_options.erase(current_option);
//   }
  
//   if(myrank==0)
//     std::cout << "Mesh File name is " << mesh_file_name << std::endl;
// }

// //==============================================================================
// //    Setup FEA Modules for the simulation
// //==============================================================================

// void Simulation_Parameters::FEA_module_setup(){
  
//   //initial buffer size for FEA module list storage
//   int buffer_size = 10 + nfea_modules;
//   FEA_Module_List.resize(buffer_size);
//   fea_module_must_read.resize(buffer_size);
//   int start_module = nfea_modules;

//   //decides which FEA modules to setup based on user decided implicit solves
//   FEA_Module_List[nfea_modules] = "Elasticity";
//   nfea_modules++;
//   FEA_Module_List[nfea_modules] = "Heat_Conduction";
//   nfea_modules++;
//   //FEA_Module_List[nfea_modules] = "Thermo_Elasticity";
//   //nfea_modules++;
//   //example for later
//   if(nfea_modules==buffer_size){
//     buffer_size += 10;
//     FEA_Module_List.resize(buffer_size);
//     fea_module_must_read.resize(buffer_size);
//   }

//   //initialize
//   for(int imodule = start_module; imodule < nfea_modules; imodule++){
//     fea_module_must_read[imodule] = false;
//   }
// }

// //==============================================================================
// //    Setup FEA Modules for the simulation from requested yaml option
// //==============================================================================

// void Simulation_Parameters::yaml_FEA_module_setup(){
//   std::string current_option;
//   //initial buffer size for FEA module list storage
//   int buffer_size = 10 + nfea_modules;
//   FEA_Module_List.resize(buffer_size);
//   fea_module_must_read.resize(buffer_size);
//   int start_module = nfea_modules;

//   std::string fea_module_base = "fea_module_";
//   std::string index;
//   std::string fea_module_name;
  
//   index = std::to_string(nfea_modules+1);
//   fea_module_name = fea_module_base + index;
  
//   // --- set of user requested FEA modules ---
//   while(set_options.find(fea_module_name+":type")!=set_options.end()){
//     //std::cout << "WHILE LOOP ENTERED " << index << std::endl;
//     if(set_options[fea_module_name+":type"]=="elasticity"){
//         //std::cout << "read FEA module " << fea_module_name+":type" << std::endl;
//         FEA_Module_List[nfea_modules] = "Elasticity";
//         set_options.erase(fea_module_name+":type");
//     }
//     else if(set_options[fea_module_name+":type"]=="steady_heat"){
//         //std::cout << "read FEA module " << fea_module_name+":type" << std::endl;
//         FEA_Module_List[nfea_modules] = "Heat_Conduction";
//         set_options.erase(fea_module_name+":type");
//     }
    
//     nfea_modules++;
    
//     //if(nfea_modules==buffer_size){
//       //buffer_size += 10;
//       //FEA_Module_List.resize(buffer_size);
//       //fea_module_must_read.resize(buffer_size);
//     //}
    
//     index = std::to_string(nfea_modules+1);
//     fea_module_name = fea_module_base + index;
//   }
  
//   /*
//   if(set_options[fea_module_name+":type"]=="elasticity"){
//         std::cout << "read FEA module " << fea_module_name+":type" << std::endl;
//         FEA_Module_List[nfea_modules] = "Elasticity";
//         set_options.erase(fea_module_name+":type");
//   }
//   nfea_modules++;
//   index = std::to_string(nfea_modules+1);
//   fea_module_name = fea_module_base + index;
//   if(set_options[fea_module_name+":type"]=="steady_heat"){
//         std::cout << "read FEA module " << fea_module_name+":type" << std::endl;
//         FEA_Module_List[nfea_modules] = "Heat_Conduction";
//         set_options.erase(fea_module_name+":type");
//   }
//   nfea_modules++;
//   */

//   if(nfea_modules==buffer_size){
//       buffer_size += 10;
//       FEA_Module_List.resize(buffer_size);
//       fea_module_must_read.resize(buffer_size);
//   }

//   //allocate request for Inertial module (needed generally so far)
//   FEA_Module_List[nfea_modules++] = "Inertial";

//   //initialize
//   for(int imodule = start_module; imodule < nfea_modules; imodule++){
//     fea_module_must_read[imodule] = false;
//   }
// }