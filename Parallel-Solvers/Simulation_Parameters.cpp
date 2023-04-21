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
  zero_index_base = false;
  nfea_modules = 0;
  element_type = "Hex8";

  //MPI info
  world = MPI_COMM_WORLD; //used for convenience to represent all the ranks in the job
  MPI_Comm_rank(world,&myrank);
  MPI_Comm_size(world,&nranks);
}

Simulation_Parameters::~Simulation_Parameters(){
}

void Simulation_Parameters::input(){
  //file input flags
  tecplot_input = false;
  restart_file = false;
  ansys_dat_input = false;
  vtk_input = false;
  zero_index_base = false;

  //simulation spatial dimension
  num_dim = 3; //simulation spatial dimension
  p_order = 0; //polynomial interpolation order
  unit_scaling = 1;
  
  //file readin parameters
  words_per_line = 1;
  vtk_words_per_line = tecplot_words_per_line = 3;
  ansys_dat_node_words_per_line = 4;
  ansys_dat_elem_words_per_line = 11;

  //default element types
  if(num_dim==3)
    element_type = "Hex8";
  else if(num_dim==2)
    element_type = "Quad4";

  if(element_type == "Hex8") 
    elem_words_per_line = 8;
  else if(element_type == "Quad4") 
    elem_words_per_line = 4;

  //debug and performance report flags
  report_runtime_flag = 1;

}


//==============================================================================
//    Read in user settings from YAML file
//==============================================================================

std::string Simulation_Parameters::yaml_input(std::string filename){
  Yaml::Node root;
  std::string current_option_outer, current_setting_outer;
  std::string current_option_inner, current_setting_inner;
  std::string current_option_center, current_setting_center;
  std::string error = "success";
  std::string colon = ":";
  bool inner_found, inner_found_nest;
  if(myrank==0){
    try
    {
        Yaml::Parse(root, filename.c_str());
    }
    catch (const Yaml::Exception e)
    {
        std::cout << "Exception " << e.Type() << ": " << e.what() << std::endl;
    }
    
    
    //std::cout << "print root Size = " << root.Size() << "\n";
    //std::cout << "print root = \n";
    //for (size_t i=0; i<root.Size(); i++){
        
        // get the outer stuff
        Yaml::Node & outer_item = root;
        
        //std::cout << "\n";
        //std::cout << "size = " << outer_item.Size() << std::endl;
        if (outer_item.Size()!=0){
            for(auto outer_it = outer_item.Begin(); outer_it != outer_item.End(); outer_it++)
            {   
                current_option_outer = (*outer_it).first;
                current_setting_outer = (*outer_it).second.As<std::string>();
                
                //find the keyword for this option out of the three multimaps of possible options with different nesting structure
                //multimap_iterator multi_outer_iterator = sgh_possible_options.find(current_option_outer);
                //multimap_iterator_nested2 multi_outer_iterator_nested2 = sgh_possible_options_nested2.find(current_option_outer);
                //multimap_iterator_nested3 multi_outer_iterator_nested3 = sgh_possible_options_nested3.find(current_option_outer);

                //check if this keyword is an option for this solver type
                /*
                if(multi_outer_iterator==possible_options.end()){
                    if(multi_outer_iterator_nested2==possible_options_nested2.end()){
                        if(multi_outer_iterator_nested3==possible_options_nested3.end()){
                            std::string message = "Unsupported option requested in YAML input file: ";
                            error = message + current_option_outer;
                        }
                    }
                  //return error;
                }
                */
                
                //std::cout << current_option_outer << " " << current_setting_outer << std::endl;
            
                Yaml::Node & inner_item = (*outer_it).second;
            
                // inner layer
                if (inner_item.Size()!=0){
                    for(auto inner_it = inner_item.Begin(); inner_it != inner_item.End();   inner_it++)
                    {
                        current_option_inner = (*inner_it).first;
                        current_setting_inner = (*inner_it).second.As<std::string>();
                        
                        //check if this option is supported
                        //std::pair<multimap_iterator_nested2,multimap_iterator_nested2> iterator_range = sgh_possible_options_nested2.equal_range(current_option_outer);
                        //std::pair<multimap_iterator_nested3,multimap_iterator_nested3> iterator_range_nested = sgh_possible_options_nested3.equal_range(current_option_outer);
                        //multimap_iterator multi_outer_iterator_nested2 = sgh_possible_options_nested2[current_option_outer].find(current_option_outer);
                        //multimap_iterator_nested2 multi_outer_iterator_nested3 = sgh_possible_options_nested3[current_option_outer].find(current_option_outer);
                        
                        inner_found = inner_found_nest = false;
                        /*
                        for(auto temp_it = iterator_range.first; temp_it != iterator_range.second; temp_it++){
                            //test if map element corresponds to a possible option
                            multimap_iterator multi_inner_iterator = temp_it->second.find(current_option_inner);
                            //std::cout << "Test print " << temp_it->first << std::endl;
                            if(multi_inner_iterator!=temp_it->second.end()){
                                inner_found = true;
                            }
                        }

                        for(auto temp_it = iterator_range_nested.first; temp_it != iterator_range_nested.second; temp_it++){
                            //test if map element corresponds to a possible option
                            multimap_iterator_nested2 multi_inner_iterator = temp_it->second.find(current_option_inner);
                            //std::cout << "Test print " << temp_it->first << std::endl;
                            if(multi_inner_iterator!=temp_it->second.end()){
                                inner_found_nest = true;
                            }
                        }
                        */
                        /*
                        if(!inner_found&&!inner_found_nest){
                            std::string message = "Unsupported option requested in YAML input file: ";
                            error = message + current_option_inner;
                        }
                        */
                        //std::cout << "    " << current_option_inner << " " << current_setting_inner << " " << inner_found << " " << inner_found_nest << std::endl;
            
                        // inner_most layer
                        Yaml::Node & center_item = (*inner_it).second;
            
                        //std::cout << "  \n";
                        if (center_item.Size()!=0){
                            for(auto center_it = center_item.Begin(); center_it !=  center_item.End();   center_it++)
                            {
                                current_option_center = (*center_it).first;
                                current_setting_center = (*center_it).second.As<std::string>();

                                //check if this option is supported
                                /*
                                for(auto temp_it = iterator_range_nested.first; temp_it != iterator_range_nested.second; temp_it++){
                                    //test if map element corresponds to a possible option
                                    std::pair<multimap_iterator_nested2,multimap_iterator_nested2> iterator_range_nested3 = temp_it->second.equal_range(current_option_inner);
                                    for(auto temp_it_inner = iterator_range_nested3.first; temp_it_inner != iterator_range_nested3.second; temp_it_inner++){
                                        //test if map element corresponds to a possible option
                                        multimap_iterator multi_center_iterator = temp_it_inner->second.find(current_option_center);
                                        if(multi_center_iterator!=temp_it_inner->second.end()){
                                            std::string message = "Unsupported option requested in YAML input file: ";
                                            error = message + current_option_center;
                                        }
                                    }
                                }
                                */
                                //std::cout << "        " << current_option_center << "   " << current_setting_center << std::endl;

                                set_options[current_option_outer + colon + current_option_inner + colon + current_option_center] = current_setting_center;
                            } // end for
                        }
                        else{
                            set_options[current_option_outer + colon + current_option_inner] = current_setting_inner;
                        }
            
                    } // end for
                } // end if inner_item.Size
                else{
                    //there were no inner items; add setting to map of user defined settings to query later
                    set_options[current_option_outer] = current_setting_outer;
                    //if(current_option_outer=="solver_type"){
                      //std::string test_string = "solver_type";
                      //std::cout << "solver type: " << set_options[test_string] << std::endl;
                    //}

                } //end else for outer item with no inner items
            
            } // end for outer_it
        } // end if outer_it
    //} // end for
  }
  size_t max_string_size = 100;
  size_t settings_map_size;
  std::string temp_option, temp_setting;
  std::map<std::string,std::string>::iterator temp_it;

  //compute maximum string size in settings
  if(myrank==0){
    for(auto temp_it = set_options.begin(); temp_it != set_options.end(); temp_it++){
          if(temp_it->first.length()>max_string_size) max_string_size = temp_it->first.length() + 100;
          if(temp_it->second.length()>max_string_size) max_string_size = temp_it->second.length() + 100;
      }
  }

  if(myrank==0) settings_map_size = set_options.size();

  MPI_Bcast(&max_string_size,1,MPI_LONG_LONG_INT,0,world);
  MPI_Bcast(&settings_map_size,1,MPI_LONG_LONG_INT,0,world);

  //std::cout << "max string length on rank " << myrank << " " << max_string_size << std::endl;

  char* read_buffer = new char[max_string_size];
  char* read_buffer2 = new char[max_string_size];
    
  //for(auto temp_it = set_options.begin(); temp_it != set_options.end(); temp_it++){
  if(myrank==0)
    temp_it = set_options.begin();
    
  for(int imap = 0; imap < settings_map_size; imap++){
    //copy map strings to read buffer for mpi
    if(myrank==0){
      strncpy(read_buffer, temp_it->first.c_str(), temp_it->first.length()+1); //includes null terminator
      strncpy(read_buffer2, temp_it->second.c_str(), temp_it->second.length()+1); //includes null terminator
      temp_it++;
    }

    //communicate map settings read in on rank 0 to other ranks
    MPI_Bcast(read_buffer,max_string_size,MPI_CHAR,0,world);
    MPI_Bcast(read_buffer2,max_string_size,MPI_CHAR,0,world);

    temp_option = read_buffer;
    temp_setting = read_buffer2;

    //insert option name into map
    set_options[temp_option] = temp_setting;
  }

  MPI_Barrier(world);

  delete[] read_buffer;
  delete[] read_buffer2;

  return error;
}

//==============================================================================
//    Communicate user settings from YAML file and apply to class members
//==============================================================================

size_t Simulation_Parameters::unapplied_settings(){

  if(myrank==0){
    //print user settings for this module
    for(auto temp_it = set_options.begin(); temp_it != set_options.end(); temp_it++){
        //print current option
        std::cout << "Unapplied user option: " << temp_it->first << "=" << temp_it->second << std::endl;
    }
  }

  return set_options.size();
}

//==============================================================================
//    Communicate user settings from YAML file and apply to class members
//==============================================================================

void Simulation_Parameters::apply_settings(){
  std::string current_option;

  current_option = "solver_type"; //string for the parameter to find
  if(set_options.find(current_option)!=set_options.end()){
      //set parameter here
      solver_type = set_options[current_option];
      set_options.erase(current_option);
  }

  if(myrank==0)
    std::cout << "Solver Type is " << solver_type << std::endl;

  current_option = "solver_options:mesh_file_name"; //string for the parameter to find
  if(set_options.find(current_option)!=set_options.end()){
      //set parameter here
      mesh_file_name = set_options[current_option];
      set_options.erase(current_option);
  }
  
  current_option = "solver_options:mesh_file_format"; //string for the parameter to find
  if(set_options.find(current_option)!=set_options.end()){
      //set parameter here
      mesh_file_format = set_options[current_option];
      set_options.erase(current_option);
  }

  current_option = "solver_options:num_dims"; //string for the parameter to find
  if(set_options.find(current_option)!=set_options.end()){
      //set parameter here
      num_dim = std::stoi(set_options[current_option]);
      set_options.erase(current_option);
  }

  current_option = "output_options:graphics_step_frequency"; //string for the parameter to find
  if(set_options.find(current_option)!=set_options.end()){
      //set parameter here
      file_output_frequency = std::stoi(set_options[current_option]);
      set_options.erase(current_option);
  }
  
  if(myrank==0)
    std::cout << "Mesh File name is " << mesh_file_name << std::endl;
}

//==============================================================================
//    Setup FEA Modules for the simulation
//==============================================================================

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
  //FEA_Module_List[nfea_modules] = "Thermo_Elasticity";
  //nfea_modules++;
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

//==============================================================================
//    Setup FEA Modules for the simulation from requested yaml option
//==============================================================================

void Simulation_Parameters::yaml_FEA_module_setup(){
  std::string current_option;
  //initial buffer size for FEA module list storage
  int buffer_size = 10 + nfea_modules;
  FEA_Module_List.resize(buffer_size);
  fea_module_must_read.resize(buffer_size);
  int start_module = nfea_modules;

  std::string fea_module_base = "fea_module_";
  std::string index;
  std::string fea_module_name;
  
  index = std::to_string(nfea_modules+1);
  fea_module_name = fea_module_base + index;
  
  // --- set of user requested FEA modules ---
  while(set_options.find(fea_module_name+":type")!=set_options.end()){
    std::cout << "WHILE LOOP ENTERED " << index << std::endl;
    if(set_options[fea_module_name+":type"]=="elasticity"){
        std::cout << "read FEA module " << fea_module_name+":type" << std::endl;
        FEA_Module_List[nfea_modules] = "Elasticity";
        set_options.erase(fea_module_name+":type");
    }
    else if(set_options[fea_module_name+":type"]=="steady_heat"){
        std::cout << "read FEA module " << fea_module_name+":type" << std::endl;
        FEA_Module_List[nfea_modules] = "Heat_Conduction";
        set_options.erase(fea_module_name+":type");
    }
    
    nfea_modules++;
    
    //if(nfea_modules==buffer_size){
      //buffer_size += 10;
      //FEA_Module_List.resize(buffer_size);
      //fea_module_must_read.resize(buffer_size);
    //}
    
    index = std::to_string(nfea_modules+1);
    fea_module_name = fea_module_base + index;
  }
  
  /*
  if(set_options[fea_module_name+":type"]=="elasticity"){
        std::cout << "read FEA module " << fea_module_name+":type" << std::endl;
        FEA_Module_List[nfea_modules] = "Elasticity";
        set_options.erase(fea_module_name+":type");
  }
  nfea_modules++;
  index = std::to_string(nfea_modules+1);
  fea_module_name = fea_module_base + index;
  if(set_options[fea_module_name+":type"]=="steady_heat"){
        std::cout << "read FEA module " << fea_module_name+":type" << std::endl;
        FEA_Module_List[nfea_modules] = "Heat_Conduction";
        set_options.erase(fea_module_name+":type");
  }
  nfea_modules++;
  */

  if(nfea_modules==buffer_size){
      buffer_size += 10;
      FEA_Module_List.resize(buffer_size);
      fea_module_must_read.resize(buffer_size);
  }

  //allocate request for Inertial module (needed generally so far)
  FEA_Module_List[nfea_modules++] = "Inertial";

  //initialize
  for(int imodule = start_module; imodule < nfea_modules; imodule++){
    fea_module_must_read[imodule] = false;
  }
}