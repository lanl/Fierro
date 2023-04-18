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
#include "Simulation_Parameters_Topology_Optimization.h"
#include <iostream>
#include "elements.h"
#include "swage.h"
#include "Implicit_Solver.h"

using namespace utils;

Simulation_Parameters_Topology_Optimization::Simulation_Parameters_Topology_Optimization(Implicit_Solver *solver_pointer){

  //initialize data and flags to defaults
  solver_pointer_ = solver_pointer;
  report_runtime_flag = false;
  nodal_density_flag = true;
  thick_condition_boundary = true;
  topology_optimization_on = shape_optimization_on = false;
  optimization_output_freq = 20;
  penalty_power = 3;
  density_epsilon = 0.0001;
  nTO_modules = 0;
  multi_objective_structure = "linear";
}

Simulation_Parameters_Topology_Optimization::~Simulation_Parameters_Topology_Optimization(){
}

void Simulation_Parameters_Topology_Optimization::input(){
  Simulation_Parameters::input();
  topology_optimization_on = true;
  //Simulation_Parameters::input();
  //initial buffer size for TO module list storage
  int buffer_size = 10;
  TO_Module_List.resize(buffer_size);
  TO_Function_Type.resize(buffer_size);
  Multi_Objective_Modules.resize(buffer_size);
  Multi_Objective_Weights.resize(buffer_size);
  Function_Arguments.resize(buffer_size);
  TO_Module_My_FEA_Module.resize(buffer_size);
  //use pushback to add arguments for each TO module
  
  //TO objectives and constraints
  /*
  TO_Module_List[nTO_modules] = "Strain_Energy_Minimize";
  TO_Function_Type[nTO_modules] = OBJECTIVE;
  nTO_modules++;
  */
  //Multi objective format
  TO_Module_List[nTO_modules] = "Multi_Objective";
  TO_Function_Type[nTO_modules] = OBJECTIVE;
  nTO_modules++;
  nmulti_objective_modules = 2;
  //add TO modules for objectives participating in the multi objective
  TO_Module_List[nTO_modules] = "Heat_Capacity_Potential_Minimize";
  TO_Function_Type[nTO_modules] = MULTI_OBJECTIVE_TERM;
  Multi_Objective_Modules[0] = nTO_modules;
  Multi_Objective_Weights[0] = 0.25;
  nTO_modules++;

  TO_Module_List[nTO_modules] = "Strain_Energy_Minimize";
  TO_Function_Type[nTO_modules] = MULTI_OBJECTIVE_TERM;
  Multi_Objective_Modules[1] = nTO_modules;
  Multi_Objective_Weights[1] = 0.75;
  nTO_modules++;

  /*
  TO_Module_List[nTO_modules] = "Thermo_Elastic_Strain_Energy_Minimize";
  TO_Function_Type[nTO_modules] = MULTI_OBJECTIVE_TERM;
  Multi_Objective_Modules[1] = nTO_modules;
  Multi_Objective_Weights[1] = 0.75;
  nTO_modules++;
  */

  //Constraints
  /*
  TO_Module_List[nTO_modules] = "Heat_Capacity_Potential_Constraint";
  TO_Function_Type[nTO_modules] = INEQUALITY_CONSTRAINT;
  Function_Arguments[nTO_modules].push_back(0);
  Function_Arguments[nTO_modules].push_back(1);
  Function_Arguments[nTO_modules].push_back(20);
  nTO_modules++;
  */
  TO_Module_List[nTO_modules] = "Mass_Constraint";
  TO_Function_Type[nTO_modules] = EQUALITY_CONSTRAINT;
  Function_Arguments[nTO_modules].push_back(0.15);
  nTO_modules++;
  
  TO_Module_List[nTO_modules] = "Moment_of_Inertia_Constraint";
  TO_Function_Type[nTO_modules] = EQUALITY_CONSTRAINT;
  Function_Arguments[nTO_modules].push_back(0.20);
  Function_Arguments[nTO_modules].push_back(1);
  nTO_modules++;
  
  TO_Module_List[nTO_modules] = "Moment_of_Inertia_Constraint";
  TO_Function_Type[nTO_modules] = EQUALITY_CONSTRAINT;
  Function_Arguments[nTO_modules].push_back(0);
  Function_Arguments[nTO_modules].push_back(3);
  nTO_modules++;


  //example for later
  if(nTO_modules==buffer_size){
    buffer_size += 10;
    TO_Module_List.resize(buffer_size);
    TO_Function_Type.resize(buffer_size);
    Multi_Objective_Modules.resize(buffer_size);
    Multi_Objective_Weights.resize(buffer_size);
    Function_Arguments.resize(buffer_size);
    TO_Module_My_FEA_Module.resize(buffer_size);
  }

  //simulation spatial dimension
  num_dim = 3;

  //polynomial interpolation order
  p_order = 0;

  //debug and performance report flags
  report_runtime_flag = true;

  //Topology Optimization flags
  nodal_density_flag = true;

  //Topology Optimization parameters
  penalty_power = 3;
  density_epsilon = 0.0001;

  // ---- boundary conditions ---- //
  NB = 6; // number of boundaries
  NBD = 2; //number of surface sets used to specify a fixed density

}

//==============================================================================
//    Communicate user settings from YAML file and apply to class members
//==============================================================================

void Simulation_Parameters_Topology_Optimization::apply_settings(){
  int buffer_size = TO_Module_List.size();
  std::string current_option;
  
  current_option = "optimization_options:optimization_process";
  if(set_options.find(current_option)!=set_options.end()){
    if(set_options[current_option]=="topology_optimization"){
      //std::cout << "FOUND TO SETTING" << std::endl;
      topology_optimization_on = true;
      set_options.erase(current_option);
    }
    else if(set_options[current_option]=="shape_optimization"){
      shape_optimization_on = true;
      set_options.erase(current_option);
    }
  }

  //read in TO parameters
  if(topology_optimization_on){
    current_option = "optimization_options:simp_penalty_power";
    if(set_options.find(current_option)!=set_options.end()){
      //std::cout << "FOUND TO SETTING" << std::endl;
      penalty_power = std::stod(set_options[current_option]);
      set_options.erase(current_option);
    }
    current_option = "optimization_options:density_epsilon";
    if(set_options.find(current_option)!=set_options.end()){
      //std::cout << "FOUND TO SETTING" << std::endl;
      density_epsilon = std::stod(set_options[current_option]);
      set_options.erase(current_option);
    }

  }
  nTO_modules = 0;
  current_option = "optimization_options:optimization_objective";
  if(set_options.find(current_option)!=set_options.end()){

    if(set_options[current_option]=="minimize_compliance"){
      TO_Module_List[nTO_modules] = "Strain_Energy_Minimize";
      TO_Function_Type[nTO_modules] = OBJECTIVE;
      nTO_modules++;
      set_options.erase(current_option);
    }

    else if(set_options[current_option]=="minimize_thermal_resistance"){
      TO_Module_List[nTO_modules] = "Heat_Capacity_Potential_Minimize";
      TO_Function_Type[nTO_modules] = OBJECTIVE;
      nTO_modules++;
      set_options.erase(current_option);
    }

    else if(set_options[current_option]=="multi-objective"){
      TO_Module_List[nTO_modules] = "Multi_Objective";
      TO_Function_Type[nTO_modules] = OBJECTIVE;
      nTO_modules++;

      if(set_options.find("optimization_options:multi-objective_structure")!=set_options.end()){
        if(set_options["optimization_options:multi-objective_structure"]=="linear"){
          multi_objective_structure = "linear";
          set_options.erase("optimization_options:multi-objective_structure");
        }
      }

      //read in multi-objective function definition and terms
      //read in multi-objective function definition and terms
      std::string optimization_module_base = "optimization_options:multi-objective_module_";
      std::string index;
      std::string optimization_module_name;
      nmulti_objective_modules = 0;

      index = std::to_string(nmulti_objective_modules+1);
      optimization_module_name = optimization_module_base + index;

      // --- set of user requested multi-objective terms---
      while(set_options.find(optimization_module_name+":type")!=set_options.end()){
  
        if(set_options[optimization_module_name+":type"]=="minimize_compliance"){
          TO_Module_List[nTO_modules] = "Strain_Energy_Minimize";
          set_options.erase(optimization_module_name+":type");
        }
        else if(set_options[optimization_module_name+":type"]=="minimize_thermal_resistance"){
          TO_Module_List[nTO_modules] = "Heat_Capacity_Potential_Minimize";
          set_options.erase(optimization_module_name+":type");
        }

        if(set_options.find(optimization_module_name+":weight_coefficient")!=set_options.end()){
          Multi_Objective_Weights[nmulti_objective_modules] = std::stod(set_options[optimization_module_name+":weight_coefficient"]);
          set_options.erase(optimization_module_name+":weight_coefficient");
        }
        Multi_Objective_Modules[nmulti_objective_modules] = nTO_modules;
        TO_Function_Type[nTO_modules] = MULTI_OBJECTIVE_TERM;
        nTO_modules++;
        if(nTO_modules==buffer_size){
          buffer_size += 10;
          TO_Module_List.resize(buffer_size);
          TO_Function_Type.resize(buffer_size);
          Multi_Objective_Modules.resize(buffer_size);
          Multi_Objective_Weights.resize(buffer_size);
          Function_Arguments.resize(buffer_size);
          TO_Module_My_FEA_Module.resize(buffer_size);
        }

        nmulti_objective_modules++;
        index = std::to_string(nmulti_objective_modules+1);
        optimization_module_name = optimization_module_base + index;
      }
      set_options.erase(current_option);
    }
  }
    
  int num_constraints;
  current_option = "optimization_options:num_optimization_constraint";
  if(set_options.find(current_option)!=set_options.end()){
      num_constraints = std::stoi(set_options[current_option]);
      set_options.erase(current_option);
  }

  //allocate constraint modules requested by the user
  std::string constraint_base = "optimization_options:constraint_";
  std::string index;
  std::string constraint_name;
  double constraint_value;

  // --- set of material specifications ---
  for(int icon = 0; icon < num_constraints; icon++){
      //readin material data
      index = std::to_string(icon+1);
      constraint_name = constraint_base + index;

      //expand storage if needed
      if(nTO_modules==buffer_size){
        buffer_size += 10;
        TO_Module_List.resize(buffer_size);
        TO_Function_Type.resize(buffer_size);
        Multi_Objective_Modules.resize(buffer_size);
        Multi_Objective_Weights.resize(buffer_size);
        Function_Arguments.resize(buffer_size);
        TO_Module_My_FEA_Module.resize(buffer_size);
      }

      //constraint request
      //function arguments for constraint; constraint value is the first
      current_option = constraint_name+":value";
      if(set_options.find(current_option)!=set_options.end()){
          constraint_value = std::stod(set_options[current_option]);
          Function_Arguments[nTO_modules].push_back(constraint_value); 
          set_options.erase(current_option);
      }

      //constraint type; some may have additional arguments that follow the value
      current_option = constraint_name+":type";
      if(set_options.find(current_option)!=set_options.end()){
          
          if(set_options[current_option]=="mass"){
            TO_Module_List[nTO_modules] = "Mass_Constraint";
            set_options.erase(current_option);
          }

          else if(set_options[current_option]=="moment_of_inertia"){
            TO_Module_List[nTO_modules] = "Moment_of_Inertia_Constraint";
            //obtain function arguments
            if(set_options.find(constraint_name+":component")!=set_options.end()){
              if(set_options[constraint_name+":component"]=="xx"){
                Function_Arguments[nTO_modules].push_back(0);
                set_options.erase(constraint_name+":component");
              }
              else if(set_options[constraint_name+":component"]=="yy"){
                Function_Arguments[nTO_modules].push_back(1);
                set_options.erase(constraint_name+":component");
              }
              else if(set_options[constraint_name+":component"]=="zz"){
                Function_Arguments[nTO_modules].push_back(2);
                set_options.erase(constraint_name+":component");
              }
              else if(set_options[constraint_name+":component"]=="xy"){
                Function_Arguments[nTO_modules].push_back(3);
                set_options.erase(constraint_name+":component");
              }
              else if(set_options[constraint_name+":component"]=="xz"){
                Function_Arguments[nTO_modules].push_back(4);
                set_options.erase(constraint_name+":component");
              }
              else if(set_options[constraint_name+":component"]=="yz"){
                Function_Arguments[nTO_modules].push_back(5);
                set_options.erase(constraint_name+":component");
              }
            }
            set_options.erase(current_option);
          }
      }

      //equality or inequality
      if(set_options.find(constraint_name+":relation")!=set_options.end()){
          if(set_options[constraint_name+":relation"]=="equality"){
            TO_Function_Type[nTO_modules] = EQUALITY_CONSTRAINT;
          }
          set_options.erase(constraint_name+":relation");
      }
      nTO_modules++;
  }
}

//==============================================================================
//    Detect pairings of TO Modules with FEA Modules
//==============================================================================

void Simulation_Parameters_Topology_Optimization::FEA_module_setup(){
  
  //initial buffer size for FEA module list storage
  bool module_found;
  int buffer_size = 10 + nfea_modules;
  FEA_Module_List.resize(buffer_size);
  fea_module_must_read.resize(buffer_size);
  int start_module = nfea_modules;
  
  for(int imodule = 0; imodule < nTO_modules; imodule++){
    module_found = false;
    //decides which FEA modules to setup based on user decided TO problem
    //automate selection list later; use std::map maybe?
    if(TO_Module_List[imodule] == "Strain_Energy_Minimize"){
      //check if module type was already allocated
      for(int ifea = 0; ifea < nfea_modules; ifea++){
        if(FEA_Module_List[ifea] == "Elasticity"){
          module_found = true;
          TO_Module_My_FEA_Module[imodule] = ifea;
        }
      }
      if(!module_found){
        *(solver_pointer_->fos) << "PROGRAM IS ENDING DUE TO ERROR; TOPOLOGY OPTIMIZATION FUNCTION (minimize_strain_energy)"
                              << "REQUIRES FEA MODULE (elasticity)" << std::endl;
        solver_pointer_->exit_solver(0);
      }
    }
    else if(TO_Module_List[imodule] == "Thermo_Elastic_Strain_Energy_Minimize"){

      //allocate required Heat Conduction Module if missing
      for(int ifea = 0; ifea < nfea_modules; ifea++){
        if(FEA_Module_List[ifea] == "Heat_Conduction"){
          module_found = true;
        }
      }
      if(!module_found){
          *(solver_pointer_->fos) << "PROGRAM IS ENDING DUE TO ERROR; TOPOLOGY OPTIMIZATION FUNCTION (minimize_strain_energy)"
                              << "REQUIRES FEA MODULE (steady_heat)" << std::endl;
          solver_pointer_->exit_solver(0);
      }
      module_found = false;
      //check if module type was already allocated
      for(int ifea = 0; ifea < nfea_modules; ifea++){
        if(FEA_Module_List[ifea] == "Thermo_Elasticity"){
          module_found = true;
          TO_Module_My_FEA_Module[imodule] = ifea;
        }
      }
      if(!module_found){
          *(solver_pointer_->fos) << "PROGRAM IS ENDING DUE TO ERROR; TOPOLOGY OPTIMIZATION FUNCTION (minimize_strain_energy)"
                              << "REQUIRES FEA MODULE (thermo_elasticity)" << std::endl;
          solver_pointer_->exit_solver(0);
      }

    }
    else if(TO_Module_List[imodule] == "Heat_Capacity_Potential_Minimize"){
      //check if module type was already allocated
      for(int ifea = 0; ifea < nfea_modules; ifea++){
        if(FEA_Module_List[ifea] == "Heat_Conduction"){
          module_found = true;
          TO_Module_My_FEA_Module[imodule] = ifea;
          
        }
      }
      if(!module_found){
        *(solver_pointer_->fos) << "PROGRAM IS ENDING DUE TO ERROR; TOPOLOGY OPTIMIZATION FUNCTION (minimize_mean_temperature)"
                              << "REQUIRES FEA MODULE (steady_heat)" << std::endl;
        solver_pointer_->exit_solver(0);
      }
    }
    else if(TO_Module_List[imodule] == "Mass_Constraint"){
      //check if module type was already allocated
      for(int ifea = 0; ifea < nfea_modules; ifea++){
        if(FEA_Module_List[ifea] == "Inertial"){
          module_found = true;
          TO_Module_My_FEA_Module[imodule] = ifea;
          
        }
      }
      if(!module_found){
        *(solver_pointer_->fos) << "PROGRAM IS ENDING DUE TO ERROR; TOPOLOGY OPTIMIZATION FUNCTION (mass)"
                              << "REQUIRES FEA MODULE (inertial)" << std::endl;
        solver_pointer_->exit_solver(0);
      }
    }
    else if(TO_Module_List[imodule] == "Moment_of_Inertia_Constraint"){
      //check if module type was already allocated
      for(int ifea = 0; ifea < nfea_modules; ifea++){
        if(FEA_Module_List[ifea] == "Inertial"){
          module_found = true;
          TO_Module_My_FEA_Module[imodule] = ifea;
        }
      }
      
      if(!module_found){
        *(solver_pointer_->fos) << "PROGRAM IS ENDING DUE TO ERROR; TOPOLOGY OPTIMIZATION FUNCTION (mass)"
                              << "REQUIRES FEA MODULE (inertial)" << std::endl;
        solver_pointer_->exit_solver(0);
      }
    }
    else if(TO_Module_List[imodule] == "Strain_Energy_Constraint"){
      //check if module type was already allocated
      for(int ifea = 0; ifea < nfea_modules; ifea++){
        if(FEA_Module_List[ifea] == "Elasticity"){
          module_found = true;
          TO_Module_My_FEA_Module[imodule] = ifea;
        }
      }
      
      if(!module_found){
        *(solver_pointer_->fos) << "PROGRAM IS ENDING DUE TO ERROR; TOPOLOGY OPTIMIZATION FUNCTION (minimize_strain_energy)"
                              << "REQUIRES FEA MODULE (elasticity)" << std::endl;
        solver_pointer_->exit_solver(0);
      }
    }
    else if(TO_Module_List[imodule] == "Heat_Capacity_Potential_Constraint"){
      //check if module type was already allocated
      for(int ifea = 0; ifea < nfea_modules; ifea++){
        if(FEA_Module_List[ifea] == "Heat_Conduction"){
          module_found = true;
          TO_Module_My_FEA_Module[imodule] = ifea;
        }
      }
      
      if(!module_found){
        *(solver_pointer_->fos) << "PROGRAM IS ENDING DUE TO ERROR; TOPOLOGY OPTIMIZATION FUNCTION (minimize_mean_temperature)"
                              << "REQUIRES FEA MODULE (steady_heat)" << std::endl;
        solver_pointer_->exit_solver(0);
      }
    }
    else if(TO_Module_List[imodule] == "Multi_Objective"){
      module_found = true;
    }
    else{
      *(solver_pointer_->fos) << "PROGRAM IS ENDING DUE TO ERROR; UNDEFINED TOPOLOGY OPTIMIZATION FUNCTION REQUESTED WITH NAME \""
                              <<TO_Module_List[imodule]<<"\" AT FEA MODULE PAIRING" << std::endl;
      solver_pointer_->exit_solver(0);
    }
    
    if(module_found){
      if(nfea_modules==buffer_size){
        buffer_size += 10;
        FEA_Module_List.resize(buffer_size);
        fea_module_must_read.resize(buffer_size);
      }
    }
  }

  //initialize
  for(int imodule = start_module; imodule < nfea_modules; imodule++){
    fea_module_must_read[imodule] = false;
  }
}
