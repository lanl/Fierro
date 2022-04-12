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
  penalty_power = 3;
  nTO_modules = 0;
}

Simulation_Parameters_Topology_Optimization::~Simulation_Parameters_Topology_Optimization(){
}

void Simulation_Parameters_Topology_Optimization::input(){
  Simulation_Parameters::input();
  //Simulation_Parameters::input();
  //initial buffer size for TO module list storage
  int buffer_size = 10;
  TO_Module_List = std::vector<std::string>(buffer_size);
  TO_Function_Type = std::vector<function_type>(buffer_size);
  Multi_Objective_Modules = std::vector<int>(buffer_size);
  Multi_Objective_Weights = std::vector<real_t>(buffer_size);
  Function_Arguments = std::vector<std::vector<real_t>>(buffer_size);
  //use pushback to add arguments for each TO module
  
  //TO objectives and constraints
  //TO_Module_List[0] = "Strain_Energy_Minimize";
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
  Function_Arguments[nTO_modules].push_back(0.12);
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

  // ---- boundary conditions ---- //
  NB = 6; // number of boundaries
  NBD = 2; //number of surface sets used to specify a fixed density

}


void Simulation_Parameters_Topology_Optimization::FEA_module_setup(){
  
  //initial buffer size for FEA module list storage
  bool module_found;
  int buffer_size = 10;
  FEA_Module_List = std::vector<std::string>(buffer_size);
  TO_Module_My_FEA_Module = std::vector<int>(buffer_size);
  fea_module_must_read = std::vector<bool>(buffer_size);
  
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
        TO_Module_My_FEA_Module[imodule] = nfea_modules;
        FEA_Module_List[nfea_modules++] = "Elasticity";
        module_found = true;
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
        TO_Module_My_FEA_Module[imodule] = nfea_modules;
        FEA_Module_List[nfea_modules++] = "Heat_Conduction";
        module_found = true;
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
        TO_Module_My_FEA_Module[imodule] = nfea_modules;
        FEA_Module_List[nfea_modules++] = "Inertial";
        module_found = true;
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
        TO_Module_My_FEA_Module[imodule] = nfea_modules;
        FEA_Module_List[nfea_modules++] = "Inertial";
        module_found = true;
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
        TO_Module_My_FEA_Module[imodule] = nfea_modules;
        FEA_Module_List[nfea_modules++] = "Elasticity";
        module_found = true;
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
        TO_Module_My_FEA_Module[imodule] = nfea_modules;
        FEA_Module_List[nfea_modules++] = "Heat_Conduction";
        module_found = true;
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
        TO_Module_My_FEA_Module.resize(buffer_size);
        fea_module_must_read.resize(buffer_size);
      }
    }
  }

  //initialize
  for(int imodule = 0; imodule < nfea_modules; imodule++){
    fea_module_must_read[imodule] = false;
}