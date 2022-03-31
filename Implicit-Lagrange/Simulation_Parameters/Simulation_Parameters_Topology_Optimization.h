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
