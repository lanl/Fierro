#ifndef SIMULATION_PARAMETERS_TOPOLOGY_OPTIMIZATION_H
#define SIMULATION_PARAMETERS_TOPOLOGY_OPTIMIZATION_H

#include "utilities.h"
#include "Simulation_Parameters.h"
using namespace utils;

class Simulation_Parameters_Topology_Optimization : public Simulation_Parameters
{
 public:
  Simulation_Parameters_Topology_Optimization();
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

  //volumes to hold density constant

  //list of TO functions needed by problem
  
  //types of TO functions
  enum function_type {OBJECTIVE, EQUALITY_CONSTRAINT, INEQUALITY_CONSTRAINT, VECTOR_EQUALITY_CONSTRAINT, VECTOR_INEQUALITY_CONSTRAINT};
  std::vector<std::string> TO_Module_List;
  std::vector<function_type> TO_Function_Type;
  std::vector<int> TO_Module_My_FEA_Module;
  std::vector<std::vector<real_t>> Function_Arguments;
  int nTO_modules;
};

#endif // end HEADER_H
