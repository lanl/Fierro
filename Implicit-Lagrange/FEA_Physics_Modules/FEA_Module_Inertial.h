#ifndef FEA_MODULE_INERTIAL_H
#define FEA_MODULE_INERTIAL_H

#include "FEA_Module.h"

class FEA_Module_Inertial: public FEA_Module{

public:
  FEA_Module_Inertial(Implicit_Solver *Solver_Pointer);
  ~FEA_Module_Inertial();

  void comm_variables(Teuchos::RCP<const MV> zp);

  void compute_element_volumes();

  void compute_element_masses(const_host_vec_array design_densities, bool max_flag);

  void compute_element_moments(const_host_vec_array design_densities, bool max_flag, int moment_component);

  void compute_element_moments_of_inertia(const_host_vec_array design_densities, bool max_flag, int inertia_component);

  void compute_nodal_gradients(const_host_vec_array design_densities, host_vec_array gradients);

  void compute_moment_gradients(const_host_vec_array design_densities, host_vec_array gradients, int moment_component);

  void compute_moment_of_inertia_gradients(const_host_vec_array design_densities, host_vec_array gradients, int intertia_component);
  
  //forward declare
  class Simulation_Parameters_Inertial *simparam;

  //Global FEA data
  Teuchos::RCP<MV> mass_gradients_distributed;
  Teuchos::RCP<MV> center_of_mass_gradients_distributed;
  Teuchos::RCP<MV> Global_Element_Volumes;
  Teuchos::RCP<MV> Global_Element_Masses;
  Teuchos::RCP<MV> Global_Element_Moments_x;
  Teuchos::RCP<MV> Global_Element_Moments_y;
  Teuchos::RCP<MV> Global_Element_Moments_z;
  Teuchos::RCP<MV> Global_Element_Moments_of_Inertia_xx;
  Teuchos::RCP<MV> Global_Element_Moments_of_Inertia_yy;
  Teuchos::RCP<MV> Global_Element_Moments_of_Inertia_zz;
  Teuchos::RCP<MV> Global_Element_Moments_of_Inertia_xy;
  Teuchos::RCP<MV> Global_Element_Moments_of_Inertia_xz;
  Teuchos::RCP<MV> Global_Element_Moments_of_Inertia_yz;

  //inertial properties
  real_t mass, center_of_mass[3], moments_of_inertia[6];

  //runtime flags
  bool mass_init, com_init[3];

  //update counters (first attempt at reducing redundant calls through ROL for Moments of Inertia and Center of Mass)
  int mass_update, com_update[3];
  int mass_gradient_update, com_gradient_update[3];

};

#endif // end HEADER_H
