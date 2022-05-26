#ifndef FEA_MODULE_HEAT_CONDUCTION_H
#define FEA_MODULE_HEAT_CONDUCTION_H

#include "FEA_Module.h"

//forward declare linear solver classes
namespace MueLu{
  template<class floattype, class local_ind, class global_ind, class nodetype> 
  class Hierarchy;
}

namespace Xpetra{
  template<class floattype, class local_ind, class global_ind, class nodetype> 
  class Operator;
  
  template<class floattype, class local_ind, class global_ind, class nodetype> 
  class MultiVector;

  template<class floattype, class local_ind, class global_ind, class nodetype> 
  class Matrix;
}

class FEA_Module_Heat_Conduction: public FEA_Module{

public:
  FEA_Module_Heat_Conduction(Implicit_Solver *Solver_Pointer);
  ~FEA_Module_Heat_Conduction();
  
  //initialize data for boundaries of the model and storage for boundary conditions and applied loads
  void init_boundaries();

  //initializes memory for arrays used in the global stiffness matrix assembly
  void init_boundary_sets(int num_boundary_sets);

  void grow_boundary_sets(int num_boundary_sets);

  void grow_temperature_condition_sets(int num_boundary_sets);

  void grow_loading_condition_sets(int num_boundary_sets);
  
  void init_assembly();

  void assemble_matrix();

  void assemble_vector();

  int solve();

  void linear_solver_parameters();

  void comm_variables(Teuchos::RCP<const MV> zp);

  void update_linear_solve(Teuchos::RCP<const MV> zp);

  void compute_adjoint_gradients(const_host_vec_array design_densities, host_vec_array gradients);

  void compute_adjoint_hessian_vec(const_host_vec_array design_densities, host_vec_array hessvec, Teuchos::RCP<const MV> direction_vec_distributed);

  void compute_nodal_heat_fluxes();

  void local_matrix(int ielem, CArrayKokkos<real_t, array_layout, device_type, memory_traits> &Local_Matrix);

  void Element_Material_Properties(size_t ielem, real_t &Element_Conductivity, real_t density);

  void Gradient_Element_Material_Properties(size_t ielem, real_t &Element_Conductivity_Derivative, real_t density);

  void Concavity_Element_Material_Properties(size_t ielem, real_t &Element_Conductivity_Derivative, real_t density);

  void Body_Term(size_t ielem, real_t density, real_t &specific_internal_energy_rate);

  void Gradient_Body_Term(size_t ielem, real_t density, real_t &gradient_specific_internal_energy_rate);

  //interfaces between user input and creating data structures for bcs
  void generate_bcs();
  
  //interfaces between user input and creating data structures for applied loads
  void generate_applied_loads();

  void Temperature_Boundary_Conditions();

  void init_output();

  void compute_output();

  void collect_output(Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > global_reduce_map);

  void node_density_constraints(host_vec_array node_densities_lower_bound);

  class Simulation_Parameters_Thermal *simparam;
  
  //Local FEA data
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Global_Conductivity_Matrix_Assembly_Map;
  RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits> Graph_Matrix; //stores global indices
  RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits> DOF_Graph_Matrix; //stores global indices
  RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout> Conductivity_Matrix;
  //CArrayKokkos<real_t, Kokkos::LayoutLeft, device_type, memory_traits> Nodal_Forces;
  CArrayKokkos<real_t, Kokkos::LayoutLeft, device_type, memory_traits> Nodal_Results; //result of linear solve; typically displacements and densities
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Conductivity_Matrix_Strides;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Graph_Matrix_Strides;
  RaggedRightArrayKokkos<real_t, array_layout, device_type, memory_traits> Original_Conductivity_Entries;
  RaggedRightArrayKokkos<LO, array_layout, device_type, memory_traits> Original_Conductivity_Entry_Indices;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Original_Conductivity_Entries_Strides;
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> Original_RHS_Entries;

  //Global FEA data
  Teuchos::RCP<MV> node_temperatures_distributed;
  Teuchos::RCP<MV> node_heat_fluxes_distributed;
  Teuchos::RCP<MV> node_temperature_gradients_distributed;
  Teuchos::RCP<MV> all_node_temperatures_distributed;
  Teuchos::RCP<MV> all_cached_node_temperatures_distributed;
  Teuchos::RCP<MV> all_node_temperature_gradients_distributed;
  Teuchos::RCP<MV> all_node_heat_fluxes_distributed;
  bool adjoints_allocated;
  Teuchos::RCP<MV> adjoint_temperatures_distributed;
  Teuchos::RCP<MV> adjoint_equation_RHS_distributed;
  Teuchos::RCP<MV> all_adjoint_temperatures_distributed;
  Teuchos::RCP<MV> Global_Nodal_Heat;
  Teuchos::RCP<MV> Global_Nodal_RHS;
  Teuchos::RCP<MAT> Global_Conductivity_Matrix;
  
  //Boundary Conditions Data
  
  enum bc_type {NONE,TEMPERATURE_CONDITION, POINT_LOADING_CONDITION, LINE_LOADING_CONDITION, SURFACE_LOADING_CONDITION};
  int max_boundary_sets, max_temp_boundary_sets, max_load_boundary_sets;
  int num_surface_temp_sets, num_surface_flux_sets;
  bool matrix_bc_reduced;
  
  //body force parameters
  bool thermal_flag, electric_flag;
  
  //stores the displacement value for the boundary condition on this nodal DOF
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> Node_Temperature_Boundary_Conditions;
  //stores applied point forces on nodal DOF
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> Node_Heat_Flux_Boundary_Conditions;
  //constant surface force densities corresponding to each boundary set (provide varying field later)
  CArrayKokkos<real_t, array_layout, HostSpace, memory_traits> Boundary_Surface_Heat_Flux;
  //constant displacement condition applied to all nodes on a boundary surface (convenient option to avoid specifying nodes)
  CArrayKokkos<real_t, array_layout, HostSpace, memory_traits> Boundary_Surface_Temperatures;

  //linear solver parameters
  Teuchos::RCP<Teuchos::ParameterList> Linear_Solve_Params;

  //multigrid solver data and functions
  Teuchos::RCP<Xpetra::Matrix<real_t,LO,GO,node_type>> xA;
  Teuchos::RCP<Xpetra::MultiVector<real_t,LO,GO,node_type>> xX;
  Teuchos::RCP<MV> X;
  Teuchos::RCP<Xpetra::MultiVector<real_t,LO,GO,node_type>> xB;
  Teuchos::RCP<MueLu::Hierarchy<real_t,LO,GO,node_type>> H;
  Teuchos::RCP<Xpetra::Operator<real_t,LO,GO,node_type>> Prec;
  bool Hierarchy_Constructed;

  //output dof data
  //Global arrays with collected data used to print
  int collected_temperature_index, collected_temperature_gradient_index, collected_heat_flux_index;
  const_host_vec_array collected_node_heat_fluxes;
  const_host_vec_array collected_node_temperatures;
  
};

#endif // end HEADER_H
