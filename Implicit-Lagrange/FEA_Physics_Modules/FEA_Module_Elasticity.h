#ifndef FEA_MODULE_ELASTICITY_H
#define FEA_MODULE_ELASTICITY_H

#include "FEA_Module.h"

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

class FEA_Module_Elasticity: public FEA_Module{

public:
  FEA_Module_Elasticity(Implicit_Solver *Solver_Pointer);
  ~FEA_Module_Elasticity();
  
  //initialize data for boundaries of the model and storage for boundary conditions and applied loads
  void init_boundaries();

  //initializes memory for arrays used in the global stiffness matrix assembly
  void init_boundary_sets(int num_boundary_sets);

  void grow_boundary_sets(int num_boundary_sets);

  void grow_displacement_condition_sets(int num_boundary_sets);

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

  void compute_nodal_strains();
  
  void local_matrix(int ielem, CArrayKokkos<real_t, array_layout, device_type, memory_traits> &Local_Matrix);

  void local_matrix_multiply(int ielem, CArrayKokkos<real_t, array_layout, device_type, memory_traits> &Local_Matrix);

  void Element_Material_Properties(size_t ielem, real_t &Element_Modulus, real_t &Poisson_Ratio, real_t density);

  void Gradient_Element_Material_Properties(size_t ielem, real_t &Element_Modulus, real_t &Poisson_Ratio, real_t density);

  void Concavity_Element_Material_Properties(size_t ielem, real_t &Element_Modulus, real_t &Poisson_Ratio, real_t density);

  void Body_Term(size_t ielem, real_t density, real_t *forces);

  void Gradient_Body_Term(size_t ielem, real_t density, real_t *forces);

  void read_conditions_ansys_dat(std::ifstream *in, std::streampos before_condition_header);

  //interfaces between user input and creating data structures for bcs
  void generate_bcs();
  
  //interfaces between user input and creating data structures for applied loads
  void generate_applied_loads();

  void Displacement_Boundary_Conditions();

  void init_output();

  void compute_output();

  void collect_output(Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > global_reduce_map);

  void node_density_constraints(host_vec_array node_densities_lower_bound);
  
  class Simulation_Parameters_Elasticity *simparam;

  //Local FEA data
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Global_Stiffness_Matrix_Assembly_Map;
  RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits> Graph_Matrix; //stores global indices
  RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits> DOF_Graph_Matrix; //stores global indices
  RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout> Stiffness_Matrix;
  //CArrayKokkos<real_t, Kokkos::LayoutLeft, device_type, memory_traits> Nodal_Forces;
  CArrayKokkos<real_t, Kokkos::LayoutLeft, device_type, memory_traits> Nodal_Results; //result of linear solve; typically displacements and densities
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Stiffness_Matrix_Strides;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Graph_Matrix_Strides;
  RaggedRightArrayKokkos<real_t, array_layout, device_type, memory_traits> Original_Stiffness_Entries;
  RaggedRightArrayKokkos<LO, array_layout, device_type, memory_traits> Original_Stiffness_Entry_Indices;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Original_Stiffness_Entries_Strides;
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> Original_RHS_Entries;

  //Global FEA data
  Teuchos::RCP<MV> node_displacements_distributed;
  Teuchos::RCP<MV> node_strains_distributed;
  Teuchos::RCP<MV> all_node_displacements_distributed;
  Teuchos::RCP<MV> all_cached_node_displacements_distributed;
  Teuchos::RCP<MV> all_node_strains_distributed;
  bool adjoints_allocated;
  Teuchos::RCP<MV> adjoint_displacements_distributed;
  Teuchos::RCP<MV> adjoint_equation_RHS_distributed;
  Teuchos::RCP<MV> all_adjoint_displacements_distributed;
  Teuchos::RCP<MAT> Global_Stiffness_Matrix;
  Teuchos::RCP<MV> Global_Nodal_RHS;
  Teuchos::RCP<MV> Global_Nodal_Forces;
  
  //Boundary Conditions Data
  
  enum bc_type {NONE,DISPLACEMENT_CONDITION, X_DISPLACEMENT_CONDITION,
   Y_DISPLACEMENT_CONDITION, Z_DISPLACEMENT_CONDITION, POINT_LOADING_CONDITION, LINE_LOADING_CONDITION, SURFACE_LOADING_CONDITION, SURFACE_PRESSURE_CONDITION};
  int max_boundary_sets, max_disp_boundary_sets, max_load_boundary_sets;
  int num_surface_disp_sets, num_surface_force_sets;
  bool matrix_bc_reduced;
  
  //body force parameters
  bool gravity_flag, thermal_flag, electric_flag;
  real_t *gravity_vector;
  
  //stores the displacement value for the boundary condition on this nodal DOF
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> Node_DOF_Displacement_Boundary_Conditions;
  //stores applied point forces on nodal DOF
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> Node_DOF_Force_Boundary_Conditions;
  //constant surface force densities corresponding to each boundary set (provide varying field later)
  CArrayKokkos<real_t, array_layout, HostSpace, memory_traits> Boundary_Surface_Force_Densities;
  //constant displacement condition applied to all nodes on a boundary surface (convenient option to avoid specifying nodes)
  CArrayKokkos<real_t, array_layout, HostSpace, memory_traits> Boundary_Surface_Displacements;

  //linear solver parameters
  Teuchos::RCP<Teuchos::ParameterList> Linear_Solve_Params;

  //multigrid solver data and functions
  Teuchos::RCP<Xpetra::Matrix<real_t,LO,GO,node_type>> xA;
  Teuchos::RCP<Xpetra::MultiVector<real_t,LO,GO,node_type>> xX;
  Teuchos::RCP<MV> X;
  Teuchos::RCP<Xpetra::MultiVector<real_t,LO,GO,node_type>> xB;
  Teuchos::RCP<MueLu::Hierarchy<real_t,LO,GO,node_type>> H;
  Teuchos::RCP<Xpetra::Operator<real_t,LO,GO,node_type>> Prec;
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > local_reduced_dof_original_map;
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > all_reduced_dof_original_map;
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > local_reduced_dof_map;
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > local_balanced_reduced_dof_map;
  CArrayKokkos<GO, array_layout, device_type, memory_traits> Free_Indices;
  bool Hierarchy_Constructed;

  //output dof data
  //Global arrays with collected data used to print
  int collected_displacement_index, collected_strain_index, collected_stress_index;
  const_host_vec_array collected_node_strains;
  const_host_vec_array collected_node_stresses;
};

#endif // end HEADER_H
