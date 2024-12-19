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

#ifndef FEA_MODULE_ELASTICITY_H
#define FEA_MODULE_ELASTICITY_H

#include "FEA_Module.h"

//Trilinos
#include <Tpetra_MultiVector.hpp>

//forward declare
namespace MueLu{
  template<class floattype, class local_ind, class global_ind, class nodetype> 
  class Hierarchy;
}

namespace Anasazi{
  template<class floattype, class vectortype> 
  class Eigensolution;
}

namespace Xpetra{
  template<class floattype, class local_ind, class global_ind, class nodetype> 
  class Operator;
  
  template<class floattype, class local_ind, class global_ind, class nodetype> 
  class MultiVector;

  template<class floattype, class local_ind, class global_ind, class nodetype> 
  class Matrix;
}

class Implicit_Solver;

class Elasticity_Parameters;

class FEA_Module_Elasticity: public FEA_Module{

public:

  typedef Tpetra::Operator<real_t>                OP;
  typedef Tpetra::Vector<real_t,LO,GO> TpetraVector;

  FEA_Module_Elasticity(Elasticity_Parameters& in_params, Solver *Solver_Pointer, const int my_fea_module_index = 0);
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

  int eigensolve();

  void linear_solver_parameters();

  void comm_variables(Teuchos::RCP<const MV> zp);

  void update_linear_solve(Teuchos::RCP<const MV> zp, int compute_step);

  void compute_adjoint_gradients(const_host_vec_array design_densities, host_vec_array gradients);

  void compute_adjoint_hessian_vec(const_host_vec_array design_densities, host_vec_array hessvec, Teuchos::RCP<const MV> direction_vec_distributed);

  void compute_nodal_strains();
  
  void local_matrix(int ielem, CArrayKokkos<real_t, array_layout, device_type, memory_traits> &Local_Matrix);

  void local_matrix_multiply(int ielem, CArrayKokkos<real_t, array_layout, device_type, memory_traits> &Local_Matrix);

  void local_mass_matrix(int ielem, CArrayKokkos<real_t, array_layout, device_type, memory_traits> &Local_Matrix);

  void Element_Material_Properties(size_t ielem, real_t &Element_Modulus, real_t &Poisson_Ratio, real_t density);

  void Gradient_Element_Material_Properties(size_t ielem, real_t &Element_Modulus, real_t &Poisson_Ratio, real_t density);

  void Concavity_Element_Material_Properties(size_t ielem, real_t &Element_Modulus, real_t &Poisson_Ratio, real_t density);

  void Element_Anisotropic_Material_Properties(size_t ielem, real_t Element_Moduli[3], real_t Poisson_Ratios[3], real_t Shear_Moduli[3], real_t density);

  void Gradient_Element_Anisotropic_Material_Properties(size_t ielem, real_t Element_Moduli[3], real_t Poisson_Ratios[3], real_t Shear_Moduli[3], real_t density);

  void Concavity_Element_Anisotropic_Material_Properties(size_t ielem, real_t Element_Moduli[3], real_t Poisson_Ratios[3], real_t Shear_Moduli[3], real_t density);

  void Body_Term(size_t ielem, real_t density, real_t *forces);

  void Gradient_Body_Term(size_t ielem, real_t density, real_t *forces);

  void read_conditions_ansys_dat(std::ifstream *in, std::streampos before_condition_header); //ANSYS .dat import of specified load and boundary conditions

  void read_conditions_abaqus_inp(std::ifstream *in, std::streampos before_condition_header); //ABAQUS .inp import of specified load and boundary conditions

  //interfaces between user input and creating data structures for bcs
  void generate_bcs();
  
  //interfaces between user input and creating data structures for applied loads
  void generate_applied_loads();

  void Displacement_Boundary_Conditions();

  void init_output();

  void compute_output();
  
  void sort_output(Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > sorted_map);

  void collect_output(Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > global_reduce_map);

  void node_density_constraints(host_vec_array node_densities_lower_bound);

  //for the displacement constraint adjoint solves
  void compute_displacement_constraint_gradients(const_host_vec_array design_densities, const_host_vec_array target_displacements, const_host_int_array active_dofs, host_vec_array gradients);

  void compute_displacement_constraint_hessian_vec(const_host_vec_array design_densities, const_host_vec_array target_displacements, const_host_int_array active_dofs, host_vec_array hessvec, Teuchos::RCP<const MV> direction_vec_distributed);
  
  Elasticity_Parameters *module_params;
  Implicit_Solver *Implicit_Solver_Pointer_;
  
  //Local FEA data
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Global_Stiffness_Matrix_Assembly_Map;
  RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits> Graph_Matrix; //stores global indices
  RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits> DOF_Graph_Matrix; //stores global indices
  RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout> Stiffness_Matrix;
  RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout> Mass_Matrix;
  //CArrayKokkos<real_t, Kokkos::LayoutLeft, device_type, memory_traits> Nodal_Forces;
  CArrayKokkos<real_t, Kokkos::LayoutLeft, device_type, memory_traits> Nodal_Results; //result of linear solve; typically displacements and densities
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Stiffness_Matrix_Strides;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Graph_Matrix_Strides;
  RaggedRightArrayKokkos<real_t, array_layout, device_type, memory_traits> Original_Stiffness_Entries;
  RaggedRightArrayKokkos<real_t, array_layout, device_type, memory_traits> Original_Mass_Entries;
  RaggedRightArrayKokkos<LO, array_layout, device_type, memory_traits> Original_Stiffness_Entry_Indices;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Original_Stiffness_Entries_Strides;
  CArrayKokkos<real_t, array_layout, device_type, memory_traits> Original_RHS_Entries;

  //Global FEA data
  Teuchos::RCP<MV> node_displacements_distributed;
  Teuchos::RCP<MV> node_strains_distributed;
  Teuchos::RCP<MV> all_node_displacements_distributed;
  Teuchos::RCP<MV> all_cached_node_displacements_distributed;
  Teuchos::RCP<MV> all_node_strains_distributed;
  bool adjoints_allocated, constraint_adjoints_allocated;
  Teuchos::RCP<MV> adjoint_displacements_distributed;
  Teuchos::RCP<MV> psi_adjoint_vector_distributed;
  Teuchos::RCP<MV> phi_adjoint_vector_distributed;
  Teuchos::RCP<MV> adjoint_equation_RHS_distributed;
  Teuchos::RCP<MV> all_adjoint_displacements_distributed;
  Teuchos::RCP<MV> all_psi_adjoint_vector_distributed;
  Teuchos::RCP<MV> all_phi_adjoint_vector_distributed;
  Teuchos::RCP<MAT> Global_Stiffness_Matrix;
  Teuchos::RCP<MAT> Global_Mass_Matrix;
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
  Teuchos::RCP<MueLu::Hierarchy<real_t,LO,GO,node_type>> eigen_H;
  Teuchos::RCP<Xpetra::Operator<real_t,LO,GO,node_type>> eigen_Prec;
  bool Hierarchy_Constructed;
  bool Eigen_Hierarchy_Constructed;

  //Eigenvalue solution data
  Teuchos::RCP<Anasazi::Eigensolution<real_t,MV>> sol;
  Teuchos::RCP<MV> evecs;
  int numev;

  //output dof data
  //Global arrays with collected data used to print
  int output_displacement_index, output_strain_index, output_stress_index;
  Teuchos::RCP<MV> sorted_node_strains_distributed;
  Teuchos::RCP<MV> sorted_node_stresses_distributed;
  Teuchos::RCP<MV> sorted_node_displacements_distributed;
  Teuchos::RCP<MV> collected_node_displacements_distributed;
  Teuchos::RCP<MV> collected_node_strains_distributed;
  Teuchos::RCP<MV> collected_node_stresses_distributed;
};

#endif // end HEADER_H
