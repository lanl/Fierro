#ifndef PARALLEL_NONLINEAR_SOLVER_H
#define PARALLEL_NONLINEAR_SOLVER_H

#include "utilities.h"
#include "../Solver.h"
#include "matar.h"
#include "elements.h"
#include "node_combination.h"
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_Parallel_Reduce.hpp>
//#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"

namespace swage{
  class mesh_t;
}

namespace elements{
  class element_selector;
  class Element3D;
  class Element2D;
  class ref_element;
}

class Parallel_Nonlinear_Solver: public Solver{

public:
  Parallel_Nonlinear_Solver();
  ~Parallel_Nonlinear_Solver();

  //Trilinos type definitions
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;

  typedef Tpetra::CrsMatrix<real_t,LO,GO> MAT;
  typedef const Tpetra::CrsMatrix<real_t,LO,GO> const_MAT;
  typedef Tpetra::MultiVector<real_t,LO,GO> MV;
  typedef Tpetra::MultiVector<GO,LO,GO> MCONN;

  typedef Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>::size_type SizeType;
  typedef Tpetra::Details::DefaultTypes::node_type node_type;
  using traits = Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>;
  
  using array_layout    = typename traits::array_layout;
  using execution_space = typename traits::execution_space;
  using device_type     = typename traits::device_type;
  using memory_traits   = typename traits::memory_traits;
  using global_size_t = Tpetra::global_size_t;
  
  typedef Kokkos::View<real_t*, Kokkos::LayoutRight, device_type, memory_traits> values_array;
  typedef Kokkos::View<GO*, array_layout, device_type, memory_traits> global_indices_array;
  typedef Kokkos::View<LO*, array_layout, device_type, memory_traits> indices_array;
  typedef Kokkos::View<SizeType*, array_layout, device_type, memory_traits> row_pointers;

  //typedef Kokkos::DualView<real_t**, Kokkos::LayoutLeft, device_type>::t_dev vec_array;
  typedef MV::dual_view_type::t_dev vec_array;
  typedef MV::dual_view_type::t_host host_vec_array;
  typedef Kokkos::View<const real_t**, array_layout, HostSpace, memory_traits> const_host_vec_array;
  typedef Kokkos::View<const real_t**, array_layout, device_type, memory_traits> const_vec_array;
  typedef MV::dual_view_type dual_vec_array;
  typedef MCONN::dual_view_type dual_elem_conn_array;
  typedef MCONN::dual_view_type::t_host host_elem_conn_array;
  typedef MCONN::dual_view_type::t_dev elem_conn_array;
  typedef Kokkos::View<const GO**, array_layout, HostSpace, memory_traits> const_host_elem_conn_array;
  typedef Kokkos::View<const GO**, array_layout, device_type, memory_traits> const_elem_conn_array;

  void run(int argc, char *argv[]);

  void read_mesh_ensight(char *MESH);

  void read_mesh_tecplot(char *MESH);
  
  //setup ghosts and element maps
  void init_maps();
  
  //initializes memory for arrays used in the global stiffness matrix assembly
  void init_assembly();

  void init_design();

  void assemble_matrix();

  void assemble_vector();

  int solve();

  void linear_solver_parameters();

  void comm_variables(Teuchos::RCP<const MV> zp);

  void update_linear_solve(Teuchos::RCP<const MV> zp);

  void collect_information();

  void compute_element_volumes();

  void compute_element_masses(const_host_vec_array design_densities, bool max_flag);

  void compute_element_moments(const_host_vec_array design_densities, bool max_flag, int moment_component);

  void compute_element_moments_of_inertia(const_host_vec_array design_densities, bool max_flag, int inertia_component);

  void compute_nodal_gradients(const_host_vec_array design_densities, host_vec_array gradients);

  void compute_moment_gradients(const_host_vec_array design_densities, host_vec_array gradients, int moment_component);

  void compute_moment_of_inertia_gradients(const_host_vec_array design_densities, host_vec_array gradients, int intertia_component);

  void compute_adjoint_gradients(const_host_vec_array design_densities, host_vec_array gradients);

  void compute_nodal_strains();

  void setup_optimization_problem();

  void local_matrix(int ielem, CArray <real_t> &Local_Matrix);

  void local_matrix_multiply(int ielem, CArray <real_t> &Local_Matrix);
  
  //interfaces between user input and creating data structures for bcs
  void generate_bcs();
  
  //finds the boundary element surfaces in this model
  void Get_Boundary_Patches();

  //void vtk_writer();

  //void ensight_writer();

  void tecplot_writer();

  void Element_Material_Properties(size_t, real_t &Element_Modulus, real_t &Poisson_Ratio, real_t density);

  void Gradient_Element_Material_Properties(size_t, real_t &Element_Modulus, real_t &Poisson_Ratio, real_t density);

  void Displacement_Boundary_Conditions();

  void init_boundary_sets(int num_boundary_sets);

  void tag_boundaries(int this_bc_tag, real_t val, int bdy_set);

  int check_boundary(Node_Combination &Patch_Nodes, int this_bc_tag, real_t val);\

  //debug and performance functions/variables
  double CPU_Time();
  void init_clock();
  double initial_CPU_time;
  int update_count;
  
  swage::mesh_t *init_mesh;
  swage::mesh_t *mesh;
  
  elements::element_selector *element_select;
  elements::Element3D *elem;
  elements::Element2D *elem2D;
  elements::ref_element  *ref_elem;
  

  class Simulation_Parameters *simparam;
  
  //Local FEA data
  size_t nlocal_nodes;
  dual_vec_array dual_node_coords; //coordinates of the nodes
  dual_vec_array dual_node_displacements; //first three indices of second dim should be positions
  dual_vec_array dual_node_densities; //topology optimization design variable
  dual_vec_array dual_nodal_forces;
  dual_elem_conn_array dual_nodes_in_elem; //dual view of element connectivity to nodes
  host_elem_conn_array nodes_in_elem; //host view of element connectivity to nodes
  CArray<elements::elem_types::elem_type> Element_Types;
  CArray<size_t> Nodes_Per_Element_Type;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Global_Stiffness_Matrix_Assembly_Map;
  RaggedRightArrayKokkos<size_t, array_layout, device_type, memory_traits> Graph_Matrix; //stores global indices
  RaggedRightArrayKokkos<GO, array_layout, device_type, memory_traits> DOF_Graph_Matrix; //stores global indices
  RaggedRightArrayKokkos<real_t, Kokkos::LayoutRight, device_type, memory_traits, array_layout> Stiffness_Matrix;
  //CArrayKokkos<real_t, Kokkos::LayoutLeft, device_type, memory_traits> Nodal_Forces;
  CArrayKokkos<real_t, Kokkos::LayoutLeft, device_type, memory_traits> Nodal_Results; //result of linear solve; typically displacements and densities
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Stiffness_Matrix_Strides;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Graph_Matrix_Strides;

  //Ghost data on this MPI rank
  size_t nghost_nodes;
  CArrayKokkos<GO, Kokkos::LayoutLeft, node_type::device_type> ghost_nodes;
  CArrayKokkos<int, array_layout, device_type, memory_traits> ghost_node_ranks;

  //Local FEA data including ghosts
  size_t nall_nodes;
  size_t rnum_elem;
  dual_vec_array dual_all_node_coords; //coordinates of the nodes including ghosts
  dual_vec_array dual_all_node_displacements; //coordinates of the nodes including ghosts
  dual_vec_array dual_all_node_densities; //includes ghost data of the topology optimization design variable

  //Global FEA data
  long long int num_nodes, num_elem;
  Teuchos::RCP<const Teuchos::Comm<int> > comm;
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > map; //map of node indices
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > ghost_node_map; //map of node indices with ghosts on each rank
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > all_node_map; //map of node indices with ghosts on each rank
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > element_map; //non overlapping map of elements owned by each rank used in reduction ops
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > all_element_map; //overlapping map of elements connected to the local nodes in each rank
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > local_dof_map; //map of local dofs (typically num_node_local*num_dim)
  Teuchos::RCP<Tpetra::Map<LO,GO,node_type> > all_dof_map; //map of local and ghost dofs (typically num_node_all*num_dim)
  Teuchos::RCP<MCONN> nodes_in_elem_distributed; //element to node connectivity table
  Teuchos::RCP<MCONN> node_nconn_distributed; //how many elements a node is connected to
  Teuchos::RCP<MV> node_coords_distributed;
  Teuchos::RCP<MV> node_displacements_distributed;
  Teuchos::RCP<MV> node_strains_distributed;
  Teuchos::RCP<MV> all_node_coords_distributed;
  Teuchos::RCP<MV> all_node_displacements_distributed;
  Teuchos::RCP<MV> all_node_strains_distributed;
  Teuchos::RCP<MV> design_node_densities_distributed;
  Teuchos::RCP<const MV> test_node_densities_distributed;
  Teuchos::RCP<MV> all_node_densities_distributed;
  Teuchos::RCP<MAT> Global_Stiffness_Matrix;
  Teuchos::RCP<MV> Global_Nodal_Forces;
  Teuchos::RCP<MV> lower_bound_node_densities_distributed;
  Teuchos::RCP<MV> upper_bound_node_densities_distributed;
  Teuchos::RCP<MV> mass_gradients_distributed;
  Teuchos::RCP<MV> Global_Element_Densities_Upper_Bound;
  Teuchos::RCP<MV> Global_Element_Densities_Lower_Bound;
  Teuchos::RCP<MV> Global_Element_Densities;
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

  //Global arrays with collected data used to print
  const_host_vec_array collected_node_coords;
  const_host_vec_array collected_node_displacements;
  const_host_vec_array collected_node_densities;
  const_host_vec_array collected_node_strains;
  const_host_elem_conn_array collected_nodes_in_elem;
  
  //Boundary Conditions Data
  //CArray <Nodal_Combination> Patch_Nodes;
  size_t nboundary_patches;
  size_t num_boundary_conditions;
  CArrayKokkos<Node_Combination, array_layout, device_type, memory_traits> Boundary_Patches;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Boundary_Condition_Patches; //set of patches corresponding to each boundary condition
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> NBoundary_Condition_Patches;
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Boundary_Condition_Patches_strides;

  //pointer to FEA solver object passed to objectives and constraints
  Teuchos::RCP<Parallel_Nonlinear_Solver> FEM_pass;

  //element selection parameters and data
  size_t max_nodes_per_element;

  //types of boundary conditions
  enum bc_type {NONE,DISPLACEMENT_CONDITION, X_DISPLACEMENT_CONDITION,
   Y_DISPLACEMENT_CONDITION, Z_DISPLACEMENT_CONDITION, LOADING_CONDITION};

  //lists what kind of boundary condition the nodal DOF is subjected to if any
  CArrayKokkos<int, array_layout, device_type, memory_traits> Node_DOF_Boundary_Condition_Type;
  //stores the displacement value for the boundary condition on this nodal DOF
  CArray<real_t> Node_DOF_Displacement_Boundary_Conditions;
  //stores applied point forces on nodal DOF
  CArray<real_t> Node_DOF_Force_Boundary_Conditions;
  //lists what kind of boundary condition each boundary set is assigned to
  CArray<int> Boundary_Condition_Type_List;
  //constant surface force densities corresponding to each boundary set (provide varying field later)
  CArray<real_t> Boundary_Surface_Force_Densities;
  //constant displacement condition applied to all nodes on a boundary surface (convenient option to avoid specifying nodes)
  CArray<real_t> Boundary_Surface_Displacements;
  
  //number of displacement boundary conditions acting on nodes; used to size the reduced global stiffness map
  size_t Number_DOF_BCS;

  //MPI data
  int myrank; //index of this mpi rank in the world communicator
  int nranks; //number of mpi ranks in the world communicator
  MPI_Comm world; //stores the default communicator object (MPI_COMM_WORLD)

  //! mapping used to get local ghost index from the global ID.
  //typedef ::Tpetra::Details::FixedHashTable<GO, LO, Kokkos::HostSpace::device_type>
    //global_to_local_table_host_type;

  //global_to_local_table_host_type global2local_map;
  //CArrayKokkos<int, Kokkos::LayoutLeft, Kokkos::HostSpace::device_type> active_ranks;

  //Pertains to local mesh information being stored as prescribed by the row map
  global_size_t local_nrows;
  global_size_t min_gid;
  global_size_t max_gid;
  global_size_t index_base;

  //allocation flags to avoid repeat MV and global matrix construction
  int Matrix_alloc;

  //file readin variables
  std::ifstream *in;
  int words_per_line, elem_words_per_line;

  //file output variables
  int file_index, nsteps_print;  //file sequence index and print frequency in # of optimization steps

  //debug flags
  int gradient_print_sync;

  //linear solver parameters
  Teuchos::RCP<Teuchos::ParameterList> Linear_Solve_Params;

  //inertial properties
  real_t mass, center_of_mass[3], moments_of_inertia[6];

  //runtime flags
  bool mass_init, com_init[3];

  //update counters (first attempt at reducing redundant calls through ROL for Moments of Inertia and Center of Mass)
  int mass_update, com_update[3];
  int mass_gradient_update, com_gradient_update[3];
  
};

#endif // end HEADER_H
