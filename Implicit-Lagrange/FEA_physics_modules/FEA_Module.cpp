#include "FEA_Module.h"
#include "Implicit_Solver.h"
using namespace utils;

FEA_Module::FEA_Module(Implicit_Solver *Solver_Pointer){

  //obtain global and local node and element counts
  num_nodes = Solver_Pointer->num_nodes;
  num_elem = Solver_Pointer->num_elem;
  nall_nodes = Solver_Pointer->nall_nodes;
  rnum_elem = Solver_Pointer->rnum_elem;

  hessvec_count = update_count = 0;
  file_index = 0;
  linear_solve_time = hessvec_time = hessvec_linear_time = 0;

  Matrix_alloc=0;
  gradient_print_sync = 0;
  //RCP pointer to *this (Parallel Nonlinear Solver Object)
  //FEM_pass = Teuchos::rcp(this);

  //Trilinos output stream
  std::ostream &out = std::cout;
  fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  (*fos).setOutputToRootOnly(0);
  
  //obtain node and element maps
  comm = Solver_Pointer->comm;
  map = Solver_Pointer->map; //map of node indices
  ghost_node_map = Solver_Pointer->ghost_node_map; //map of node indices with ghosts on each rank
  all_node_map = Solver_Pointer->all_node_map; //map of node indices with ghosts on each rank
  element_map = Solver_Pointer->element_map; //non overlapping map of elements owned by each rank used in reduction ops
  all_element_map = Solver_Pointer->all_element_map; //overlapping map of elements connected to the local nodes in each rank
  local_dof_map = Solver_Pointer->local_dof_map; //map of local dofs (typically num_node_local*num_dim)
  all_dof_map = Solver_Pointer->all_dof_map //map of local and ghost dofs (typically num_node_all*num_dim)

  //obtain mesh coordinates, densities, and element connectivity
  nodes_in_elem_distributed = Solver_Pointer->nodes_in_elem_distributed; //element to node connectivity table
  node_nconn_distributed = Solver_Pointer->node_nconn_distributed; //how many elements a node is connected to
  node_coords_distributed = Solver_Pointer->node_coords_distributed;
  all_node_coords_distributed = Solver_Pointer->all_node_coords_distributed;
  design_node_densities_distributed = Solver_Pointer->design_node_densities_distributed;
  test_node_densities_distributed = Solver_Pointer->test_node_densities_distributed;
  all_node_densities_distributed = Solver_Pointer->all_node_densities_distributed;
  Global_Element_Densities = Solver_Pointer->Global_Element_Densities;

  //obtain boundary condition and loading data

}

FEA_Module::~FEA_Module() {}