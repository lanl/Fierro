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

#ifndef FEA_MODULE_H
#define FEA_MODULE_H

#include <string>
#include <vector>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Kokkos_Core.hpp>
// #include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"
#include "utilities.h"
#include "matar.h"
#include "elements.h"
#include "node_combination.h"
// #include "Simulation_Parameters/Simulation_Parameters.h"

using namespace mtr;

// forward declare
class Solver;

// forward declarations
namespace ROL
{
template<class datatype>
class Problem;
} // namespace ROL

class Simulation_Parameters;
class FierroOptimizationObjective;

enum class FEA_MODULE_TYPE;
enum class BOUNDARY_TYPE;

class FEA_Module
{
public:
    FEA_Module(Solver* Solver_Pointer);
    virtual ~FEA_Module();

    // Trilinos type definitions
    typedef Tpetra::Map<>::local_ordinal_type LO;
    typedef Tpetra::Map<>::global_ordinal_type GO;

    typedef Tpetra::CrsMatrix<real_t, LO, GO> MAT;
    typedef const Tpetra::CrsMatrix<real_t, LO, GO> const_MAT;
    typedef Tpetra::MultiVector<real_t, LO, GO> MV;
    typedef Tpetra::MultiVector<GO, LO, GO> MCONN;

    typedef Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>::size_type SizeType;
    typedef Tpetra::Details::DefaultTypes::node_type node_type;
    using traits = Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>;

    using array_layout    = typename traits::array_layout;
    using execution_space = typename traits::execution_space;
    using device_type     = typename traits::device_type;
    using memory_traits   = typename traits::memory_traits;
    using global_size_t   = Tpetra::global_size_t;

    typedef Kokkos::View<real_t*, Kokkos::LayoutRight, device_type, memory_traits> values_array;
    typedef Kokkos::View<GO*, array_layout, device_type, memory_traits> global_indices_array;
    typedef Kokkos::View<LO*, array_layout, device_type, memory_traits> indices_array;
    // typedef Kokkos::View<SizeType*, array_layout, device_type, memory_traits> row_pointers;
    typedef MAT::local_graph_device_type::row_map_type::non_const_type row_pointers;
    // typedef Kokkos::DualView<real_t**, Kokkos::LayoutLeft, device_type>::t_dev vec_array;
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
    typedef Kokkos::View<const int**, array_layout, HostSpace, memory_traits> const_host_int_array;
    typedef Kokkos::View<const bool**, array_layout, HostSpace, memory_traits> const_host_bool_array;

    // initializes memory for arrays used in the global stiffness matrix assembly
    // initialize data for boundaries of the model and storage for boundary conditions and applied loads
    virtual void init_boundaries() {}

    // initializes memory for arrays used in the global stiffness matrix assembly
    virtual void init_boundary_sets(int num_boundary_sets) {}

    virtual void init_assembly() {}

    virtual void assemble_matrix() {}

    virtual void assemble_vector() {}

    // interfaces between user input and creating data structures for bcs
    virtual void generate_bcs() {}

    // interfaces between user input and creating data structures for applied loads
    virtual void generate_applied_loads() {}

    virtual void setup() {}

    virtual void module_cleanup() {}

    virtual int solve() { return 0; }

    virtual int eigensolve() { return 0; }

    virtual void read_conditions_ansys_dat(std::ifstream* in, std::streampos before_condition_header) {}

    virtual void read_conditions_abaqus_inp(std::ifstream* in, std::streampos before_condition_header) {}

    virtual void linear_solver_parameters() {}

    virtual void comm_variables(Teuchos::RCP<const MV> zp) {}

    virtual void comm_densities(Teuchos::RCP<const MV> zp);

    virtual void comm_filtered_densities();

    virtual void update_linear_solve(Teuchos::RCP<const MV> zp, int compute_step) {}

    virtual void update_forward_solve(Teuchos::RCP<const MV> zp) {}

    virtual void local_matrix(int ielem, CArrayKokkos<real_t, array_layout, device_type, memory_traits>& Local_Matrix) {}

    virtual void local_matrix_multiply(int ielem, CArrayKokkos<real_t, array_layout, device_type, memory_traits>& Local_Matrix) {}

    virtual void Element_Material_Properties(size_t ielem, real_t& Element_Modulus, real_t& Poisson_Ratio, real_t density) {}

    virtual void Gradient_Element_Material_Properties(size_t ielem, real_t& Element_Modulus, real_t& Poisson_Ratio, real_t density) {}

    virtual void Concavity_Element_Material_Properties(size_t ielem, real_t& Element_Modulus, real_t& Poisson_Ratio, real_t density) {}

    virtual void Body_Term(size_t ielem, real_t density, real_t* forces) {}

    virtual void Gradient_Body_Term(size_t ielem, real_t density, real_t* forces) {}

    virtual void tag_boundaries(int this_bc_tag, real_t val, int bdy_set, real_t* patch_limits = NULL);

    virtual int check_boundary(Node_Combination& Patch_Nodes, int this_bc_tag, real_t val, real_t* patch_limits);

    virtual void compute_output() {}

    virtual void output_control() {}

    virtual void write_data(std::map<std::string, const double*>&                    point_data_scalars_double,
                            std::map<std::string, const double*>&                    point_data_vectors_double,
                            std::map<std::string, const double*>&                    cell_data_scalars_double,
                            std::map<std::string, const int*>&                       cell_data_scalars_int,
                            std::map<std::string, std::pair<const double*, size_t>>& cell_data_fields_double) {}

    virtual void sort_output(Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> sorted_map) {} // node data outputs

    virtual void sort_element_output(Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> sorted_element_map) {} // element data outputs

    virtual void collect_output(Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> global_reduce_map) {}

    virtual void node_density_constraints(host_vec_array node_densities_lower_bound) {}

    virtual void compute_topology_optimization_adjoint() {} // Force does not depend on node coords and velocity

    virtual void compute_topology_optimization_adjoint_full() {} // Force depends on node coords and velocity

    virtual void compute_topology_optimization_gradient(const_vec_array design_densities, vec_array gradients) {}

    virtual void compute_topology_optimization_gradient_full(Teuchos::RCP<const MV> design_densities_distributed, Teuchos::RCP<MV> design_gradients_distributed) {}

    // interfacing information
    FEA_MODULE_TYPE Module_Type;
    int last_compute_step;

    // output stream
    Teuchos::RCP<Teuchos::FancyOStream> fos;

    std::shared_ptr<elements::element_selector> element_select;
    elements::Element3D* elem;
    elements::Element2D* elem2D;

    Simulation_Parameters* simparam;
    Solver* Solver_Pointer_;
    int     num_dim;
    int     num_gauss_points;
    int     my_fea_module_index_;

    // Local FEA data
    size_t nlocal_nodes;
    dual_vec_array       dual_node_coords; // coordinates of the nodes
    dual_vec_array       dual_node_densities; // topology optimization design variable
    dual_elem_conn_array dual_nodes_in_elem; // dual view of element connectivity to nodes
    host_elem_conn_array nodes_in_elem; // host view of element connectivity to nodes
    CArrayKokkos<elements::elem_types::elem_type, array_layout, HostSpace, memory_traits> Element_Types;
    CArrayKokkos<size_t, array_layout, HostSpace, memory_traits>   Nodes_Per_Element_Type;
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> corner_value_storage;
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> corner_vector_storage;
    CArrayKokkos<real_t, array_layout, device_type, memory_traits> corner_gradient_storage;

    // Ghost data on this MPI rank
    size_t nghost_nodes;
    CArrayKokkos<GO, Kokkos::LayoutLeft, node_type::device_type> ghost_nodes;
    CArrayKokkos<int, array_layout, device_type, memory_traits>  ghost_node_ranks;

    // Local FEA data including ghosts
    size_t nall_nodes;
    size_t rnum_elem, num_corners;

    // Global FEA data
    long long int num_nodes, num_elem;
    Teuchos::RCP<const Teuchos::Comm<int>>       comm;
    Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> map; // map of node indices
    Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> ghost_node_map; // map of node indices with ghosts on each rank
    Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> all_node_map; // map of node indices with ghosts on each rank
    Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> element_map; // non overlapping map of elements owned by each rank used in reduction ops
    Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> all_element_map; // overlapping map of elements connected to the local nodes in each rank
    Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> local_dof_map; // map of local dofs (typically num_node_local*num_dim)
    Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> all_dof_map; // map of local and ghost dofs (typically num_node_all*num_dim)
    Teuchos::RCP<MCONN>    global_nodes_in_elem_distributed; // element to node connectivity table
    Teuchos::RCP<MCONN>    node_nconn_distributed; // how many elements a node is connected to
    Teuchos::RCP<MV>       node_coords_distributed;
    Teuchos::RCP<MV>       all_node_coords_distributed;
    Teuchos::RCP<MV>       initial_node_coords_distributed;
    Teuchos::RCP<MV>       all_initial_node_coords_distributed;
    Teuchos::RCP<MV>       design_node_densities_distributed;
    Teuchos::RCP<MV>       design_node_coords_distributed;
    Teuchos::RCP<MV>       filtered_node_densities_distributed;
    Teuchos::RCP<const MV> test_node_densities_distributed;
    Teuchos::RCP<MV>       all_node_densities_distributed;
    Teuchos::RCP<MV>       all_design_node_coords_distributed;
    Teuchos::RCP<MV>       all_filtered_node_densities_distributed;
    Teuchos::RCP<MV>       Global_Element_Densities;

    // Boundary Conditions Data
    // CArray <Nodal_Combination> Patch_Nodes;
    size_t nboundary_patches;
    size_t num_boundary_conditions;
    int    current_bdy_id;
    bool   bcs_initialized;
    CArrayKokkos<Node_Combination, array_layout, HostSpace, memory_traits> Boundary_Patches;
    CArrayKokkos<size_t, array_layout, device_type, memory_traits> Boundary_Condition_Patches; // set of patches corresponding to each boundary condition
    CArrayKokkos<size_t, array_layout, device_type, memory_traits> NBoundary_Condition_Patches;
    CArrayKokkos<size_t, array_layout, device_type, memory_traits> Boundary_Condition_Patches_strides;

    // element selection parameters and data
    size_t max_nodes_per_element;

    // determines if rhs gets a contribution from bcs
    bool nonzero_bc_flag;

    // body force parameters
    bool    body_term_flag;
    real_t* gravity_vector;

    // lists what kind of boundary condition the nodal DOF is subjected to if any
    CArrayKokkos<int, array_layout, device_type, memory_traits> Node_DOF_Boundary_Condition_Type;
    // lists what kind of boundary condition each boundary set is assigned to
    CArrayKokkos<int, array_layout, HostSpace, memory_traits> Boundary_Condition_Type_List;

    // number of displacement boundary conditions acting on nodes; used to size the reduced global stiffness map
    size_t Number_DOF_BCS;

    // MPI data
    int      myrank; // index of this mpi rank in the world communicator
    int      nranks; // number of mpi ranks in the world communicator
    MPI_Comm world; // stores the default communicator object (MPI_COMM_WORLD)
    Teuchos::RCP<Tpetra::Import<LO, GO>> importer; // all node comms
    Teuchos::RCP<Tpetra::Import<LO, GO>> ghost_importer; // ghost node comms
    Teuchos::RCP<Tpetra::Import<LO, GO>> node_sorting_importer; // sorted node comms
    Teuchos::RCP<Tpetra::Import<LO, GO>> element_sorting_importer; // sorted element comms
    Teuchos::RCP<Tpetra::Import<LO, GO>> dof_importer; // ghost dof comms

    // ! mapping used to get local ghost index from the global ID.
    // typedef ::Tpetra::Details::FixedHashTable<GO, LO, Kokkos::HostSpace::device_type>
    // global_to_local_table_host_type;

    // global_to_local_table_host_type global2local_map;
    // CArrayKokkos<int, Kokkos::LayoutLeft, Kokkos::HostSpace::device_type> active_ranks;

    // Pertains to local mesh information being stored as prescribed by the row map
    global_size_t min_gid;
    global_size_t max_gid;
    global_size_t index_base;

    // allocation flags to avoid repeat MV and global matrix construction
    int Matrix_alloc;

    // debug flags
    int gradient_print_sync;

    // Topology Optimization parameter
    int  penalty_power;
    bool nodal_density_flag;

    // runtime and counters for performance output
    double linear_solve_time, hessvec_time, hessvec_linear_time;
    int    update_count, hessvec_count;

    // nodal DOF output data
    enum vector_styles { NODAL, DOF }; // multivector can store as ndof by 1 or nnode by vector_size
    int noutput;
    int displacement_index;
    std::vector<std::vector<std::string>> output_dof_names;
    std::vector<const_host_vec_array>     module_outputs;
    std::vector<vector_styles> vector_style;
    std::vector<int> output_vector_sizes;

    // Pointer to ROL Problem for optimization solves
    Teuchos::RCP<ROL::Problem<real_t>> problem;

    // Explicit BC data kept for now until bc data is refactored/consolidated
    // node ids in bdy_patch set
    RaggedRightArrayKokkos<size_t> bdy_nodes_in_set;
    DCArrayKokkos<size_t> num_bdy_nodes_in_set;
    bool node_specified_bcs; // currently happens with ansys import

    // patch ids in bdy set
    size_t num_bdy_sets;
    DynamicRaggedRightArrayKokkos<size_t> bdy_patches_in_set;

    //optimization data
    FierroOptimizationObjective* objective_function;
};

#endif // end HEADER_H
