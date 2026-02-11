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

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>  // fmin, fmax, abs note: fminl is long
#include <sys/stat.h>
#include <mpi.h>
#include <chrono>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"
#include "Tpetra_Details_FixedHashTable.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Import_Util2.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include <set>

#include "elements.h"
#include "swage.h"
#include "matar.h"
#include "utilities.h"
#include "node_combination.h"
#include "Simulation_Parameters/FEA_Module/Elasticity_Parameters.h"
#include "Simulation_Parameters/FEA_Module/Eulerian_Parameters.h"
#include "FEA_Module_Eulerian.h"
#include "Explicit_Solver_Eulerian.h"

// optimization
#include "ROL_Algorithm.hpp"
#include "ROL_Solver.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Types.hpp"
#include "ROL_Elementwise_Reduce.hpp"
#include "ROL_Stream.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_StdBoundConstraint.hpp"
#include "ROL_ParameterList.hpp"
#include <ROL_TpetraMultiVector.hpp>
#include "Kinetic_Energy_Minimize.h"

#define MAX_ELEM_NODES 8
#define STRAIN_EPSILON 0.000000001
#define DENSITY_EPSILON 0.0001
#define BC_EPSILON 1.0e-6
#define BUFFER_GROW 100

using namespace utils;

FEA_Module_Eulerian::FEA_Module_Eulerian(Solver* Solver_Pointer, mesh_t& mesh, const int my_fea_module_index) : FEA_Module(Solver_Pointer), mesh(mesh)
{
    // assign interfacing index
    my_fea_module_index_ = my_fea_module_index;

    // recast solver pointer for non-base class access
    Explicit_Solver_Pointer_ = dynamic_cast<Explicit_Solver_Eulerian*>(Solver_Pointer);

    // create parameter object
    simparam = &Explicit_Solver_Pointer_->simparam;
    // ---- Read input file, define state and boundary conditions ---- //
    // simparam->input();

    // TO parameters
    simparam_dynamic_opt = Explicit_Solver_Pointer_->simparam_dynamic_opt;

    // create ref element object
    // ref_elem = new elements::ref_element();
    // create mesh objects
    // init_mesh = new swage::mesh_t(simparam);
    // mesh = new swage::mesh_t(simparam);

    // boundary condition data
    max_boundary_sets = 0;
    Local_Index_Boundary_Patches = Explicit_Solver_Pointer_->Local_Index_Boundary_Patches;

    // set Tpetra vector pointers
    initial_node_velocities_distributed = Explicit_Solver_Pointer_->initial_node_velocities_distributed;
    initial_node_coords_distributed     = Explicit_Solver_Pointer_->initial_node_coords_distributed;
    all_initial_node_coords_distributed = Explicit_Solver_Pointer_->all_initial_node_coords_distributed;
    node_coords_distributed         = Explicit_Solver_Pointer_->node_coords_distributed;
    node_velocities_distributed     = Explicit_Solver_Pointer_->node_velocities_distributed;
    all_node_velocities_distributed = Explicit_Solver_Pointer_->all_node_velocities_distributed;
    if (simparam_dynamic_opt.topology_optimization_on || simparam_dynamic_opt.shape_optimization_on) {
        all_cached_node_velocities_distributed = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
        force_gradient_velocity    = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
        force_gradient_position    = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
        force_gradient_design      = Teuchos::rcp(new MV(all_node_map, 1));
        corner_value_storage       = Solver_Pointer->corner_value_storage;
        corner_vector_storage      = Solver_Pointer->corner_vector_storage;
        relative_element_densities = DCArrayKokkos<double>(rnum_elem, "relative_element_densities");
    }

    if (simparam_dynamic_opt.topology_optimization_on || simparam_dynamic_opt.shape_optimization_on || simparam->num_dims == 2) {
        node_masses_distributed = Teuchos::rcp(new MV(map, 1));
        ghost_node_masses_distributed  = Teuchos::rcp(new MV(ghost_node_map, 1));
        adjoint_vector_distributed     = Teuchos::rcp(new MV(map, simparam->num_dims));
        phi_adjoint_vector_distributed = Teuchos::rcp(new MV(map, simparam->num_dims));
    }

    // setup output
    noutput = 0;
    init_output();

    // optimization flags
    kinetic_energy_objective = false;

    // set parameters
    Time_Variables tv = simparam->time_variables;
    time_value        = simparam->time_value;
    time_final        = tv.time_final;
    dt_max            = tv.dt_max;
    dt_min            = tv.dt_min;
    dt_cfl            = tv.dt_cfl;
    graphics_time     = simparam->graphics_options.graphics_time;
    graphics_cyc_ival = simparam->graphics_options.graphics_cyc_ival;
    graphics_dt_ival  = simparam->graphics_options.graphics_dt_ival;
    cycle_stop        = tv.cycle_stop;
    rk_num_stages     = simparam->rk_num_stages;
    dt    = tv.dt;
    fuzz  = tv.fuzz;
    tiny  = tv.tiny;
    small = tv.small;
    graphics_times = simparam->graphics_options.graphics_times;
    graphics_id    = simparam->graphics_options.graphics_id;

    if (simparam_dynamic_opt.topology_optimization_on) {
        max_time_steps = BUFFER_GROW;
        forward_solve_velocity_data = Teuchos::rcp(new std::vector<Teuchos::RCP<MV>>(max_time_steps + 1));
        time_data.resize(max_time_steps + 1);
        forward_solve_coordinate_data = Teuchos::rcp(new std::vector<Teuchos::RCP<MV>>(max_time_steps + 1));
        adjoint_vector_data     = Teuchos::rcp(new std::vector<Teuchos::RCP<MV>>(max_time_steps + 1));
        phi_adjoint_vector_data = Teuchos::rcp(new std::vector<Teuchos::RCP<MV>>(max_time_steps + 1));
        // assign a multivector of corresponding size to each new timestep in the buffer
        for (int istep = 0; istep < max_time_steps + 1; istep++) {
            (*forward_solve_velocity_data)[istep]   = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
            (*forward_solve_coordinate_data)[istep] = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
            (*adjoint_vector_data)[istep]     = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
            (*phi_adjoint_vector_data)[istep] = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
        }
    }
}

FEA_Module_Eulerian::~FEA_Module_Eulerian()
{
    // delete simparam;
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn read_conditions_ansys_dat
///
/// \brief Read ANSYS dat format mesh file
///
/// \param Input file path
/// \param Before conditions
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Eulerian::read_conditions_ansys_dat(std::ifstream* in, std::streampos before_condition_header)
{
    char ch;
    std::string skip_line, read_line, substring, token;
    std::stringstream line_parse, line_parse2;

    int num_dim = simparam->num_dims;
    int buffer_lines = 1000;
    int max_word     = 30;
    int p_order = simparam->p_order;
    int local_node_index, current_column_index;
    int buffer_loop, buffer_iteration, buffer_iterations, scan_loop, nodes_per_element, words_per_line;

    real_t unit_scaling = simparam->get_unit_scaling();

    size_t strain_count;
    size_t read_index_start, node_rid, elem_gid;

    CArrayKokkos<char, array_layout, HostSpace, memory_traits> read_buffer;
    CArrayKokkos<long long int, array_layout, HostSpace, memory_traits> read_buffer_indices;

    LO     local_dof_id;
    GO     node_gid;
    real_t dof_value;
    host_vec_array node_densities;
} // end read_conditions_ansys_dat

/////////////////////////////////////////////////////////////////////////////
///
/// \fn init_boundaries
///
/// \brief Initialize sets of element boundary surfaces and arrays for input conditions
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Eulerian::init_boundaries()
{
    max_boundary_sets = simparam->NB;
    int num_dim = simparam->num_dims;

    // set the number of boundary sets
    if (myrank == 0) {
        std::cout << "building boundary sets " << std::endl;
    }

    // initialize to 1 since there must be at least 1 boundary set anyway; read in may occure later
    if (max_boundary_sets == 0) {
        max_boundary_sets = 1;
    }
    // std::cout << "NUM BOUNDARY CONDITIONS ON RANK " << myrank << " FOR INIT " << num_boundary_conditions <<std::endl;
    init_boundary_sets(max_boundary_sets);

    // allocate nodal data
    Node_DOF_Boundary_Condition_Type = CArrayKokkos<int, array_layout, device_type, memory_traits>(nall_nodes * num_dim, "Node_DOF_Boundary_Condition_Type");

    // initialize
    for (int init = 0; init < nall_nodes * num_dim; init++) {
        Node_DOF_Boundary_Condition_Type(init) = NONE;
    }

    Number_DOF_BCS = 0;
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn init_boundary_sets
///
/// \brief initialize storage for element boundary surfaces corresponding
///        to user BCs
///
/// \param Number of boundary sets
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Eulerian::init_boundary_sets(int num_sets)
{
    if (num_sets == 0) {
        std::cout << " Warning: number of boundary conditions = 0";
        return;
    }
    // initialize maximum
    max_boundary_sets = num_sets;
    // std::cout << " DEBUG PRINT "<<num_sets << " " << nboundary_patches << std::endl;
    Boundary_Condition_Type_List = CArrayKokkos<int, array_layout, HostSpace, memory_traits>(num_sets, "Boundary_Condition_Type_List");
    NBoundary_Condition_Patches  = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_sets, "NBoundary_Condition_Patches");
    // std::cout << "NBOUNDARY PATCHES ON RANK " << myrank << " FOR INIT IS " << nboundary_patches <<std::endl;
    Boundary_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(num_sets, nboundary_patches, "Boundary_Condition_Patches");

    // initialize data
    for (int iset = 0; iset < num_sets; iset++) {
        NBoundary_Condition_Patches(iset) = 0;
    }

    // initialize
    for (int ibdy = 0; ibdy < num_sets; ibdy++) {
        Boundary_Condition_Type_List(ibdy) = NONE;
    }
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn grow_boundary_sets
///
/// \brief Grow boundary conditions sets of element boundary surfaces
///
/// \param Number of boundary sets
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Eulerian::grow_boundary_sets(int num_sets)
{
    int num_dim = simparam->num_dims;

    if (num_sets == 0) {
        std::cout << " Warning: number of boundary conditions being set to 0";
        return;
    }

    // std::cout << " DEBUG PRINT "<<num_sets << " " << nboundary_patches << std::endl;
    if (num_sets > max_boundary_sets) {
        // temporary storage for previous data
        CArrayKokkos<int, array_layout, HostSpace, memory_traits> Temp_Boundary_Condition_Type_List     = Boundary_Condition_Type_List;
        CArrayKokkos<size_t, array_layout, device_type, memory_traits> Temp_NBoundary_Condition_Patches = NBoundary_Condition_Patches;
        CArrayKokkos<size_t, array_layout, device_type, memory_traits> Temp_Boundary_Condition_Patches  = Boundary_Condition_Patches;

        max_boundary_sets = num_sets + 5; // 5 is an arbitrary buffer
        Boundary_Condition_Type_List = CArrayKokkos<int, array_layout, HostSpace, memory_traits>(max_boundary_sets, "Boundary_Condition_Type_List");
        NBoundary_Condition_Patches  = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(max_boundary_sets, "NBoundary_Condition_Patches");
        // std::cout << "NBOUNDARY PATCHES ON RANK " << myrank << " FOR GROW " << nboundary_patches <<std::endl;
        Boundary_Condition_Patches = CArrayKokkos<size_t, array_layout, device_type, memory_traits>(max_boundary_sets, nboundary_patches, "Boundary_Condition_Patches");

        // copy previous data back over
        // std::cout << "NUM BOUNDARY CONDITIONS ON RANK " << myrank << " FOR COPY " << max_boundary_sets <<std::endl;
        for (int iset = 0; iset < num_boundary_conditions; iset++) {
            Boundary_Condition_Type_List(iset) = Temp_Boundary_Condition_Type_List(iset);
            NBoundary_Condition_Patches(iset)  = Temp_NBoundary_Condition_Patches(iset);
            for (int ipatch = 0; ipatch < nboundary_patches; ipatch++) {
                Boundary_Condition_Patches(iset, ipatch) = Temp_Boundary_Condition_Patches(iset, ipatch);
            }
        }

        // initialize data
        for (int iset = num_boundary_conditions; iset < max_boundary_sets; iset++) {
            NBoundary_Condition_Patches(iset) = 0;
        }

        // initialize
        for (int ibdy = num_boundary_conditions; ibdy < max_boundary_sets; ibdy++) {
            Boundary_Condition_Type_List(ibdy) = NONE;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn generate_bcs
///
/// \brief Assign sets of element boundary surfaces corresponding to user BCs
///
/// Unused in this context
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Eulerian::generate_bcs()
{
} // end generate_bcs

/////////////////////////////////////////////////////////////////////////////
///
/// \fn generate_bcs
///
/// \brief Loop through applied boundary conditions and tag node ids to remove
///        necessary rows and columns from the assembled linear system
///
/// Unused in this context
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Eulerian::Displacement_Boundary_Conditions()
{
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn init_output
///
/// \brief Initialize output data structures
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Eulerian::init_output()
{
    // check user parameters for output
    bool output_velocity_flag = simparam->output_options.output_velocity;
    bool output_strain_flag   = simparam->output_options.output_strain;
    bool output_stress_flag   = simparam->output_options.output_stress;
    int  num_dim = simparam->num_dims;
    int  Brows;
    if (num_dim == 3) {
        Brows = 6;
    }
    else{
        Brows = 3;
    }

    if (output_velocity_flag) {
        // displacement_index is accessed by writers at the solver level for deformed output
        output_velocity_index = noutput;
        noutput += 1;
        module_outputs.resize(noutput);

        vector_style.resize(noutput);
        vector_style[noutput - 1] = DOF;

        output_vector_sizes.resize(noutput);
        output_vector_sizes[noutput - 1] = num_dim;

        output_dof_names.resize(noutput);
        output_dof_names[noutput - 1].resize(num_dim);
        output_dof_names[noutput - 1][0] = "vx";
        output_dof_names[noutput - 1][1] = "vy";
        if (num_dim == 3) {
            output_dof_names[noutput - 1][2] = "vz";
        }
    }
    if (output_strain_flag) {
        output_strain_index = noutput;
        noutput += 1;
        module_outputs.resize(noutput);

        vector_style.resize(noutput);
        vector_style[noutput - 1] = NODAL;

        output_vector_sizes.resize(noutput);
        output_vector_sizes[noutput - 1] = Brows;

        output_dof_names.resize(noutput);
        output_dof_names[noutput - 1].resize(Brows);
        if (num_dim == 2) {
            output_dof_names[noutput - 1][0] = "strain_xx";
            output_dof_names[noutput - 1][1] = "strain_yy";
            output_dof_names[noutput - 1][2] = "strain_xy";
        }
        if (num_dim == 3) {
            output_dof_names[noutput - 1][0] = "strain_xx";
            output_dof_names[noutput - 1][1] = "strain_yy";
            output_dof_names[noutput - 1][2] = "strain_zz";
            output_dof_names[noutput - 1][3] = "strain_xy";
            output_dof_names[noutput - 1][4] = "strain_xz";
            output_dof_names[noutput - 1][5] = "strain_yz";
        }
    }
    if (output_stress_flag) {
        output_stress_index = noutput;
        noutput += 1;
        module_outputs.resize(noutput);

        vector_style.resize(noutput);
        vector_style[noutput - 1] = NODAL;

        output_vector_sizes.resize(noutput);
        output_vector_sizes[noutput - 1] = Brows;

        output_dof_names.resize(noutput);
        output_dof_names[noutput - 1].resize(Brows);
        if (num_dim == 2) {
            output_dof_names[noutput - 1][0] = "stress_xx";
            output_dof_names[noutput - 1][1] = "stress_yy";
            output_dof_names[noutput - 1][3] = "stress_xy";
        }
        if (num_dim == 3) {
            output_dof_names[noutput - 1][0] = "stress_xx";
            output_dof_names[noutput - 1][1] = "stress_yy";
            output_dof_names[noutput - 1][2] = "stress_zz";
            output_dof_names[noutput - 1][3] = "stress_xy";
            output_dof_names[noutput - 1][4] = "stress_xz";
            output_dof_names[noutput - 1][5] = "stress_yz";
        }
    }
}

/* -------------------------------------------------------------------------------------------
   Prompts sorting for elastic response output data. For now, nodal strains.
---------------------------------------------------------------------------------------------- */

void FEA_Module_Eulerian::sort_output(Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> sorted_map)
{
}

/* -------------------------------------------------------------------------------------------
   Prompts computation of elastic response output data. For now, nodal strains.
---------------------------------------------------------------------------------------------- */

void FEA_Module_Eulerian::collect_output(Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> global_reduce_map)
{
}

/* -------------------------------------------------------------------------------------------
   Prompts computation of elastic response output data. For now, nodal strains.
---------------------------------------------------------------------------------------------- */

void FEA_Module_Eulerian::compute_output()
{
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn comm_node_masses
///
/// \brief Communicate updated nodal velocities to ghost nodes
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Eulerian::comm_node_masses()
{
    // debug print of design vector
    // std::ostream &out = std::cout;
    // Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    // if(myrank==0)
    // *fos << "Density data :" << std::endl;
    // node_densities_distributed->describe(*fos,Teuchos::VERB_EXTREME);
    // *fos << std::endl;
    // std::fflush(stdout);

    // communicate design densities
    // create import object using local node indices map and all indices map
    // Tpetra::Import<LO, GO> importer(map, ghost_node_map);

    // comms to get ghosts
    ghost_node_masses_distributed->doImport(*node_masses_distributed, *ghost_importer, Tpetra::INSERT);
    // all_node_map->describe(*fos,Teuchos::VERB_EXTREME);
    // all_node_velocities_distributed->describe(*fos,Teuchos::VERB_EXTREME);

    // update_count++;
    // if(update_count==1){
    // MPI_Barrier(world);
    // MPI_Abort(world,4);
    // }
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn comm_adjoint_vectors
///
/// \brief  Communicate updated nodal adjoint vectors to ghost nodes
///
/// \param Optimization iteration cycle count
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Eulerian::comm_adjoint_vectors(int cycle)
{
    // debug print of design vector
    // std::ostream &out = std::cout;
    // Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    // if(myrank==0)
    // *fos << "Density data :" << std::endl;
    // node_densities_distributed->describe(*fos,Teuchos::VERB_EXTREME);
    // *fos << std::endl;
    // std::fflush(stdout);

    // communicate design densities
    // create import object using local node indices map and all indices map
    // Tpetra::Import<LO, GO> importer(map, all_node_map);

    // comms to get ghosts
    (*adjoint_vector_data)[cycle]->doImport(*adjoint_vector_distributed, *importer, Tpetra::INSERT);
    (*phi_adjoint_vector_data)[cycle]->doImport(*phi_adjoint_vector_distributed, *importer, Tpetra::INSERT);
    // all_node_map->describe(*fos,Teuchos::VERB_EXTREME);
    // all_node_velocities_distributed->describe(*fos,Teuchos::VERB_EXTREME);

    // update_count++;
    // if(update_count==1){
    // MPI_Barrier(world);
    // MPI_Abort(world,4);
    // }
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn comm_variables
///
/// \brief Communicate ghosts using the current optimization design data
///
/// \param Design variables
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Eulerian::comm_variables(Teuchos::RCP<const MV> zp)
{
    if (simparam_dynamic_opt.topology_optimization_on) {
        // set density vector to the current value chosen by the optimizer
        test_node_densities_distributed = zp;

        // debug print of design vector
        // std::ostream &out = std::cout;
        // Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
        // if(myrank==0)
        // *fos << "Density data :" << std::endl;
        // node_densities_distributed->describe(*fos,Teuchos::VERB_EXTREME);
        // *fos << std::endl;
        // std::fflush(stdout);

        // communicate design densities
        // create import object using local node indices map and all indices map
        // Tpetra::Import<LO, GO> importer(map, all_node_map);

        // comms to get ghosts
        all_node_densities_distributed->doImport(*test_node_densities_distributed, *importer, Tpetra::INSERT);
    }
    else if (simparam_dynamic_opt.shape_optimization_on) {
        // clause to communicate boundary node data if the boundary nodes are ghosts on this rank
    }
}

/* -------------------------------------------------------------------------------------------
   enforce constraints on nodes due to BCS
---------------------------------------------------------------------------------------------- */

void FEA_Module_Eulerian::node_density_constraints(host_vec_array node_densities_lower_bound)
{
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn setup
///
/// \brief Setup Eulerian solver data
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Eulerian::setup()
{
    const size_t rk_level      = simparam->rk_num_bins - 1;
    const size_t num_fills     = simparam->region_options.size();
    const size_t rk_num_bins   = simparam->rk_num_bins;
    const size_t num_bcs       = simparam->boundary_conditions.size();
    const size_t num_materials = simparam->material_options.size();
    const int    num_dim       = simparam->num_dims;

    const DCArrayKokkos<mat_fill_t> mat_fill = simparam->mat_fill;
    const DCArrayKokkos<boundary_t> boundary = simparam->boundary;
    const DCArrayKokkos<material_t> material = simparam->material;
    eos_global_vars = simparam->eos_global_vars;
    strength_global_vars = simparam->strength_global_vars;
    eos_state_vars = DCArrayKokkos<double>(rnum_elem, simparam->max_num_eos_state_vars);
    strength_state_vars = DCArrayKokkos<double>(rnum_elem, simparam->max_num_strength_state_vars);

    // --- calculate bdy sets ---//
    mesh.num_nodes_in_patch  = 2 * (num_dim - 1); // 2 (2D) or 4 (3D)
    mesh.num_patches_in_elem = 2 * num_dim; // 4 (2D) or 6 (3D)

    mesh.init_bdy_sets(num_bcs);
    num_bdy_sets = mesh.num_bdy_sets;
    printf("Num BC's = %lu\n", num_bcs);

    // patch ids in bdy set
    bdy_patches_in_set = mesh.bdy_patches_in_set;
    if (num_dim == 2) {
        bdy_nodes = mesh.bdy_nodes;
    }

    // node ids in bdy_patch set
    bdy_nodes_in_set     = mesh.bdy_nodes_in_set;
    num_bdy_nodes_in_set = mesh.num_bdy_nodes_in_set;

    // assign mesh views needed by the FEA module

    // elem ids in elem
    elems_in_elem     = mesh.elems_in_elem;
    num_elems_in_elem = mesh.num_elems_in_elem;

    // corners
    num_corners_in_node = mesh.num_corners_in_node;
    corners_in_node     = mesh.corners_in_node;
    corners_in_elem     = mesh.corners_in_elem;

    // elem-node conn & node-node conn
    elems_in_node = mesh.elems_in_node;
    if (num_dim == 2) {
        nodes_in_node     = mesh.nodes_in_node;
        num_nodes_in_node = mesh.num_nodes_in_node;
        // patch conn

        patches_in_elem = mesh.patches_in_elem;
        nodes_in_patch  = mesh.nodes_in_patch;
        elems_in_patch  = mesh.elems_in_patch;
    }

    // initialize if topology optimization is used
    if (simparam_dynamic_opt.topology_optimization_on || simparam_dynamic_opt.shape_optimization_on) {
        // create parameter object
        simparam_elasticity = Simulation_Parameters_Elasticity();
        init_assembly();
        assemble_matrix();
    }

    return;
} // end of setup

/////////////////////////////////////////////////////////////////////////////
///
/// \fn euler_solve
///
/// \brief Eulerian solver loop
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Eulerian::euler_solve()
{
    Time_Variables tv = simparam->time_variables;

    const size_t rk_level = simparam->rk_num_bins - 1;
    time_value        = simparam->time_value;
    time_final        = tv.time_final;
    dt_max            = tv.dt_max;
    dt_min            = tv.dt_min;
    dt_cfl            = tv.dt_cfl;
    graphics_time     = simparam->graphics_options.graphics_time;
    graphics_cyc_ival = simparam->graphics_options.graphics_cyc_ival;
    graphics_dt_ival  = simparam->graphics_options.graphics_dt_ival;
    cycle_stop        = tv.cycle_stop;
    rk_num_stages     = simparam->rk_num_stages;
    dt    = tv.dt;
    fuzz  = tv.fuzz;
    tiny  = tv.tiny;
    small = tv.small;
    graphics_times = simparam->graphics_options.graphics_times;
    graphics_id    = simparam->graphics_options.graphics_id;
    size_t num_bdy_nodes = mesh.num_bdy_nodes;
    const DCArrayKokkos<boundary_t> boundary = simparam->boundary;
    const DCArrayKokkos<material_t> material = simparam->material;
    int       nTO_modules;
    int       old_max_forward_buffer;
    size_t    cycle;
    const int num_dim = simparam->num_dims;
    real_t    objective_accumulation, global_objective_accumulation;
    std::vector<std::vector<int>> FEA_Module_My_TO_Modules = simparam_dynamic_opt.FEA_Module_My_TO_Modules;
    problem = Explicit_Solver_Pointer_->problem; // Pointer to ROL optimization problem object
    ROL::Ptr<ROL::Objective<real_t>> obj_pointer;

    // reset time accumulating objective and constraints
    /*
    for(int imodule = 0 ; imodule < FEA_Module_My_TO_Modules[my_fea_module_index_].size(); imodule++){
    current_module_index = FEA_Module_My_TO_Modules[my_fea_module_index_][imodule];
    //test if module needs reset
    if(){

    }
    }
    */
    // simple setup to just request KE for now; above loop to be expanded and used later for scanning modules
    if (simparam_dynamic_opt.topology_optimization_on) {
        obj_pointer = problem->getObjective();
        KineticEnergyMinimize_TopOpt& kinetic_energy_minimize_function = dynamic_cast<KineticEnergyMinimize_TopOpt&>(*obj_pointer);
        kinetic_energy_minimize_function.objective_accumulation = 0;
        global_objective_accumulation = objective_accumulation = 0;
        kinetic_energy_objective = true;
        if (max_time_steps + 1 > forward_solve_velocity_data->size()) {
            old_max_forward_buffer = forward_solve_velocity_data->size();
            time_data.resize(max_time_steps + 1);
            forward_solve_velocity_data->resize(max_time_steps + 1);
            forward_solve_coordinate_data->resize(max_time_steps + 1);
            adjoint_vector_data->resize(max_time_steps + 1);
            phi_adjoint_vector_data->resize(max_time_steps + 1);
            // assign a multivector of corresponding size to each new timestep in the buffer
            for (int istep = old_max_forward_buffer; istep < max_time_steps + 1; istep++) {
                (*forward_solve_velocity_data)[istep]   = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
                (*forward_solve_coordinate_data)[istep] = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
                (*adjoint_vector_data)[istep]     = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
                (*phi_adjoint_vector_data)[istep] = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
            }
        }
    }

    if (simparam_dynamic_opt.topology_optimization_on) {
        nTO_modules = simparam_dynamic_opt.TO_Module_List.size();
    }

    int myrank = Explicit_Solver_Pointer_->myrank;
    if (simparam->output_options.output_file_format == OUTPUT_FORMAT::vtk) {
        if (myrank == 0) {
            printf("Writing outputs to file at %f \n", time_value);
        }

        double comm_time1 = Explicit_Solver_Pointer_->CPU_Time();
        // Explicit_Solver_Pointer_->write_outputs_new();
        double comm_time2 = Explicit_Solver_Pointer_->CPU_Time();
        Explicit_Solver_Pointer_->output_time += comm_time2 - comm_time1;
    }
    /*
    write_outputs(mesh,
                  Explicit_Solver_Pointer_,
                  node_coords,
                  node_vel,
                  node_mass,
                  elem_den,
                  elem_pres,
                  elem_stress,
                  elem_sspd,
                  elem_sie,
                  elem_vol,
                  elem_mass,
                  elem_mat_id,
                  graphics_times,
                  graphics_id,
                  time_value);
      */

    CArrayKokkos<double> node_extensive_mass(nall_nodes, "node_extensive_mass");

    // extensive energy tallies over the mesh elements local to this MPI rank
    double IE_t0 = 0.0;
    double KE_t0 = 0.0;
    double TE_t0 = 0.0;

    double IE_sum = 0.0;
    double KE_sum = 0.0;

    double IE_loc_sum = 0.0;
    double KE_loc_sum = 0.0;

    // extensive energy tallies over the entire mesh
    double global_IE_t0 = 0.0;
    double global_KE_t0 = 0.0;
    double global_TE_t0 = 0.0;

    // ---- Calculate energy tallies ----
    double IE_tend = 0.0;
    double KE_tend = 0.0;
    double TE_tend = 0.0;

    double global_IE_tend = 0.0;
    double global_KE_tend = 0.0;
    double global_TE_tend = 0.0;

    int nlocal_elem_non_overlapping = Explicit_Solver_Pointer_->nlocal_elem_non_overlapping;

    // extensive IE
    FOR_REDUCE_SUM_CLASS(elem_gid, 0, nlocal_elem_non_overlapping, IE_loc_sum, {
        IE_loc_sum += elem_mass(elem_gid) * elem_sie(rk_level, elem_gid);
    }, IE_sum);
    IE_t0 = IE_sum;

    MPI_Allreduce(&IE_t0, &global_IE_t0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // extensive KE
    FOR_REDUCE_SUM_CLASS(node_gid, 0, nlocal_nodes, KE_loc_sum, {
        double ke = 0;
        for (size_t dim = 0; dim < num_dim; dim++) {
            ke += node_vel(rk_level, node_gid, dim) * node_vel(rk_level, node_gid, dim); // 1/2 at end
        } // end for

        if (num_dim == 2) {
            KE_loc_sum += node_mass(node_gid) * node_coords(rk_level, node_gid, 1) * ke;
        }
        else{
            KE_loc_sum += node_mass(node_gid) * ke;
        }
    }, KE_sum);
    Kokkos::fence();
    KE_t0 = 0.5 * KE_sum;

    MPI_Allreduce(&KE_t0, &global_KE_t0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // extensive TE
    global_TE_t0 = global_IE_t0 + global_KE_t0;
    TE_t0 = global_TE_t0;
    KE_t0 = global_KE_t0;
    IE_t0 = global_IE_t0;

    // save the nodal mass
    FOR_ALL_CLASS(node_gid, 0, nall_nodes, {
        double radius = 1.0;
        if (num_dim == 2) {
            radius = node_coords(rk_level, node_gid, 1);
        }
        node_extensive_mass(node_gid) = node_mass(node_gid) * radius;
    }); // end parallel for

    // a flag to exit the calculation
    size_t stop_calc = 0;

    auto time_1 = std::chrono::high_resolution_clock::now();

    // save initial data
    if (simparam_dynamic_opt.topology_optimization_on || simparam_dynamic_opt.shape_optimization_on) {
        time_data[0] = 0;
        // assign current velocity data to multivector
        // view scope
        {
            vec_array node_velocities_interface = Explicit_Solver_Pointer_->node_velocities_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
            vec_array node_coords_interface     = Explicit_Solver_Pointer_->node_coords_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
            FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
                for (int idim = 0; idim < num_dim; idim++) {
                    node_velocities_interface(node_gid, idim) = node_vel(rk_level, node_gid, idim);
                    node_coords_interface(node_gid, idim)     = node_coords(rk_level, node_gid, idim);
                }
      });
        } // end view scope
        Kokkos::fence();

        // communicate ghosts
        double comm_time1 = Explicit_Solver_Pointer_->CPU_Time();

        // active view scope; triggers host comms from updated data on device
        {
            const_host_vec_array node_velocities_host = Explicit_Solver_Pointer_->node_velocities_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
            const_host_vec_array node_coords_host     = Explicit_Solver_Pointer_->node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
        }
        double comm_time2 = Explicit_Solver_Pointer_->CPU_Time();
        Explicit_Solver_Pointer_->dev2host_time += comm_time2 - comm_time1;

        // communicate ghost velocities
        Explicit_Solver_Pointer_->comm_velocities();
        Explicit_Solver_Pointer_->comm_coordinates();

        double comm_time3 = Explicit_Solver_Pointer_->CPU_Time();

        // view scope
        {
            const_vec_array node_velocities_interface = Explicit_Solver_Pointer_->node_velocities_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            const_vec_array ghost_node_velocities_interface = Explicit_Solver_Pointer_->ghost_node_velocities_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            const_vec_array node_coords_interface = Explicit_Solver_Pointer_->node_coords_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            const_vec_array ghost_node_coords_interface = Explicit_Solver_Pointer_->ghost_node_coords_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);

            vec_array all_node_coords_interface     = Explicit_Solver_Pointer_->all_node_coords_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
            vec_array all_node_velocities_interface = Explicit_Solver_Pointer_->all_node_velocities_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);

            FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
                for (int idim = 0; idim < num_dim; idim++) {
                    all_node_velocities_interface(node_gid, idim) = node_velocities_interface(node_gid, idim);
                    all_node_coords_interface(node_gid, idim)     = node_coords_interface(node_gid, idim);
                }
            }); // end parallel for
            Kokkos::fence();

            FOR_ALL_CLASS(node_gid, nlocal_nodes, nlocal_nodes + nghost_nodes, {
                for (int idim = 0; idim < num_dim; idim++) {
                    all_node_velocities_interface(node_gid, idim) = ghost_node_velocities_interface(node_gid - nlocal_nodes, idim);
                    all_node_coords_interface(node_gid, idim)     = ghost_node_coords_interface(node_gid - nlocal_nodes, idim);
                }
            }); // end parallel for
            Kokkos::fence();
        } // end view scope

        (*forward_solve_velocity_data)[0]->assign(*Explicit_Solver_Pointer_->all_node_velocities_distributed);
        (*forward_solve_coordinate_data)[0]->assign(*Explicit_Solver_Pointer_->all_node_coords_distributed);
    }

    // loop over the max number of time integration cycles
    for (cycle = 0; cycle < cycle_stop; cycle++) {
        // get the step
        if (num_dim == 2) {
            get_timestep2D(mesh,
                           node_coords,
                           node_vel,
                           elem_sspd,
                           elem_vol);
        }
        else{
            get_timestep(mesh,
                         node_coords,
                         node_vel,
                         elem_sspd,
                         elem_vol);
        } // end if 2D

        double global_dt;
        MPI_Allreduce(&dt, &global_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        dt = global_dt;

        // stop calculation if flag
        // if (stop_calc == 1) break;

        if (cycle == 0) {
            if (myrank == 0) {
                printf("cycle = %lu, time = %12.5e, time step = %12.5e \n", cycle, time_value, dt);
            }
        }
        // print time step every 10 cycles
        else if (cycle % 20 == 0) {
            if (myrank == 0) {
                printf("cycle = %lu, time = %12.5e, time step = %12.5e \n", cycle, time_value, dt);
            }
        } // end if

        // ---------------------------------------------------------------------
        //  integrate the solution forward to t(n+1) via Runge Kutta (RK) method
        // ---------------------------------------------------------------------

        // save the values at t_n
        rk_init(node_coords,
                node_vel,
                elem_sie,
                elem_stress,
                rnum_elem,
                nall_nodes);

        // integrate solution forward in time
        for (size_t rk_stage = 0; rk_stage < rk_num_stages; rk_stage++) {
            // ---- RK coefficient ----
            double rk_alpha = 1.0 / ((double)rk_num_stages - (double)rk_stage);

            // current interface has differing velocity arrays; this equates them until we unify memory
            // first comm time interval point
            double comm_time1 = Explicit_Solver_Pointer_->CPU_Time();
            // view scope
            {
                vec_array node_velocities_interface = Explicit_Solver_Pointer_->node_velocities_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
                FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
                    for (int idim = 0; idim < num_dim; idim++) {
                        node_velocities_interface(node_gid, idim) = node_vel(rk_level, node_gid, idim);
                    }
              }); // end parallel for
            } // end view scope
            Kokkos::fence();

            // active view scope
            {
                const_host_vec_array node_velocities_host = Explicit_Solver_Pointer_->node_velocities_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
            }
            double comm_time2 = Explicit_Solver_Pointer_->CPU_Time();
            Explicit_Solver_Pointer_->dev2host_time += comm_time2 - comm_time1;
            // communicate ghost velocities
            Explicit_Solver_Pointer_->comm_velocities();

            double comm_time3 = Explicit_Solver_Pointer_->CPU_Time();
            // this is forcing a copy to the device
            // view scope
            {
                vec_array ghost_node_velocities_interface = Explicit_Solver_Pointer_->ghost_node_velocities_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);

                FOR_ALL_CLASS(node_gid, nlocal_nodes, nall_nodes, {
                    for (int idim = 0; idim < num_dim; idim++) {
                        node_vel(rk_level, node_gid, idim) = ghost_node_velocities_interface(node_gid - nlocal_nodes, idim);
                    }
                }); // end parallel for
            } // end view scope
            Kokkos::fence();

            double comm_time4 = Explicit_Solver_Pointer_->CPU_Time();
            Explicit_Solver_Pointer_->host2dev_time += comm_time4 - comm_time3;
            Explicit_Solver_Pointer_->communication_time += comm_time4 - comm_time1;
            // debug print vector values on a rank
            /*
            if(myrank==0)
             for(int i = 0; i < nall_nodes; i++){
               std::cout << Explicit_Solver_Pointer_->all_node_map->getGlobalElement(i) << " " << node_vel(rk_level,i,0) << " " << node_vel(rk_level,i,1) << " " << node_vel(rk_level,i,2) << std::endl;
             }
            */

            // calculate the new corner masses if 2D
            if (num_dim == 2) {
                // calculate the nodal areal mass
                FOR_ALL_CLASS(node_gid, 0, nall_nodes, {
                    node_mass(node_gid) = 0.0;

                    if (node_coords(rk_level, node_gid, 1) > tiny) {
                        node_mass(node_gid) = node_extensive_mass(node_gid) / node_coords(rk_level, node_gid, 1);
                    }
                    // if(cycle==0&&node_gid==1&&myrank==0)
                    // std::cout << "index " << node_gid << " on rank " << myrank << " node vel " << node_vel(rk_level,node_gid,0) << "  " << node_mass(node_gid) << std::endl << std::flush;
                }); // end parallel for over node_gid
                Kokkos::fence();

                // current interface has differing density arrays; this equates them until we unify memory
                // view scope
                {
                    vec_array node_mass_interface = node_masses_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
                    FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
                        node_mass_interface(node_gid, 0) = node_mass(node_gid);
                  }); // end parallel for
                } // end view scope
                Kokkos::fence();
                // communicate ghost densities
                comm_node_masses();

                // this is forcing a copy to the device
                // view scope
                {
                    vec_array ghost_node_mass_interface = ghost_node_masses_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);

                    FOR_ALL_CLASS(node_gid, nlocal_nodes, nall_nodes, {
                        node_mass(node_gid) = ghost_node_mass_interface(node_gid - nlocal_nodes, 0);
                  }); // end parallel for
                } // end view scope
                Kokkos::fence();

                // -----------------------------------------------
                // Calcualte the areal mass for nodes on the axis
                // -----------------------------------------------
                // The node order of the 2D element is
                //
                //   J
                //   |
                // 3---2
                // |   |  -- I
                // 0---1
                /*
                FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {

                    // loop over the corners of the element and calculate the mass
                    for (size_t node_lid=0; node_lid<4; node_lid++){

                        size_t node_gid = nodes_in_elem(elem_gid, node_lid);
                        size_t node_minus_gid;
                        size_t node_plus_gid;


                        if (node_coords(rk_level,node_gid,1) < tiny){
                            // node is on the axis

                            // minus node
                            if (node_lid==0){
                                node_minus_gid = nodes_in_elem(elem_gid, 3);
                            } else {
                                node_minus_gid = nodes_in_elem(elem_gid, node_lid-1);
                            }

                            // plus node
                            if (node_lid==3){
                                node_plus_gid = nodes_in_elem(elem_gid, 0);
                            } else {
                                node_plus_gid = nodes_in_elem(elem_gid, node_lid+1);
                            }

                            node_mass(node_gid) = fmax(node_mass(node_plus_gid), node_mass(node_minus_gid))/2.0;

                        } // end if

                    } // end for over corners

                }); // end parallel for over elem_gid
                Kokkos::fence();
                 */

                FOR_ALL_CLASS(node_bdy_gid, 0, num_bdy_nodes, {
                    // FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
                    size_t node_gid = bdy_nodes(node_bdy_gid);

                    if (node_coords(rk_level, node_gid, 1) < tiny) {
                        // node is on the axis

                        for (size_t node_lid = 0; node_lid < num_nodes_in_node(node_gid); node_lid++) {
                            size_t node_neighbor_gid = nodes_in_node(node_gid, node_lid);

                            // if the node is off the axis, use it's areal mass on the boundary
                            if (node_coords(rk_level, node_neighbor_gid, 1) > tiny) {
                                node_mass(node_gid) = fmax(node_mass(node_gid), node_mass(node_neighbor_gid) / 2.0);
                            }
                        } // end for over neighboring nodes
                    } // end if
                }); // end parallel for over elem_gid
            } // end of if 2D-RZ
        } // end of RK loop

        // increment the time
        time_value += dt;
        simparam->time_value = time_value;

        if (simparam_dynamic_opt.topology_optimization_on || simparam_dynamic_opt.shape_optimization_on) {
            if (cycle >= max_time_steps) {
                max_time_steps = cycle + 1;
            }

            if (max_time_steps + 1 > forward_solve_velocity_data->size()) {
                old_max_forward_buffer = forward_solve_velocity_data->size();
                time_data.resize(max_time_steps + BUFFER_GROW + 1);
                forward_solve_velocity_data->resize(max_time_steps + BUFFER_GROW + 1);
                forward_solve_coordinate_data->resize(max_time_steps + BUFFER_GROW + 1);
                adjoint_vector_data->resize(max_time_steps + BUFFER_GROW + 1);
                phi_adjoint_vector_data->resize(max_time_steps + BUFFER_GROW + 1);
                // assign a multivector of corresponding size to each new timestep in the buffer
                for (int istep = old_max_forward_buffer; istep < max_time_steps + BUFFER_GROW + 1; istep++) {
                    (*forward_solve_velocity_data)[istep]   = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
                    (*forward_solve_coordinate_data)[istep] = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
                    (*adjoint_vector_data)[istep]     = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
                    (*phi_adjoint_vector_data)[istep] = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
                }
            }

            time_data[cycle + 1] = dt + time_data[cycle];

            // assign current velocity data to multivector
            // view scope
            {
                vec_array node_velocities_interface = Explicit_Solver_Pointer_->node_velocities_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
                vec_array node_coords_interface     = Explicit_Solver_Pointer_->node_coords_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
                FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
                    for (int idim = 0; idim < num_dim; idim++) {
                        node_velocities_interface(node_gid, idim) = node_vel(rk_level, node_gid, idim);
                        node_coords_interface(node_gid, idim)     = node_coords(rk_level, node_gid, idim);
                    }
                });
            } // end view scope
            Kokkos::fence();

            // communicate ghosts
            double comm_time1 = Explicit_Solver_Pointer_->CPU_Time();

            // active view scope; triggers host comms from updated data on device
            {
                const_host_vec_array node_velocities_host = Explicit_Solver_Pointer_->node_velocities_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
                const_host_vec_array node_coords_host     = Explicit_Solver_Pointer_->node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
            }
            double comm_time2 = Explicit_Solver_Pointer_->CPU_Time();
            Explicit_Solver_Pointer_->dev2host_time += comm_time2 - comm_time1;

            // communicate ghost velocities
            Explicit_Solver_Pointer_->comm_velocities();
            Explicit_Solver_Pointer_->comm_coordinates();

            double comm_time3 = Explicit_Solver_Pointer_->CPU_Time();

            // view scope
            {
                const_vec_array node_velocities_interface = Explicit_Solver_Pointer_->node_velocities_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                const_vec_array ghost_node_velocities_interface = Explicit_Solver_Pointer_->ghost_node_velocities_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                vec_array all_node_velocities_interface     = Explicit_Solver_Pointer_->all_node_velocities_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
                const_vec_array node_coords_interface       = Explicit_Solver_Pointer_->node_coords_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                const_vec_array ghost_node_coords_interface = Explicit_Solver_Pointer_->ghost_node_coords_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                vec_array all_node_coords_interface = Explicit_Solver_Pointer_->all_node_coords_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
                FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
                    for (int idim = 0; idim < num_dim; idim++) {
                        all_node_velocities_interface(node_gid, idim) = node_velocities_interface(node_gid, idim);
                        all_node_coords_interface(node_gid, idim)     = node_coords_interface(node_gid, idim);
                    }
                }); // end parallel for
                Kokkos::fence();

                FOR_ALL_CLASS(node_gid, nlocal_nodes, nlocal_nodes + nghost_nodes, {
                    for (int idim = 0; idim < num_dim; idim++) {
                        all_node_velocities_interface(node_gid, idim) = ghost_node_velocities_interface(node_gid - nlocal_nodes, idim);
                        all_node_coords_interface(node_gid, idim)     = ghost_node_coords_interface(node_gid - nlocal_nodes, idim);
                    }
                }); // end parallel for
                Kokkos::fence();
            } // end view scope

            double comm_time4 = Explicit_Solver_Pointer_->CPU_Time();
            Explicit_Solver_Pointer_->host2dev_time += comm_time4 - comm_time3;
            Explicit_Solver_Pointer_->communication_time += comm_time4 - comm_time1;

            (*forward_solve_velocity_data)[cycle + 1]->assign(*Explicit_Solver_Pointer_->all_node_velocities_distributed);
            (*forward_solve_coordinate_data)[cycle + 1]->assign(*Explicit_Solver_Pointer_->all_node_coords_distributed);

            // kinetic energy accumulation
            if (kinetic_energy_objective) {
                const_vec_array node_velocities_interface = (*forward_solve_velocity_data)[cycle + 1]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                const_vec_array previous_node_velocities_interface = (*forward_solve_velocity_data)[cycle]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
                KE_loc_sum = 0.0;
                KE_sum     = 0.0;
                // extensive KE
                FOR_REDUCE_SUM_CLASS(node_gid, 0, nlocal_nodes, KE_loc_sum, {
                    double ke = 0;
                    for (size_t dim = 0; dim < num_dim; dim++) {
                        // midpoint integration approximation
                        ke += (node_velocities_interface(node_gid, dim) + node_velocities_interface(node_gid, dim)) * (node_velocities_interface(node_gid, dim) + node_velocities_interface(node_gid,
                        dim)) / 4;                                                                                                                                                       // 1/2 at end
                    } // end for

                    if (num_dim == 2) {
                        KE_loc_sum += node_mass(node_gid) * node_coords(rk_level, node_gid, 1) * ke;
                    }
                    else{
                        KE_loc_sum += node_mass(node_gid) * ke;
                    }
                }, KE_sum);
                Kokkos::fence();
                KE_sum = 0.5 * KE_sum;
                objective_accumulation += KE_sum * dt;
            }
        } // end topology optimization if check

        size_t write = 0;
        if ((cycle + 1) % graphics_cyc_ival == 0 && cycle > 0) {
            write = 1;
        }
        else if (cycle == cycle_stop) {
            write = 1;
        }
        else if (time_value >= time_final) {
            write = 1;
        }
        else if (time_value >= graphics_time) {
            write = 1;
        }

        // write outputs
        if (write == 1) {
            // interface nodal coordinate data (note: this is not needed if using write_outputs_new())
            // view scope
            {
                vec_array node_coords_interface = Explicit_Solver_Pointer_->node_coords_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
                FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
                    for (int idim = 0; idim < num_dim; idim++) {
                        node_coords_interface(node_gid, idim) = node_coords(rk_level, node_gid, idim);
                    }
              }); // end parallel for
            } // end view scope
            if (simparam->output_options.output_file_format == OUTPUT_FORMAT::vtk) {
                if (myrank == 0) {
                    printf("Writing outputs to file at %f \n", graphics_time);
                }

                double comm_time1 = Explicit_Solver_Pointer_->CPU_Time();
                // Explicit_Solver_Pointer_->write_outputs_new();

                double comm_time2 = Explicit_Solver_Pointer_->CPU_Time();
                Explicit_Solver_Pointer_->output_time += comm_time2 - comm_time1;
            }
            // Explicit_Solver_Pointer_->parallel_vtk_writer();
            // Explicit_Solver_Pointer_->parallel_vtk_writer_new();
            // Explicit_Solver_Pointer_->parallel_tecplot_writer();
            /*
          write_outputs(mesh,
                        Explicit_Solver_Pointer_,
                        node_coords,
                        node_vel,
                        node_mass,
                        elem_den,
                        elem_pres,
                        elem_stress,
                        elem_sspd,
                        elem_sie,
                        elem_vol,
                        elem_mass,
                        elem_mat_id,
                        graphics_times,
                        graphics_id,
                        time_value);
          */
            graphics_time = time_value + graphics_dt_ival;
        } // end if

        // end of calculation
        if (time_value >= time_final) {
            break;
        }
    } // end for cycle loop

    last_time_step = cycle;

    // simple setup to just calculate KE minimize objective for now
    if (simparam_dynamic_opt.topology_optimization_on) {
        KineticEnergyMinimize_TopOpt& kinetic_energy_minimize_function = dynamic_cast<KineticEnergyMinimize_TopOpt&>(*obj_pointer);

        // collect local objective values
        MPI_Allreduce(&objective_accumulation, &global_objective_accumulation, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        kinetic_energy_minimize_function.objective_accumulation = global_objective_accumulation;

        if (myrank == 0) {
            std::cout << "CURRENT TIME INTEGRAL OF KINETIC ENERGY " << global_objective_accumulation << std::endl;
        }
    }

    auto time_2 = std::chrono::high_resolution_clock::now();
    auto time_difference = time_2 - time_1;
    // double calc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(diff).count();
    double calc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(time_difference).count();
    if (myrank == 0) {
        printf("\nCalculation time in seconds: %f \n", calc_time * 1e-09);
    }

    IE_loc_sum = 0.0;
    KE_loc_sum = 0.0;
    IE_sum     = 0.0;
    KE_sum     = 0.0;

    // extensive IE
    FOR_REDUCE_SUM_CLASS(elem_gid, 0, nlocal_elem_non_overlapping, IE_loc_sum, {
        IE_loc_sum += elem_mass(elem_gid) * elem_sie(rk_level, elem_gid);
    }, IE_sum);
    IE_tend = IE_sum;

    // reduce over MPI ranks
    MPI_Allreduce(&IE_tend, &global_IE_tend, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // extensive KE
    FOR_REDUCE_SUM_CLASS(node_gid, 0, nlocal_nodes, KE_loc_sum, {
        double ke = 0;
        for (size_t dim = 0; dim < num_dim; dim++) {
            ke += node_vel(rk_level, node_gid, dim) * node_vel(rk_level, node_gid, dim); // 1/2 at end
        } // end for

        if (num_dim == 2) {
            KE_loc_sum += node_mass(node_gid) * node_coords(rk_level, node_gid, 1) * ke;
        }
        else{
            KE_loc_sum += node_mass(node_gid) * ke;
        }
    }, KE_sum);
    Kokkos::fence();
    KE_tend = 0.5 * KE_sum;

    // reduce over MPI ranks
    MPI_Allreduce(&KE_tend, &global_KE_tend, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // extensive TE
    TE_tend = IE_tend + KE_tend;
    KE_tend = global_KE_tend;
    IE_tend = global_IE_tend;

    // extensive TE
    TE_tend = IE_tend + KE_tend;

    // reduce over MPI ranks

    if (myrank == 0) {
        printf("Time=0:   KE = %20.15f, IE = %20.15f, TE = %20.15f \n", KE_t0, IE_t0, TE_t0);
    }
    if (myrank == 0) {
        printf("Time=End: KE = %20.15f, IE = %20.15f, TE = %20.15f \n", KE_tend, IE_tend, TE_tend);
    }
    if (myrank == 0) {
        printf("total energy conservation error %= %e \n\n", 100 * (TE_tend - TE_t0) / TE_t0);
    }

    return;
} // end of Eulerian solve
