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

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include "Tpetra_Import.hpp"
#include "Tpetra_Import_Util2.hpp"

#include "elements.h"
#include "swage.h"
#include "matar.h"
#include "utilities.h"
#include "node_combination.h"
#include "FEA_Module_Dynamic_Elasticity.h"
#include "Explicit_Solver.h"
#include "Simulation_Parameters/FEA_Module/Dynamic_Elasticity_Parameters.h"
#include "Simulation_Parameters/Simulation_Parameters_Explicit.h"

// optimization
#include "ROL_Solver.hpp"
#include "Kinetic_Energy_Minimize.h"

#define MAX_ELEM_NODES 8
#define STRAIN_EPSILON 0.000000001
#define DENSITY_EPSILON 0.0001
#define BC_EPSILON 1.0e-6
#define BUFFER_GROW 100

using namespace utils;

FEA_Module_Dynamic_Elasticity::FEA_Module_Dynamic_Elasticity(
    Dynamic_Elasticity_Parameters& params, Solver* Solver_Pointer,
    std::shared_ptr<mesh_t> mesh_in, const int my_fea_module_index)
    : FEA_Module(Solver_Pointer)
{
    // assign interfacing index
    my_fea_module_index_ = my_fea_module_index;
    Module_Type = FEA_MODULE_TYPE::Dynamic_Elasticity;

    // recast solver pointer for non-base class access
    Explicit_Solver_Pointer_ = dynamic_cast<Explicit_Solver*>(Solver_Pointer);
    module_params = &params;
    simparam = &(Explicit_Solver_Pointer_->simparam);

    mesh = mesh_in;

    // boundary condition data
    max_boundary_sets = 0;
    Local_Index_Boundary_Patches = Explicit_Solver_Pointer_->Local_Index_Boundary_Patches;

    // set Tpetra vector pointers
    initial_node_velocities_distributed = Explicit_Solver_Pointer_->initial_node_velocities_distributed;
    node_velocities_distributed     = Explicit_Solver_Pointer_->node_velocities_distributed;
    all_node_velocities_distributed = Explicit_Solver_Pointer_->all_node_velocities_distributed;

    if (simparam->topology_optimization_on || simparam->shape_optimization_on) {
        all_cached_node_velocities_distributed = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
        force_gradient_velocity    = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
        force_gradient_position    = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
        force_gradient_design      = Teuchos::rcp(new MV(all_node_map, 1));
        corner_value_storage       = Solver_Pointer->corner_value_storage;
        corner_vector_storage      = Solver_Pointer->corner_vector_storage;
        relative_element_densities = DCArrayKokkos<double>(rnum_elem, "relative_element_densities");
    }

    if (simparam->topology_optimization_on || simparam->shape_optimization_on || simparam->num_dims == 2) {
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
    Dynamic_Options dynamic_options = simparam->dynamic_options;

    dt = dynamic_options.dt;
    time_value  = dynamic_options.time_value;
    time_final  = dynamic_options.time_final;
    dt_max      = dynamic_options.dt_max;
    dt_min      = dynamic_options.dt_min;
    dt_cfl      = dynamic_options.dt_cfl;
    rk_num_bins = simparam->dynamic_options.rk_num_bins;

    graphics_time     = simparam->output_options.graphics_time;
    graphics_dt_ival  = simparam->output_options.graphics_dt_ival;
    graphics_cyc_ival = simparam->output_options.graphics_cyc_ival;
    graphics_times    = simparam->output_options.graphics_times;
    graphics_id = simparam->output_options.graphics_id;

    cycle_stop    = dynamic_options.cycle_stop;
    rk_num_stages = dynamic_options.rk_num_stages;

    fuzz  = dynamic_options.fuzz;
    tiny  = dynamic_options.tiny;
    small = dynamic_options.small;

    if (simparam->topology_optimization_on) {
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

FEA_Module_Dynamic_Elasticity::~FEA_Module_Dynamic_Elasticity()
{
    // delete simparam;
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn read_conditions_ansys_dat
///
/// \brief Read ANSYS dat format mesh file
///
/// \param Initial conditions header
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::read_conditions_ansys_dat(std::ifstream* in, std::streampos before_condition_header)
{
    Input_Options input_options = simparam->input_options.value();

    char ch;
    std::string skip_line, read_line, substring, token;
    std::stringstream line_parse, line_parse2;

    int num_dim = simparam->num_dims;
    int buffer_lines = 1000;
    int max_word     = 30;
    int local_node_index, current_column_index;
    int p_order = input_options.p_order;
    int buffer_loop, buffer_iteration, buffer_iterations, scan_loop, nodes_per_element, words_per_line;

    size_t strain_count;
    size_t read_index_start;
    size_t node_rid;
    size_t elem_gid;

    real_t unit_scaling = input_options.unit_scaling;

    CArrayKokkos<char, array_layout, HostSpace, memory_traits> read_buffer;
    CArrayKokkos<long long int, array_layout, HostSpace, memory_traits> read_buffer_indices;

    LO     local_dof_id;
    GO     node_gid;
    real_t dof_value;
    host_vec_array node_densities;
} // end read_conditions_ansys_dat

/////////////////////////////////////////////////////////////////////////////
///
/// \fn elastic_interface_setup
///
/// \brief Interfaces read in data with the elastic solver data; currently a hack to streamline
///
/// \param Nodal state data
/// \param Element state data
/// \param Corner state data
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::elastic_interface_setup(node_t& node,
    elem_t&   elem,
    corner_t& corner)
{
    const size_t num_dim     = simparam->num_dims;
    const size_t rk_num_bins = simparam->dynamic_options.rk_num_bins;

    num_nodes_in_elem = 1;
    for (int dim = 0; dim < num_dim; dim++) {
        num_nodes_in_elem *= 2;
    }

    // --- Read in the nodes in the mesh ---

    nall_nodes = Explicit_Solver_Pointer_->nall_nodes;
    int myrank = Explicit_Solver_Pointer_->myrank;
    int nranks = Explicit_Solver_Pointer_->nranks;
    // printf("Num nodes assigned to MPI rank %lu is %lu\n" , myrank, nall_nodes);

    // intialize node variables
    mesh->initialize_nodes(nall_nodes);
    mesh->initialize_local_nodes(Explicit_Solver_Pointer_->nlocal_nodes);
    node.initialize(rk_num_bins, nall_nodes, num_dim);
    // std::cout << "Bin counts " << rk_num_bins << " Node counts " << nall_nodes << " Num dim " << num_dim << std::endl;

    // view scope
    {
        host_vec_array interface_node_coords = Explicit_Solver_Pointer_->all_node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
        // save node data to node.coords
        // std::cout << "NODE DATA ON RANK " << myrank << std::endl;
        if (num_dim == 2) {
            for (int inode = 0; inode < nall_nodes; inode++) {
                // std::cout << "Node index " << inode+1 << " ";
                node.coords.host(0, inode, 0) = interface_node_coords(inode, 0);
                // std::cout << host_node_coords_state(0,inode,0)+1<< " ";
                node.coords.host(0, inode, 1) = interface_node_coords(inode, 1);
                // std::cout << host_node_coords_state(0,inode,1)+1<< " ";
            }
        }
        else if (num_dim == 3) {
            for (int inode = 0; inode < nall_nodes; inode++) {
                // std::cout << "Node index " << inode+1 << " ";
                node.coords.host(0, inode, 0) = interface_node_coords(inode, 0);
                // std::cout << host_node_coords_state(0,inode,0)+1<< " ";
                node.coords.host(0, inode, 1) = interface_node_coords(inode, 1);
                // std::cout << host_node_coords_state(0,inode,1)+1<< " ";

                node.coords.host(0, inode, 2) = interface_node_coords(inode, 2);
                // std::cout << host_node_coords_state(0,inode,2)+1<< std::endl;
            }
        }
    } // end view scope
      // --- read in the elements in the mesh ---

    rnum_elem = Explicit_Solver_Pointer_->rnum_elem;
    // printf("Num elems assigned to MPI rank %lu is %lu\n" , myrank, rnum_elem);

    // intialize elem variables
    mesh->initialize_elems(rnum_elem, num_dim);
    elem.initialize(rk_num_bins, nall_nodes, 3); // always 3D here, even for 2D
    nodes_in_elem = mesh->nodes_in_elem;
    // save data to nodes_in_elem.host
    // CArrayKokkos<size_t, DefaultLayout, HostSpace> host_mesh_nodes_in_elem(rnum_elem, num_nodes_in_elem);
    // view scope
    {
        host_elem_conn_array interface_nodes_in_elem = Explicit_Solver_Pointer_->global_nodes_in_elem_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
        // save node data to node.coords
        // std::cout << "ELEMENT CONNECTIVITY ON RANK " << myrank << std::endl;
        for (int ielem = 0; ielem < rnum_elem; ielem++) {
            // std::cout << "Element index " << ielem+1 << " ";
            for (int inode = 0; inode < num_nodes_in_elem; inode++) {
                nodes_in_elem.host(ielem, inode) = Explicit_Solver_Pointer_->all_node_map->getLocalElement(interface_nodes_in_elem(ielem, inode));
                // debug print
                // std::cout << nodes_in_elem.get_kokkos_dual_view().h_view(ielem*num_nodes_in_elem + inode)+1<< " ";
            }
            // std::cout << std::endl;
        }
    }
    // update device side
    nodes_in_elem.update_device();

    // debug print

    // CArrayKokkos<size_t> device_mesh_nodes_in_elem(rnum_elem, num_nodes_in_elem);
    // device_mesh_nodes_in_elem.get_kokkos_view() = nodes_in_elem.get_kokkos_dual_view().d_view;
    // host_mesh_nodes_in_elem.get_kokkos_view() = nodes_in_elem.get_kokkos_dual_view().view_host();
    /*
    if(myrank==1){
    std::cout << "ELEMENT CONNECTIVITY ON RANK 1 in LOCAL INDICES" << myrank << std::endl;
    for(int ielem = 0; ielem < rnum_elem; ielem++){
        std::cout << "Element index " << ielem+1 << " ";
        for(int inode = 0; inode < num_nodes_in_elem; inode++){
            //debug print
            //device_mesh_nodes_in_elem(ielem,inode) = Explicit_Solver_Pointer_->all_node_map->getLocalElement(interface_nodes_in_elem(ielem,inode));
            std::cout << nodes_in_elem(ielem, inode)+1<< " ";
        }
        std::cout << std::endl;
    }
    }
    */
    /*
    std::cout.flush();
    if(myrank==1){
    std::cout << "ELEMENT CONNECTIVITY ON RANK 1 in GLOBAL INDICES" << myrank << std::endl;
    std::cout << "local node index of global index 275 on rank 1 " << Explicit_Solver_Pointer_->all_node_map->getLocalElement(275) << std::endl;
    for(int ielem = 0; ielem < rnum_elem; ielem++){
        std::cout << ielem << " ";
        for(int inode = 0; inode < num_nodes_in_elem; inode++){
            //debug print
            //device_mesh_nodes_in_elem(ielem,inode) = Explicit_Solver_Pointer_->all_node_map->getLocalElement(interface_nodes_in_elem(ielem,inode));
            std::cout << Explicit_Solver_Pointer_->all_node_map->getGlobalElement(nodes_in_elem(ielem, inode))<< " ";
        }
        std::cout << std::endl;
    }
    }
    std::cout.flush();
    */
    /*
    size_t nall_nodes = Explicit_Solver_Pointer_->nall_nodes;
    node.all_coords = DCArrayKokkos <double> (rk_num_bins, nall_nodes, num_dim);
    node.all_vel    = DCArrayKokkos <double> (rk_num_bins, nall_nodes, num_dim);
    node.all_mass   = DCArrayKokkos <double> (nall_nodes);

    //save all data (nlocal +nghost)
    CArrayKokkos<double, DefaultLayout, HostSpace> host_all_node_coords_state(rk_num_bins, nall_nodes, num_dim);
    host_vec_array interface_all_node_coords = Explicit_Solver_Pointer_->all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);
    host_all_node_coords_state.get_kokkos_view() = node.all_coords.get_kokkos_dual_view().view_host();
    //host_node_coords_state = CArrayKokkos<double, DefaultLayout, HostSpace>(rk_num_bins, nall_nodes, num_dim);
    //host_all_node_coords_state.get_kokkos_view() = Kokkos::View<double*,DefaultLayout, HostSpace>("debug", rk_num_bins*nall_nodes*num_dim);
    //save node data to node.coords

    //std::cout << "ALL NODE DATA ON RANK " << myrank << std::endl;
    for(int inode = 0; inode < nall_nodes; inode++){
        //std::cout << "Node index " << inode+1 << " ";
        node.all_coords.host(0,inode,0) = interface_all_node_coords(inode,0);
        //std::cout << host_all_node_coords_state(0,inode,0)+1<< " ";
        node.all_coords.host(0,inode,1) = interface_all_node_coords(inode,1);
        //std::cout << host_all_node_coords_state(0,inode,1)+1<< " ";
        node.all_coords.host(0,inode,2) = interface_all_node_coords(inode,2);
        //std::cout << host_all_node_coords_state(0,inode,2)+1<< std::endl;
    }
    */

    // save the node coords to the current RK value
    for (size_t node_gid = 0; node_gid < nall_nodes; node_gid++) {
        for (int rk = 1; rk < rk_num_bins; rk++) {
            for (int dim = 0; dim < num_dim; dim++) {
                node.coords.host(rk, node_gid, dim) = node.coords.host(0, node_gid, dim);
            } // end for dim
        } // end for rk
    } // end parallel for

    /*
    // save the node coords to the current RK value
    for (size_t node_gid=0; node_gid<nall_nodes; node_gid++){

        for(int rk=1; rk<rk_num_bins; rk++){
            for (int dim = 0; dim < num_dim; dim++){
                node.all_coords.host(rk, node_gid, dim) = node.all_coords.host(0, node_gid, dim);
            } // end for dim
        } // end for rk

    } // end parallel for
    */

    node.coords.update_device();
    // node.all_coords.update_device();

    // intialize corner variables
    int num_corners = rnum_elem * num_nodes_in_elem;
    mesh->initialize_corners(num_corners);
    corner.initialize(num_corners, num_dim);

    /*
    for(int inode = 0; inode < nall_nodes; inode++){
        std::cout << "Node index " << inode+1 << " ";
        for(int rk=0; rk<rk_num_bins; rk++){
          std::cout << "rk index " << rk+1 << " ";
          std::cout << node.coords(rk,inode,0)+1<< " ";
          std::cout << node.coords(rk,inode,1)+1<< " ";
          std::cout << node.coords(rk,inode,2)+1<< std::endl;
        }
    }
    */
    // Close mesh input file
    // fclose(in);

    return;
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn init_boundaries
///
/// \brief Initialize sets of element boundary surfaces and arrays for input conditions
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::init_boundaries()
{
    max_boundary_sets = module_params->boundary_conditions.size();
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
/// \brief Initialize storage for element boundary surfaces corresponding to user BCs
///
/// \param Number of boundary sets
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::init_boundary_sets(int num_sets)
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
void FEA_Module_Dynamic_Elasticity::grow_boundary_sets(int num_sets)
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
void FEA_Module_Dynamic_Elasticity::generate_bcs()
{
} // end generate_bcs

/////////////////////////////////////////////////////////////////////////////
///
/// \fn Displacement_Boundary_Conditions
///
/// \brief Loop through applied boundary conditions and tag node ids to remove
///        necessary rows and columns from the assembled linear system
///
/// Unused in this context
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::Displacement_Boundary_Conditions()
{
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn init_output
///
/// \brief Initialize output data structures
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::init_output()
{
    // check user parameters for output
    bool output_velocity_flag = simparam->output(FIELD::velocity);
    bool output_strain_flag   = simparam->output(FIELD::strain);
    bool output_stress_flag   = simparam->output(FIELD::stress);

    int num_dim = simparam->num_dims;
    int Brows;

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

/////////////////////////////////////////////////////////////////////////////
///
/// \fn sort_output
///
/// \brief Prompts sorting for elastic response output data. For now, nodal strains.
///
/// Unused in this context
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::sort_output(Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> sorted_map)
{
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn write_data
///
/// \brief Populate requests this module makes for output data
///
/// \param Map of scalar point data
/// \param Map of vector point data
/// \param Map of scalar cell data
/// \param Map of int cell data
/// \param Map of cell field data
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::write_data(std::map<std::string, const double*>& point_data_scalars_double,
    std::map<std::string, const double*>& point_data_vectors_double,
    std::map<std::string, const double*>& cell_data_scalars_double,
    std::map<std::string, const int*>&    cell_data_scalars_int,
    std::map<std::string, std::pair<const double*, size_t>>& cell_data_fields_double)
{
    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;

    for (const FIELD& field_name : simparam->output_options.output_fields) {
        switch (field_name) {
        case FIELD::velocity:
            // node "velocity"
            node_vel.update_host();
            point_data_vectors_double["velocity"] = &node_vel.host(rk_level, 0, 0);
            break;

        case FIELD::element_density:
            // element "density"
            elem_den.update_host();
            cell_data_scalars_double["element_density"] = elem_den.host_pointer();
            break;

        case FIELD::pressure:
            // element "pressure"
            elem_pres.update_host();
            cell_data_scalars_double["pressure"] = elem_pres.host_pointer();
            break;

        case FIELD::volume:
            // element "volume"
            elem_vol.update_host();
            cell_data_scalars_double["volume"] = elem_vol.host_pointer();
            break;

        case FIELD::mass:
            // element "mass"
            elem_mass.update_host();
            cell_data_scalars_double["mass"] = elem_mass.host_pointer();
            break;

        case FIELD::material_id:
            // element "material_id"
            elem_mat_id.update_host();
            cell_data_scalars_int["material_id"] = reinterpret_cast<int*>(elem_mat_id.host_pointer());
            break;

        case FIELD::user_vars:
            // element "user_vars"
            elem_user_output_vars.update_host();
            cell_data_fields_double["user_vars"] = std::make_pair(elem_user_output_vars.host_pointer(),
                                                                   elem_user_output_vars.dims(1));
        case FIELD::stress:
            // element "stress"
            elem_stress.update_host();
            cell_data_fields_double["stress"] = std::make_pair(&elem_stress.host(rk_level, 0, 0, 0), 9);
            break;

        default:
            break;
        } // end switch
    } // end if

    // element "mat_id" //uncomment if needed (works fine)
    // sgh_module->elem_mat_id.update_host();
    // cell_data_scalars_int["mat_id"] = reinterpret_cast<int*>(&sgh_module->elem_mat_id.host(0));

    // element "user_output_vars" //uncomment if needed (works fine)
    // sgh_module->elem_user_output_vars.update_host();
    // cell_data_fields_double["user_output_vars"] = std::make_pair(&sgh_module->elem_user_output_vars.host_pointer(),
    //                                                             sgh_module->elem_user_output_vars.dims(1));

    // element "stress" //uncomment if needed (works fine)
    // sgh_module->elem_stress.update_host();
    // cell_data_fields_double["stress"] = std::make_pair(&sgh_module->elem_stress.host(rk_level,0,0,0), 9);
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn sort_element_output
///
/// \brief  Prompts sorting for elastic response output data. For now, element densities
///
/// \param Sorted element data map
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::sort_element_output(Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> sorted_map)
{
    // interface element density data
    {
        host_vec_array Element_Densities = Global_Element_Densities->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
        elem_den.update_host();
        for (int ielem = 0; ielem < rnum_elem; ielem++) {
            Element_Densities(ielem, 0) = elem_den.host(ielem);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn collect_output
///
/// \brief Prompts computation of elastic response output data.
///
/// Unused in this context
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::collect_output(Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> global_reduce_map)
{
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn compute_output
///
/// \brief Prompts computation of elastic response output data.
///
/// Unused in this context
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::compute_output()
{
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn comm_node_masses
///
/// \brief Communicate updated nodal velocities to ghost nodes
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::comm_node_masses()
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
/// \fn comm_variables
///
/// \brief Communicate ghosts using the current optimization design data
///
/// \param Design variables
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::comm_variables(Teuchos::RCP<const MV> zp)
{
    if (simparam->topology_optimization_on) {
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
    else if (simparam->shape_optimization_on) {
        // clause to communicate boundary node data if the boundary nodes are ghosts on this rank
    }
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn node_density_constraints
///
/// \brief Enforce constraints on nodes due to BCS
///
/// \param Lower bound of nodal densities
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::node_density_constraints(host_vec_array& node_densities_lower_bound)
{
    const_vec_array all_initial_node_coords = all_initial_node_coords_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);

    const size_t num_dim = mesh->num_dims;
    const size_t num_lcs = module_params->loading_conditions.size();

    const DCArrayKokkos<mat_fill_t> mat_fill = simparam->mat_fill;
    const DCArrayKokkos<loading_t>  loading  = module_params->loading;

    // walk over the nodes to update the velocity
    FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
        double current_node_coords[3];
        double radius;
        for (size_t dim = 0; dim < num_dim; dim++) {
            current_node_coords[dim] = all_initial_node_coords(node_gid, dim);
        } // end for dim
        radius = sqrt(current_node_coords[0] * current_node_coords[0] + current_node_coords[1] * current_node_coords[1] + current_node_coords[2] * current_node_coords[2]);
        for (size_t ilc = 0; ilc < num_lcs; ilc++) {
            // debug check
            // std::cout << "LOADING CONDITION VOLUME TYPE: " << to_string(loading(ilc).volume) << std::endl;

            bool fill_this = loading(ilc).volume.contains(current_node_coords);
            if (fill_this) {
                node_densities_lower_bound(node_gid, 0) = 1;
            }
        }
    }); // end for parallel for over nodes
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn setup
///
/// \brief Setup elastic solver data
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::setup()
{
    Dynamic_Options dynamic_options = simparam->dynamic_options;

    const size_t rk_level      = dynamic_options.rk_num_bins - 1;
    const size_t num_fills     = simparam->regions.size();
    const size_t rk_num_bins   = dynamic_options.rk_num_bins;
    const size_t num_bcs       = module_params->boundary_conditions.size();
    const size_t num_materials = simparam->materials.size();

    const int num_dim = simparam->num_dims;

    // ---------------------------------------------------------------------
    //    obtain mesh data
    // ---------------------------------------------------------------------
    elastic_interface_setup(node_interface, elem_interface, corner_interface);
    mesh->build_corner_connectivity();
    mesh->build_elem_elem_connectivity();
    mesh->num_bdy_patches = nboundary_patches;
    if (num_dim == 2) {
        mesh->build_patch_connectivity();
        mesh->build_node_node_connectivity();
    }

    // ---------------------------------------------------------------------
    //    allocate memory
    // ---------------------------------------------------------------------

    // shorthand names
    const size_t num_nodes   = mesh->num_nodes;
    const size_t num_elems   = mesh->num_elems;
    const size_t num_corners = mesh->num_corners;

    // --- make dual views of data on CPU and GPU ---
    //  Notes:
    //     Instead of using a struct of dual types like the mesh type,
    //     individual dual views will be made for all the state
    //     variables.  The motivation is to reduce memory movement
    //     when passing state into a function.  Passing a struct by
    //     reference will copy the meta data and pointers for the
    //     variables held inside the struct.  Since all the mesh
    //     variables are typically used by most functions, a single
    //     mesh struct or passing the arrays will be roughly equivalent
    //     for memory movement.

    // create Dual Views of the individual node struct variables
    node_coords = DViewCArrayKokkos<double>(node_interface.coords.get_kokkos_dual_view().view_host().data(), rk_num_bins, num_nodes, num_dim);
    node_vel    = DViewCArrayKokkos<double>(node_interface.vel.get_kokkos_dual_view().view_host().data(), rk_num_bins, num_nodes, num_dim);
    node_mass   = DViewCArrayKokkos<double>(node_interface.mass.get_kokkos_dual_view().view_host().data(), num_nodes);

    // create Dual Views of the individual elem struct variables
    elem_den    = DViewCArrayKokkos<double>(&elem_interface.den(0), num_elems);
    elem_pres   = DViewCArrayKokkos<double>(&elem_interface.pres(0), num_elems);
    elem_stress = DViewCArrayKokkos<double>(&elem_interface.stress(0, 0, 0, 0), rk_num_bins, num_elems, 3, 3); // always 3D even in 2D-RZ
    elem_sspd   = DViewCArrayKokkos<double>(&elem_interface.sspd(0), num_elems);
    elem_sie    = DViewCArrayKokkos<double>(&elem_interface.sie(0, 0), rk_num_bins, num_elems);
    elem_vol    = DViewCArrayKokkos<double>(&elem_interface.vol(0), num_elems);
    elem_div    = DViewCArrayKokkos<double>(&elem_interface.div(0), num_elems);
    elem_mass   = DViewCArrayKokkos<double>(&elem_interface.mass(0), num_elems);
    elem_mat_id = DViewCArrayKokkos<size_t>(&elem_interface.mat_id(0), num_elems);

    // create Dual Views of the corner struct variables
    corner_force = DViewCArrayKokkos<double>(&corner_interface.force(0, 0), num_corners, num_dim);
    corner_mass  = DViewCArrayKokkos<double>(&corner_interface.mass(0), num_corners);

    // allocate elem_vel_grad
    elem_vel_grad = DCArrayKokkos<double>(num_elems, 3, 3);

    // allocate material models
    elem_eos = DCArrayKokkos<eos_t>(num_elems);
    elem_strength = DCArrayKokkos<strength_t>(num_elems);

    // ---------------------------------------------------------------------
    //   calculate geometry
    // ---------------------------------------------------------------------
    node_coords.update_device();
    Kokkos::fence();

    get_vol();

    // FEA_Module bc variable
    num_boundary_conditions = num_bcs;

    const DCArrayKokkos<boundary_t> boundary = module_params->boundary;
    const DCArrayKokkos<mat_fill_t> mat_fill = simparam->mat_fill;
    const DCArrayKokkos<material_t> material = simparam->material;
    eos_global_vars = simparam->eos_global_vars;
    strength_global_vars = simparam->strength_global_vars;
    eos_state_vars = DCArrayKokkos<double>(rnum_elem, simparam->max_num_eos_state_vars);
    strength_state_vars   = DCArrayKokkos<double>(rnum_elem, simparam->max_num_strength_state_vars);
    elem_user_output_vars = DCArrayKokkos<double>(rnum_elem, simparam->output_options.max_num_user_output_vars);

    // --- calculate bdy sets ---//
    mesh->num_nodes_in_patch  = 2 * (num_dim - 1); // 2 (2D) or 4 (3D)
    mesh->num_patches_in_elem = 2 * num_dim; // 4 (2D) or 6 (3D)

    mesh->init_bdy_sets(num_bcs);
    num_bdy_sets = mesh->num_bdy_sets;
    printf("Num BC's = %lu\n", num_bcs);

    // patch ids in bdy set
    bdy_patches_in_set = mesh->bdy_patches_in_set;
    if (num_dim == 2) {
        bdy_nodes = mesh->bdy_nodes;
    }

    // tag boundary patches in the set
    tag_bdys(boundary, *mesh, node_coords);

    build_boundry_node_sets(*mesh);

    // node ids in bdy_patch set
    bdy_nodes_in_set     = mesh->bdy_nodes_in_set;
    num_bdy_nodes_in_set = mesh->num_bdy_nodes_in_set;

    // assign mesh views needed by the FEA module

    // elem ids in elem
    elems_in_elem     = mesh->elems_in_elem;
    num_elems_in_elem = mesh->num_elems_in_elem;

    // corners
    num_corners_in_node = mesh->num_corners_in_node;
    corners_in_node     = mesh->corners_in_node;
    corners_in_elem     = mesh->corners_in_elem;

    // elem-node conn & node-node conn
    elems_in_node = mesh->elems_in_node;
    if (num_dim == 2) {
        nodes_in_node     = mesh->nodes_in_node;
        num_nodes_in_node = mesh->num_nodes_in_node;
        // patch conn

        patches_in_elem = mesh->patches_in_elem;
        nodes_in_patch  = mesh->nodes_in_patch;
        elems_in_patch  = mesh->elems_in_patch;
    }

    // loop over BCs
    for (size_t this_bdy = 0; this_bdy < num_bcs; this_bdy++) {
        RUN_CLASS({
            printf("Boundary Condition number %lu \n", this_bdy);
            printf("  Num bdy patches in this set = %lu \n", bdy_patches_in_set.stride(this_bdy));
            printf("  Num bdy nodes in this set = %lu \n", bdy_nodes_in_set.stride(this_bdy));
        });
        Kokkos::fence();
    } // end for

    // elem_mat_id needs to be initialized before initialization of material models
    for (int f_id = 0; f_id < num_fills; f_id++) {
        FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
            elem_mat_id(elem_gid) = mat_fill(f_id).material_id;
        });
    }
    elem_mat_id.update_host();

    // function for initializing state_vars
    init_state_vars(material,
                    elem_mat_id,
                    eos_state_vars,
                    strength_state_vars,
                    eos_global_vars,
                    strength_global_vars,
                    elem_user_output_vars,
                    rnum_elem);

    // initialize strength model
    init_strength_model(elem_strength,
                        material,
                        elem_mat_id,
                        eos_state_vars,
                        strength_state_vars,
                        eos_global_vars,
                        strength_global_vars,
                        elem_user_output_vars,
                        rnum_elem);

    // initialize eos model
    init_eos_model(elem_eos,
                   material,
                   elem_mat_id,
                   eos_state_vars,
                   strength_state_vars,
                   eos_global_vars,
                   strength_global_vars,
                   elem_user_output_vars,
                   rnum_elem);

    // --- apply the fill instructions over each of the Elements---//

    // initialize if topology optimization is used
    if (simparam->topology_optimization_on) {
        for (int elem_id = 0; elem_id < rnum_elem; elem_id++) {
            relative_element_densities.host(elem_id) = 1;
        } // for
        relative_element_densities.update_device();
    }

    // loop over the fill instructures
    for (int f_id = 0; f_id < num_fills; f_id++) {
        // parallel loop over elements in mesh
        FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
            // calculate the coordinates and radius of the element
            double elem_coords[3]; // note:initialization with a list won't work
            elem_coords[0] = 0.0;
            elem_coords[1] = 0.0;
            elem_coords[2] = 0.0;

            // get the coordinates of the element center
            for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
                elem_coords[0] += node_coords(rk_level, nodes_in_elem(elem_gid, node_lid), 0);
                elem_coords[1] += node_coords(rk_level, nodes_in_elem(elem_gid, node_lid), 1);
                if (num_dim == 3) {
                    elem_coords[2] += node_coords(rk_level, nodes_in_elem(elem_gid, node_lid), 2);
                }
                else{
                    elem_coords[2] = 0.0;
                }
            } // end loop over nodes in element
            elem_coords[0] = elem_coords[0] / num_nodes_in_elem;
            elem_coords[1] = elem_coords[1] / num_nodes_in_elem;
            elem_coords[2] = elem_coords[2] / num_nodes_in_elem;

            // default is not to fill the element
            bool fill_this = mat_fill(f_id).volume.contains(elem_coords);

            // paint the material state on the element
            if (fill_this) {
                // density
                elem_den(elem_gid) = mat_fill(f_id).den;

                // mass
                elem_mass(elem_gid) = elem_den(elem_gid) * elem_vol(elem_gid);

                // specific internal energy
                elem_sie(rk_level, elem_gid) = mat_fill(f_id).sie;

                size_t mat_id = elem_mat_id(elem_gid); // short name

                // --- stress tensor ---
                // always 3D even for 2D-RZ
                for (size_t i = 0; i < 3; i++) {
                    for (size_t j = 0; j < 3; j++) {
                        elem_stress(rk_level, elem_gid, i, j) = 0.0;
                    }
                }  // end for

                // --- Pressure ---
                elem_eos(elem_gid).calc_pressure(elem_pres,
                                         elem_stress,
                                         elem_gid,
                                         elem_mat_id(elem_gid),
                                         eos_state_vars,
                                         strength_state_vars,
                                         eos_global_vars,
                                         strength_global_vars,
                                         elem_user_output_vars,
                                         elem_sspd,
                                         elem_den(elem_gid),
                                         elem_sie(rk_level, elem_gid));

                // --- Sound speed ---
                elem_eos(elem_gid).calc_sound_speed(elem_pres,
                                            elem_stress,
                                            elem_gid,
                                            elem_mat_id(elem_gid),
                                            eos_state_vars,
                                            strength_state_vars,
                                            eos_global_vars,
                                            strength_global_vars,
                                            elem_user_output_vars,
                                            elem_sspd,
                                            elem_den(elem_gid),
                                            elem_sie(rk_level, elem_gid));

                // loop over the nodes of this element and apply velocity
                for (size_t node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
                    // get the mesh node index
                    size_t node_gid = nodes_in_elem(elem_gid, node_lid);

                    // --- Velocity ---
                    switch (mat_fill(f_id).velocity) {
                        case VELOCITY_TYPE::cartesian:
                            {
                                node_vel(rk_level, node_gid, 0) = mat_fill(f_id).u;
                                node_vel(rk_level, node_gid, 1) = mat_fill(f_id).v;
                                if (num_dim == 3) {
                                    node_vel(rk_level, node_gid, 2) = mat_fill(f_id).w;
                                }

                                break;
                            }
                        case VELOCITY_TYPE::radial:
                            {
                                // Setting up cylindrical
                                double dir[2];
                                dir[0] = 0.0;
                                dir[1] = 0.0;
                                double radius_val = 0.0;

                                for (int dim = 0; dim < 2; dim++) {
                                    dir[dim]    = node_coords(rk_level, node_gid, dim);
                                    radius_val += node_coords(rk_level, node_gid, dim) * node_coords(rk_level, node_gid, dim);
                                } // end for
                                radius_val = sqrt(radius_val);

                                for (int dim = 0; dim < 2; dim++) {
                                    if (radius_val > 1.0e-14) {
                                        dir[dim] /= (radius_val);
                                    }
                                    else{
                                        dir[dim] = 0.0;
                                    }
                                } // end for

                                node_vel(rk_level, node_gid, 0) = mat_fill(f_id).speed * dir[0];
                                node_vel(rk_level, node_gid, 1) = mat_fill(f_id).speed * dir[1];
                                if (num_dim == 3) {
                                    node_vel(rk_level, node_gid, 2) = 0.0;
                                }

                                break;
                            }
                        case VELOCITY_TYPE::spherical:
                            {
                                // Setting up spherical
                                double dir[3];
                                dir[0] = 0.0;
                                dir[1] = 0.0;
                                dir[2] = 0.0;
                                double radius_val = 0.0;

                                for (int dim = 0; dim < 3; dim++) {
                                    dir[dim]    = node_coords(rk_level, node_gid, dim);
                                    radius_val += node_coords(rk_level, node_gid, dim) * node_coords(rk_level, node_gid, dim);
                                } // end for
                                radius_val = sqrt(radius_val);

                                for (int dim = 0; dim < 3; dim++) {
                                    if (radius_val > 1.0e-14) {
                                        dir[dim] /= (radius_val);
                                    }
                                    else{
                                        dir[dim] = 0.0;
                                    }
                                } // end for

                                node_vel(rk_level, node_gid, 0) = mat_fill(f_id).speed * dir[0];
                                node_vel(rk_level, node_gid, 1) = mat_fill(f_id).speed * dir[1];
                                if (num_dim == 3) {
                                    node_vel(rk_level, node_gid, 2) = mat_fill(f_id).speed * dir[2];
                                }

                                break;
                            }
                        case VELOCITY_TYPE::radial_linear:
                            {
                                break;
                            }
                        case VELOCITY_TYPE::spherical_linear:
                            {
                                break;
                            }
                        case VELOCITY_TYPE::tg_vortex:
                            {
                                node_vel(rk_level, node_gid, 0) = sin(PI * node_coords(rk_level, node_gid, 0)) * cos(PI * node_coords(rk_level, node_gid, 1));
                                node_vel(rk_level, node_gid, 1) =  -1.0 * cos(PI * node_coords(rk_level, node_gid, 0)) * sin(PI * node_coords(rk_level, node_gid, 1));
                                if (num_dim == 3) {
                                    node_vel(rk_level, node_gid, 2) = 0.0;
                                }

                                break;
                            }
                    } // end of switch
                } // end loop over nodes of element

                if (mat_fill(f_id).velocity == VELOCITY_TYPE::tg_vortex) {
                    elem_pres(elem_gid) = 0.25 * (cos(2.0 * PI * elem_coords[0]) + cos(2.0 * PI * elem_coords[1]) ) + 1.0;

                    // p = rho*ie*(gamma - 1)
                    size_t mat_id = f_id;

                    double gamma = eos_global_vars(mat_id, 0); // gamma value

                    elem_sie(rk_level, elem_gid) =
                        elem_pres(elem_gid) / (mat_fill(f_id).den * (gamma - 1.0));
                } // end if
            } // end if fill
        }); // end FOR_ALL_CLASS element loop
        Kokkos::fence();
    } // end for loop over fills

    // apply BC's to velocity
    FEA_Module_Dynamic_Elasticity::boundary_velocity(*mesh, boundary, node_vel);

    // calculate the corner massess if 2D
    if (num_dim == 2) {
        FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
            // facial area of the corners
            double corner_areas_array[4];

            ViewCArrayKokkos<double> corner_areas(&corner_areas_array[0], 4);
            ViewCArrayKokkos<size_t> elem_node_gids(&nodes_in_elem(elem_gid, 0), 4);

            get_area_weights2D(corner_areas,
                               elem_gid,
                               node_coords,
                               elem_node_gids,
                               rk_level);

            // loop over the corners of the element and calculate the mass
            for (size_t corner_lid = 0; corner_lid < 4; corner_lid++) {
                size_t corner_gid = corners_in_elem(elem_gid, corner_lid);
                corner_mass(corner_gid) = corner_areas(corner_lid) * elem_den(elem_gid); // node radius is added later
            } // end for over corners
        }); // end of FOR_ALL_CLASS
    } // end of

    // calculate the nodal mass
    FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
        node_mass(node_gid) = 0.0;

        if (num_dim == 3) {
            for (size_t elem_lid = 0; elem_lid < num_corners_in_node(node_gid); elem_lid++) {
                size_t elem_gid      = elems_in_node(node_gid, elem_lid);
                node_mass(node_gid) += 1.0 / 8.0 * elem_mass(elem_gid);
            } // end for elem_lid
        } // end if dims=3
        else{
            // 2D-RZ
            for (size_t corner_lid = 0; corner_lid < num_corners_in_node(node_gid); corner_lid++) {
                size_t corner_gid    = corners_in_node(node_gid, corner_lid);
                node_mass(node_gid) += corner_mass(corner_gid);  // sans the radius so it is areal node mass

                corner_mass(corner_gid) *= node_coords(rk_level, node_gid, 1); // true corner mass now
            } // end for elem_lid
        } // end else
    }); // end FOR_ALL_CLASS
    Kokkos::fence();

    // current interface has differing mass arrays; this equates them until we unify memory
    // view scope
    if (simparam->topology_optimization_on || simparam->shape_optimization_on || simparam->num_dims == 2) {
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
    } // endif

    // initialize if topology optimization is used
    if (simparam->topology_optimization_on || simparam->shape_optimization_on) {
        init_assembly();
        assemble_matrix();
    }

    // update host copies of arrays modified in this function
    elem_den.update_host();
    elem_mass.update_host();
    elem_sie.update_host();
    elem_stress.update_host();
    elem_pres.update_host();
    elem_sspd.update_host();

    return;
} // end of setup

/////////////////////////////////////////////////////////////////////////////
///
/// \fn module_cleanup
///
/// \brief Deallocate memory associated with the module
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::module_cleanup()
{
    cleanup_material_models();
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn cleanup_material_models
///
/// \brief Deallocate memory used for  material models
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::cleanup_material_models()
{
    const DCArrayKokkos<material_t> material = simparam->material;

    // destroy strength model
    destroy_strength_model(elem_strength,
                           material,
                           elem_mat_id,
                           eos_state_vars,
                           strength_state_vars,
                           eos_global_vars,
                           strength_global_vars,
                           elem_user_output_vars,
                           rnum_elem);

    // destroy eos model
    destroy_eos_model(elem_eos,
                      material,
                      elem_mat_id,
                      eos_state_vars,
                      strength_state_vars,
                      eos_global_vars,
                      strength_global_vars,
                      elem_user_output_vars,
                      rnum_elem);
    return;
} // end cleanup_user_strength_model;

/////////////////////////////////////////////////////////////////////////////
///
/// \fn tag_bdys
///
/// \brief Determines which of the boundary patches are associated with which boundary.
///        Modifies: bdy_patches_in_set
///
/// \param Array of boundaries
/// \param The simulation mesh
/// \param Nodal coordinates
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::tag_bdys(const DCArrayKokkos<boundary_t>& boundary,
    mesh_t& mesh,
    const DViewCArrayKokkos<double>& node_coords)
{
    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;
    size_t num_dim = simparam->num_dims;
    int    nboundary_patches  = Explicit_Solver_Pointer_->nboundary_patches;
    int    num_nodes_in_patch = mesh.num_nodes_in_patch;

    // if (bdy_set == mesh.num_bdy_sets){
    //    printf(" ERROR: number of boundary sets must be increased by %zu",
    //              bdy_set-mesh.num_bdy_sets+1);
    //    exit(0);
    // } // end if

    // error and debug flag
    // DCArrayKokkos<bool> print_flag(1, "print_flag");
    // print_flag.host(0) = false;
    // print_flag.update_device();

    FOR_ALL_CLASS(bdy_set, 0, num_bdy_sets, {
        // for (size_t bdy_set = 0; bdy_set < num_bdy_sets; bdy_set++) {

        // tag boundaries
        BOUNDARY_TYPE bc_type = boundary(bdy_set).surface.type;
        double val = boundary(bdy_set).surface.plane_position;

        // save the boundary patches to this set that are on the plane, spheres, etc.
        for (size_t bdy_patch_lid = 0; bdy_patch_lid < nboundary_patches; bdy_patch_lid++) {
            // save the patch index
            size_t bdy_patch_gid = bdy_patch_lid;

            // check to see if this patch is on the specified plane
            bool is_on_bdy = check_bdy(bdy_patch_gid,
                                         num_dim,
                                         num_nodes_in_patch,
                                         bc_type,
                                         val,
                                         node_coords,
                                         rk_level); // no=0, yes=1 WARNING: POSSIBLE BUG
            if (is_on_bdy) {
                size_t index = bdy_patches_in_set.stride(bdy_set);

                // increment the number of boundary patches saved
                bdy_patches_in_set.stride(bdy_set)++;

                bdy_patches_in_set(bdy_set, index) = bdy_patch_gid;
            } // end if
        } // end for bdy_patch
          // }
    });  // end FOR_ALL_CLASS bdy_sets

    // debug check
    // print_flag.update_host();
    // if(print_flag.host(0)) std::cout << "found boundary node with id 549412" << std::endl;

    return;
} // end tag

/////////////////////////////////////////////////////////////////////////////
///
/// \fn check_bdy
///
/// \brief Determines which of the boundary patches are associated with which boundary.
///
/// \param Global index of the patch
/// \param Number of spatial dimensions
/// \param Number of the nodes in the patch
/// \param Type of boundary condition
/// \param Boundary value
/// \param Nodal coordinates
/// \param Runge Kutta integration level
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_INLINE_FUNCTION
bool FEA_Module_Dynamic_Elasticity::check_bdy(const size_t patch_gid,
    const int num_dim,
    const int num_nodes_in_patch,
    const BOUNDARY_TYPE bc_type,
    const double val,
    const DViewCArrayKokkos<double>& node_coords,
    const size_t rk_level) const
{
    // default bool is not on the boundary
    size_t is_on_bdy = 0;

    // the patch coordinates
    double these_patch_coords[3];  // Note: cannot allocated array with num_dim

    // loop over the nodes on the patch
    for (size_t patch_node_lid = 0; patch_node_lid < num_nodes_in_patch; patch_node_lid++) {
        // get the nodal_gid for this node in the patch
        // size_t node_gid = mesh.nodes_in_patch(patch_gid, patch_node_lid);
        size_t node_gid = Local_Index_Boundary_Patches(patch_gid, patch_node_lid);

        for (size_t dim = 0; dim < num_dim; dim++) {
            these_patch_coords[dim] = node_coords(rk_level, node_gid, dim);  // (rk, node_gid, dim)
        }

        if (bc_type == BOUNDARY_TYPE::x_plane) {
            if (fabs(these_patch_coords[0] - val) <= 1.0e-7) {
                is_on_bdy += 1;
            }
        }
        else if (bc_type == BOUNDARY_TYPE::y_plane) {
            if (fabs(these_patch_coords[1] - val) <= 1.0e-7) {
                is_on_bdy += 1;
            }
        }
        else if (bc_type == BOUNDARY_TYPE::z_plane) {
            if (fabs(these_patch_coords[2] - val) <= 1.0e-7) {
                is_on_bdy += 1;
            }
        }
        else if (bc_type == BOUNDARY_TYPE::cylinder) {
            real_t R = sqrt(these_patch_coords[0] * these_patch_coords[0] +
                          these_patch_coords[1] * these_patch_coords[1]);

            if (fabs(R - val) <= 1.0e-7) {
                is_on_bdy += 1;
            }
        }
        else if (bc_type == BOUNDARY_TYPE::sphere) {
            real_t R = sqrt(these_patch_coords[0] * these_patch_coords[0] +
                          these_patch_coords[1] * these_patch_coords[1] +
                          these_patch_coords[2] * these_patch_coords[2]);

            if (fabs(R - val) <= 1.0e-7) {
                is_on_bdy += 1;
            }
        }
    }

    // if all nodes in the patch are on the surface
    return is_on_bdy == num_nodes_in_patch;
} // end method to check bdy

/////////////////////////////////////////////////////////////////////////////
///
/// \fn build_boundry_node_sets
///
/// \brief Build set of nodes assigned to each boundary condition
///
/// \param The simulation mesh
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::build_boundry_node_sets(mesh_t& mesh)
{
    // build boundary nodes in each boundary set
    int nboundary_patches  = Explicit_Solver_Pointer_->nboundary_patches;
    int num_nodes_in_patch = mesh.num_nodes_in_patch;
    num_bdy_nodes_in_set = mesh.num_bdy_nodes_in_set = DCArrayKokkos<size_t>(num_bdy_sets, "num_bdy_nodes_in_set");
    CArrayKokkos<long long int> temp_count_num_bdy_nodes_in_set(num_bdy_sets, nall_nodes, "temp_count_num_bdy_nodes_in_set");

    DynamicRaggedRightArrayKokkos<size_t> temp_nodes_in_set(mesh.num_bdy_sets, nboundary_patches * mesh.num_nodes_in_patch, "temp_nodes_in_set");

    // Parallel loop over boundary sets on device
    FOR_ALL_CLASS(bdy_set, 0, num_bdy_sets, {
        // finde the number of patches_in_set
        size_t num_bdy_patches_in_set = bdy_patches_in_set.stride(bdy_set);

        num_bdy_nodes_in_set(bdy_set) = 0;

        // Loop over boundary patches in boundary set
        for (size_t bdy_patch_gid = 0; bdy_patch_gid < num_bdy_patches_in_set; bdy_patch_gid++) {
            // get the global id for this boundary patch
            size_t patch_gid = bdy_patches_in_set(bdy_set, bdy_patch_gid);

            // apply boundary condition at nodes on boundary
            for (size_t node_lid = 0; node_lid < num_nodes_in_patch; node_lid++) {
                size_t node_gid = Local_Index_Boundary_Patches(patch_gid, node_lid);

                temp_count_num_bdy_nodes_in_set(bdy_set, node_gid) = -1;
            }     // end for node_lid
        } // end for bdy_patch_gid

        // Loop over boundary patches in boundary set
        for (size_t bdy_patch_gid = 0; bdy_patch_gid < num_bdy_patches_in_set; bdy_patch_gid++) {
            // get the global id for this boundary patch
            size_t patch_gid = bdy_patches_in_set(bdy_set, bdy_patch_gid);

            // apply boundary condition at nodes on boundary
            for (size_t node_lid = 0; node_lid < num_nodes_in_patch; node_lid++) {
                size_t node_gid = Local_Index_Boundary_Patches(patch_gid, node_lid);

                if (temp_count_num_bdy_nodes_in_set(bdy_set, node_gid) == -1) {
                    size_t num_saved = num_bdy_nodes_in_set(bdy_set);

                    num_bdy_nodes_in_set(bdy_set)++;

                    // replace -1 with node_gid to denote the node was already saved
                    temp_count_num_bdy_nodes_in_set(bdy_set, node_gid) = node_gid;

                    // increment the number of saved nodes, create memory
                    temp_nodes_in_set.stride(bdy_set)++;
                    temp_nodes_in_set(bdy_set, num_saved) = node_gid;
                }     // end if
            }     // end for node_lid
        } // end for bdy_patch_gid
    }); // end FOR_ALL_CLASS bdy_set
    Kokkos::fence();

    // allocate the RaggedRight bdy_nodes_in_set array
    bdy_nodes_in_set = mesh.bdy_nodes_in_set = RaggedRightArrayKokkos<size_t>(mesh.num_bdy_nodes_in_set, "bdy_nodes_in_set");

    FOR_ALL_CLASS(bdy_set, 0, num_bdy_sets, {
        // Loop over boundary patches in boundary set
        for (size_t bdy_node_lid = 0; bdy_node_lid < num_bdy_nodes_in_set(bdy_set); bdy_node_lid++) {
            // save the bdy_node_gid
            bdy_nodes_in_set(bdy_set, bdy_node_lid) = temp_nodes_in_set(bdy_set, bdy_node_lid);
        } // end for
    }); // end FOR_ALL_CLASS bdy_set

    // update the host side for the number nodes in a bdy_set
    num_bdy_nodes_in_set.update_host();

    return;
} // end method to build boundary nodes

/////////////////////////////////////////////////////////////////////////////
///
/// \fn solve
///
/// \brief Solve function called by solver
///
/// Calls elastic_solve
///
/// \return 0 if the solver finishes
///
/////////////////////////////////////////////////////////////////////////////
int FEA_Module_Dynamic_Elasticity::solve()
{
    elastic_solve();

    return 0;
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn elastic_solve
///
/// \brief Elastic solver loop
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_Dynamic_Elasticity::elastic_solve()
{
    Dynamic_Options dynamic_options = simparam->dynamic_options;

    const size_t rk_level = dynamic_options.rk_num_bins - 1;
    const int    num_dim  = simparam->num_dims;

    time_value       = dynamic_options.time_initial;
    time_final       = dynamic_options.time_final;
    dt_max           = dynamic_options.dt_max;
    dt_min           = dynamic_options.dt_min;
    dt_cfl           = dynamic_options.dt_cfl;
    graphics_time    = simparam->output_options.graphics_step;
    graphics_dt_ival = simparam->output_options.graphics_step;
    cycle_stop       = dynamic_options.cycle_stop;
    rk_num_stages    = dynamic_options.rk_num_stages;
    dt    = dynamic_options.dt;
    fuzz  = dynamic_options.fuzz;
    tiny  = dynamic_options.tiny;
    small = dynamic_options.small;
    graphics_times = simparam->output_options.graphics_times;
    graphics_id    = simparam->output_options.graphics_id;
    size_t num_bdy_nodes = mesh->num_bdy_nodes;

    const DCArrayKokkos<boundary_t> boundary = module_params->boundary;
    const DCArrayKokkos<material_t> material = simparam->material;
    int old_max_forward_buffer;

    unsigned long cycle;

    real_t objective_accumulation, global_objective_accumulation;

    problem = Explicit_Solver_Pointer_->problem; // Pointer to ROL optimization problem object
    ROL::Ptr<ROL::Objective<real_t>> obj_pointer;

    // reset time accumulating objective and constraints
    if (simparam->topology_optimization_on) {
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

    int myrank = Explicit_Solver_Pointer_->myrank;
    if (simparam->output_options.output_file_format == OUTPUT_FORMAT::vtk && simparam->output_options.write_initial) {
        if (myrank == 0) {
            printf("Writing outputs to file at %f \n", time_value);
        }

        double comm_time1 = Explicit_Solver_Pointer_->CPU_Time();
        Explicit_Solver_Pointer_->write_outputs();
        double comm_time2 = Explicit_Solver_Pointer_->CPU_Time();
        Explicit_Solver_Pointer_->output_time += comm_time2 - comm_time1;
    }

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
    if (simparam->topology_optimization_on || simparam->shape_optimization_on) {
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

        (*forward_solve_velocity_data)[0]->assign(*Explicit_Solver_Pointer_->all_node_velocities_distributed);
        (*forward_solve_coordinate_data)[0]->assign(*Explicit_Solver_Pointer_->all_node_coords_distributed);
    }

    // loop over the max number of time integration cycles
    for (cycle = 0; cycle < cycle_stop; cycle++) {
        // get the step
        if (num_dim == 2) {
            get_timestep2D(*mesh,
                           node_coords,
                           node_vel,
                           elem_sspd,
                           elem_vol);
        }
        else{
            get_timestep(*mesh,
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

        if (simparam->dynamic_options.output_time_sequence_level >= TIME_OUTPUT_LEVEL::high) {
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
        }

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

            // ---- Calculate velocity diveregence for the element ----
            get_divergence(elem_div,
                               *mesh,
                               node_coords,
                               node_vel,
                               elem_vol);

            // ---- Update nodal velocities ---- //
            get_force_elastic(material,
                            *mesh,
                            node_coords,
                            node_vel,
                            node_mass,
                            elem_den,
                            elem_vol,
                            elem_div,
                            elem_mat_id,
                            corner_force,
                            rk_alpha,
                            cycle);
            applied_forces(material,
                            *mesh,
                            node_coords,
                            node_vel,
                            node_mass,
                            elem_den,
                            elem_vol,
                            elem_div,
                            elem_mat_id,
                            corner_force,
                            rk_alpha,
                            cycle);

            // ---- apply force boundary conditions to the boundary patches----
            boundary_velocity(*mesh, boundary, node_vel);

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

            // ---- Update nodal positions ----
            update_position_elastic(rk_alpha,
                                nall_nodes,
                                node_coords,
                                node_vel);

            // ---- Calculate cell volume for next time step ----
            get_vol();

            // ---- Calculate elem state (den, pres, sound speed, stress) for next time step ----
            /*
            if(num_dim==2){
                update_state2D(material,
                               *mesh,
                               node_coords,
                               node_vel,
                               elem_den,
                               elem_pres,
                               elem_stress,
                               elem_sspd,
                               elem_sie,
                               elem_vol,
                               elem_mass,
                               elem_mat_id,
                               rk_alpha,
                               cycle);
            }
            else{
                update_state(material,
                             *mesh,
                             node_coords,
                             node_vel,
                             elem_den,
                             elem_pres,
                             elem_stress,
                             elem_sspd,
                             elem_sie,
                             elem_vol,
                             elem_mass,
                             elem_mat_id,
                             rk_alpha,
                             cycle);
            }
            */
            // ----
            // Notes on strength:
            //    1) hyper-elastic strength models are called in update_state
            //    2) hypo-elastic strength models are called in get_force
            //    3) strength models must be added by the user in user_mat.cpp

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
        Explicit_Solver_Pointer_->time_value = simparam->dynamic_options.time_value = time_value += dt;

        if (simparam->topology_optimization_on || simparam->shape_optimization_on) {
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
                        ke += (node_velocities_interface(node_gid, dim) + previous_node_velocities_interface(node_gid, dim)) * (node_velocities_interface(node_gid,
                        dim) + previous_node_velocities_interface(node_gid, dim)) / 4;                                                                                                                     // 1/2 at end
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
        }

        size_t write = 0;
        if ((cycle + 1) % graphics_cyc_ival == 0 && cycle > 0) {
            write = 1;
        }
        else if (cycle == cycle_stop) {
            write = 1;
        }
        else if (time_value >= time_final && simparam->output_options.write_final) {
            write = 1;
        }
        else if (time_value >= graphics_time) {
            write = 1;
        }

        // write outputs
        if (write == 1) {
            // interface nodal coordinate data (note: this is not needed if using write_outputs())
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
                Explicit_Solver_Pointer_->write_outputs();

                double comm_time2 = Explicit_Solver_Pointer_->CPU_Time();
                Explicit_Solver_Pointer_->output_time += comm_time2 - comm_time1;
            }

            graphics_time = time_value + graphics_dt_ival;
        } // end if

        // end of calculation
        if (time_value >= time_final) {
            break;
        }
    } // end for cycle loop

    last_time_step = cycle;

    // simple setup to just calculate KE minimize objective for now
    if (simparam->topology_optimization_on) {
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
        printf("total energy conservation error = %e \n\n", 100 * (TE_tend - TE_t0) / TE_t0);
    }

    return;
} // end of elastic solve
