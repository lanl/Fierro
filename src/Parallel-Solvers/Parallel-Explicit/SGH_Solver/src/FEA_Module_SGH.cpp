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
#include "dynamic_checkpoint.h"
#include "Simulation_Parameters/Simulation_Parameters_Explicit.h"
#include "Simulation_Parameters/FEA_Module/SGH_Parameters.h"
#include "FEA_Module_SGH.h"
#include "Explicit_Solver.h"

// optimization
#include "ROL_Solver.hpp"
#include "Fierro_Optimization_Objective.hpp"

#define MAX_ELEM_NODES 8
#define STRAIN_EPSILON 0.000000001
#define DENSITY_EPSILON 0.0001
#define BC_EPSILON 1.0e-6
#define BUFFER_GROW 100

// #define DEBUG

using namespace utils;

FEA_Module_SGH::FEA_Module_SGH(
    SGH_Parameters& params, Solver* Solver_Pointer,
    std::shared_ptr<mesh_t> mesh_in, const int my_fea_module_index)
    : FEA_Module(Solver_Pointer)
{
    // assign interfacing index
    my_fea_module_index_ = my_fea_module_index;
    Module_Type = FEA_MODULE_TYPE::SGH;

    // recast solver pointer for non-base class access
    Explicit_Solver_Pointer_ = dynamic_cast<Explicit_Solver*>(Solver_Pointer);
    simparam = &(Explicit_Solver_Pointer_->simparam);
    module_params = &params;

    // create ref element object
    // ref_elem = new elements::ref_element();
    // create mesh objects
    // init_mesh = new swage::mesh_t(simparam);
    // mesh = new swage::mesh_t(simparam);

    mesh = mesh_in;

    // boundary condition data
    max_boundary_sets = 0;
    Local_Index_Boundary_Patches = Explicit_Solver_Pointer_->Local_Index_Boundary_Patches;

    // set Tpetra vector pointers
    initial_node_velocities_distributed = Explicit_Solver_Pointer_->initial_node_velocities_distributed;
    node_velocities_distributed     = Explicit_Solver_Pointer_->node_velocities_distributed;
    all_node_velocities_distributed = Explicit_Solver_Pointer_->all_node_velocities_distributed;

    // Switch for optimization solver
    if (simparam->topology_optimization_on || simparam->shape_optimization_on) {
        cached_design_gradients_distributed    = Teuchos::rcp(new MV(map, 1));
        all_cached_node_velocities_distributed = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
        force_gradient_velocity                = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
        force_gradient_position                = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
        force_gradient_design                  = Teuchos::rcp(new MV(all_node_map, 1));
        corner_value_storage                   = Solver_Pointer->corner_value_storage;
        corner_vector_storage                  = Solver_Pointer->corner_vector_storage;
        corner_gradient_storage                = Solver_Pointer->corner_gradient_storage;
        relative_element_densities             = DCArrayKokkos<double>(rnum_elem, "relative_element_densities");
        all_adjoint_vector_distributed         = Teuchos::rcp(new MV(all_node_map, num_dim));
        adjoint_vector_distributed             = Teuchos::rcp(new MV(*all_adjoint_vector_distributed, map));
        all_phi_adjoint_vector_distributed     = Teuchos::rcp(new MV(all_node_map, num_dim));
        phi_adjoint_vector_distributed         = Teuchos::rcp(new MV(*all_phi_adjoint_vector_distributed, map));
        psi_adjoint_vector_distributed         = Teuchos::rcp(new MV(all_element_map, 1));
    }

    if (simparam->topology_optimization_on || simparam->shape_optimization_on || simparam->num_dims == 2) {
        node_masses_distributed = Teuchos::rcp(new MV(map, 1));
        ghost_node_masses_distributed = Teuchos::rcp(new MV(ghost_node_map, 1));
    }

    // setup output
    noutput = 0;
    init_output();

    // set parameters
    Dynamic_Options dynamic_options = simparam->dynamic_options;
    time_value        = dynamic_options.time_value;
    time_final        = dynamic_options.time_final;
    dt_max            = dynamic_options.dt_max;
    dt_min            = dynamic_options.dt_min;
    dt_cfl            = dynamic_options.dt_cfl;
    graphics_time     = simparam->output_options.graphics_time;
    graphics_dt_ival  = simparam->output_options.graphics_dt_ival;
    graphics_cyc_ival = simparam->output_options.graphics_cyc_ival;
    cycle_stop        = dynamic_options.cycle_stop;
    rk_num_stages     = dynamic_options.rk_num_stages;
    dt    = dynamic_options.dt;
    fuzz  = dynamic_options.fuzz;
    tiny  = dynamic_options.tiny;
    small = dynamic_options.small;
    graphics_times = simparam->output_options.graphics_times;
    graphics_id    = simparam->output_options.graphics_id;
    rk_num_bins    = simparam->dynamic_options.rk_num_bins;

    if (simparam->topology_optimization_on) {
        if(simparam->optimization_options.use_solve_checkpoints){
            max_time_steps                               = simparam->optimization_options.num_solve_checkpoints;
            dynamic_checkpoint_set                       = Teuchos::rcp(new std::set<Dynamic_Checkpoint>());
            cached_dynamic_checkpoints                   = Teuchos::rcp(new std::vector<Dynamic_Checkpoint>());
            previous_node_velocities_distributed         = Teuchos::rcp(new MV(all_node_map, num_dim));
            previous_node_coords_distributed             = Teuchos::rcp(new MV(all_node_map, num_dim));
            previous_element_internal_energy_distributed = Teuchos::rcp(new MV(all_element_map, 1));
            previous_adjoint_vector_distributed          = Teuchos::rcp(new MV(all_node_map, num_dim));
            previous_phi_adjoint_vector_distributed      = Teuchos::rcp(new MV(all_node_map, num_dim));
            previous_psi_adjoint_vector_distributed      = Teuchos::rcp(new MV(all_element_map, 1));
            midpoint_adjoint_vector_distributed          = Teuchos::rcp(new MV(all_node_map, num_dim));
            midpoint_phi_adjoint_vector_distributed      = Teuchos::rcp(new MV(all_node_map, num_dim));
            midpoint_psi_adjoint_vector_distributed      = Teuchos::rcp(new MV(all_element_map, 1));
        }
        else{
            max_time_steps = BUFFER_GROW;
        }
        elem_power_dgradients = DCArrayKokkos<real_t>(rnum_elem);
        element_internal_energy_distributed = Teuchos::rcp(new MV(all_element_map, 1));
    }

    if (simparam->topology_optimization_on) {
        forward_solve_velocity_data   = Teuchos::rcp(new std::vector<Teuchos::RCP<MV>>(max_time_steps + 1));
        forward_solve_coordinate_data = Teuchos::rcp(new std::vector<Teuchos::RCP<MV>>(max_time_steps + 1));
        forward_solve_internal_energy_data = Teuchos::rcp(new std::vector<Teuchos::RCP<MV>>(max_time_steps + 1));
        if(!simparam->optimization_options.use_solve_checkpoints){
            time_data.resize(max_time_steps + 1);
            adjoint_vector_data     = Teuchos::rcp(new std::vector<Teuchos::RCP<MV>>(max_time_steps + 1));
            phi_adjoint_vector_data = Teuchos::rcp(new std::vector<Teuchos::RCP<MV>>(max_time_steps + 1));
            psi_adjoint_vector_data = Teuchos::rcp(new std::vector<Teuchos::RCP<MV>>(max_time_steps + 1));
        }

        // assign a multivector of corresponding size to each new timestep in the buffer
        for (int istep = 0; istep < max_time_steps + 1; istep++) {
            (*forward_solve_velocity_data)[istep]   = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
            (*forward_solve_coordinate_data)[istep] = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
            (*forward_solve_internal_energy_data)[istep] = Teuchos::rcp(new MV(all_element_map, 1));
            
            if(!simparam->optimization_options.use_solve_checkpoints){
                (*adjoint_vector_data)[istep]     = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
                (*phi_adjoint_vector_data)[istep] = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
                (*psi_adjoint_vector_data)[istep] = Teuchos::rcp(new MV(all_element_map, 1));
            }
        }
    }
    
    num_active_checkpoints = 0;
    have_loading_conditions = false;
    // if(module_params->matar_mpi_test){
    //     //test partition map wrappers
    //     DCArrayKokkos<long long int, array_layout, device_type, memory_traits> input_indices = DCArrayKokkos<long long int, array_layout, device_type, memory_traits>(all_node_map->getLocalNumElements(), "map_indices_test");
    //     auto map_indices = all_node_map->getMyGlobalIndices();
    //     for(int i =0; i < all_node_map->getLocalNumElements(); i++){
    //         input_indices(i) = map_indices(i);
    //     }
    //     mtr_map = TpetraPartitionMap<long long int, array_layout, device_type, memory_traits>(input_indices);

    //     DCArrayKokkos<long long int, array_layout, device_type, memory_traits> input_local_indices = DCArrayKokkos<long long int, array_layout, device_type, memory_traits>(map->getLocalNumElements(), "local_map_indices_test");
    //     map_indices = map->getMyGlobalIndices();
    //     for(int i =0; i < map->getLocalNumElements(); i++){
    //         input_local_indices(i) = map_indices(i);
    //     }
    //     mtr_local_map = TpetraPartitionMap<long long int, array_layout, device_type, memory_traits>(input_local_indices);

    //     //mtr_map.tpetra_map->describe(*fos, Teuchos::VERB_EXTREME);
    //     //mtr_node_velocities_distributed = TpetraMVArray<real_t, array_layout, device_type, memory_traits>(all_node_map->getLocalNumElements(), num_dim, all_node_map, "mtr_node_velocities_distributed");
    //     mtr_node_velocities_distributed = TpetraMVArray<real_t, array_layout, device_type, memory_traits>(all_node_map->getLocalNumElements(), num_dim, mtr_map, "mtr_node_velocities_distributed");
    //     //all_node_map->describe(*fos, Teuchos::VERB_EXTREME);
    //     // for(int i = 0; i < all_node_map->getLocalNumElements(); i++){
    //     //     mtr_node_velocities_distributed(i) = 1;
    //     //     //std::cout << mtr_node_velocities_distributed(i);
    //     // }
    //     mtr_node_velocities_distributed.own_comm_setup(mtr_local_map);
    // }
}

FEA_Module_SGH::~FEA_Module_SGH()
{
    // delete simparam;
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn read_conditions_ansys_dat
///
/// \brief Read ANSYS dat format mesh file
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::read_conditions_ansys_dat(std::ifstream* in, std::streampos before_condition_header)
{
    auto input_options = simparam->input_options.value();

    char ch;
    std::string skip_line, read_line, substring, token;
    std::stringstream line_parse, line_parse2;

    int num_dim = simparam->num_dims;
    int buffer_lines = 1000;
    int max_word     = 30;
    int local_node_index, current_column_index;
    int p_order = input_options.p_order;
    int buffer_loop, buffer_iteration, buffer_iterations, scan_loop, nodes_per_element, words_per_line;

    size_t read_index_start, node_rid, elem_gid;
    size_t strain_count;

    real_t unit_scaling = input_options.unit_scaling;
    real_t dof_value;

    CArrayKokkos<char, array_layout, HostSpace, memory_traits> read_buffer;
    CArrayKokkos<long long int, array_layout, HostSpace, memory_traits> read_buffer_indices;

    LO local_dof_id;
    GO node_gid;

    host_vec_array node_densities;
} // end read_conditions_ansys_dat

/////////////////////////////////////////////////////////////////////////////
///
/// \fn generate_bcs
///
/// \brief Assign sets of element boundary surfaces corresponding to user BCs
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::generate_bcs()
{
} // end generate_bcs

/////////////////////////////////////////////////////////////////////////////
///
/// \fn Displacement_Boundary_Conditions
///
/// \brief Apply displacement boundary conditions
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::Displacement_Boundary_Conditions()
{
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn output_control
///
/// \brief Output field settings and file settings
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::output_control()
{
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn init_output
///
/// \brief  Initialize output data structures
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::init_output()
{
    // check user parameters for output
    bool output_velocity_flag = simparam->output(FIELD::velocity);
    bool output_strain_flag   = simparam->output(FIELD::strain);
    bool output_stress_flag   = simparam->output(FIELD::stress);
    int  num_dim = simparam->num_dims;
    int  Brows;
    if (num_dim == 3) {
        Brows = 6;
    }
    else{
        Brows = 3;
    }

    // Implicit compliant code
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
///        Inactive in this context
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::sort_output(Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> sorted_map)
{
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn write_data
///
/// \brief Populate requests this module makes for output data
///
/// \param Scalar point data
/// \param Vector point data
/// \param Scalar cell data (double)
/// \param Scalar cell data (int)
/// \param Cell field data
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::write_data(std::map<std::string, const double*>& point_data_scalars_double,
    std::map<std::string, const double*>& point_data_vectors_double,
    std::map<std::string, const double*>& cell_data_scalars_double,
    std::map<std::string, const int*>&    cell_data_scalars_int,
    std::map<std::string, std::pair<const double*, size_t>>& cell_data_fields_double)
{
    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;

    for (const auto& field_name : simparam->output_options.output_fields) {
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

        case FIELD::SIE:
            // element "SIE"
            elem_sie.update_host();
            cell_data_scalars_double["SIE"] = &elem_sie.host(rk_level, 0);
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

        case FIELD::sound_speed:
            // element "sspd"
            elem_sspd.update_host();
            cell_data_scalars_double["sound_speed"] = elem_sspd.host_pointer();
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
/// \brief Prompts sorting for elementoutput data. For now, density.
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::sort_element_output(Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> sorted_map)
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
/// \brief Inactive in this context
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::collect_output(Teuchos::RCP<Tpetra::Map<LO, GO, node_type>> global_reduce_map)
{
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn compute_output
///
/// \brief Inactive in this context
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::compute_output()
{
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn comm_node_masses
///
/// \brief Communicate updated nodal mass to ghost nodes
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::comm_node_masses()
{
    // debug print of design vector
#ifdef DEBUG
    std::ostream& out = std::cout;
    Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    if (myrank == 0) {
        *fos << "Density data :" << std::endl;
    }
    node_densities_distributed->describe(*fos, Teuchos::VERB_EXTREME);
    *fos << std::endl;
    std::fflush(stdout);
    communicate design densities
    create import object using local node indices map and all indices map
    Tpetra::Import<LO, GO> importer(map, ghost_node_map);
#endif

    // comms to get ghosts
    ghost_node_masses_distributed->doImport(*node_masses_distributed, *ghost_importer, Tpetra::INSERT);

#ifdef DEBUG
    all_node_map->describe(*fos, Teuchos::VERB_EXTREME);
    all_node_velocities_distributed->describe(*fos, Teuchos::VERB_EXTREME);

    update_count++;
    if (update_count == 1) {
        MPI_Barrier(world);
        MPI_Abort(world, 4);
    }
#endif
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn module_cleanup
///
/// \brief Cleanup function called by solver
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::module_cleanup()
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
void FEA_Module_SGH::cleanup_material_models()
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
///
/// Modifies: bdy_patches_in_set
///
/// \param Array of boundary sets
/// \param Simulation mesh
/// \param View into nodal position data
///
/// \return <return type and definition description if not void>
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::tag_bdys(const DCArrayKokkos<boundary_t>& boundary,
    mesh_t& mesh,
    const DViewCArrayKokkos<double>& node_coords)
{
    const size_t rk_level = simparam->dynamic_options.rk_num_bins - 1;
    size_t num_dim = simparam->num_dims;
    int    nboundary_patches  = Explicit_Solver_Pointer_->nboundary_patches;
    int    num_nodes_in_patch = mesh.num_nodes_in_patch;
#ifdef DEBUG
    if (bdy_set == mesh.num_bdy_sets) {
        printf(" ERROR: number of boundary sets must be increased by %zu",
                 bdy_set - mesh.num_bdy_sets + 1);
        exit(0);
    } // end if

    // error and debug flag
    DCArrayKokkos<bool> print_flag(1, "print_flag");
    print_flag.host(0) = false;
    print_flag.update_device();
#endif

    FOR_ALL_CLASS(bdy_set, 0, num_bdy_sets, {
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
                                         rk_level); // no=0, yes=1

            // debug check
#ifdef DEBUG
            for (size_t patch_node_lid = 0; patch_node_lid < mesh.num_nodes_in_patch; patch_node_lid++) {
                size_t node_gid = mesh.nodes_in_patch(bdy_patch_gid, patch_node_lid);
                if (bdy_node_gid == 549412) {
                    print_flag(0) = true;
                }
            }
#endif

            if (is_on_bdy) {
                size_t index = bdy_patches_in_set.stride(bdy_set);

                // increment the number of boundary patches saved
                bdy_patches_in_set.stride(bdy_set)++;

                bdy_patches_in_set(bdy_set, index) = bdy_patch_gid;
            } // end if
        } // end for bdy_patch
    });  // end FOR_ALL_CLASS bdy_sets

#ifdef DEBUG
    print_flag.update_host();
    if (print_flag.host(0)) {
        std::cout << "found boundary node with id 549412" << std::endl;
    }
#endif

    return;
} // end tag

/////////////////////////////////////////////////////////////////////////////
///
/// \fn check_bdy
///
/// \brief Checks if a given patch is on a defined boundary
///
/// The test is done by checking if the centroid of the patch is within
/// 1.0e-7 of the boundary
///
/// \param Global index of a patch
/// \param Number of spatial dimensions
/// \param Number of nodes in a patch
/// \param Type of boundary condition applied
/// \param View of nodal position data
/// \param Current Runge Kutta time integration step
///
/// \return True if patch is on the surface of a defined boundary, else false
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_INLINE_FUNCTION
bool FEA_Module_SGH::check_bdy(const size_t patch_gid,
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
/// \fn solve
///
/// \brief Solve function called by solver
///
/////////////////////////////////////////////////////////////////////////////
int FEA_Module_SGH::solve()
{
    sgh_solve();

    return 0;
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn sgh_solve
///
/// \brief SGH solver loop
///
/////////////////////////////////////////////////////////////////////////////
void FEA_Module_SGH::sgh_solve()
{
    Dynamic_Options dynamic_options = simparam->dynamic_options;

    const int    num_dim  = simparam->num_dims;
    const size_t rk_level = dynamic_options.rk_num_bins - 1;

    const DCArrayKokkos<boundary_t> boundary = module_params->boundary;
    const DCArrayKokkos<material_t> material = simparam->material;

    time_value = dynamic_options.time_initial;
    time_final = dynamic_options.time_final;
    dt_max     = dynamic_options.dt_max;
    dt_min     = dynamic_options.dt_min;
    dt     = dynamic_options.dt;
    dt_cfl = dynamic_options.dt_cfl;
    int print_cycle = dynamic_options.print_cycle;

    graphics_time    = simparam->output_options.graphics_step;
    graphics_dt_ival = simparam->output_options.graphics_step;
    cycle_stop     = dynamic_options.cycle_stop;
    rk_num_stages  = dynamic_options.rk_num_stages;
    graphics_times = simparam->output_options.graphics_times;
    graphics_id    = simparam->output_options.graphics_id;

    fuzz  = dynamic_options.fuzz;
    tiny  = dynamic_options.tiny;
    small = dynamic_options.small;

    double cached_pregraphics_dt = fuzz;

    size_t num_bdy_nodes = mesh->num_bdy_nodes;
    size_t cycle;
    real_t objective_accumulation, global_objective_accumulation;
    int old_max_forward_buffer;

    problem = Explicit_Solver_Pointer_->problem; // Pointer to ROL optimization problem object
    ROL::Ptr<ROL::Objective<real_t>> obj_pointer;
    
    bool topology_optimization_on = simparam->topology_optimization_on;
    bool shape_optimization_on    = simparam->shape_optimization_on;
    bool use_solve_checkpoints    = simparam->optimization_options.use_solve_checkpoints;
    int  num_solve_checkpoints    = simparam->optimization_options.num_solve_checkpoints;
    std::set<Dynamic_Checkpoint>::iterator current_checkpoint, last_raised_checkpoint, dispensable_checkpoint, search_end;
    int  last_raised_level = 0;
    bool dispensable_found = false;
    bool optimization_on = simparam->topology_optimization_on||simparam->shape_optimization_on;
    num_active_checkpoints = 0;
    bool time_accumulation;

    TpetraDFArray<double> test_nodes(num_nodes, num_dim);

    
    if(simparam->optimization_options.disable_forward_solve_output){
        //sets the value large enough to not write during the sgh loop
        graphics_time = time_final + graphics_time;
    }

    // reset time accumulating objective and constraints
    if (topology_optimization_on||shape_optimization_on) {
        obj_pointer = problem->getObjective();
        objective_function = dynamic_cast<FierroOptimizationObjective*>(obj_pointer.getRawPtr());
        time_accumulation = objective_function->time_accumulation;
        if(!use_solve_checkpoints){
            if (max_time_steps + 1 > forward_solve_velocity_data->size()) {
                old_max_forward_buffer = forward_solve_velocity_data->size();
                time_data.resize(max_time_steps + 1);
                forward_solve_velocity_data->resize(max_time_steps + 1);
                forward_solve_coordinate_data->resize(max_time_steps + 1);
                forward_solve_internal_energy_data->resize(max_time_steps + 1);
                adjoint_vector_data->resize(max_time_steps + 1);
                phi_adjoint_vector_data->resize(max_time_steps + 1);
                psi_adjoint_vector_data->resize(max_time_steps + 1);
                // assign a multivector of corresponding size to each new timestep in the buffer
                for (int istep = old_max_forward_buffer; istep < max_time_steps + 1; istep++) {
                    (*forward_solve_velocity_data)[istep]   = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
                    (*forward_solve_coordinate_data)[istep] = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
                    (*forward_solve_internal_energy_data)[istep] = Teuchos::rcp(new MV(all_element_map, 1));
                    (*adjoint_vector_data)[istep]     = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
                    (*phi_adjoint_vector_data)[istep] = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
                    (*psi_adjoint_vector_data)[istep] = Teuchos::rcp(new MV(all_element_map, 1));
                }
            }
        }
    }

    int myrank = Explicit_Solver_Pointer_->myrank;
    if (simparam->output_options.write_initial&&!simparam->optimization_options.disable_forward_solve_output) {
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
    if (topology_optimization_on || shape_optimization_on) {
        if(!use_solve_checkpoints){
            time_data[0] = 0;
        }
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
            vec_array all_node_velocities_interface = Explicit_Solver_Pointer_->all_node_velocities_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
            vec_array element_internal_energy     = element_internal_energy_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
            const_vec_array node_coords_interface = Explicit_Solver_Pointer_->node_coords_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
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

            // interface for element internal energies
            FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
                element_internal_energy(elem_gid, 0) = elem_sie(rk_level, elem_gid);
            }); // end parallel for
            Kokkos::fence();
        } // end view scope

        if(use_solve_checkpoints){
            //reset containers
            dynamic_checkpoint_set->clear();
            cached_dynamic_checkpoints->clear();
            //always assign t=0 as a checkpoint
            Dynamic_Checkpoint temp(3,0,time_value, dt, std::numeric_limits<int>::max());
            temp.change_vector(U_DATA, (*forward_solve_coordinate_data)[0]);
            temp.change_vector(V_DATA, (*forward_solve_velocity_data)[0]);
            temp.change_vector(SIE_DATA, (*forward_solve_internal_energy_data)[0]);
            temp.assign_vector(U_DATA,all_node_coords_distributed);
            temp.assign_vector(V_DATA,all_node_velocities_distributed);
            temp.assign_vector(SIE_DATA,element_internal_energy_distributed);
            dynamic_checkpoint_set->insert(temp);
        }
        if(!use_solve_checkpoints){
            (*forward_solve_internal_energy_data)[0]->assign(*element_internal_energy_distributed);
            (*forward_solve_velocity_data)[0]->assign(*Explicit_Solver_Pointer_->all_node_velocities_distributed);
            (*forward_solve_coordinate_data)[0]->assign(*Explicit_Solver_Pointer_->all_node_coords_distributed);
        }
    }

    // loop over the max number of time integration cycles
    for (cycle = 0; cycle < cycle_stop; cycle++) {

        //save timestep from before graphics output contraction
        cached_pregraphics_dt = dt;

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
            else if (cycle % print_cycle == 0) {
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

        if(use_solve_checkpoints&&optimization_on){
                previous_node_velocities_distributed->assign(*all_node_velocities_distributed);
                previous_node_coords_distributed->assign(*all_node_coords_distributed);
                previous_element_internal_energy_distributed->assign(*element_internal_energy_distributed);
        }

        // integrate solution forward in time
        for (size_t rk_stage = 0; rk_stage < rk_num_stages; rk_stage++) {
            // ---- RK coefficient ----
            double rk_alpha = 1.0 / ((double)rk_num_stages - (double)rk_stage);

            // ---- Calculate velocity diveregence for the element ----
            if (num_dim == 2) {
                get_divergence2D(elem_div,
                                 node_coords,
                                 node_vel,
                                 elem_vol);
            }
            else{
                get_divergence(elem_div,
                               node_coords,
                               node_vel,
                               elem_vol);
            } // end if 2D

            // ---- calculate the forces on the vertices and evolve stress (hypo model) ----
            if (num_dim == 2) {
                get_force_sgh2D(material,
                                *mesh,
                                node_coords,
                                node_vel,
                                elem_den,
                                elem_sie,
                                elem_pres,
                                elem_stress,
                                elem_sspd,
                                elem_vol,
                                elem_div,
                                elem_mat_id,
                                corner_force,
                                rk_alpha,
                                cycle);
            }
            else{
                get_force_sgh(material,
                              *mesh,
                              node_coords,
                              node_vel,
                              elem_den,
                              elem_sie,
                              elem_pres,
                              elem_stress,
                              elem_sspd,
                              elem_vol,
                              elem_div,
                              elem_mat_id,
                              corner_force,
                              rk_alpha,
                              cycle);
            }

            if (have_loading_conditions) {
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
            }

#ifdef DEBUG
            if (myrank == 1) {
                std::cout << "rk_alpha = " << rk_alpha << ", dt = " << dt << std::endl;
                for (int i = 0; i < nall_nodes; i++) {
                    double node_force[3];
                    for (size_t dim = 0; dim < num_dim; dim++) {
                        node_force[dim] = 0.0;
                    } // end for dim

                    // loop over all corners around the node and calculate the nodal force
                    for (size_t corner_lid = 0; corner_lid < mesh.num_corners_in_node(i); corner_lid++) {
                        // Get corner gid
                        size_t corner_gid = mesh.corners_in_node(i, corner_lid);
                        std::cout << Explicit_Solver_Pointer_->all_node_map->getGlobalElement(i) << " " << corner_gid << " " << corner_force(corner_gid, 0) << " " << corner_force(corner_gid,
                        1) << " " << corner_force(corner_gid, 2) << std::endl;
                        // loop over dimension
                        for (size_t dim = 0; dim < num_dim; dim++) {
                            node_force[dim] += corner_force(corner_gid, dim);
                        } // end for dim
                    } // end for corner_lid
                }
            }

            // debug print vector values on a rank

            if (myrank == 0) {
                for (int i = 0; i < nall_nodes; i++) {
                    std::cout << Explicit_Solver_Pointer_->all_node_map->getGlobalElement(i) << " " << node_vel(rk_level, i, 0) << " " << node_vel(rk_level, i, 1) << " " << node_vel(rk_level, i,
                    2) << std::endl;
                }
            }
#endif

            // ---- Update nodal velocities ---- //
            update_velocity_sgh(rk_alpha,
                              node_vel,
                              node_mass,
                              corner_force);


            // ---- apply force boundary conditions to the boundary patches----
            boundary_velocity(*mesh, boundary, node_vel);

            // current interface has differing velocity arrays; this equates them until we unify memory
            // first comm time interval point
            double comm_time1 = Explicit_Solver_Pointer_->CPU_Time();
            
            // if(module_params->matar_mpi_test){
            //     FOR_ALL_CLASS(node_gid, 0, nlocal_nodes, {
            //         for (int idim = 0; idim < num_dim; idim++) {
            //             mtr_node_velocities_distributed(node_gid, idim) = node_vel(rk_level, node_gid, idim);
            //         }
            //     }); // end parallel for
            //     Kokkos::fence();
            //     mtr_node_velocities_distributed.update_host();
            // }

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
            // if(module_params->matar_mpi_test){
            //     mtr_node_velocities_distributed.perform_comms();
            //     mtr_node_velocities_distributed.update_device();
            // }

            Explicit_Solver_Pointer_->comm_velocities();

            double comm_time3 = Explicit_Solver_Pointer_->CPU_Time();
            // if(module_params->matar_mpi_test){
            //     FOR_ALL_CLASS(node_gid, nlocal_nodes, nall_nodes, {
            //         for (int idim = 0; idim < num_dim; idim++) {
            //             node_vel(rk_level, node_gid, idim) = mtr_node_velocities_distributed(node_gid, idim);
            //         }
            //     }); // end parallel for
            //     Kokkos::fence();
            // }

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

#ifdef DEBUG
            // debug print vector values on a rank
            if (myrank == 0) {
                for (int i = 0; i < nall_nodes; i++) {
                    std::cout << Explicit_Solver_Pointer_->all_node_map->getGlobalElement(i) << " " << node_vel(rk_level, i, 0) << " " << node_vel(rk_level, i, 1) << " " << node_vel(rk_level, i,
                    2) << std::endl;
                }
            }
#endif
            // ---- Update specific internal energy in the elements ----
            update_energy_sgh(rk_alpha,
                              *mesh,
                              node_vel,
                              node_coords,
                              elem_sie,
                              elem_mass,
                              corner_force);

            // ---- Update nodal positions ----
            update_position_sgh(rk_alpha,
                                nall_nodes,
                                node_coords,
                                node_vel);

            // ---- Calculate cell volume for next time step ----
            get_vol();

            // ---- Calculate elem state (den, pres, sound speed, stress) for next time step ----
            if (num_dim == 2) {
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

        if (topology_optimization_on || shape_optimization_on) {
            if (cycle >= max_time_steps) {
                max_time_steps = cycle + 1;
            }
            
            if(!use_solve_checkpoints){
                if (max_time_steps + 1 > forward_solve_velocity_data->size()) {
                    old_max_forward_buffer = forward_solve_velocity_data->size();
                    time_data.resize(max_time_steps + BUFFER_GROW + 1);
                    forward_solve_velocity_data->resize(max_time_steps + BUFFER_GROW + 1);
                    forward_solve_coordinate_data->resize(max_time_steps + BUFFER_GROW + 1);
                    forward_solve_internal_energy_data->resize(max_time_steps + BUFFER_GROW + 1);
                    adjoint_vector_data->resize(max_time_steps + BUFFER_GROW + 1);
                    phi_adjoint_vector_data->resize(max_time_steps + BUFFER_GROW + 1);
                    psi_adjoint_vector_data->resize(max_time_steps + BUFFER_GROW + 1);
                    // assign a multivector of corresponding size to each new timestep in the buffer
                    for (int istep = old_max_forward_buffer; istep < max_time_steps + BUFFER_GROW + 1; istep++) {
                        (*forward_solve_velocity_data)[istep]   = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
                        (*forward_solve_coordinate_data)[istep] = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
                        (*forward_solve_internal_energy_data)[istep] = Teuchos::rcp(new MV(all_element_map, 1));
                        (*adjoint_vector_data)[istep]     = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
                        (*phi_adjoint_vector_data)[istep] = Teuchos::rcp(new MV(all_node_map, simparam->num_dims));
                        (*psi_adjoint_vector_data)[istep] = Teuchos::rcp(new MV(all_element_map, 1));
                    }
                }

                time_data[cycle + 1] = dt + time_data[cycle];
            }
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
                vec_array all_node_velocities_interface = Explicit_Solver_Pointer_->all_node_velocities_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
                vec_array element_internal_energy     = element_internal_energy_distributed->getLocalView<device_type>(Tpetra::Access::ReadWrite);
                const_vec_array node_coords_interface = Explicit_Solver_Pointer_->node_coords_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
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

                // interface for element internal energies
                FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
                    element_internal_energy(elem_gid, 0) = elem_sie(rk_level, elem_gid);
                }); // end parallel for
                Kokkos::fence();
            } // end view scope

            double comm_time4 = Explicit_Solver_Pointer_->CPU_Time();
            Explicit_Solver_Pointer_->host2dev_time += comm_time4 - comm_time3;
            Explicit_Solver_Pointer_->communication_time += comm_time4 - comm_time1;

            if(use_solve_checkpoints){
                //add level 0 checkpoints sequentially until requested total limit is reached
                if(num_active_checkpoints < num_solve_checkpoints){
                    Dynamic_Checkpoint temp(3,cycle+1,time_value, dt);
                    temp.change_vector(U_DATA, (*forward_solve_coordinate_data)[num_active_checkpoints + 1]);
                    temp.change_vector(V_DATA, (*forward_solve_velocity_data)[num_active_checkpoints + 1]);
                    temp.change_vector(SIE_DATA, (*forward_solve_internal_energy_data)[num_active_checkpoints + 1]);
                    temp.assign_vector(U_DATA,all_node_coords_distributed);
                    temp.assign_vector(V_DATA,all_node_velocities_distributed);
                    temp.assign_vector(SIE_DATA,element_internal_energy_distributed);
                    dynamic_checkpoint_set->insert(temp);
                    num_active_checkpoints++;
                    //initializes to the end of the set until allowed total number of checkpoints is reached
                    last_raised_checkpoint = dynamic_checkpoint_set->end();
                    --last_raised_checkpoint;
                }
                //if limit of checkpoints was reached; search for dispensable checkpoints or raise level of recent checkpoint
                else{
                    //a dispensable checkpoint has a lower level than another checkpoint located later in time
                    //find if there is a dispenable checkpoint to remove
                    current_checkpoint = last_raised_checkpoint;
                    --current_checkpoint; // dont need to check against itself
                    search_end = dynamic_checkpoint_set->begin();
                    dispensable_found = false;
                    while(current_checkpoint!=search_end){
                        if(current_checkpoint->level<last_raised_level){
                            dispensable_checkpoint = current_checkpoint;
                            dispensable_found = true;
                            break;
                        }
                        --current_checkpoint;
                    }
                    if(dispensable_found){
                        //add replacement checkpoint
                        Dynamic_Checkpoint temp(3,cycle+1,time_value, dt);
                        //get pointers to vector buffers from the checkpoint we're about to delete
                        temp.copy_vectors(*dispensable_checkpoint);
                        //remove checkpoint at timestep = cycle
                        dynamic_checkpoint_set->erase(dispensable_checkpoint);
                        //assign current phase data to vector buffers
                        temp.assign_vector(U_DATA,all_node_coords_distributed);
                        temp.assign_vector(V_DATA,all_node_velocities_distributed);
                        temp.assign_vector(SIE_DATA,element_internal_energy_distributed);
                        dynamic_checkpoint_set->insert(temp);

                    }
                    else{
                        //since no dispensable checkpoints were found, raise level of new one to be one higher than previous checkpoint
                        current_checkpoint = dynamic_checkpoint_set->end();
                        --current_checkpoint; //reduce iterator by 1 so it doesnt point to the sentinel past the last element
                        last_raised_level = current_checkpoint->level;

                        //add replacement checkpoint
                        Dynamic_Checkpoint temp(3,cycle+1,time_value, dt, ++last_raised_level);
                        //get pointers to vector buffers from the checkpoint we're about to delete
                        temp.copy_vectors(*current_checkpoint);
                        //remove checkpoint at timestep = cycle
                        dynamic_checkpoint_set->erase(current_checkpoint);
                        //assign current phase data to vector buffers
                        temp.assign_vector(U_DATA,all_node_coords_distributed);
                        temp.assign_vector(V_DATA,all_node_velocities_distributed);
                        temp.assign_vector(SIE_DATA,element_internal_energy_distributed);
                        dynamic_checkpoint_set->insert(temp);
                        last_raised_checkpoint = dynamic_checkpoint_set->end();
                        --last_raised_checkpoint; //save iterator for this checkpoint to expedite dispensable search
                    }
                }
            }
            else{
                (*forward_solve_internal_energy_data)[cycle + 1]->assign(*element_internal_energy_distributed);
                (*forward_solve_velocity_data)[cycle + 1]->assign(*Explicit_Solver_Pointer_->all_node_velocities_distributed);
                (*forward_solve_coordinate_data)[cycle + 1]->assign(*Explicit_Solver_Pointer_->all_node_coords_distributed);
            }

            // kinetic energy accumulation
            if (time_accumulation) {
                objective_function->step_accumulation(dt, cycle, rk_level);
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

            if (myrank == 0) {
                printf("Writing outputs to file at %f \n", graphics_time);
            }

            double comm_time1 = Explicit_Solver_Pointer_->CPU_Time();
            Explicit_Solver_Pointer_->write_outputs();

            double comm_time2 = Explicit_Solver_Pointer_->CPU_Time();
            Explicit_Solver_Pointer_->output_time += comm_time2 - comm_time1;

            graphics_time = time_value + graphics_dt_ival;
            dt = cached_pregraphics_dt;
        } // end if

        // end of calculation
        if (time_value >= time_final) {
            break;
        }
    } // end for cycle loop

    last_time_step = cycle;

    //debug print of checkpoint timesteps
    // if(use_solve_checkpoints){
    //     int count = 0;
    //     for(auto it = dynamic_checkpoint_set->begin(); it != dynamic_checkpoint_set->end(); it++){
    //         *fos << "Checkpoint # " << count++ << " is at timestep " << (*it).saved_timestep << " with timestep "
    //          << (*it).saved_dt << " at time " << (*it).saved_time << std::endl;
    //     }
    // }

    // simple setup to just calculate KE minimize objective for now
    if (topology_optimization_on) {
        objective_function->global_reduction();
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
    if (simparam->dynamic_options.output_time_sequence_level >= TIME_OUTPUT_LEVEL::low) {
        if (myrank == 0) {
            printf("Time=0:   KE = %20.15f, IE = %20.15f, TE = %20.15f \n", KE_t0, IE_t0, TE_t0);
        }
        if (myrank == 0) {
            printf("Time=End: KE = %20.15f, IE = %20.15f, TE = %20.15f \n", KE_tend, IE_tend, TE_tend);
        }
        if (myrank == 0) {
            printf("total energy conservation error = %e \n\n", 100 * (TE_tend - TE_t0) / TE_t0);
        }
    }

    return;
} // end of SGH solve
