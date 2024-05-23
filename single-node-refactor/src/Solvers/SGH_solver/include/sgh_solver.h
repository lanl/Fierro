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

#ifndef SGH_SOLVER_H
#define SGH_SOLVER_H

#include "matar.h"
#include "solver.h"
#include "geometry_new.h"
#include "contact.h"
// #include "io_utils.h"

#include "simulation_parameters.h"

using namespace mtr; // matar namespace

/////////////////////////////////////////////////////////////////////////////
///
/// \class SGH
///
/// \brief Class for containing functions required to perform SGH
///
/// This class contains the requisite functions requited to perform
/// staggered grid hydrodynamics (SGH) which is equivalent to a lumped
/// mass finite element (FE) scheme.
///
/////////////////////////////////////////////////////////////////////////////
class SGH : public Solver
{
public:

    double dt = 0.0;
    double time_value = 0.0;

    double time_initial = 0.0;  // Starting time
    double time_final   = 1.0;  // Final simulation time
    double dt_min   = 1e-8;     // Minimum time step
    double dt_max   = 1e-2;     // Maximum time step
    double dt_start = 1e-5;     // Starting time step
    double dt_cfl   = 0.4;      // CFL multiplier for time step calculation

    double graphics_dt_ival  = 1.0; // time increment for graphics output
    int    graphics_cyc_ival = 2000000; // Cycle count for graphics output

    int rk_num_stages = 2;
    int cycle_stop    = 1000000000;

    contact_patches_t contact_bank;  // keeps track of contact patches
    bool doing_contact = false;  // Condition used in SGH::execute

    SGH()  : Solver()
    {
    }

    ~SGH() = default;

    // Initialize data specific to the SGH solver
    void initialize(simulation_parameters_t& sim_param) override
    {
        // Dimensions
        num_dims = 3;

        graphics_times = CArray<double>(20000);

        // NOTE: Possible remove this and pass directly
        fuzz  = sim_param.dynamic_options.fuzz;
        tiny  = sim_param.dynamic_options.tiny;
        small = sim_param.dynamic_options.small;

        time_initial = sim_param.dynamic_options.time_initial;
        time_final   = sim_param.dynamic_options.time_final;
        dt_min   = sim_param.dynamic_options.dt_min;
        dt_max   = sim_param.dynamic_options.dt_max;
        dt_start = sim_param.dynamic_options.dt_start;
        dt_cfl   = sim_param.dynamic_options.dt_cfl;

        graphics_dt_ival  = sim_param.output_options.graphics_time_step;
        graphics_cyc_ival = sim_param.output_options.graphics_iteration_step;

        rk_num_stages = sim_param.dynamic_options.rk_num_stages;

        cycle_stop = sim_param.dynamic_options.cycle_stop;

        // initialize time, time_step, and cycles
        time_value = 0.0;
        dt = sim_param.dynamic_options.dt_start;

        graphics_id = 0;
        graphics_times(0) = sim_param.output_options.graphics_time_step;
        graphics_time     = sim_param.output_options.graphics_time_step;
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn setup
    ///
    /// \brief Calls setup_sgh, which initializes state, and material data
    ///
    /////////////////////////////////////////////////////////////////////////////
    void setup(simulation_parameters_t& sim_param, mesh_t& mesh, node_t& node, elem_t& elem, corner_t& corner) override
    {
        std::cout << "INSIDE SETUP FOR SGH SOLVER" << std::endl;

        std::cout << "Applying initial boundary conditions" << std::endl;
        boundary_velocity(mesh, sim_param.boundary_conditions, node.vel, time_value);

        // Setting up contact
        for (size_t i = 0; i < mesh.num_bdy_sets; i++) {
            boundary_condition_t bound = sim_param.boundary_conditions(i);
            if (bound.type == boundary_conds::contact && bound.geometry == boundary_conds::global) {
                std::cout << "Setting up global contact" << std::endl;
                doing_contact = true;

                contact_bank.initialize(mesh, mesh.bdy_patches, node, corner);
                // run_contact_tests(contact_bank, mesh, node, corner, sim_param);
                break;
            } else if (bound.type == boundary_conds:: contact && bound.geometry != boundary_conds::global) {
                doing_contact = true;
                std::cerr << "Contact boundary conditions are only supported for global at the moment." << std::endl;
                exit(1);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn execute
    ///
    /// \brief Calls the solve function which evolves the state
    ///
    ///
    /////////////////////////////////////////////////////////////////////////////
    void execute(simulation_parameters_t& sim_param, mesh_t& mesh, node_t& node, elem_t& elem, corner_t& corner) override;


    void finalize(simulation_parameters_t& sim_param) override
    {
        // Any finalize goes here, remove allocated memory, etc
    }

    // **** Functions defined in boundary.cpp **** //
    void boundary_velocity(
        const mesh_t& mesh,
        const CArrayKokkos<boundary_condition_t>& boundary,
        DCArrayKokkos<double>& node_vel,
        const double time_value);

    void boundary_contact(const double &del_t);

    // **** Functions defined in energy_sgh.cpp **** //
    void update_energy(
        double rk_alpha,
        double dt,
        const mesh_t& mesh,
        const DCArrayKokkos<double>& node_vel,
        const DCArrayKokkos<double>& node_coords,
        DCArrayKokkos<double>& elem_sie,
        const DCArrayKokkos<double>& elem_mass,
        const DCArrayKokkos<double>& corner_force);

    // **** Functions defined in force_sgh.cpp **** //
    void get_force(
        const CArrayKokkos<material_t>& material,
        const mesh_t& mesh,
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& node_vel,
        const DCArrayKokkos<double>& elem_den,
        const DCArrayKokkos<double>& elem_sie,
        const DCArrayKokkos<double>& elem_pres,
        const DCArrayKokkos<double>& elem_stress,
        const DCArrayKokkos<double>& elem_sspd,
        const DCArrayKokkos<double>& elem_vol,
        const DCArrayKokkos<double>& elem_div,
        const DCArrayKokkos<size_t>& elem_mat_id,
        DCArrayKokkos<double>& corner_force,
        const double fuzz,
        const double small,
        const DCArrayKokkos<double>& elem_statev,
        const double dt,
        const double rk_alpha);

    void get_force_2D(
        const CArrayKokkos<material_t>& material,
        const mesh_t& mesh,
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& node_vel,
        const DCArrayKokkos<double>& elem_den,
        const DCArrayKokkos<double>& elem_sie,
        const DCArrayKokkos<double>& elem_pres,
        const DCArrayKokkos<double>& elem_stress,
        const DCArrayKokkos<double>& elem_sspd,
        const DCArrayKokkos<double>& elem_vol,
        const DCArrayKokkos<double>& elem_div,
        const DCArrayKokkos<size_t>& elem_mat_id,
        DCArrayKokkos<double>& corner_force,
        const double fuzz,
        const double small,
        const DCArrayKokkos<double>& elem_statev,
        const double dt,
        const double rk_alpha);

    // **** Functions defined in geometry.cpp **** //
    void update_position(
        double rk_alpha,
        double dt,
        const size_t num_dims,
        const size_t num_nodes,
        DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& node_vel);

    // **** Functions defined in momentum.cpp **** //
    void update_velocity(
        double rk_alpha,
        double dt,
        const mesh_t& mesh,
        DCArrayKokkos<double>& node_vel,
        const DCArrayKokkos<double>& node_mass,
        const DCArrayKokkos<double>& corner_force,
        const CArrayKokkos<contact_node_t> &contact_nodes);

    KOKKOS_FUNCTION
    void get_velgrad(
        ViewCArrayKokkos<double>& vel_grad,
        const ViewCArrayKokkos<size_t>& elem_node_gids,
        const DCArrayKokkos<double>&    node_vel,
        const ViewCArrayKokkos<double>& b_matrix,
        const double elem_vol,
        const size_t elem_gid);

    KOKKOS_FUNCTION
    void get_velgrad2D(
        ViewCArrayKokkos<double>& vel_grad,
        const ViewCArrayKokkos<size_t>& elem_node_gids,
        const DCArrayKokkos<double>&    node_vel,
        const ViewCArrayKokkos<double>& b_matrix,
        const double elem_vol,
        const double elem_area,
        const size_t elem_gid);

    void get_divergence(
        DCArrayKokkos<double>& elem_div,
        const mesh_t mesh,
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& node_vel,
        const DCArrayKokkos<double>& elem_vol);

    void get_divergence2D(
        DCArrayKokkos<double>& elem_div,
        const mesh_t mesh,
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& node_vel,
        const DCArrayKokkos<double>& elem_vol);

    KOKKOS_FUNCTION
    void decompose_vel_grad(
        ViewCArrayKokkos<double>& D_tensor,
        ViewCArrayKokkos<double>& W_tensor,
        const ViewCArrayKokkos<double>& vel_grad,
        const ViewCArrayKokkos<size_t>& elem_node_gids,
        const size_t elem_gid,
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& node_vel,
        const double vol);

    // **** Functions defined in properties.cpp **** //
    void update_state(
        const CArrayKokkos<material_t>& material,
        const mesh_t& mesh,
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& node_vel,
        DCArrayKokkos<double>& elem_den,
        DCArrayKokkos<double>& elem_pres,
        DCArrayKokkos<double>& elem_stress,
        DCArrayKokkos<double>& elem_sspd,
        const DCArrayKokkos<double>& elem_sie,
        const DCArrayKokkos<double>& elem_vol,
        const DCArrayKokkos<double>& elem_mass,
        const DCArrayKokkos<size_t>& elem_mat_id,
        const DCArrayKokkos<double>& elem_statev,
        const double dt,
        const double rk_alpha);

    void update_state2D(
        const CArrayKokkos<material_t>& material,
        const mesh_t& mesh,
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& node_vel,
        DCArrayKokkos<double>& elem_den,
        DCArrayKokkos<double>& elem_pres,
        DCArrayKokkos<double>& elem_stress,
        DCArrayKokkos<double>& elem_sspd,
        const DCArrayKokkos<double>& elem_sie,
        const DCArrayKokkos<double>& elem_vol,
        const DCArrayKokkos<double>& elem_mass,
        const DCArrayKokkos<size_t>& elem_mat_id,
        const DCArrayKokkos<double>& elem_statev,
        const double dt,
        const double rk_alpha);

    // **** Functions defined in time_integration.cpp **** //
    // NOTE: Consider pulling up
    void rk_init(
        DCArrayKokkos<double>& node_coords,
        DCArrayKokkos<double>& node_vel,
        DCArrayKokkos<double>& elem_sie,
        DCArrayKokkos<double>& elem_stress,
        const size_t num_dims,
        const size_t num_elems,
        const size_t num_nodes);

    void get_timestep(
        mesh_t& mesh,
        DCArrayKokkos<double>& node_coords,
        DCArrayKokkos<double>& node_vel,
        DCArrayKokkos<double>& elem_sspd,
        DCArrayKokkos<double>& elem_vol,
        double time_value,
        const double graphics_time,
        const double time_final,
        const double dt_max,
        const double dt_min,
        const double dt_cfl,
        double&      dt,
        const double fuzz);

    void get_timestep2D(
        mesh_t& mesh,
        DCArrayKokkos<double>& node_coords,
        DCArrayKokkos<double>& node_vel,
        DCArrayKokkos<double>& elem_sspd,
        DCArrayKokkos<double>& elem_vol,
        double time_value,
        const double graphics_time,
        const double time_final,
        const double dt_max,
        const double dt_min,
        const double dt_cfl,
        double&      dt,
        const double fuzz);

};

#endif // end HEADER_H
