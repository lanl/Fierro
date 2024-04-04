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
// #include "mesh.h"
// #include "state.h"
// #include "material.h"
#include "geometry_new.h"
#include "io_utils.h"

#include "simulation_parameters.h"

using namespace mtr; // matar namespace

/////////////////////////////////////////////////////////////////////////////
///
/// \class SGH
///
/// \brief Class for containing functions required to perform SGH
///
/// This class containts the requisite functions requited to perform
/// staggered grid hydrodynamics (SGH) which is equivalent to a lumped
/// mass finite element (FE) scheme.
///
/////////////////////////////////////////////////////////////////////////////
class SGH : public Solver
{
public:

    char* mesh_file;

    MeshReader* reader;

    simulation_parameters_t sim_param;
    mesh_t mesh;
    node_t node;
    elem_t elem;
    corner_t corner;


    SGH(simulation_parameters_t& io, mesh_t &mesh_in, node_t &node_in, elem_t &elem_in, corner_t &corner_in)  : Solver()
    {
        sim_param = io;
        mesh = mesh_in;

        node = node_in;
        elem = elem_in;
        corner = corner_in;
    }

    ~SGH() = default;

    // Initialize data for the SGH solver
    // This will be where we take in parsed data from YAML
    void initialize()
    {
        std::cout << "In initialize function in sgh solver" << std::endl;

        // Dimensions
        num_dims = 3;

        graphics_times    = CArray<double>(20000);

        // // ---- time varaibles and cycle info ----
        // time_final = 1.0;  // 1.0 for Sedov
        // dt_min     = 1.e-8;
        // dt_max     = 1.e-2;
        // dt_start   = 1.e-5;
        // cycle_stop = 100000;

        // // ---- graphics information ----
        // graphics_times    = CArray<double>(20000);
        // graphics_cyc_ival = 1000000;
        // graphics_dt_ival  = 0.25;

        // // --- number of material regions ---
        // num_materials = 1;
        // material = CArrayKokkos<material_t>(num_materials);      // create material

        // // --- declare model state variable array size ---
        // state_vars = CArrayKokkos<double>(num_materials, max_num_state_vars); // init values

        // // --- number of fill regions ---
        // num_fills = 2;  // =2 for Sedov
        // region_fill  = CArrayKokkos<reg_fill_t>(num_fills);  // create fills

        // // --- number of boundary conditions ---
        // num_bcs  = 6; // =6 for Sedov
        // boundary = CArrayKokkos<boundary_condition_t>(num_bcs);  // create boundaries

        // RUN({
        //     // gamma law model
        //     // statev(0) = gamma
        //     // statev(1) = minimum sound speed
        //     // statev(2) = specific heat
        //     // statev(3) = ref temperature
        //     // statev(4) = ref density
        //     // statev(5) = ref specific internal energy

        //     material(0).eos_model = ideal_gas; // EOS model is required

        //     material(0).strength_type  = model::none;
        //     material(0).strength_setup = model_init::input; // not need, the input is the default
        //     material(0).strength_model = NULL;  // not needed, but illistrates the syntax

        //     material(0).q1   = 1.0;            // accoustic coefficient
        //     material(0).q2   = 1.3333;         // linear slope of UsUp for Riemann solver
        //     material(0).q1ex = 1.0;            // accoustic coefficient in expansion
        //     material(0).q2ex = 0.0;            // linear slope of UsUp in expansion

        //     material(0).num_state_vars = 3;  // actual num_state_vars
        //     state_vars(0, 0) = 5.0 / 3.0; // gamma value
        //     state_vars(0, 1) = 1.0E-14; // minimum sound speed
        //     state_vars(0, 2) = 1.0;     // specific heat

        //     // global initial conditions
        //     region_fill(0).volume = region::global; // fill everywhere
        //     region_fill(0).mat_id = 0;              // material id
        //     region_fill(0).den    = 1.0;            // intial density
        //     region_fill(0).sie    = 1.e-10;         // intial specific internal energy

        //     region_fill(0).velocity = init_conds::cartesian;
        //     region_fill(0).u = 0.0;   // initial x-dir velocity
        //     region_fill(0).v = 0.0;   // initial y-dir velocity
        //     region_fill(0).w = 0.0;   // initial z-dir velocity

        //     // energy source initial conditions
        //     region_fill(1).volume  = region::sphere; // fill a sphere
        //     region_fill(1).mat_id  = 0;             // material id
        //     region_fill(1).radius1 = 0.0;           // inner radius of fill region
        //     region_fill(1).radius2 = 1.2 / 8.0;       // outer radius of fill region
        //     region_fill(1).den     = 1.0;           // initial density
        //     region_fill(1).sie     = (963.652344 *
        //                            pow((1.2 / 30.0), 3)) / pow((region_fill(1).radius2), 3);

        //     region_fill(1).velocity = init_conds::cartesian;
        //     region_fill(1).u = 0.0;   // initial x-dir velocity
        //     region_fill(1).v = 0.0;   // initial y-dir velocity
        //     region_fill(1).w = 0.0;   // initial z-dir velocity

        //     // ---- boundary conditions ---- //

        //     // Tag X plane
        //     boundary(0).surface  = bdy::x_plane; // planes, cylinder, spheres, or a files
        //     boundary(0).value    = 0.0;
        //     boundary(0).hydro_bc = bdy::reflected;

        //     // Tag Y plane
        //     boundary(1).surface  = bdy::y_plane;
        //     boundary(1).value    = 0.0;
        //     boundary(1).hydro_bc = bdy::reflected;

        //     // Tag Z plane
        //     boundary(2).surface  = bdy::z_plane;
        //     boundary(2).value    = 0.0;
        //     boundary(2).hydro_bc = bdy::reflected;

        //     // Tag X plane
        //     boundary(3).surface  = bdy::x_plane; // planes, cylinder, spheres, or a files
        //     boundary(3).value    = 1.2;
        //     boundary(3).hydro_bc = bdy::reflected;

        //     // Tag Y plane
        //     boundary(4).surface  = bdy::y_plane;
        //     boundary(4).value    = 1.2;
        //     boundary(4).hydro_bc = bdy::reflected;

        //     // Tag Z plane
        //     boundary(5).surface  = bdy::z_plane;
        //     boundary(5).value    = 1.2;
        //     boundary(5).hydro_bc = bdy::reflected;
        // });  // end RUN



        // ---------------------------------------------------------------------
        //    allocate memory
        // ---------------------------------------------------------------------

        std::cout << "Num elements = " << mesh.num_elems << std::endl;
        std::cout << "Num nodes = " << mesh.num_nodes << std::endl;


        // ---------------------------------------------------------------------
        //   calculate geometry
        // ---------------------------------------------------------------------

        std::cout << "Before update device call"  << std::endl;
        node.coords.update_device();
        std::cout << "After update device call"  << std::endl;

        Kokkos::fence();

        std::cout << "Before get volume call"  << std::endl;
        geometry::get_vol(elem.vol, node.coords, mesh);


        std::cout << "After get volume call"  << std::endl;
        // intialize time, time_step, and cycles
        time_value = 0.0;
        dt = 0.005;
        graphics_id = 0;
        graphics_times(0) = 0.0;
        graphics_time     = graphics_dt_ival; // the times for writing graphics dump
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn setup
    ///
    /// \brief Calls setup_sgh, which initializes mesh, state, and material data
    ///
    /////////////////////////////////////////////////////////////////////////////
    void setup()
    {
        std::cout<<"INSIDE SETUP FOR SGH SOLVER"<<std::endl;


        printf("Writing outputs to file at %f \n", time_value);
        write_outputs(mesh,
                      node.coords,
                      node.vel,
                      node.mass,
                      elem.den,
                      elem.pres,
                      elem.stress,
                      elem.sspd,
                      elem.sie,
                      elem.vol,
                      elem.mass,
                      elem.mat_id,
                      graphics_times,
                      graphics_id,
                      time_value);

        // setup_sgh(
        //     sim_param.materials,
        //     sim_param.region_fills,
        //     sim_param.boundary_conditions,
        //     mesh,
        //     node.coords,
        //     node.vel,
        //     node.mass,
        //     elem.den,
        //     elem.pres,
        //     elem.stress,
        //     elem.sspd,
        //     elem.sie,
        //     elem.vol,
        //     elem.mass,
        //     elem.mat_id,
        //     elem.statev,
        //     state_vars,
        //     corner.mass,
        //     num_fills,
        //     rk_num_bins,
        //     num_bcs,
        //     num_materials,
        //     2); // num_state_vars
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn run
    ///
    /// \brief Calls the solve function which evolves the state
    ///
    ///
    /////////////////////////////////////////////////////////////////////////////
    void execute()
    {
        std::cout << "In execute function in sgh solver" << std::endl;

        std::cout << "Num elements = " << mesh.num_elems << std::endl;
        std::cout << "Num nodes = " << mesh.num_nodes << std::endl;

        std::cout << "Elem 0 pressure =  = " << elem.pres(0) << std::endl;
        std::cout << "Elem 0 sie rk0  = " << elem.sie(0,0) << std::endl;
        std::cout << "Elem 0 sie rk1  = " << elem.sie(1,0) << std::endl;

        solve(sim_param.materials,
              sim_param.boundary_conditions,
              mesh,
              node.coords,
              node.vel,
              node.mass,
              elem.den,
              elem.pres,
              elem.stress,
              elem.sspd,
              elem.sie,
              elem.vol,
              elem.div,
              elem.mass,
              elem.mat_id,
              elem.statev,
              corner.force,
              corner.mass,
              time_value,
              time_final,
              dt_max,
              dt_min,
              dt_cfl,
              graphics_time,
              graphics_cyc_ival,
              graphics_dt_ival,
              5, //sim_param.dynamic_options.cycle_stop,
              rk_num_stages,
              dt,
              fuzz,
              tiny,
              small,
              graphics_times,
              graphics_id);
    }

    void setup_sgh(
        const CArrayKokkos<material_t>& material,
        const CArrayKokkos<reg_fill_t>& region_fill,
        const CArrayKokkos<boundary_condition_t>& boundary,
        mesh_t& mesh,
        const DCArrayKokkos<double>& node_coords,
        DCArrayKokkos<double>& node_vel,
        DCArrayKokkos<double>& node_mass,
        const DCArrayKokkos<double>& elem_den,
        const DCArrayKokkos<double>& elem_pres,
        const DCArrayKokkos<double>& elem_stress,
        const DCArrayKokkos<double>& elem_sspd,
        const DCArrayKokkos<double>& elem_sie,
        const DCArrayKokkos<double>& elem_vol,
        const DCArrayKokkos<double>& elem_mass,
        const DCArrayKokkos<size_t>& elem_mat_id,
        const DCArrayKokkos<double>& elem_statev,
        const CArrayKokkos<double>&      state_vars,
        const DCArrayKokkos<double>& corner_mass,
        const size_t num_fills,
        const size_t rk_num_bins,
        const size_t num_bcs,
        const size_t num_materials,
        const size_t num_state_vars);

    void write_outputs(
        const mesh_t& mesh,
        DCArrayKokkos<double>& node_coords,
        DCArrayKokkos<double>& node_vel,
        DCArrayKokkos<double>& node_mass,
        DCArrayKokkos<double>& elem_den,
        DCArrayKokkos<double>& elem_pres,
        DCArrayKokkos<double>& elem_stress,
        DCArrayKokkos<double>& elem_sspd,
        DCArrayKokkos<double>& elem_sie,
        DCArrayKokkos<double>& elem_vol,
        DCArrayKokkos<double>& elem_mass,
        DCArrayKokkos<size_t>& elem_mat_id,
        CArray<double>& graphics_times,
        size_t&      graphics_id,
        const double time_value);

    void solve(CArrayKokkos<material_t>& material,
               CArrayKokkos<boundary_condition_t>& boundary,
               mesh_t& mesh,
               DCArrayKokkos<double>& node_coords,
               DCArrayKokkos<double>& node_vel,
               DCArrayKokkos<double>& node_mass,
               DCArrayKokkos<double>& elem_den,
               DCArrayKokkos<double>& elem_pres,
               DCArrayKokkos<double>& elem_stress,
               DCArrayKokkos<double>& elem_sspd,
               DCArrayKokkos<double>& elem_sie,
               DCArrayKokkos<double>& elem_vol,
               DCArrayKokkos<double>& elem_div,
               DCArrayKokkos<double>& elem_mass,
               DCArrayKokkos<size_t>& elem_mat_id,
               DCArrayKokkos<double>& elem_statev,
               DCArrayKokkos<double>& corner_force,
               DCArrayKokkos<double>& corner_mass,
               double&      time_value,
               const double time_final,
               const double dt_max,
               const double dt_min,
               const double dt_cfl,
               double&      graphics_time,
               size_t graphics_cyc_ival,
               double graphics_dt_ival,
               const size_t cycle_stop,
               const size_t rk_num_stages,
               double dt,
               const double    fuzz,
               const double    tiny,
               const double    small,
               CArray<double>& graphics_times,
               size_t& graphics_id);

    // **** Functions defined in boundary.cpp **** //
    void boundary_velocity(
        const mesh_t& mesh,
        const CArrayKokkos<boundary_condition_t>& boundary,
        DCArrayKokkos<double>&      node_vel,
        const double time_value);

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
        const DCArrayKokkos<double>& corner_force);

    KOKKOS_FUNCTION
    void get_velgrad(
        ViewCArrayKokkos<double>& vel_grad,
        const ViewCArrayKokkos<size_t>&  elem_node_gids,
        const DCArrayKokkos<double>& node_vel,
        const ViewCArrayKokkos<double>&  b_matrix,
        const double elem_vol,
        const size_t elem_gid);

    KOKKOS_FUNCTION
    void get_velgrad2D(
        ViewCArrayKokkos<double>& vel_grad,
        const ViewCArrayKokkos<size_t>&  elem_node_gids,
        const DCArrayKokkos<double>& node_vel,
        const ViewCArrayKokkos<double>&  b_matrix,
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

    // **** Functions defined in user_mat.cpp **** //
    // NOTE: Pull up into high level
    KOKKOS_FUNCTION
    void user_eos_model(
        const DCArrayKokkos<double>& elem_pres,
        const DCArrayKokkos<double>& elem_stress,
        const size_t elem_gid,
        const size_t mat_id,
        const DCArrayKokkos<double>& elem_state_vars,
        const DCArrayKokkos<double>& elem_sspd,
        const double den,
        const double sie);

    KOKKOS_FUNCTION
    void user_strength_model(
        const DCArrayKokkos<double>& elem_pres,
        const DCArrayKokkos<double>& elem_stress,
        const size_t elem_gid,
        const size_t mat_id,
        const DCArrayKokkos<double>& elem_state_vars,
        const DCArrayKokkos<double>& elem_sspd,
        const double den,
        const double sie,
        const ViewCArrayKokkos<double>&  vel_grad,
        const ViewCArrayKokkos<size_t>&  elem_node_gids,
        const DCArrayKokkos<double>& node_coords,
        const DCArrayKokkos<double>& node_vel,
        const double vol,
        const double dt,
        const double rk_alpha);
};

#endif // end HEADER_H
