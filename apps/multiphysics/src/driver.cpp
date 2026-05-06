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

#include "driver.hpp"

// Resolve relative mesh file paths relative to the YAML input directory.
#include <filesystem>

// Headers for solver classes
#include "sgh_solver_3D.hpp"
#include "sgh_solver_rz.hpp"
#include "sgtm_solver_3D.hpp"
#include "level_set_solver.hpp"
#include "tlqs_solver_3D.hpp"

#include "region_fill.hpp"


// Initialize driver data.  Solver type, number of solvers
// Will be parsed from YAML input
void Driver::initialize()
{
    log_.info("Initializing Driver\n");
    Yaml::Node root;
    try
    {
        Yaml::Parse(root, yaml_file);
    }
    catch (const Yaml::Exception e)
    {
        log_.error("YAML parse exception %d: %s\n", (int)e.Type(), e.what());
        log_.flush();
        exit(0);
    }

    // Read the YAML input file
    parse_yaml(root, SimulationParamaters, Materials, BoundaryConditions);
    log_.info("Finished parsing YAML file\n");


    // Create initial mesh on rank 0
    swage::Mesh initial_mesh;
    MPICArrayKokkos<double> initial_node_coords;
    MPICArrayKokkos<double> final_node_coords;

    // Initialize communication plans
    // These are set up by elements::partition_mesh
    element_communication_plan.initialize(MPI_COMM_WORLD);
    node_communication_plan.initialize(MPI_COMM_WORLD);

    if (SimulationParamaters.mesh_input.source == mesh_input::file) {

        // Create and/or read mesh on rank 0
        if (rank == 0) {
            // Make mesh paths relative to the YAML location (so `file_path: meshes/abaqus.inp`
            // works regardless of where you run `./app/Fierro`).
            try {
                std::filesystem::path yaml_path(yaml_file ? yaml_file : "");
                std::filesystem::path mesh_path(SimulationParamaters.mesh_input.file_path);
                if (!yaml_path.empty() && !mesh_path.empty() && !mesh_path.is_absolute()) {
                    mesh_path = yaml_path.parent_path() / mesh_path;
                    mesh_path = mesh_path.lexically_normal();
                    SimulationParamaters.mesh_input.file_path = mesh_path.string();
                }
            } catch (...) {
                // Fall back to the original path if path resolution fails.
            }

            // Create and/or read mesh
            log_.info("Mesh file path: %s\n", SimulationParamaters.mesh_input.file_path.c_str());
            mesh_reader.set_mesh_file(SimulationParamaters.mesh_input.file_path.data());
            mesh_reader.read_mesh(initial_mesh, 
                                  initial_node_coords,
                                  SimulationParamaters.mesh_input, 
                                  num_dims);
        } // end if rank == 0
    }
    else if (SimulationParamaters.mesh_input.source == mesh_input::generate) {
        // Build mesh on rank 0
        if (rank == 0) {
            mesh_builder.build_mesh(initial_mesh, 
                                    initial_node_coords,
                                    SimulationParamaters);
        } // end if rank == 0
    }
    else{
        throw std::runtime_error("**** NO MESH INPUT OPTIONS PROVIDED IN YAML ****");
        return;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
    // Partition the mesh to all ranks
    if(world_size != 1) { // pass through the partitioning function if not a single rank
        elements::partition_mesh(initial_mesh, mesh, initial_node_coords, final_node_coords, element_communication_plan, node_communication_plan, world_size, rank);   
        // Verify communication plans (matches ELEMENTS decomp_example pattern)
        // element_communication_plan.verify_graph_communicator();
        // node_communication_plan.verify_graph_communicator();
        MPI_Barrier(MPI_COMM_WORLD);
    } else {
        mesh = initial_mesh;
        mesh.num_elems = initial_mesh.num_elems;
        mesh.num_nodes = initial_mesh.num_nodes;
        mesh.num_owned_elems = initial_mesh.num_elems;
        mesh.num_owned_nodes = initial_mesh.num_nodes;
        final_node_coords = initial_node_coords;
    }

    // partition_mesh() builds corner connectivity (corners_in_elem, corners_in_node)
    // but does not set the scalar counter mesh.num_corners. Many solvers size
    // State.corner.* using mesh.num_corners, so if it is left at 0 any
    // corner_mass(corner_gid) write becomes an out-of-bounds / heap-corrupting
    // store. There is one corner per node-in-element pair.
    mesh.num_corners = mesh.num_elems * mesh.num_nodes_in_elem;
    MPI_Barrier(MPI_COMM_WORLD);

    // Report per-rank mesh information. With the Logger we no longer need a
    // serialized-by-barriers print loop: every rank appends to its own buffer
    // (non-collective), and the collective flush() below dumps them to stdout
    // in rank order AND to each rank's per-rank file (Fierro_log_<rank>).
    {
        size_t num_ghost_nodes = mesh.num_nodes > mesh.num_owned_nodes ? mesh.num_nodes - mesh.num_owned_nodes : 0;
        size_t num_ghost_elems = mesh.num_elems > mesh.num_owned_elems ? mesh.num_elems - mesh.num_owned_elems : 0;
        log_.info("num_nodes=%zu (owned=%zu, ghost=%zu), num_elems=%zu (owned=%zu, ghost=%zu)\n",
                  mesh.num_nodes, mesh.num_owned_nodes, num_ghost_nodes,
                  mesh.num_elems, mesh.num_owned_elems, num_ghost_elems);
    }
    
    // For totals, root rank gathers and prints
    size_t local_nodes = mesh.num_nodes;
    size_t local_elems = mesh.num_elems;
    size_t local_owned_nodes = mesh.num_owned_nodes;
    size_t local_owned_elems = mesh.num_owned_elems;

    size_t total_nodes = 0, total_elems = 0, total_owned_nodes = 0, total_owned_elems = 0;

    MPI_Reduce(&local_nodes, &total_nodes, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_elems, &total_elems, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_owned_nodes, &total_owned_nodes, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_owned_elems, &total_owned_elems, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        size_t total_ghost_nodes = total_nodes > total_owned_nodes ? total_nodes - total_owned_nodes : 0;
        size_t total_ghost_elems = total_elems > total_owned_elems ? total_elems - total_owned_elems : 0;
        log_.info("Totals across all ranks: num_nodes=%zu (owned=%zu, ghost=%zu), num_elems=%zu (owned=%zu, ghost=%zu)\n",
                  total_nodes, total_owned_nodes, total_ghost_nodes,
                  total_elems, total_owned_elems, total_ghost_elems);
    }

    // Build boundary conditions
    const int num_bcs = BoundaryConditions.num_bcs;

    

    // Initialize node state
    std::vector<node_state> required_node_state = { node_state::coords };
    State.node.initialize(mesh.num_nodes, num_dims, required_node_state, node_communication_plan);

    // Copy the partitioned node coordinates to the state
    FOR_ALL(node_gid, 0, mesh.num_nodes, {
        for (size_t dim = 0; dim < num_dims; dim++) {
            State.node.coords(node_gid, dim) = final_node_coords(node_gid, dim);
        }
    });

    // --- calculate bdy sets ---//
    mesh.init_bdy_sets(num_bcs);
    tag_bdys(BoundaryConditions, mesh, State.node.coords);
    build_boundry_node_sets(mesh);


    // Setup the Solvers
    double time_final = SimulationParamaters.dynamic_options.time_final;
    for (size_t solver_id = 0; solver_id < SimulationParamaters.solver_inputs.size(); solver_id++) {

        if (SimulationParamaters.solver_inputs[solver_id].method == solver_input::SGH3D) {

            log_.info("Initializing dynx_FE solver\n");
            SGH3D* sgh_solver = new SGH3D(); 
            sgh_solver->set_logger(log_);

            sgh_solver->initialize(SimulationParamaters, 
                                   Materials, 
                                   mesh, 
                                   BoundaryConditions,
                                   State);

            // set the variables in the solver class
            setup_solver_vars(sgh_solver, 
                              solver_id);

            solvers.push_back(sgh_solver);


        } // end if SGH solver
        else if (SimulationParamaters.solver_inputs[solver_id].method == solver_input::SGHRZ) {

            log_.info("Initializing dynx_FE_RZ solver\n");
            SGHRZ* sgh_solver_rz = new SGHRZ(); 
            sgh_solver_rz->set_logger(log_);

            sgh_solver_rz->initialize(SimulationParamaters, 
                                   Materials, 
                                   mesh, 
                                   BoundaryConditions,
                                   State);

            // set the variables in the solver class
            setup_solver_vars(sgh_solver_rz, 
                              solver_id);

            solvers.push_back(sgh_solver_rz);
        } // end if SGHRZ solver
        else if (SimulationParamaters.solver_inputs[solver_id].method == solver_input::SGTM3D) {

            log_.info("Initializing thrmex_FE solver\n");
            SGTM3D* sgtm_solver_3d = new SGTM3D(); 
            sgtm_solver_3d->set_logger(log_);
        
            sgtm_solver_3d->initialize(SimulationParamaters, 
                                       Materials, 
                                       mesh, 
                                       BoundaryConditions,
                                       State);

            // set the variables in the solver class
            setup_solver_vars(sgtm_solver_3d, 
                              solver_id);

            solvers.push_back(sgtm_solver_3d);

        } // end if SGTM solver
        else if (SimulationParamaters.solver_inputs[solver_id].method == solver_input::levelSet) {

            log_.info("Initializing level set solver\n");
            LevelSet* level_set_solver = new LevelSet(); 
            level_set_solver->set_logger(log_);
        
            level_set_solver->initialize(SimulationParamaters, 
                                         Materials, 
                                         mesh, 
                                         BoundaryConditions,
                                         State);

            // set the variables in the solver class
            setup_solver_vars(level_set_solver, 
                              solver_id);

            solvers.push_back(level_set_solver);

        } // end if level set solver
        else if (SimulationParamaters.solver_inputs[solver_id].method == solver_input::TLQS3D) {
            log_.info("Initializing TLQS solver\n");
            TLQS3D* tlqs_solver = new TLQS3D(); 
            tlqs_solver->set_logger(log_);

            tlqs_solver->initialize(SimulationParamaters, 
                                    Materials, 
                                    mesh, 
                                    BoundaryConditions, State);

            // set the variables in the solver class
            setup_solver_vars(tlqs_solver, 
                              solver_id);

            solvers.push_back(tlqs_solver);

            log_.info("TLQS solver added to solvers vector\n");

        } // end if TLQS solver
        else {
            throw std::runtime_error("**** NO SOLVER INPUT OPTIONS PROVIDED IN YAML, OR OPTION NOT UNDERSTOOD ****");
            return;
        }

    } // end for loop over solvers

    // Solvers allocate nodal velocity with node_t::initialize(..., states) without a comm plan.
    // MPICArrayKokkos::communicate() requires initialize_comm_plan(node_communication_plan).
    if (State.node.vel.size() > 0) {
        State.node.vel.initialize_comm_plan(node_communication_plan);
    }
    if (State.node.vel_n0.size() > 0) {
        State.node.vel_n0.initialize_comm_plan(node_communication_plan);
    }

    // ----
    // setup the simulation by applying all the fills to the mesh

    fillGaussState_t fillGaussState;
    fillElemState_t  fillElemState;

    log_.info("Applying fills to the mesh\n");

    simulation_setup(SimulationParamaters, 
                     Materials, 
                     mesh, 
                     BoundaryConditions,
                     State,
                     fillGaussState,
                     fillElemState);

    log_.info("Fills applied to the mesh\n");
    // Allocate material state
    for (auto& solver : solvers) {
        solver->initialize_material_state(SimulationParamaters, 
                      Materials, 
                      mesh, 
                      BoundaryConditions,
                      State);
    } // end for over solvers


    // populate the material point state
    log_.info("Populating the material point state\n");
    material_state_setup(SimulationParamaters, 
                         Materials, 
                         mesh, 
                         BoundaryConditions,
                         State,
                         fillGaussState,
                         fillElemState);
    log_.info("Material point state populated\n");

    log_.info("Driver initialization complete\n");

    // Collective dump of the full initialization transcript to stdout (in rank
    // order) and to per-rank log files (Fierro_log_<rank>). Safe: every rank
    // reaches this point exactly once.
    log_.flush();
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn setup
///
/// \brief Calls the setup function for each of the created solvers
///
/////////////////////////////////////////////////////////////////////////////
void Driver::setup()
{
    log_.info("Inside driver setup\n");

    // allocate state, setup models, and apply fill instructions
    for (auto& solver : solvers) {
        solver->setup(SimulationParamaters, 
                      Materials, 
                      mesh, 
                      BoundaryConditions,
                      State);
    } // end for over solvers

    // Collective: dump buffered output. Every rank reaches this after all
    // solvers have returned from setup().
    log_.flush();

} // end setup function of driver


/////////////////////////////////////////////////////////////////////////////
///
/// \fn execute
///
/// \brief Calls the execute function for each of the created solvers
///
/////////////////////////////////////////////////////////////////////////////
void Driver::execute()
{
    log_.info("Inside driver execute\n");
    for (auto& solver : solvers) {
        solver->execute(SimulationParamaters, 
                        Materials, 
                        BoundaryConditions, 
                        mesh, 
                        State);
    } // loop over solvers

    // Collective: dump any remaining buffered output from the run. Solvers
    // are free to call log_.flush() at their own cycle boundaries; this
    // catches anything they didn't.
    log_.flush();
}


/////////////////////////////////////////////////////////////////////////////
///
/// \fn finalize
///
/// \brief Calls the finalize function of each of the solvers assuming the 
///        finalize function exists and deletes the solver
///
///
/////////////////////////////////////////////////////////////////////////////
void Driver::finalize()
{
    log_.info("Inside driver finalize\n");
    for (auto& solver : solvers) {
        if (solver->finalize_flag) {
            solver->finalize(SimulationParamaters, 
                             Materials, 
                             BoundaryConditions);
        }
    }
    // destroy FEA modules
    for (auto& solver : solvers) {
        log_.info("Deleting solver\n");
        delete solver;
    }

    // Last collective flush before Driver destruction. ~Logger will also
    // flush if MPI is still initialized, but doing it here keeps the
    // observable output ordered next to the finalize phase.
    log_.flush();
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn setup_solver_vars
///
/// \brief a function to set the solver variables based on user yaml inputs
///
///
/////////////////////////////////////////////////////////////////////////////
template <typename T>
void Driver::setup_solver_vars(T& a_solver, 
                               const size_t solver_id){


    // the final time of the simulation
    double time_final = this->SimulationParamaters.dynamic_options.time_final;

    // save the solver_id
    a_solver->solver_id = solver_id;

    // add other solver vars here

    // Check if the solver is of type SGTM3D and set use_moving_heat_source if specified
    if (this->SimulationParamaters.solver_inputs[solver_id].method == solver_input::SGTM3D) {
        // By default, do not use moving heat source unless specified in input
        if (this->SimulationParamaters.solver_inputs[solver_id].use_moving_heat_source) {
            a_solver->use_moving_heat_source = true;
        }

        else {
            a_solver->use_moving_heat_source = false;
        }
    }
    


    // -------
    // setting the ending times are tricky, requiring logic
            
    // set the start and ending times
    double t_end = this->SimulationParamaters.solver_inputs[solver_id].time_end;  // default is t=0
    if(solver_id==0){
        a_solver->time_start = 0.0;

        if(t_end <= 1.0e-14){
            // time wasn't set so end the solver at the final time of the calculation
            a_solver->time_end = fmax(0.0, time_final);
        }
        else {
            // use the specified time in the input file
            a_solver->time_end = t_end;
        } // end if time was set
    }
    else {
        // this is not the first solver run
        double time_prior_solver_ends = this->solvers[solver_id-1]->time_end;

        a_solver->time_start = time_prior_solver_ends;

        if (t_end <= 1.0e-14){
            // time wasn't set so end the solver at the final time of the calculation
            a_solver->time_end = time_final;
        }
        else {
            // use the specified time in the input file
            a_solver->time_end = fmin(t_end, time_final); // ensure t_end is bounded by final time
        } // end if time was set
        
    } // end if solver=0

    log_.info("Solver %zu start time = %g, ending time = %g\n",
              solver_id, a_solver->time_start, a_solver->time_end);

    return;
} // end setup solver vars function
