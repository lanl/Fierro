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

#include "sgh_solver_3D.hpp"
#include "simulation_parameters.hpp"
#include "material.hpp"
#include "boundary_conditions.hpp"
//#include "mesh.hpp""
#include "state.hpp"
#include "geometry_new.hpp"
#include "mesh_io.hpp"
#include "tipton_equilibration.hpp"
#include "fracture_stress_bc.hpp"
#include "reorientation_kinematics.hpp"
#include "user_defined_velocity_bc.hpp"
#include "mach_shock.hpp"


/////////////////////////////////////////////////////////////////////////////
///
/// \fn solve
///
/// Evolve the state according to the SGH method
///
/////////////////////////////////////////////////////////////////////////////

void SGH3D::execute(SimulationParameters_t& SimulationParamaters, 
                    Material_t& Materials, 
                    BoundaryCondition_t& BoundaryConditions, 
                    swage::Mesh_t& mesh, 
                    State_t& State)
{

    // Get MPI ranks and num ranks
    int rank;
    int num_ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

    // if (log) log->set_level(fierro::LogLevel::Off);

    contact_state_t Contact_State; // keeps track of contact variables
    if (doing_contact) {
        if (num_ranks > 1) {
            Kokkos::abort("MPI is not supported with contact.");
        }
        Contact_State.initialize(mesh.num_dims, mesh.num_nodes_in_patch, mesh.bdy_patches, mesh.num_bdy_nodes, mesh.num_bdy_patches,
                                 mesh.patches_in_elem, mesh.elems_in_patch, mesh.nodes_in_elem, mesh.nodes_in_patch,
                                 mesh.bdy_nodes, mesh.num_patches, mesh.num_nodes, State.node.coords,
                                 BoundaryConditions.contact_max_local_iter, BoundaryConditions.contact_max_global_iter);
        printf("MAX ITERATIONS PER CONTACT PAIR: %zu\n", Contact_State.max_local_iter);
        printf("MAX ITERATIONS FOR GLOBAL CONTACT: %zu\n", Contact_State.max_global_iter);
    }
    

    // cohesive zone initialization for fracture
    if (doing_fracture) {
        if (num_ranks > 1) {
            Kokkos::abort("MPI is not supported with cohesive zone fracture.");
        }
        const int fracture_bdy_set = BoundaryConditions.fracture_bc_id;
        if (fracture_bdy_set >= 0) {
            cohesive_zones_bank.initialize_fracture_bc(
            mesh.num_nodes,
            BoundaryConditions,
            fracture_bdy_set
            );
        }
    }

    // fracture reorientation validation mode initialization
    // scans user-defined velocity BC parameters for the reorientation mode flag
    // if enabled, stores the reorientation parameters and allocated validation-mode data needed
    if (doing_fracture) {
        cohesive_zones_bank.initialize_reorientation_mode(
            mesh.nodes_in_elem,
            mesh.num_nodes,
            mesh.num_elems,
            mesh.num_nodes_in_elem,
            State.node.coords,
            BoundaryConditions,
            doing_fracture
        );
    }
 
    double fuzz  = SimulationParamaters.DynamicOptions.fuzz;  // 1.e-16
    double tiny  = SimulationParamaters.DynamicOptions.tiny;  // 1.e-12
    double small = SimulationParamaters.DynamicOptions.small; // 1.e-8

    double graphics_dt_ival  = SimulationParamaters.OutputOptions.graphics_time_step;
    int    graphics_cyc_ival = SimulationParamaters.OutputOptions.graphics_iteration_step;

    // double time_initial = SimulationParamaters.DynamicOptions.time_initial;
    double time_final   = this->time_end; //SimulationParamaters.DynamicOptions.time_final;
    double dt_min   = SimulationParamaters.DynamicOptions.dt_min;
    double dt_max   = SimulationParamaters.DynamicOptions.dt_max;
    double dt_start = SimulationParamaters.DynamicOptions.dt_start;
    double dt_cfl   = SimulationParamaters.DynamicOptions.dt_cfl;

    int rk_num_stages = SimulationParamaters.DynamicOptions.rk_num_stages;
    int cycle_stop    = SimulationParamaters.DynamicOptions.cycle_stop;

    // initialize time, time_step, and cycles
    double time_value = this->time_start;  // was 0.0
    double dt = dt_start;

    // local memory for this solver
    CArrayKokkos <double> GaussPoint_pres(mesh.num_elems*mesh.num_gauss_in_elem);
    CArrayKokkos <double> GaussPoint_pres_denominator(mesh.num_elems*mesh.num_gauss_in_elem);
    CArrayKokkos <double> GaussPoint_volfrac_min(mesh.num_elems*mesh.num_gauss_in_elem);
    CArrayKokkos <double> GaussPoint_volfrac_limiter(mesh.num_elems*mesh.num_gauss_in_elem);
    
    // Create mesh writer
    MeshWriter mesh_writer; // Note: Pull to driver after refactoring evolution

    // --- Graphics vars ----
    CArray<double> graphics_times = CArray<double>(20000);
    graphics_times(0) = this->time_start; // was zero
    double graphics_time = this->time_start; // the times for writing graphics dump, was started at 0.0
    size_t output_id = 0; // the id for the outputs written

    if (log) log->info("Applying initial boundary conditions\n");
    boundary_velocity(mesh, BoundaryConditions, State.node.vel, time_value); // Time value = 0.0;

    // extensive energy tallies over the entire mesh
    double IE_t0 = 0.0;
    double KE_t0 = 0.0;
    double TE_t0 = 0.0;

    double cached_pregraphics_dt = fuzz;

    // the number of materials specified by the user input
    const size_t num_mats = Materials.num_mats;
    // Computing initial conditions for conservation checks
    // extensive IE
    IE_t0 = sum_domain_internal_energy(mesh,
                                       State.MeshtoMaterialMaps,
                                       State.MaterialPoints.mass,
                                       State.MaterialPoints.sie);

    if (log) log->info("Extensive IE = %f \n", IE_t0);
    // extensive KE
    KE_t0 = sum_domain_kinetic_energy(mesh,
                                      State.node.vel,
                                      State.node.mass);
    // extensive TE
    TE_t0 = IE_t0 + KE_t0;

    // domain mass for each material (they are at material points)
    double mass_domain_all_mats_t0 = sum_domain_material_mass(
        mesh,
        State.MeshtoMaterialMaps,
        State.MaterialPoints.mass);

    // node mass of the domain. 
    double mass_domain_nodes_t0 = 0.0;
    mass_domain_nodes_t0 = sum_domain_node_mass(mesh,
                                                State.node.coords,
                                                State.node.mass);

    if (log) log->info("nodal mass domain = %f \n", mass_domain_nodes_t0);

    // a flag to exit the calculation
    size_t stop_calc = 0;

    int64_t comm_time_ns_total = 0;

    auto time_1 = std::chrono::high_resolution_clock::now();

    // Write initial state at t=0
    if (log) log->info("Writing outputs to file at %f \n", graphics_time);
    if (log) log->flush();
    
    
    mesh_writer.write_mesh(
        mesh, 
        State, 
        SimulationParamaters,
        dt, 
        time_value, 
        graphics_times,
        SGH3D_State::required_node_state,
        SGH3D_State::required_gauss_pt_state,
        SGH3D_State::required_material_pt_state,
        this->solver_id);

    output_id++; // saved an output file

    graphics_time = time_value + graphics_dt_ival;


    // loop over the max number of time integration cycles
    for (size_t cycle = 0; cycle < cycle_stop; cycle++) {
        // stop calculation if flag
        if (stop_calc == 1) {
            break;
        }

        cached_pregraphics_dt = dt;

        // the smallest time step across all materials
        double min_dt_calc = dt_max;

        // calculating time step per material
        for(size_t mat_id = 0; mat_id < num_mats; mat_id++){

            // initialize the material dt
            double dt_mat = dt;

            // get the stable time step
            get_timestep(mesh,
                         State.node.coords,
                         State.node.vel,
                         State.GaussPoints.vol,
                         State.MaterialPoints.sspd,
                         State.MaterialPoints.eroded,
                         State.MaterialToMeshMaps.elem_in_mat_elem,
                         State.MaterialToMeshMaps.num_mat_elems.host(mat_id),
                         time_value,
                         graphics_time,
                         time_final,
                         dt_max,
                         dt_min,
                         dt_cfl,
                         dt_mat,
                         fuzz,
                         tiny,
                         mat_id);

            // NOTE: We need to use this same sort of pattern to calculate the max sounds speed (and minimum phi) for the 

            // save the smallest dt of all materials
            min_dt_calc = fmin(dt_mat, min_dt_calc);
        } // end for loop over all mats

        dt = min_dt_calc;  // save this dt time step

        // Global minimum dt across MPI ranks (each rank's CFL limit is local to its partition).
        {
            int init = 0;
            if (MPI_Initialized(&init) == MPI_SUCCESS && init) {
                MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            }
        }

        if (cycle == 0) {
            if (log) log->info("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
            if (log) log->flush();
        }
        // print time step every 10 cycles
        else if (cycle % 20 == 0) {
            if (log) log->info("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
            if (log) log->flush();
        } // end if


        // ---------------------------------------------------------------------
        //  integrate the solution forward to t(n+1) via Runge Kutta (RK) method
        // ---------------------------------------------------------------------
        for(size_t mat_id = 0; mat_id < num_mats; mat_id++){
            // save the values at t_n
            rk_init(State.node.coords,
                    State.node.coords_n0,
                    State.node.vel,
                    State.node.vel_n0,
                    State.MaterialPoints.sie,
                    State.MaterialPoints.sie_n0,
                    State.MaterialPoints.stress,
                    State.MaterialPoints.stress_n0,
                    mesh.num_dims,
                    mesh.num_elems,
                    mesh.num_nodes,
                    State.MaterialPoints.num_material_points.host(mat_id),
                    mat_id);
        } // end for mat_id


        // integrate solution forward in time
        for (size_t rk_stage = 0; rk_stage < rk_num_stages; rk_stage++) {
            // ---- RK coefficient ----
            double rk_alpha = 1.0 / ((double)rk_num_stages - (double)rk_stage);

            // dt_stage is the effective stage time increment used inside the cohesive
            // zone predictor/update so that the Forward-Euler-style incrementilization
            // stage increment is consistent with the RK stage 
            double dt_stage = rk_alpha * dt;
        
            // check if reorientation validation mode testing is on
            const bool reorient_mode = (doing_fracture && cohesive_zones_bank.reorientation_validation_mode);

            if (reorient_mode) {
                // prescribe reorientation kinematics to all nodes (skip over normal SGH evolution)
                ReorientationKinematics::prescribe_reorientation_kinematics(
                    mesh, State,
                    cohesive_zones_bank.initial_coords,
                    cohesive_zones_bank.cz_b_side_flag,
                    time_value,
                    dt_stage,
                    cohesive_zones_bank.omega_y, 
                    cohesive_zones_bank.omega_z,
                    cohesive_zones_bank.cz_opening_rate
                );
            }

            // if not in reorientation mode, proceed with normal SGH evolution
            if (!reorient_mode){

                // ---- Calculate velocity gradient for the element ----

                get_velgrad(State.GaussPoints.vel_grad,
                            mesh,
                            State.node.coords,
                            State.node.vel,
                            State.GaussPoints.vol);

                
                // compute the shock detector at each gauss point for each material
                State.GaussPoints.shock_detector.set_values(0.0);
                MATAR_FENCE();

                MachShockDetector::detect_shock(mesh,
                    State.GaussPoints.vel_grad,
                    State.GaussPoints.shock_detector,
                    State.GaussPoints.vol,
                    State.MaterialPoints.sspd,
                    State.MaterialToMeshMaps.elem_in_mat_elem,
                    State.MaterialToMeshMaps.num_mat_elems,
                    fuzz,
                    small,
                    num_mats);
                
                // Communicate shock detector based on the element communication plan
                State.GaussPoints.shock_detector.communicate();

                set_corner_force_zero(mesh, State.corner.force);


                // ---- calculate the forces on the vertices and evolve stress (hypo model) ----
                for(size_t mat_id = 0; mat_id < num_mats; mat_id++){

                    get_force(Materials,
                            mesh,
                            State.GaussPoints.vol,
                            State.GaussPoints.vel_grad,
                            State.GaussPoints.shock_detector,
                            State.MaterialPoints.eroded,
                            State.corner.force,
                            State.node.coords,
                            State.node.vel,
                            State.MaterialPoints.den,
                            State.MaterialPoints.sie,
                            State.MaterialPoints.pres,
                            State.MaterialPoints.stress,
                            State.MaterialPoints.sspd,
                            State.MaterialCorners.force,
                            State.MaterialPoints.mat_volfrac,
                            State.MaterialPoints.geo_volfrac,
                            State.corners_in_mat_elem,
                            State.MaterialToMeshMaps.elem_in_mat_elem,
                            State.MaterialToMeshMaps.num_mat_elems.host(mat_id),
                            mat_id,
                            fuzz,
                            small,
                            dt,
                            rk_alpha);

                    if (Materials.MaterialEnums.host(mat_id).StrengthType == model::incrementBased) {
                        update_stress(Materials,
                                    mesh,
                                    State.GaussPoints.vol,
                                    State.node.coords,
                                    State.node.vel,
                                    State.GaussPoints.vel_grad,
                                    State.MaterialPoints.den,
                                    State.MaterialPoints.sie,
                                    State.MaterialPoints.pres,
                                    State.MaterialPoints.stress,
                                    State.MaterialPoints.stress_n0,
                                    State.MaterialPoints.sspd,
                                    State.MaterialPoints.eos_state_vars,
                                    State.MaterialPoints.strength_state_vars,
                                    State.MaterialPoints.shear_modulii,
                                    State.MaterialToMeshMaps.elem_in_mat_elem,
                                    State.MaterialToMeshMaps.num_mat_elems.host(mat_id),
                                    mat_id,
                                    fuzz,
                                    small,
                                    time_value,
                                    dt,
                                    rk_alpha,
                                    cycle);
                    } // end if on increment

                } // end for mat_id


                // ---- Calculate boundary and body forces ---- //
                // setting nodal force to zero here, 
                // the node force stores the BCs supplied forces and body forces
                State.node.force.set_values(0.0);

                // call stress BC's routine
                boundary_stress(mesh, 
                                BoundaryConditions, 
                                State.node.force, 
                                State.node.coords,
                                time_value);

            } else {
                // kinematics-only: start forces at zero (you can still add cohesive forces below)
                State.node.force.set_values(0.0);
                set_corner_force_zero(mesh, State.corner.force);
            }

            // apply cohesive zone nodal forces (fracture)
            if (doing_fracture) {
                boundary_fracture_force(State, mesh, dt_stage, cohesive_zones_bank,
                                        time_value, cycle, rk_stage, rk_num_stages);
            } // end if doing_fracture

            // SKIP SGH SOLVER EVOLUTION (REORIENTATION TESTING MODE)
            if (reorient_mode) {
                continue; // skip update velocity, boundary vel, contact, energy, position, etc.
            }

            // apply contact forces to boundary patches
            if (doing_contact && rk_stage == rk_num_stages-1) 
            {
                //contact_bank.update_nodes(mesh, State);
                Contact_State.contact_forces.set_values(0);
                Contact_State.contact_force.set_values(0);
                if (doing_preload) {
                    double preload_time = (time_final-time_value)/2;
                    //preload_time = 1;
                    if (time_value < preload_time) {
                        boundary_contact_force(State, mesh, preload_time, Contact_State);
                    } else {
                        boundary_contact_force(State, mesh, 5*dt*rk_alpha, Contact_State);
                    }
                } else {
                    //boundary_contact_force(State, mesh, 5*dt*rk_alpha, Contact_State);
                    boundary_contact_force(State, mesh, dt, Contact_State);
                }
            }

            // ---- Update nodal velocities ---- //
            update_velocity(rk_alpha,
                            dt,
                            mesh,
                            State.node.vel,
                            State.node.vel_n0,
                            State.node.mass,
                            State.node.force,
                            State.corner.force,
                            Contact_State.contact_force,
                            doing_contact);

            // ---- apply velocity boundary conditions to the boundary patches----
            boundary_velocity(mesh, BoundaryConditions, State.node.vel, time_value);

            // ---- apply contact boundary conditions to the boundary patches----
            // boundary_contact(mesh, BoundaryConditions, State.node.vel, time_value);

            // ----- Communication of the nodal velocity -----
            MPI_Barrier(MPI_COMM_WORLD);
            auto comm_t0 = std::chrono::high_resolution_clock::now();

            State.node.vel.communicate();
            State.node.vel_n0.communicate();
            auto comm_t1 = std::chrono::high_resolution_clock::now();
            comm_time_ns_total += std::chrono::duration_cast<std::chrono::nanoseconds>(comm_t1 - comm_t0).count();

            for (size_t mat_id = 0; mat_id < num_mats; mat_id++) {
                // ---- Update specific internal energy in the elements ----
                update_energy(rk_alpha,
                              dt,
                              mesh,
                              State.node.vel,
                              State.node.vel_n0,
                              State.MaterialPoints.sie,
                              State.MaterialPoints.sie_n0,
                              State.MaterialPoints.mass,
                              State.MaterialCorners.force,
                              State.corners_in_mat_elem,
                              State.MaterialToMeshMaps.elem_in_mat_elem,
                              State.MaterialToMeshMaps.num_mat_elems.host(mat_id),
                              mat_id);
            } // end for mat_id

            // ---- Update nodal positions ----
            update_position(rk_alpha,
                            dt,
                            mesh.num_dims,
                            mesh.num_nodes,
                            State.node.coords,
                            State.node.coords_n0,
                            State.node.vel,
                            State.node.vel_n0);

            // ---- Calculate cell volume for next time step ----
            geometry::get_vol(State.GaussPoints.vol, State.node.coords, mesh);



            // ---- Calculate MaterialPoints state (den, pres, sound speed, stress) for next time step ----
            for(size_t mat_id = 0; mat_id < num_mats; mat_id++){

                update_state(Materials,
                             mesh,
                             State.node.coords,
                             State.node.vel,
                             State.GaussPoints.vel_grad,
                             State.MaterialPoints.den,
                             State.MaterialPoints.pres,
                             State.MaterialPoints.stress,
                             State.MaterialPoints.stress_n0,
                             State.MaterialPoints.sspd,
                             State.MaterialPoints.sie,
                             State.MaterialPoints.mat_volfrac,
                             State.MaterialPoints.geo_volfrac,
                             State.GaussPoints.vol,
                             State.MaterialPoints.mass,
                             State.MaterialPoints.eos_state_vars,
                             State.MaterialPoints.strength_state_vars,
                             State.MaterialPoints.eroded,
                             State.MaterialPoints.shear_modulii,
                             State.MaterialToMeshMaps.elem_in_mat_elem,
                             time_value,
                             dt,
                             rk_alpha,
                             cycle,
                             State.MaterialToMeshMaps.num_mat_elems.host(mat_id),
                             mat_id);
            } // end for mat_id

            // ----
            // Notes on strength:
            //    1) hyper-elastic strength models are called in update_state
            //    2) hypo-elastic strength models are called in get_force
            //    3) strength models must be added by the user in user_mat.cpp


            // apply pressure relaxation on material volume fractions
            if(Materials.EquilibrationModels != model::noEquilibration){
                TiptonEquilibrationModel::mat_equilibration(
                    Materials, 
                    mesh, 
                    State,
                    GaussPoint_pres,
                    GaussPoint_pres_denominator,
                    GaussPoint_volfrac_min,
                    GaussPoint_volfrac_limiter,
                    dt,
                    rk_alpha,
                    fuzz,
                    small);
            } // end if on applying equilibration

            // apply pressure relaxation on geometric volume fractions
            if(Materials.GeoEquilibrationModels != model::noEquilibration){
                TiptonEquilibrationModel::geo_equilibration(
                    Materials, 
                    mesh, 
                    State,
                    GaussPoint_pres,
                    GaussPoint_pres_denominator,
                    GaussPoint_volfrac_min,
                    GaussPoint_volfrac_limiter,
                    dt,
                    rk_alpha,
                    fuzz,
                    small);
            } // end if on applying geometric equilibration

        } // end of RK loop

        // Apply ALE per material with ALE active

        // 1. Compute the mesh velocity. Final location minus current (from Lagrange)/(dt (d\tau = 1).  (maybe make DT = 1, pseudotime).
        // Communicate the mesh velocity from the lagrangian step 
        std::cout << "Computing mesh velocity" << std::endl;
        FOR_ALL(node_gid, 0, mesh.num_nodes, {
            for(size_t dim = 0; dim < mesh.num_dims; dim++){
                this->mesh_node_velocity(node_gid, dim) = (this->mesh_node_target_coords(node_gid, dim) - State.node.coords(node_gid, dim)) / 1.0; // dt = 1.0 for now
            }
        });
        

        // ================================================================
        // Step 1: build the volume matrix for nodal DG at tau=0
        // ================================================================
        build_element_geometry(mesh, tables, State.node.coords, this->elem_det_jac, this->inv_jac_ijq,
            mesh.num_elems, Quad.num_qpts_in_elem, mesh.num_nodes_in_elem);
        build_lumped_volume(FERefElem, Quad, tables, elem_det_jac, State.corner.volume,
            mesh.num_elems, Quad.num_qpts_in_elem, mesh.num_nodes_in_elem);
        Kokkos::fence();
            

        // 2. Compute advection CFL, and pseudo DT. Solver embedded inside of the solver.
        // -----------------------------------------------------
        const double max_vel = 1.0; // the CFL velocity used for calculating d_tau
        double h_cfl = 1.e-6;       // the CFL length scale for calculating d_tau
        double dt_tau = 1.e-6;          // dt_tau from CFL at start, this time is pseudo time
        double tau = 0.0; 
        double tau_final = 1.0;


        // Useful tmp vairables
        const size_t num_elems = mesh.num_elems;
        const size_t num_nodes = mesh.num_nodes;
        const size_t num_nodes_in_elem = mesh.num_nodes_in_elem;
        const size_t num_qpts_in_elem = Quad.num_qpts_in_elem;
        const size_t num_surfs = mesh.num_surfs;
        const size_t num_qpts_in_surf = SurfQuad.num_qpts_in_surf;
        const size_t num_surfs_in_elem = mesh.num_surfs_in_elem;
        const size_t elem_dims = mesh.num_dims;

        // 3. Do Time integrator with multiple RK stages until a time of 1.

        // --------------------------------------------------
        // Time integration loop
        // WARNING: Be careful to not use intermediate state from remapped fields in not yet remapped fields.
        for(size_t cycle = 0; cycle < max_cycles; cycle++){
            
            if(cycle%10 == 0) printf(" time = %.4f \n", tau);


            // --------------------------------------------------
            // Step 1a: Store time level n state

            FOR_ALL(idx, 0, num_elems*num_nodes_in_elem, {

                const size_t elem_gid = idx / num_nodes_in_elem;
                const size_t node_lid = idx % num_nodes_in_elem;

                elem_corner_vol_n(elem_gid, node_lid) = elem_corner_vol(elem_gid, node_lid);

                const size_t corner_gid = Mesh.corners_in_elem(elem_gid, node_lid);
                corner_field_n(corner_gid) = corner_field(corner_gid);
            });

            FOR_ALL(node_gid, 0, num_nodes, {
                for(size_t dim=0; dim<elem_dims; dim++){
                    node_coords_n(node_gid, dim)   = node_coords(node_gid, dim);
                    node_velocity_n(node_gid, dim) = node_velocity(node_gid, dim);
                }
            });


            // ------------------------------------------------------
            // Step 1b: get CFL time step for moving mesh

            // x -> x^(1/3) is monotone, so the smallest length comes from the
            // smallest quadrature volume: one transcendental instead of one per qpt.
            double min_vol_loc;
            double min_vol_qpt;
            FOR_REDUCE_MIN(idx, 0, num_elems*num_qpts_in_elem,
                        min_vol_loc, {
                const size_t elem_gid = idx / num_qpts_in_elem;
                const size_t qpt_lid  = idx % num_qpts_in_elem;

                const double vol_qpt = Quad.qpt_weights(qpt_lid)*elem_det_jac(elem_gid, qpt_lid);
                if(vol_qpt < min_vol_loc) min_vol_loc = vol_qpt;
            }, min_vol_qpt);
            h_cfl = pow(min_vol_qpt, 0.3333333);

            dt = 0.1*h_cfl/max_vel; // pseudo time step used for the remap

            // A tangled element gives a non-positive quadrature volume, so the CFL
            // length becomes NaN. NaN then defeats both the max_time test and the
            // mass conservation check below, so the run has to stop here.
            if(!(dt > 0.0)){
                printf("\n STOPPING at time = %.6f, cycle %zu: CFL length is %g,"
                    " the mesh has tangled.\n", time, cycle, h_cfl);
                break;
            }


            // Runge Kutta time integration levels
            for(size_t rk_stage=0; rk_stage<rk_num_stages; rk_stage++){

                // RK coefficient
                const double rk_alpha = 1.0 / ((double)rk_num_stages - (double)rk_stage);


                // ------------------------------------------------------
                // Step 2: calculate the mesh velocity at time level k

                FOR_ALL(node_gid, 0, num_nodes,{
                    // new velocity, it is Taylor-Green vortex
                    // PI is defined in mesh class
                    node_velocity(node_gid, 0) =  sin(PI*node_coords(node_gid, 0))*cos(PI*node_coords(node_gid, 1));
                    node_velocity(node_gid, 1) = -cos(PI*node_coords(node_gid, 0))*sin(PI*node_coords(node_gid, 1));
                    node_velocity(node_gid, 2) = 0.0;
                });


                // ----------------------------------------------------------
                // Step 3: Calculate the surface fluxes at quadrature points

                build_surface_flux(Mesh, RefSurf, SurfQuad, tables, node_coords, node_velocity, corner_field,
                                surf_qpt_qpt_map, RHS_surf_flux,
                                num_surfs, num_qpts_in_surf, num_nodes_in_elem);


                // -------------------------------------------------
                // Step 4: Build RHS of DG equations in the element

                assemble_rhs(Mesh, RefSurf, Quad, tables, elem_det_jac, inv_jac_ijq,
                            corner_field, corner_field_n, node_velocity, elem_corner_vol_n,
                            RHS_surf_flux, RHS_elem, qpt_vol_flux,
                            rk_alpha, dt,
                            num_elems, num_qpts_in_elem, num_nodes_in_elem,
                            num_surfs_in_elem, num_qpts_in_surf, elem_dims);


                // ================================================================
                // Step 5: Move the mesh to the new location
                FOR_ALL(node_gid, 0, num_nodes,{
                    // new position of the mesh
                    node_coords(node_gid, 0) = node_coords_n(node_gid, 0) + 0.5*(node_velocity(node_gid, 0)+node_velocity_n(node_gid, 0)) * rk_alpha * dt; 
                    node_coords(node_gid, 1) = node_coords_n(node_gid, 1) + 0.5*(node_velocity(node_gid, 1)+node_velocity_n(node_gid, 1)) * rk_alpha * dt;
                    // z-coords never change
                });


                // ================================================================
                // Step 6: build the diagonal volume matrix for nodal DG after the mesh moved
                build_element_geometry(Mesh, tables, node_coords, elem_det_jac, inv_jac_ijq,
                                    num_elems, num_qpts_in_elem, num_nodes_in_elem);
                build_lumped_volume(FERefElem, Quad, tables, elem_det_jac, elem_corner_vol,
                                    num_elems, num_qpts_in_elem, num_nodes_in_elem);


                // -----------------------------------------------------
                // 7. Solve M * u^{n+1} = RHS where M is diagonal

                FOR_ALL(idx, 0, num_elems*num_nodes_in_elem, {

                    const size_t elem_gid = idx / num_nodes_in_elem;
                    const size_t dof_lid  = idx % num_nodes_in_elem;

                    const size_t corner_gid = Mesh.corners_in_elem(elem_gid, dof_lid);
                    corner_field(corner_gid) = RHS_elem(elem_gid, dof_lid)/elem_corner_vol(elem_gid, dof_lid);
                });


                // -----------------------------------------------------
                // 8. A slope/bound limiter would be applied to corner_field here;
                //    see limit_corner_field in remap_dg_lumped_test.cpp.

            } // end Runge Kutta time level loop


            // ================================================================
            // Step 7: update time
            time += dt;

            
        } // end loop over cycle

        

        
        for(size_t mat_id = 0; mat_id < num_mats; mat_id++){

            if(Materials.MaterialEnums.host(mat_id).ALEType == model::ALE){
                std::cout << "Applying ALE for material " << mat_id << std::endl;
                // CALL ALE HERE
            } // end if on applying ALE
        } // end for mat_id

        // increment the time
        time_value += dt;

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
            if (log) log->info("Writing outputs to file at %f \n", graphics_time);
            if (log) log->flush();
            mesh_writer.write_mesh(mesh,
                                   State,
                                   SimulationParamaters,
                                   dt,
                                   time_value,
                                   graphics_times,
                                   SGH3D_State::required_node_state,
                                   SGH3D_State::required_gauss_pt_state,
                                   SGH3D_State::required_material_pt_state,
                                   this->solver_id);
            output_id++;
            graphics_time = (double)(output_id) * graphics_dt_ival;

            dt = cached_pregraphics_dt;
        } // end if

        // end of calculation
        if (time_value >= time_final) {
            break;
        }
    } // end for cycle loop

    auto time_2    = std::chrono::high_resolution_clock::now();
    auto calc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(time_2 - time_1).count();

    if (log) log->info("Calculation time in seconds: %f \n", calc_time * 1e-9);
    
    
    
    // Print communication time and node info in rank order
    for (int print_rank = 0; print_rank < num_ranks; ++print_rank) {
        if (rank == print_rank) {
            if (log) log->info("SGH3D execute: rank %d / %d total nodal velocity communication time: %f s\n",
                   rank, num_ranks, static_cast<double>(comm_time_ns_total) * 1e-9);

            if (log) log->info("rank %d / %d owned nodes: %zu, ghost nodes: %zu\n",
                   rank, num_ranks, mesh.num_owned_nodes, mesh.num_ghost_nodes);
            
            // Also, print the communication volume on each rank for the node communication plan

            // Print communication volume for node communication plan (total send and recv counts)
            if (State.node.vel.comm_plan_ != nullptr) {
                CommunicationPlan* node_comm = State.node.vel.comm_plan_;
                int total_send = node_comm->total_send_count;
                int total_recv = node_comm->total_recv_count;
                if (log) log->info("rank %d node communication volume - send: %d, recv: %d, total: %d\n",
                       rank, total_send, total_recv, total_send + total_recv);
            }
 
            fflush(stdout);
        }
        // Synchronize: ensure only one rank prints at a time
 
        MPI_Barrier(MPI_COMM_WORLD);

    }



        
    // ---- Calculate energy tallies ----
    double IE_tend = 0.0;
    double KE_tend = 0.0;
    double TE_tend = 0.0;

    // // extensive IE
    IE_tend = sum_domain_internal_energy(mesh,
                                       State.MeshtoMaterialMaps,
                                       State.MaterialPoints.mass,
                                       State.MaterialPoints.sie);

    // extensive KE
    KE_tend = sum_domain_kinetic_energy(mesh,
                                        State.node.vel,
                                        State.node.mass);
    // extensive TE
    TE_tend = IE_tend + KE_tend;

    if (log) log->info("Time=0:   KE = %.14f, IE = %.14f, TE = %.14f \n", KE_t0, IE_t0, TE_t0);
    if (log) log->info("Time=End: KE = %.14f, IE = %.14f, TE = %.14f \n", KE_tend, IE_tend, TE_tend);
    if (log) log->info("total energy change = %.15e \n", TE_tend - TE_t0);

    // domain mass for each material (they are at material points)
    double mass_domain_all_mats_tend = 0.0;
    mass_domain_all_mats_tend = sum_domain_material_mass(
        mesh,
        State.MeshtoMaterialMaps,
        State.MaterialPoints.mass);

    // node mass of the domain
    double mass_domain_nodes_tend = 0.0;
    mass_domain_nodes_tend = sum_domain_node_mass(mesh,
                                                  State.node.coords,
                                                  State.node.mass);

    if (log) log->info("material mass conservation error = %f \n", mass_domain_all_mats_tend - mass_domain_all_mats_t0);
    if (log) log->info("nodal mass conservation error = %f \n", mass_domain_nodes_tend - mass_domain_nodes_t0);
    if (log) log->info("nodal and material mass error = %f \n\n", mass_domain_nodes_tend - mass_domain_all_mats_tend);
} // end of SGH execute

/////////////////////////////////////////////////////////////////////////////
///
/// \fn max_Eigen3D
///
/// \brief Get the maximum eigenvalues of a given tensor
///
/// \param Input tensor
///
/////////////////////////////////////////////////////////////////////////////
KOKKOS_FUNCTION
double max_Eigen3D(const ViewCArrayKokkos<double> tensor)
{
    // Compute largest eigenvalue of a 3x3 tensor
    // Algorithm only works if tensor is symmetric
    double pi  = 3.141592653589793;
    size_t dim = tensor.dims(0);
    double trace, det;
    trace  = tensor(0, 0) + tensor(1, 1) + tensor(2, 2);
    det    = tensor(0, 0) * (tensor(1, 1) * tensor(2, 2) - tensor(1, 2) * tensor(2, 1));
    det   -= tensor(0, 1) * (tensor(1, 0) * tensor(2, 2) - tensor(1, 2) * tensor(2, 0));
    det   += tensor(0, 2) * (tensor(1, 0) * tensor(2, 1) - tensor(1, 1) * tensor(2, 0));
    trace /= 3.; // easier for computation
    double p2 = pow((tensor(0, 0) - trace), 2) + pow((tensor(1, 1) - trace), 2) +
        pow((tensor(2, 2) - trace), 2);
    p2 += 2. * (pow(tensor(0, 1), 2) + pow(tensor(0, 2), 2) + pow(tensor(1, 2), 2));

    double p = sqrt(p2 / 6.);

    // check for nan
    if (det != det) {
        return 0;
    }

    if (det == 0) {
        return 0;
    }

    double B_array[9];
    ViewCArrayKokkos<double> B(B_array, 3, 3);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i == j) {
                B(i, i) = (1 / p) * (tensor(i, i) - trace);
            }
            else{
                B(i, j) = (1 / p) * tensor(i, j);
            } // end if
        } // end for j
    } // end for i

    double r, phi;
    r  =  B(0, 0) * (B(1, 1) * B(2, 2) - B(1, 2) * B(2, 1));
    r -= B(0, 1) * (B(1, 0) * B(2, 2) - B(1, 2) * B(2, 0));
    r += B(0, 2) * (B(1, 0) * B(2, 1) - B(1, 1) * B(2, 0));
    r /= 2;
    // first two cases are to handle numerical difficulties.
    // In exact math -1 <= r <= 1
    if (r <= -1) {
        phi = pi / 3;
    }
    else if (r >= 1) {
        phi = 0;
    }
    else{
        phi = acos(r) / 3.;
    } // end if

    double eig1, eig2, eig3;
    eig1 = trace + 2 * p * cos(phi);
    eig2 = trace + 2 * p * cos(phi + (2 * pi / 3));
    eig3 = 3 * trace - (eig1 + eig2);

    double abs_max_val = fmax(fmax(fabs(eig1), fabs(eig2)), fabs(eig3));

    return abs_max_val;
}

// a function to tally the internal energy
double sum_domain_internal_energy(
    const swage::Mesh_t& mesh,
    const MeshtoMaterialMap_t& MeshtoMaterialMaps,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_mass,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_sie)
{
    double IE_loc_sum = 0.0;
    double sum = 0.0;
    FOR_REDUCE_SUM(elem_gid, 0, mesh.num_owned_elems, sum, {
        size_t num_materials = MeshtoMaterialMaps.num_mats_in_elem(elem_gid);

        for(size_t mat_lid = 0; mat_lid < num_materials; mat_lid++){

            size_t mat_id = MeshtoMaterialMaps.mats_in_elem(elem_gid, mat_lid);
            size_t mat_elem_gid = MeshtoMaterialMaps.mat_elems_in_elem(elem_gid, mat_lid);

            sum += MaterialPoints_mass(mat_id, mat_elem_gid) * MaterialPoints_sie(mat_id, mat_elem_gid);
        }
    }, IE_loc_sum);
    Kokkos::fence();
    double IE_sum = 0.0;
    MPI_Allreduce(&IE_loc_sum, &IE_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return IE_sum;
} // end function

/////////////////////////////////////////////////////////////////////////////
///
/// \fn sum_domain_kinetic_energy
///
/// \brief <insert brief description>
///
/// <Insert longer more detailed description which
/// can span multiple lines if needed>
///
/// \param <function parameter description>
/// \param <function parameter description>
/// \param <function parameter description>
///
/// \return <return type and definition description if not void>
///
/////////////////////////////////////////////////////////////////////////////
double sum_domain_kinetic_energy(
    const swage::Mesh_t& mesh,
    const MPICArrayKokkos<double>& node_vel,
    const DCArrayKokkos<double>& node_mass)
{
    // extensive KE
    double KE_sum = 0.0;
    double KE_loc_sum = 0.0;
    double KE_global_sum = 0.0;

    FOR_REDUCE_SUM(node_gid, 0, mesh.num_owned_nodes, KE_loc_sum, {
        if(mesh.shared_tally_owned_nodes(node_gid)){
            double ke = 0;

            for (size_t dim = 0; dim < mesh.num_dims; dim++) {
                ke += node_vel(node_gid, dim) * node_vel(node_gid, dim); // 1/2 at end
            } // end for

            KE_loc_sum += node_mass(node_gid) * ke;
        }
    }, KE_sum);
    Kokkos::fence();

    MPI_Allreduce(&KE_sum, &KE_global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return 0.5 * KE_global_sum;
} // end function

// // a function to tally the material point masses
// double sum_domain_material_mass(
//     const DRaggedRightArrayKokkos<double>& MaterialPoints_mass,
//     const size_t num_mat_points,
//     const size_t mat_id)
// {
//     double mass_domain = 0.0;
//     double mass_loc_domain;

//     FOR_REDUCE_SUM(matpt_lid, 0, num_mat_points, mass_loc_domain, {
//         mass_loc_domain += MaterialPoints_mass(mat_id,matpt_lid);
//     }, mass_domain);
//     Kokkos::fence();

//     return mass_domain;
// } // end function

// a function to tally the material point masses
double sum_domain_material_mass(
    const swage::Mesh_t& mesh,
    const MeshtoMaterialMap_t& MeshtoMaterialMaps,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_mass)
{
    double mass_loc_sum = 0.0;
    double sum = 0.0;
    FOR_REDUCE_SUM(elem_gid, 0, mesh.num_owned_elems, sum, {
        size_t num_materials = MeshtoMaterialMaps.num_mats_in_elem(elem_gid);

        for(size_t mat_lid = 0; mat_lid < num_materials; mat_lid++){

            size_t mat_id = MeshtoMaterialMaps.mats_in_elem(elem_gid, mat_lid);
            size_t mat_elem_gid = MeshtoMaterialMaps.mat_elems_in_elem(elem_gid, mat_lid);

            sum += MaterialPoints_mass(mat_id, mat_elem_gid);
        }
    }, mass_loc_sum);
    Kokkos::fence();
    double mass_sum = 0.0;
    MPI_Allreduce(&mass_loc_sum, &mass_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return mass_sum;
} // end function

/////////////////////////////////////////////////////////////////////////////
///
/// \fn sum_domain_node_mass
///
/// \brief <insert brief description>
///
/// <Insert longer more detailed description which
/// can span multiple lines if needed>
///
/// \param <function parameter description>
/// \param <function parameter description>
/// \param <function parameter description>
///
/// \return <return type and definition description if not void>
///
/////////////////////////////////////////////////////////////////////////////
double sum_domain_node_mass(const swage::Mesh_t& mesh,
    const MPICArrayKokkos<double>& node_coords,
    const DCArrayKokkos<double>& node_mass)
{
    double mass_domain = 0.0;
    double mass_loc_domain;
    double total_mass = 0.0;

    FOR_REDUCE_SUM(node_gid, 0, mesh.num_owned_nodes, mass_loc_domain, {
        if(mesh.shared_tally_owned_nodes(node_gid)){
            if (mesh.num_dims == 2) {
                mass_loc_domain += node_mass(node_gid) * node_coords(node_gid, 1);
            }
            else{
                mass_loc_domain += node_mass(node_gid);
            }
        }
    }, mass_domain);
    Kokkos::fence();

    MPI_Allreduce(&mass_domain, &total_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return total_mass;
} // end function


// set the corner forces to zero
void set_corner_force_zero(const swage::Mesh_t& mesh,
    const DCArrayKokkos<double>& corner_force)
{
    // set corner force to zero
    FOR_ALL(corner_gid, 0, mesh.num_corners, {
        for (size_t dim = 0; dim < mesh.num_dims; dim++) {
            corner_force(corner_gid, dim) = 0.0;
        }
    }); // end parallel for corners
} // end function


