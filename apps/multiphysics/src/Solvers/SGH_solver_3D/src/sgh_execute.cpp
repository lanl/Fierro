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
#include "fracture_stress_bc.h"
#include "reorientation_kinematics.hpp"
#include "user_defined_velocity_bc.h"


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
                    swage::Mesh& mesh, 
                    State_t& State)
{

    contact_state_t Contact_State; // keeps track of contact variables
    if (doing_contact) {
        Contact_State.initialize(mesh.num_dims, mesh.num_nodes_in_patch, mesh.bdy_patches, mesh.num_bdy_nodes, mesh.num_bdy_patches,
                                 mesh.patches_in_elem, mesh.elems_in_patch, mesh.nodes_in_elem, mesh.nodes_in_patch,
                                 mesh.bdy_nodes, mesh.num_patches, mesh.num_nodes, State.node.coords);
    }

// THROTTLE OUTPUT COMMENT OUT FOR DESCALING PRINTS 
    // writing output file for time, alpha, lambda, traction test
    FILE* cz_fp = fopen("time_alpha_lambda_traction_outputs.txt", "w");
    if (!cz_fp) { perror("fopen"); }

    // no buffering (writes show up immediately)
    setvbuf(cz_fp, nullptr, _IONBF, 0);

    fprintf(cz_fp, "time,alpha,lambda,Tn\n");
    fflush(cz_fp);

    // output every x us
    constexpr double CZ_OUT_DT = 1e-3;

    // stride variables (cycle-based)
    static bool   cz_stride_init = false;
    static size_t cz_stride = 0;
    static size_t cz_next_cycle = 0;
// THROTTLE OUTPUT COMMENT OUT FOR DESCALING PRINTS 

    double fuzz  = SimulationParamaters.dynamic_options.fuzz;  // 1.e-16
    double tiny  = SimulationParamaters.dynamic_options.tiny;  // 1.e-12
    double small = SimulationParamaters.dynamic_options.small; // 1.e-8

    double graphics_dt_ival  = SimulationParamaters.output_options.graphics_time_step;
    int    graphics_cyc_ival = SimulationParamaters.output_options.graphics_iteration_step;

    // double time_initial = SimulationParamaters.dynamic_options.time_initial;
    double time_final   = this->time_end; //SimulationParamaters.dynamic_options.time_final;
    double dt_min   = SimulationParamaters.dynamic_options.dt_min;
    double dt_max   = SimulationParamaters.dynamic_options.dt_max;
    double dt_start = SimulationParamaters.dynamic_options.dt_start;
    double dt_cfl   = SimulationParamaters.dynamic_options.dt_cfl;

    int rk_num_stages = SimulationParamaters.dynamic_options.rk_num_stages;
    int cycle_stop    = SimulationParamaters.dynamic_options.cycle_stop;

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
    size_t output_id=0; // the id for the outputs written

    std::cout << "Applying initial boundary conditions" << std::endl;
    boundary_velocity(mesh, BoundaryConditions, State.node.vel, time_value); // Time value = 0.0;

    // extensive energy tallies over the entire mesh
    double IE_t0 = 0.0;
    double KE_t0 = 0.0;
    double TE_t0 = 0.0;

    double cached_pregraphics_dt = fuzz;

    // the number of materials specified by the user input
    const size_t num_mats = Materials.num_mats;

// REORIENTATION TESTING
// =======================================================================================
// DETECT REORIENTATION VALIDATION MODE FROM USER DEFINED VELOCITY BOUNDARY CONDITION
// reorientation valdation mode is a kinematics-driven verification case for fracture
// it bypasses normal SGH motion updates so that cohesive zone geometry and kinematics can be 
// presribed exactly for code verification/validation
// =======================================================================================
    bool reorientation_validation_mode = false; // reorientation mode flag
    double omega_y = 0.0; // rotation about x2
    double omega_z = 0.0; // rotation about x3
    double cz_opening_rate = 0.0; // constant opening rate for cohesive zone nodes
    double x_interface = 0.0; // x location of interface for b-side flagging of cohesive zone nodes

    // device arrays to hold reorientation parameters
    DCArrayKokkos<double> reorient_params(5, "reorient_params");
    DCArrayKokkos<int> found_reorient(1, "found_reorient");
    found_reorient.set_values(0); // initialize to zero
    reorient_params.set_values(0.0); // initialize to zero
    //found_reorient.host(0) = 0; 
    found_reorient.update_device();

    // local references for device access, so RUN kernel can read BC data
    auto bc_enums = BoundaryConditions.BoundaryConditionEnums;
    auto vel_bc_vars = BoundaryConditions.velocity_bc_global_vars;
    size_t num_bcs = BoundaryConditions.num_bcs;

    RUN({
        // search boundary conditions for user-defined velocity BC with reorientation !reo
        for (size_t bdy_set = 0; bdy_set < num_bcs; bdy_set++) {
            // check if this BC has reorientation parameters set (omega_y != 0 or omega_z != 0)
            if (bc_enums(bdy_set).BCVelocityModel
                != boundary_conditions::userDefinedVelocityBC) {
                continue;  // skip non-user-defined BCs
            }
    
            // read validation mode flag from BC parameter array
            double mode_flag = vel_bc_vars(
                bdy_set, UserDefinedVelocityBC::BCVars::reorientation_mode);
    
            // treat values > 0.5 as true
            if (mode_flag > 0.5) {
                found_reorient(0) = 1;

                // store relevant BC parameters into compact array for copying back to host
                reorient_params(0) = mode_flag;
                reorient_params(1) = vel_bc_vars(bdy_set, UserDefinedVelocityBC::BCVars::omega_y);
                reorient_params(2) = vel_bc_vars(bdy_set, UserDefinedVelocityBC::BCVars::omega_z);
                reorient_params(3) = vel_bc_vars(bdy_set, UserDefinedVelocityBC::BCVars::cz_opening_rate);
                reorient_params(4) = vel_bc_vars(bdy_set, UserDefinedVelocityBC::BCVars::x_interface);
                break;
                }  
        }
    }); // end RUN
    Kokkos::fence();

    // copying reults back to host
    found_reorient.update_host();
    reorient_params.update_host();

    // set host variables from the copied data
    reorientation_validation_mode = (found_reorient.host(0) == 1);
    if (reorientation_validation_mode) {
        omega_y = reorient_params.host(1);
        omega_z = reorient_params.host(2);
        cz_opening_rate = reorient_params.host(3);
        x_interface = reorient_params.host(4);

        // temporary prints for validation
        printf("=== REORIENTATION VALIDATION MODE ENABLED ===\n");
        printf("  omega_y         = %.10f rad/us\n", omega_y);
        printf("  omega_z         = %.10f rad/us\n", omega_z);
        printf("  cz_opening_rate = %.10e cm/us\n", cz_opening_rate);
        printf("  x_interface     = %.4f cm\n", x_interface);
        printf("==============================================\n");
    }    
   
    // allocate storage for reorientation test
    CArrayKokkos<double> initial_coords(mesh.num_nodes, 3, "initial_coords");
    CArrayKokkos<int>    cz_b_side_flag(mesh.num_nodes, "cz_b_side_flag");
    cz_b_side_flag.set_values(0);

    if (reorientation_validation_mode && doing_fracture) {
        // store initial coordinates for all nodes
        FOR_ALL(n, 0, mesh.num_nodes, {
            initial_coords(n,0) = State.node.coords(n,0);
            initial_coords(n,1) = State.node.coords(n,1);
            initial_coords(n,2) = State.node.coords(n,2);
        });
        Kokkos::fence();

        // initialize b-side flags using x_interface parameter from .yaml
        const size_t nne = mesh.num_nodes_in_elem;
        FOR_ALL(e, 0, mesh.num_elems, {
            double xc = 0.0;
            // compute element centroid x coordinate
            for (size_t a = 0; a < nne; ++a) {
                const size_t gid = mesh.nodes_in_elem(e,a);
                xc += initial_coords(gid,0);
            }
            xc /= (double)nne;

            // if the element centroid is on the B side, flag all of its nodes
            // as B side nodes for cohesive zone opening in the kinematics prescription
            if (xc > x_interface) {
                for (size_t a = 0; a < nne; ++a) {
                    const size_t gid = mesh.nodes_in_elem(e,a);
                    cz_b_side_flag(gid) = 1;
                }
            }
        });
        Kokkos::fence();
    }
// REORIENTATION TESTING

    // extensive IE
    for (size_t mat_id = 0; mat_id < num_mats; mat_id++) {

        IE_t0 += sum_domain_internal_energy(State.MaterialPoints.mass,
                                            State.MaterialPoints.sie,
                                            State.MaterialPoints.num_material_points.host(mat_id),
                                            mat_id);
    } // end loop over mat_id

    // extensive KE
    KE_t0 = sum_domain_kinetic_energy(mesh,
                                      State.node.vel,
                                      State.node.coords,
                                      State.node.mass);
    // extensive TE
    TE_t0 = IE_t0 + KE_t0;

    // domain mass for each material (they are at material points)
    double mass_domain_all_mats_t0 = 0.0;
    double mass_domain_nodes_t0    = 0.0;

    for (size_t mat_id = 0; mat_id < num_mats; mat_id++) {

        double mass_domain_mat = sum_domain_material_mass(State.MaterialPoints.mass,
                                                          State.MaterialPoints.num_material_points.host(mat_id),
                                                          mat_id);

        mass_domain_all_mats_t0 += mass_domain_mat;
        printf("material %zu mass in domain = %f \n", mat_id, mass_domain_mat);
    } // end for

    // node mass of the domain
    mass_domain_nodes_t0 = sum_domain_node_mass(mesh,
                                                State.node.coords,
                                                State.node.mass);

    printf("nodal mass domain = %f \n", mass_domain_nodes_t0);

    // a flag to exit the calculation
    size_t stop_calc = 0;

    auto time_1 = std::chrono::high_resolution_clock::now();

    // Write initial state at t=0
    printf("Writing outputs to file at %f \n", graphics_time);
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

            // save the smallest dt of all materials
            min_dt_calc = fmin(dt_mat, min_dt_calc);
        } // end for loop over all mats

        dt = min_dt_calc;  // save this dt time step

// THROTTLE OUTPUT COMMENT OUT DESCALING PRINTS
// initialize stride once using the actual dt the solver is taking
    if (!cz_stride_init) {
        cz_stride = (size_t) llround(CZ_OUT_DT / dt);  // e.g. 0.5 / 1e-8 = 50,000,000
        if (cz_stride < 1) cz_stride = 1;
        cz_next_cycle = 0;
        cz_stride_init = true;
    }
// THROTTLE OUTPUT COMMENT OUT DESCALING PRINTS

        if (cycle == 0) {
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        }
        // print time step every 10 cycles
        else if (cycle % 20 == 0) {
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
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
            
// REORIENTATIONT TESTING:
        // check if reorientation validation mode testing is on
        const bool reorient_mode = (reorientation_validation_mode && doing_fracture);

        if (reorient_mode) {
            // prescribe reorientation kinematics to all nodes (skip over normal SGH evolution)
            ReorientationKinematics::prescribe_reorientation_kinematics(
                mesh, State,
                initial_coords,
                cz_b_side_flag,
                time_value,
                dt_stage,
                omega_y, omega_z,
                cz_opening_rate
            );
        }

        // if not in reorientation mode, proceed with normal SGH evolution
        if (!reorient_mode){         
            // ---- Calculate velocity gradient for the element ---- 
// REORIENTATIONT TESTING:
            get_velgrad(State.GaussPoints.vel_grad,
                        mesh,
                        State.node.coords,
                        State.node.vel,
                        State.GaussPoints.vol);

            set_corner_force_zero(mesh, State.corner.force);


            // ---- calculate the forces on the vertices and evolve stress (hypo model) ----
            for(size_t mat_id = 0; mat_id < num_mats; mat_id++){

                get_force(Materials,
                          mesh,
                          State.GaussPoints.vol,
                          State.GaussPoints.vel_grad,
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
                          State.MaterialPoints.volfrac,
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
        

            // call body forces routine
            if (doing_fracture) {

                // index of the fracture stress BC set
                const int fracture_bdy_set = BoundaryConditions.fracture_bc_id;
                    // proceed only if a fracture BC was actually defined
                    if (fracture_bdy_set >= 0) {

                    const auto &BC     = BoundaryConditions;
                    // number of overlapping node pairs (cohesive zone node pairs)
                    const size_t npairs = cohesive_zones_bank.overlapping_node_gids.dims(0);

                    if (npairs > 0) {
                        // sync BC data to host
                        DCArrayKokkos<double> bc_params(30, "bc_params");

                        RUN({
                            // store the cohsive zone consituitive parameters ina compact temp array
                            // cohesive zone parameters 0 through 5
                            // E_inf, a1, n_exp, u_n_star, u_t_star, num_prony_terms
                            bc_params(0) = BC.stress_bc_global_vars(fracture_bdy_set,
                                           fractureStressBC::BCVars::E_inf);
                            bc_params(1) = BC.stress_bc_global_vars(fracture_bdy_set,
                                           fractureStressBC::BCVars::a1);
                            bc_params(2) = BC.stress_bc_global_vars(fracture_bdy_set,
                                           fractureStressBC::BCVars::n_exp);
                            bc_params(3) = BC.stress_bc_global_vars(fracture_bdy_set,
                                           fractureStressBC::BCVars::u_n_star);
                            bc_params(4) = BC.stress_bc_global_vars(fracture_bdy_set,
                                           fractureStressBC::BCVars::u_t_star);
                            bc_params(5) = BC.stress_bc_global_vars(fracture_bdy_set,
                                           fractureStressBC::BCVars::num_prony_terms);
                        
                            // copy prony coefficients
                            // prony terms (E_j, tau_j pairs starting at index 6)
                            int nprony = (int)(bc_params(5) + 0.5);
                            for (int j = 0; j < nprony && j < 10; ++j) {
                                int prony_base = fractureStressBC::BCVars::prony_base + 2*j;
                                bc_params(6 + 2*j)     = BC.stress_bc_global_vars(fracture_bdy_set, prony_base);     // E_j
                                bc_params(6 + 2*j + 1) = BC.stress_bc_global_vars(fracture_bdy_set, prony_base + 1); // tau_j
                            }
                        }); // end RUN
                        Kokkos::fence();
                        bc_params.update_host();
            
                        // read in cohesive zone parameters on host
                        const double E_inf_host      = bc_params.host(0);
                        const double a1_host         = bc_params.host(1);
                        const double n_exp_host      = bc_params.host(2);
                        const double u_n_star_host   = bc_params.host(3);
                        const double u_t_star_host   = bc_params.host(4);
                        const int num_prony_terms    = static_cast<int>(bc_params.host(5) + 0.5);
                        const int width = 4 + num_prony_terms;

                        // device-accessible Prony parameters (E_j, tau_j pairs)
                        DCArrayKokkos<double> prony_params(2 * num_prony_terms, "cz_prony_params");
                        for (int j = 0; j < num_prony_terms; ++j) {
                            prony_params.host(2*j)     = bc_params.host(6 + 2*j);     // E_j
                            prony_params.host(2*j + 1) = bc_params.host(6 + 2*j + 1); // tau_j
                        }
                        prony_params.update_device();


                        // ensure persistent storage for cohesive internal vars
                        // internal_vars = persistent history carried across full tme steps
                        if (cohesive_zones_bank.internal_vars.dims(0) != npairs ||
                            cohesive_zones_bank.internal_vars.dims(1) != width) {
                            cohesive_zones_bank.internal_vars =
                                DCArrayKokkos<double>(npairs, width, "cz_internal_vars");
                            cohesive_zones_bank.internal_vars.set_values(0.0);
                        }             

                        // ensure persistent storage for cohesive delta internal vars
                        // delta_internal_vars = stage local predicted updates for current RK stage
                        if (cohesive_zones_bank.delta_internal_vars.dims(0) != npairs ||
                            cohesive_zones_bank.delta_internal_vars.dims(1) != width) {
                            cohesive_zones_bank.delta_internal_vars =
                                DCArrayKokkos<double>(npairs, width, "cz_delta_internal_vars");
                            cohesive_zones_bank.delta_internal_vars.set_values(0.0);
                        }

                        // reset delta internal vars to zero each RK stage
                        cohesive_zones_bank.delta_internal_vars.set_values(0.0);

                    
                        // update mesh.nodes_in_elem for device-side kernels
                        mesh.nodes_in_elem.update_device();

                        // cohsive zone RK-stage pipeline:
                        // 1) compute interface orientation
                        // 2) map global nodal motion to local openings
                        // 3) evaluate consituitive response for this stage
                        // 4) asssemble equivalent nodal cohesive forces
                        // 5) commit constituitive updates to internal variables at the end of the last RK stage

                        // 1) orientation (normal at t and t+dt)
                        const double tol = 1.0e-8;
                        DCArrayKokkos<double> cz_orientation(
                            npairs, 6, "cz_orientation");
                        cz_orientation.set_values(0.0);

                        cohesive_zones_bank.oriented(
                            mesh.nodes_in_elem,
                            State.node.coords, // current config
                            cohesive_zones_bank.overlapping_node_gids,
                            cohesive_zones_bank.cz_info,
                            cohesive_zones_bank.max_elem_in_cohesive_zone,
                            tol,
                            cz_orientation);

                        // 2) local openings (un_t, utan_t, un_tdt, utan_tdt)
                        DCArrayKokkos<double> local_opening(
                            npairs, 4, "cz_local_opening");
                        local_opening.set_values(0.0);

                        cohesive_zones_bank.ucmap(
                            State.node.coords,
                            State.node.vel,
                            cz_orientation,
                            cohesive_zones_bank.overlapping_node_gids,
                            dt_stage, 
                            local_opening);

                        // calling cohesive_zone_var_update with extracted parameters
                        // 3) cohesive law: update internal_vars + compute increments
                        cohesive_zones_bank.cohesive_zone_var_update(
                            local_opening,
                            dt_stage, 
                            time_value,
                            cohesive_zones_bank.overlapping_node_gids,
                            E_inf_host,
                            a1_host,
                            n_exp_host,
                            u_n_star_host,
                            u_t_star_host,
                            num_prony_terms,
                            prony_params,
                            cohesive_zones_bank.internal_vars,
                            cohesive_zones_bank.delta_internal_vars);
                        
                        // 4) nodal cohesive forces
                        CArrayKokkos<double> pair_area(
                            npairs, "cz_pair_area");
                        pair_area.set_values(0.0);

                        const size_t num_nodes = mesh.num_nodes;
                        CArrayKokkos<double> F_cz(
                            3 * num_nodes, "cz_nodal_forces");
                        F_cz.set_values(0.0);

                        cohesive_zones_bank.cohesive_zone_loads(
                            mesh.nodes_in_elem,
                            State.node.coords,
                            cohesive_zones_bank.overlapping_node_gids,
                            cz_orientation,
                            cohesive_zones_bank.cz_info,
                            cohesive_zones_bank.max_elem_in_cohesive_zone,
                            cohesive_zones_bank.internal_vars,
                            cohesive_zones_bank.delta_internal_vars,
                            pair_area,
                            F_cz
                        );
                        Kokkos::fence();
                        
// REORIENTATION TESTING
if (cz_fp &&
    rk_stage == rk_num_stages - 1 &&
    cycle == cz_next_cycle)
{
    //const double u_n_star = BC.stress_bc_global_vars(fracture_bdy_set, fractureStressBC::BCVars::u_n_star);
    //const double u_t_star = BC.stress_bc_global_vars(fracture_bdy_set, fractureStressBC::BCVars::u_t_star);
    const double u_n_star   = bc_params.host(3);
    const double u_t_star   = bc_params.host(4);

    local_opening.update_host();
    cohesive_zones_bank.overlapping_node_gids.update_host();
    cohesive_zones_bank.internal_vars.update_host();

    // time,pair,gidA,gidB,un_t,ut_t,un_tdt,ut_tdt,lambda,alpha,Tn,n0x,n0y,n0z,n1x,n1y,n1z

        size_t p = 0; // printing first cohesive zone node pair values
    //for (size_t p = 0; p < npairs; ++p) {
        const size_t gidA = cohesive_zones_bank.overlapping_node_gids.host(p,0);
        const size_t gidB = cohesive_zones_bank.overlapping_node_gids.host(p,1);

        const double un_t   = local_opening.host(p,0);
        const double ut_t   = local_opening.host(p,1);
        const double un_tdt = local_opening.host(p,2);
        const double ut_tdt = local_opening.host(p,3);

        const double lambda_t = sqrt(
            (un_t/u_n_star)*(un_t/u_n_star) + (ut_t/u_t_star)*(ut_t/u_t_star));

        const double alpha_t = cohesive_zones_bank.internal_vars.host(p,1);
        const double Tn_t    = cohesive_zones_bank.internal_vars.host(p,2);
        // fprintf(cz_fp,
        // "time,pair,gidA,gidB,un_t,ut_t,un_tdt,ut_tdt,lambda,alpha,Tn,n0x,n0y,n0z,n1x,n1y,n1z\n");
        // fprintf(cz_fp,
        //     "%.9e,%zu,%zu,%zu,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,"
        //     "%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n",
        //     time_value, p, gidA, gidB,
        //     un_t, ut_t, un_tdt, ut_tdt,
        //     lambda_t, alpha_t, std::fabs(Tn_t),
        //     cz_orientation(p,0), cz_orientation(p,1), cz_orientation(p,2),
        //     cz_orientation(p,3), cz_orientation(p,4), cz_orientation(p,5));

            fprintf(cz_fp, "%.9e, %.9e, %.9e, %.9e\n",
                    time_value, alpha_t, lambda_t, fabs(Tn_t));
    //}

    cz_next_cycle += cz_stride;
    fflush(cz_fp);
}

                        // 5) update global state: internal vars and nodal forces
                        // ensuring the internal vars are updated only at the last RK stage
                        if (rk_stage == rk_num_stages - 1){

                            // local aliases to the cohesive zone internal/delta_internal variable arrays for cleaner code in the RUN loop
                            auto cz_internal_vars = cohesive_zones_bank.internal_vars;
                            auto cz_delta_internal_vars = cohesive_zones_bank.delta_internal_vars;
                            
                            // cache number of prony terms on host so it can be captured into device kernel
                            const int npr_host = num_prony_terms;

                            // launch device kernel to update internal variables with stage increments
                            RUN({

                            // use the smaller of the two row counts as a safety gaurd in case mismatched numbers of cohesive zone node pairs
                            const size_t npairs_use = (cz_internal_vars.dims(0) < cz_delta_internal_vars.dims(0))
                                ? cz_internal_vars.dims(0) : cz_delta_internal_vars.dims(0);
                            
                            // use the smaller of the two column counts as a safety gaurd in case mismatched widths
                            const size_t width_use = (cz_internal_vars.dims(1) < cz_delta_internal_vars.dims(1))
                                ? cz_internal_vars.dims(1) : cz_delta_internal_vars.dims(1);

                            // columns 0..3 reserved for:
                            // 0: lambda_dot_t
                            // 1: alpha
                            // 2: Tn (normal traction)
                            // 3: Tt (tangential traction)
                            
                            // any columns beyond 4 are Prony-history terms
                            // max number of Prony terms supported by this array width is width_use - 4
                            const int npr_max = (width_use > 4) ? static_cast<int>(width_use - 4) : 0;

                            // use smaller of:
                            // number of Prony terms specified by the BC parameters, or
                            // the number of Prony columns that actually fit in the array
                            const int npr_use = (npr_host < npr_max) ? npr_host : npr_max;

                            // loop over all cohesive zone node pairs 
                            for (size_t i = 0; i < npairs_use; ++i) {

                                // 0: lambda_dot_t (store current rate)
                                if (width_use > 0) {
                                    cz_internal_vars(i, 0) = cz_delta_internal_vars(i, 0);
                                }

                                // 1: alpha (accumulate damage)
                                // accumulated damage increment
                                if (width_use > 1) {
                                    cz_internal_vars(i, 1) += cz_delta_internal_vars(i, 1);
                                }

                                // 2, 3: tractions at t+dt become the “current” tractions for next step
                                // update committted tractions using the stage increment
                                if (width_use > 2) {
                                    cz_internal_vars(i, 2) += cz_delta_internal_vars(i, 2);
                                }
                                if (width_use > 3) {
                                    cz_internal_vars(i, 3) += cz_delta_internal_vars(i, 3);
                                }

                                // 4..(4+num_prony_terms-1): Prony stresses at t+dt
                                for (int j = 0; j < npr_use; ++j) {
                                    const int col = 4 + j;
                                    // commit updated Prony stress/history
                                    cz_internal_vars(i, col) = cz_delta_internal_vars(i, col);
                                }
                            }
                            });
                            Kokkos::fence();
                        } // end rk_stage == rk_num_stages - 1

                            // 5) add F_cz into solver's global nodal force array State.node.force

                            // local alias to global nodal force array for cleaner code in the RUN loop
                            auto node_force = State.node.force;
                            auto F_cz_local = F_cz;

                            // launch device kernel 
                            RUN({
                                // total number of mesh node to loop over
                                const size_t nmax = num_nodes;

                                // total length of the cohesive force array
                                // size = 3*num_nodes
                                const size_t fzlen = F_cz_local.size();
                                
                                // loop over all mesh nodes
                                for (size_t n = 0; n < nmax; ++n) {

                                    // index of z component of node n in the cohesive force array
                                    const size_t idx2 = 3*n + 2;

                                    // safety gaurd: if cohesive force array smaller than expected,
                                    // skip this node to avoid out-of-bounds memory access
                                    if (idx2 >= fzlen) {
                                        continue;
                                    }

                                    // add cohesive force x component into global nodal force array
                                    node_force(n,0) += F_cz_local(3*n    );

                                    // add cohesive force y component into global nodal force array
                                    node_force(n,1) += F_cz_local(3*n + 1);

                                    // add cohesive force z component into global nodal force array
                                    node_force(n,2) += F_cz_local(3*n + 2);
                                }
                            });
                            Kokkos::fence();
                } // end for loop for npairs
            } // end gaurd if fracture_bdy_set > 0
        } // end if doing_fracture
    
        
// REORIENTATION TESTING
            // SKIP SGH SOLVER EVOLUTION (REORIENTATION TESTING MODE)
            if (reorient_mode) {
                continue; // skip update velocity, boundary vel, contact, energy, position, etc.
            }
// REORIENTATION TESTING           
            // apply contact forces to boundary patches
            if (doing_contact) 
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
                    boundary_contact_force(State, mesh, 5*dt*rk_alpha, Contact_State);
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
            boundary_contact(mesh, BoundaryConditions, State.node.vel, time_value);

            // mpi_coms();

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
                             State.MaterialPoints.volfrac,
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
            printf("Writing outputs to file at %f \n", graphics_time);
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
    } 

// end for cycle loop
    
// THROTTLE OUTPUT COMMENT OUT FOR DESCALING PRINTS 
                     
            if (cz_fp) {
                fclose(cz_fp);
                cz_fp = nullptr;
            }
// THROTTLE OUTPUT COMMENT OUT FOR DESCALING PRINTS

    auto time_2    = std::chrono::high_resolution_clock::now();
    auto calc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(time_2 - time_1).count();

    printf("\nCalculation time in seconds: %f \n", calc_time * 1e-9);

    // ---- Calculate energy tallies ----
    double IE_tend = 0.0;
    double KE_tend = 0.0;
    double TE_tend = 0.0;

    // extensive IE
    for(size_t mat_id = 0; mat_id < num_mats; mat_id++){

        IE_tend += sum_domain_internal_energy(State.MaterialPoints.mass,
                                              State.MaterialPoints.sie,
                                              State.MaterialPoints.num_material_points.host(mat_id),
                                              mat_id);
    } // end loop over mat_id

    // extensive KE
    KE_tend = sum_domain_kinetic_energy(mesh,
                                        State.node.vel,
                                        State.node.coords,
                                        State.node.mass);
    // extensive TE
    TE_tend = IE_tend + KE_tend;

    printf("Time=0:   KE = %.14f, IE = %.14f, TE = %.14f \n", KE_t0, IE_t0, TE_t0);
    printf("Time=End: KE = %.14f, IE = %.14f, TE = %.14f \n", KE_tend, IE_tend, TE_tend);
    printf("total energy change = %.15e \n\n", TE_tend - TE_t0);

    // domain mass for each material (they are at material points)
    double mass_domain_all_mats_tend = 0.0;
    double mass_domain_nodes_tend    = 0.0;

    for(size_t mat_id = 0; mat_id < num_mats; mat_id++){

        double mass_domain_mat = sum_domain_material_mass(State.MaterialPoints.mass,
                                                          State.MaterialPoints.num_material_points.host(mat_id),
                                                          mat_id);

        mass_domain_all_mats_tend += mass_domain_mat;
    } // end for

    // node mass of the domain
    mass_domain_nodes_tend = sum_domain_node_mass(mesh,
                                                  State.node.coords,
                                                  State.node.mass);

    printf("material mass conservation error = %f \n", mass_domain_all_mats_tend - mass_domain_all_mats_t0);
    printf("nodal mass conservation error = %f \n", mass_domain_nodes_tend - mass_domain_nodes_t0);
    printf("nodal and material mass error = %f \n\n", mass_domain_nodes_tend - mass_domain_all_mats_tend);
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
    const DRaggedRightArrayKokkos<double>& MaterialPoints_mass,
    const DRaggedRightArrayKokkos<double>& MaterialPoints_sie,
    const size_t num_mat_points,
    const size_t mat_id)
{
    double IE_sum = 0.0;
    double IE_loc_sum;

    // loop over the material points and tally IE
    FOR_REDUCE_SUM(matpt_lid, 0, num_mat_points, IE_loc_sum, {
        IE_loc_sum += MaterialPoints_mass(mat_id,matpt_lid) * MaterialPoints_sie(mat_id,matpt_lid);
    }, IE_sum);
    Kokkos::fence();

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
    const swage::Mesh& mesh,
    const DCArrayKokkos<double>& node_vel,
    const DCArrayKokkos<double>& node_coords,
    const DCArrayKokkos<double>& node_mass)
{
    // extensive KE
    double KE_sum = 0.0;
    double KE_loc_sum;

    FOR_REDUCE_SUM(node_gid, 0, mesh.num_nodes, KE_loc_sum, {
        double ke = 0;

        for (size_t dim = 0; dim < mesh.num_dims; dim++) {
            ke += node_vel(node_gid, dim) * node_vel(node_gid, dim); // 1/2 at end
        } // end for

        KE_loc_sum += node_mass(node_gid) * ke;

    }, KE_sum);
    Kokkos::fence();

    return 0.5 * KE_sum;
} // end function

// a function to tally the material point masses
double sum_domain_material_mass(
    const DRaggedRightArrayKokkos<double>& MaterialPoints_mass,
    const size_t num_mat_points,
    const size_t mat_id)
{
    double mass_domain = 0.0;
    double mass_loc_domain;

    FOR_REDUCE_SUM(matpt_lid, 0, num_mat_points, mass_loc_domain, {
        mass_loc_domain += MaterialPoints_mass(mat_id,matpt_lid);
    }, mass_domain);
    Kokkos::fence();

    return mass_domain;
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
double sum_domain_node_mass(const swage::Mesh& mesh,
    const DCArrayKokkos<double>& node_coords,
    const DCArrayKokkos<double>& node_mass)
{
    double mass_domain = 0.0;
    double mass_loc_domain;

    FOR_REDUCE_SUM(node_gid, 0, mesh.num_nodes, mass_loc_domain, {
        if (mesh.num_dims == 2) {
            mass_loc_domain += node_mass(node_gid) * node_coords(node_gid, 1);
        }
        else{
            mass_loc_domain += node_mass(node_gid);
        }
    }, mass_domain);
    Kokkos::fence();

    return mass_domain;
} // end function


// set the corner forces to zero
void set_corner_force_zero(const swage::Mesh& mesh,
    const DCArrayKokkos<double>& corner_force)
{
    // set corner force to zero
    FOR_ALL(corner_gid, 0, mesh.num_corners, {
        for (size_t dim = 0; dim < mesh.num_dims; dim++) {
            corner_force(corner_gid, dim) = 0.0;
        }
    }); // end parallel for corners
} // end function
