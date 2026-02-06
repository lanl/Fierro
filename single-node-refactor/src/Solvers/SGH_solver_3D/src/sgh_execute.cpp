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

#include "sgh_solver_3D.h"

#include "simulation_parameters.h"
#include "material.h"
#include "boundary_conditions.h"
#include "mesh.h"
#include "state.h"
#include "geometry_new.h"
#include "mesh_io.h"
#include "tipton_equilibration.hpp"
#include "fracture_stress_bc.h"


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
                    Mesh_t& mesh, 
                    State_t& State)
{

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
    CArrayKokkos <double> GaussPoint_pres(mesh.num_elems*mesh.num_leg_gauss_in_elem);
    CArrayKokkos <double> GaussPoint_pres_denominator(mesh.num_elems*mesh.num_leg_gauss_in_elem);
    CArrayKokkos <double> GaussPoint_volfrac_min(mesh.num_elems*mesh.num_leg_gauss_in_elem);
    CArrayKokkos <double> GaussPoint_volfrac_limiter(mesh.num_elems*mesh.num_leg_gauss_in_elem);
    

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
                         State.MaterialToMeshMaps.elem,
                         State.MaterialToMeshMaps.num_material_elems.host(mat_id),
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
// // THROTTLE OUTPUT COMMENT OUT
//     printf("cycle = %zu, time = %.9e, dt = %.9e\n",
//            cycle, time_value, dt);
// // THROTTLE OUTPUT COMMENT OUT
// THROTTLE OUTPUT COMMENT OUT //
//         if (cycle == 0) {
//             printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
//         }
//         // print time step every 10 cycles
//         else if (cycle % 20 == 0) {
//             printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
//         } // end if
// // THROTTLE OUTPUT COMMENT OUT //

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

// // THROTTLE OUTPUT COMMENT OUT
// // printing only every .1 us to downsample data load
//         static double next_cz_print_time = 0.0;
//         constexpr double CZ_PRINT_DT = 1.0; // amount of us (microseconds) between prints
// // THROTTLE OUTPUT COMMENT OUT

        // integrate solution forward in time
        for (size_t rk_stage = 0; rk_stage < rk_num_stages; rk_stage++) {
            // ---- RK coefficient ----
            double rk_alpha = 1.0 / ((double)rk_num_stages - (double)rk_stage);

            // 2/2 add
            double dt_stage = rk_alpha * dt;
            // 2/2 add
            // ---- Calculate velocity gradient for the element ----

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
                          State.MaterialToMeshMaps.elem,
                          State.MaterialToMeshMaps.num_material_elems.host(mat_id),
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
                                  State.MaterialToMeshMaps.elem,
                                  State.MaterialToMeshMaps.num_material_elems.host(mat_id),
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

            // call body forces routine
            
            if (doing_fracture) {

                const int fracture_bdy_set = BoundaryConditions.fracture_bc_id;
                    if (fracture_bdy_set >= 0) {

                    const auto &BC     = BoundaryConditions;
                    const size_t npairs = cohesive_zones_bank.overlapping_node_gids.dims(0);

                    if (npairs > 0) {
                        const int num_prony_terms =
                            static_cast<int>(BC.stress_bc_global_vars(
                                    fracture_bdy_set,
                                    fractureStressBC::BCVars::num_prony_terms) + 0.5);
                        const int width = 4 + num_prony_terms;

                        // ensure persistent storage for cohesive internal vars
                        if (cohesive_zones_bank.internal_vars.dims(0) != npairs ||
                            cohesive_zones_bank.internal_vars.dims(1) != width) {
                            cohesive_zones_bank.internal_vars =
                                CArrayKokkos<double>(npairs, width, "cz_internal_vars");
                            cohesive_zones_bank.internal_vars.set_values(0.0);
                        }             

                        if (cohesive_zones_bank.delta_internal_vars.dims(0) != npairs ||
                            cohesive_zones_bank.delta_internal_vars.dims(1) != width) {
                            cohesive_zones_bank.delta_internal_vars =
                                CArrayKokkos<double>(npairs, width, "cz_delta_internal_vars");
                            cohesive_zones_bank.delta_internal_vars.set_values(0.0);
                        }

                        ViewCArrayKokkos<double> cz_internal_vars_view(
                            &cohesive_zones_bank.internal_vars(0,0), npairs, width);

                        ViewCArrayKokkos<double> cz_delta_internal_vars_view(
                            &cohesive_zones_bank.delta_internal_vars(0,0), npairs, width);

                        // reset delta internal vars to zero each RK stage
                        cohesive_zones_bank.delta_internal_vars.set_values(0.0);

                        // 1) orientation (normal at t and t+dt)
                        const double tol = 1.0e-8;
                        CArrayKokkos<double> cz_orientation(
                            npairs, 6, "cz_orientation");
                        cz_orientation.set_values(0.0);

                        cohesive_zones_bank.oriented(
                            mesh,
                            State.node.coords, // current config
                            cohesive_zones_bank.overlapping_node_gids,
                            cohesive_zones_bank.cz_info,
                            cohesive_zones_bank.max_elem_in_cohesive_zone,
                            tol,
                            cz_orientation);

// THROTTLE OUTPUT COMMENT OUT FOR DESCALING PRINTS            
                        // // calling debug_oriented
                        // cohesive_zones_bank.debug_oriented(
                        //     mesh,
                        //     State, 
                        //     cohesive_zones_bank.overlapping_node_gids,
                        //     cohesive_zones_bank.cz_info,
                        //     cohesive_zones_bank.max_elem_in_cohesive_zone,
                        //     tol);
// THROTTLE OUTPUT COMMENT OUT FOR DESCALING PRINTS 

                        // 2) local openings (un_t, utan_t, un_tdt, utan_tdt)
                        CArrayKokkos<double> local_opening(
                            npairs, 4, "cz_local_opening");
                        local_opening.set_values(0.0);

                        cohesive_zones_bank.ucmap(
                            State.node.coords,
                            State.node.vel,
                            cz_orientation,
                            cohesive_zones_bank.overlapping_node_gids,
                            //dt, // 2/2 comment out
                            dt_stage, // 2/2 add
                            local_opening);

// THROTTLE OUTPUT COMMENT OUT FOR DESCALING PRINTS                             
                        // // calling debug_ucmap
                        // cohesive_zones_bank.debug_ucmap(
                        //     State.node.coords,
                        //     State.node.vel,
                        //     dt,
                        //     cz_orientation,
                        //     cohesive_zones_bank.overlapping_node_gids,
                        //     local_opening);     
// THROTTLE OUTPUT COMMENT OUT FOR DESCALING PRINTS 

                        // 3) cohesive law: update internal_vars + compute increments
                        cohesive_zones_bank.cohesive_zone_var_update(
                            local_opening,
                            //dt, // 2/2 comment out
                            dt_stage, // 2/2 add 
                            time_value,
                            cohesive_zones_bank.overlapping_node_gids,
                            BC.stress_bc_global_vars,
                            fracture_bdy_set,
                            cz_internal_vars_view,
                            cz_delta_internal_vars_view);

// // THROTTLE OUTPUT COMMENT OUT FOR DESCALING PRINTS                             
//                         // comment out 1-12-2026, for downsampling prints of time, lambda, alpha, traction
//                         // remove comment after done downsampling
//                         // calling debug_cohesive_zone_var_update
//                         cohesive_zones_bank.debug_cohesive_zone_var_update(
//                            local_opening,
//                            dt,
//                            time_value,
//                            cohesive_zones_bank.overlapping_node_gids,
//                            BC.stress_bc_global_vars,
//                            fracture_bdy_set,
//                            cz_internal_vars_view,
//                            cz_delta_internal_vars_view);      
// // THROTTLE OUTPUT COMMENT OUT FOR DESCALING PRINTS 

// // THROTTLE OUTPUT COMMENT OUT
// // Print only once per timestep (last RK stage) and only every CZ_PRINT_DT in time
//                         if (rk_stage == rk_num_stages - 1 &&
//                             (time_value + 0.5*dt >= next_cz_print_time || cycle == 0)) {

//                             cohesive_zones_bank.debug_cohesive_zone_var_update(
//                                 local_opening,
//                                 dt,
//                                 time_value,
//                                 cohesive_zones_bank.overlapping_node_gids,
//                                 BC.stress_bc_global_vars,
//                                 fracture_bdy_set,
//                                 cz_internal_vars_view,
//                                 cz_delta_internal_vars_view);

//                             next_cz_print_time += CZ_PRINT_DT;
//                         }
// // THROTTLE OUTPUT COMMENT OUT

                        // 4) nodal cohesive forces
                        CArrayKokkos<double> pair_area(
                            npairs, "cz_pair_area");
                        pair_area.set_values(0.0);

                        const size_t num_nodes = mesh.num_nodes;
                        CArrayKokkos<double> F_cz(
                            3 * num_nodes, "cz_nodal_forces");
                        F_cz.set_values(0.0);
                        ViewCArrayKokkos<double> F_cz_view(&F_cz(0), 3 * num_nodes);

                        cohesive_zones_bank.cohesive_zone_loads(
                            mesh,
                            State.node.coords,
                            cohesive_zones_bank.overlapping_node_gids,
                            cz_orientation,
                            cohesive_zones_bank.cz_info,
                            cohesive_zones_bank.max_elem_in_cohesive_zone,
                            cz_internal_vars_view,
                            cz_delta_internal_vars_view,
                            pair_area,
                            F_cz_view
                        );

// THROTTLE OUTPUT COMMENT OUT FOR DESCALING PRINTS  //2/3 add
                        if (cz_fp && rk_stage == rk_num_stages - 1 && cycle == cz_next_cycle) {
                            const double u_n_star = BC.stress_bc_global_vars(fracture_bdy_set, fractureStressBC::BCVars::u_n_star);
                            const double u_t_star = BC.stress_bc_global_vars(fracture_bdy_set, fractureStressBC::BCVars::u_t_star);
                            const double E_inf = BC.stress_bc_global_vars(fracture_bdy_set, fractureStressBC::BCVars::E_inf);
                            
                            // fprintf(cz_fp, "\n=== TIME %.9e, cycle %zu ===\n", time_value, cycle);
                            // fprintf(cz_fp, "dt=%.9e dt_stage=%.9e\n", dt, dt_stage);
                            
                            // Loop over ALL overlapping node pairs
                            // for (size_t p = 0; p < npairs; p++) { // printing all cohesive zone node pairs
                                const size_t p = 0; // printing one cohesive zone node pair (first)

                                const size_t gidA = cohesive_zones_bank.overlapping_node_gids(p, 0);
                                const size_t gidB = cohesive_zones_bank.overlapping_node_gids(p, 1);
                                
                                // Local openings (raw)
                                const double un_t = local_opening(p, 0);
                                const double ut_t = local_opening(p, 1);
                                const double un_tdt = local_opening(p, 2);
                                const double ut_tdt = local_opening(p, 3);
                              
                                // // Lambda values
                                const double lambda_t = sqrt((un_t/u_n_star)*(un_t/u_n_star) + (ut_t/u_t_star)*(ut_t/u_t_star));
                                const double lambda_tdt = sqrt((un_tdt/u_n_star)*(un_tdt/u_n_star) + (ut_tdt/u_t_star)*(ut_tdt/u_t_star));
                                const double lambda_dot_t = (lambda_tdt - lambda_t) / dt_stage;
                                
                                // Internal vars
                                const double alpha_t = cz_internal_vars_view(p, 1);
                                const double Tn_t = cz_internal_vars_view(p, 2);
                                // const double Tt_t = cz_internal_vars_view(p, 3);
                                
                                // // Delta internal vars
                                // const double delta_lambda_dot = cz_delta_internal_vars_view(p, 0);
                                // const double delta_alpha = cz_delta_internal_vars_view(p, 1);
                                // const double delta_Tn = cz_delta_internal_vars_view(p, 2);
                                // const double delta_Tt = cz_delta_internal_vars_view(p, 3);
                                
                                // // Positions and velocities of the pair
                                // const double velA_x = State.node.vel(gidA, 0);
                                // const double velA_y = State.node.vel(gidA, 1);
                                // const double velA_z = State.node.vel(gidA, 2);
                                // const double velB_x = State.node.vel(gidB, 0);
                                // const double velB_y = State.node.vel(gidB, 1);
                                // const double velB_z = State.node.vel(gidB, 2);
                                
                                // // Relative velocity magnitude
                                // const double v_rel_x = velB_x - velA_x;
                                // const double v_rel_y = velB_y - velA_y;
                                // const double v_rel_z = velB_z - velA_z;
                                // const double v_rel_mag = sqrt(v_rel_x*v_rel_x + v_rel_y*v_rel_y + v_rel_z*v_rel_z);
                                
                                // // Cohesive force on this pair
                                // const double Fcz_Ax = F_cz(3*gidA);
                                // const double Fcz_Ay = F_cz(3*gidA + 1);
                                // const double Fcz_Az = F_cz(3*gidA + 2);
                                // const double Fcz_mag = sqrt(Fcz_Ax*Fcz_Ax + Fcz_Ay*Fcz_Ay + Fcz_Az*Fcz_Az);
                                
                                // // Effective area
                                // const double area = pair_area(p);

                                // // Cohesive normals/orientation (print both sets if stored as 6 components)
                                // const double n0x = cz_orientation(p, 0);
                                // const double n0y = cz_orientation(p, 1);
                                // const double n0z = cz_orientation(p, 2);

                                // const double n1x = cz_orientation(p, 3);
                                // const double n1y = cz_orientation(p, 4);
                                // const double n1z = cz_orientation(p, 5);

                                // const double n0mag = sqrt(n0x*n0x + n0y*n0y + n0z*n0z);
                                // const double n1mag = sqrt(n1x*n1x + n1y*n1y + n1z*n1z);

                                // const double dx = State.node.coords(gidB,0) - State.node.coords(gidA,0);
                                // const double dy = State.node.coords(gidB,1) - State.node.coords(gidA,1);
                                // const double dz = State.node.coords(gidB,2) - State.node.coords(gidA,2);

                                // const double udotn0 = dx*n0x + dy*n0y + dz*n0z;
                                // const double Tn_tdt = Tn_t + delta_Tn;
                                // const double Tt_tdt = Tt_t + delta_Tt;
                                
                            //     fprintf(cz_fp, "Pair %zu: gidA=%zu gidB=%zu\n", p, gidA, gidB);
                            //     fprintf(cz_fp, "  sep·n(0-2)=%.6e   sep·n(3-5)=%.6e\n",
                            //             dx*n0x + dy*n0y + dz*n0z,
                            //             dx*n1x + dy*n1y + dz*n1z);
                            //     fprintf(cz_fp, "  n(0-2)=[%.6e %.6e %.6e] |n|=%.6e\n", n0x, n0y, n0z, n0mag);
                            //     fprintf(cz_fp, "  n(3-5)=[%.6e %.6e %.6e] |n|=%.6e\n", n1x, n1y, n1z, n1mag);
                            //     fprintf(cz_fp, "  Openings: un_t=%.6e ut_t=%.6e un_tdt=%.6e ut_tdt=%.6e\n", 
                            //             un_t, ut_t, un_tdt, ut_tdt);
                            //     fprintf(cz_fp, "  Lambda     t=%.6e tdt=%.6e dot=%.6e\n", 
                            //             lambda_t, lambda_tdt, lambda_dot_t);
                            //     fprintf(cz_fp, "  Alpha: t=%.6e delta=%.6e\n", alpha_t, delta_alpha);
                            //     fprintf(cz_fp, "  Traction: Tn_t=%.6e Tt_t=%.6e dTn=%.6e dTt=%.6e\n", 
                            //             Tn_t, Tt_t, delta_Tn, delta_Tt);
                            //     fprintf(cz_fp, "  Traction NEW: Tn_tdt=%.6e Tt_tdt=%.6e\n", Tn_tdt, Tt_tdt);
                            //     fprintf(cz_fp, "  Rel vel: vx=%.6e vy=%.6e vz=%.6e |v|=%.6e\n", 
                            //             v_rel_x, v_rel_y, v_rel_z, v_rel_mag);
                            //     fprintf(cz_fp, "  F_cz on A: Fx=%.6e Fy=%.6e Fz=%.6e |F|=%.6e\n", 
                            //             Fcz_Ax, Fcz_Ay, Fcz_Az, Fcz_mag);
                            //     fprintf(cz_fp, "  Area=%.6e\n", area);
                            //     // Flag compression (negative normal opening)
                            //     if (un_t < 0.0 || un_tdt < 0.0) {
                            //         fprintf(cz_fp, "  *** COMPRESSION DETECTED: un_t=%s un_tdt=%s ***\n",
                            //                 (un_t < 0.0) ? "NEGATIVE" : "ok",
                            //                 (un_tdt < 0.0) ? "NEGATIVE" : "ok");     
                            //     }
                            //     if (udotn0 <= 0.0) fprintf(cz_fp,"  *** CLOSED/COMPRESSION (sep·n<=0) ***\n");                         
                            //     // Flag suspicious values
                            //     if (std::fabs(lambda_dot_t) > 1e3 || std::fabs(Tn_t) > 1e6 || 
                            //         std::fabs(v_rel_mag) > 1e3 || std::fabs(delta_Tn) > 1e6) {
                            //         fprintf(cz_fp, "  *** WARNING: SUSPICIOUS VALUES ***\n");
                            //     }
                            // }
                            // fprintf(cz_fp, "=== END TIME %.9e ===\n\n", time_value);
        
                            

                            fprintf(cz_fp, "%.9e, %.9e, %.9e, %.9e\n",
                                    time_value, alpha_t, lambda_t, std::fabs(Tn_t));
                            
                            //fprintf(cz_fp, "%.9e, %.9e, %.9e, %.9e\n",
                            //        un_t, ut_t, u_n_star, u_t_star);
                            
                            //fprintf(cz_fp, "RK info: rk_num_stages=%zu rk_stage=%zu rk_alpha=%.17g dt=%.9e\n",
                            //rk_num_stages, rk_stage, rk_alpha, dt);

                            // watch it live; can remove for speed
                            fflush(cz_fp);

                            cz_next_cycle += cz_stride;
                        }
// THROTTLE OUTPUT COMMENT OUT FOR DESCALING PRINTS  2/3 add

// THROTTLE OUTPUT COMMENT OUT FOR DESCALING PRINTS                              
                        // // calling debug_cohesive_zone_loads
                        // cohesive_zones_bank.debug_cohesive_zone_loads(
                        //     mesh,
                        //     State.node.coords,
                        //     cohesive_zones_bank.overlapping_node_gids,
                        //     cz_orientation,
                        //     //cohesive_zones_bank.cz_info,
                        //     //cohesive_zones_bank.max_elem_in_cohesive_zone,
                        //     cz_internal_vars_view,
                        //     cz_delta_internal_vars_view,
                        //     pair_area,
                        //     F_cz_view); 
// THROTTLE OUTPUT COMMENT OUT FOR DESCALING PRINTS  

                            // 5) update global state: internal vars and nodal forces
                        // ensuring the internal vars are updated only at the last RK stage
                        if (rk_stage == rk_num_stages - 1){

                            for (size_t i = 0; i < npairs; ++i) {
                                // 0: lambda_dot_t (store current rate)
                                cz_internal_vars_view(i, 0) = cz_delta_internal_vars_view(i, 0);

                                // 1: alpha (accumulate damage)
                                cz_internal_vars_view(i, 1) += cz_delta_internal_vars_view(i, 1);

                                // 2, 3: tractions at t+dt become the “current” tractions for next step
                                cz_internal_vars_view(i, 2) += cz_delta_internal_vars_view(i, 2);
                                cz_internal_vars_view(i, 3) += cz_delta_internal_vars_view(i, 3);

                                // 4..(4+num_prony_terms-1): Prony stresses at t+dt
                                for (int j = 0; j < num_prony_terms; ++j) {
                                    const int col = 4 + j;
                                    cz_internal_vars_view(i, col) = cz_delta_internal_vars_view(i, col);
                                }
                            }
                        } // end rk_stage == rk_num_stages - 1

                                // 5) add F_cz into global cd binnodal force vector
                            for (size_t n = 0; n < num_nodes; ++n) {
                                State.node.force(n,0) += F_cz(3*n    );
                                State.node.force(n,1) += F_cz(3*n + 1);
                                State.node.force(n,2) += F_cz(3*n + 2);
                                }
                    } // end for loop for npairs
                } // end gaurd if fracture_bdy_set > 0
            } // end if doing_fracture
    
            // // 2/3 add
            // auto dump_node = [&](size_t gid){
            //     fprintf(cz_fp, "NODE %zu: mass=%.6e  F=(%.6e,%.6e,%.6e)  v=(%.6e,%.6e,%.6e)  x=(%.6e,%.6e,%.6e)\n",
            //         gid,
            //         State.node.mass(gid),
            //         State.node.force(gid,0), State.node.force(gid,1), State.node.force(gid,2),
            //         State.node.vel(gid,0),   State.node.vel(gid,1),   State.node.vel(gid,2),
            //         State.node.coords(gid,0),State.node.coords(gid,1),State.node.coords(gid,2)
            //     );
            // };

            // dump_node(1);  dump_node(8);
            // dump_node(5);  dump_node(12);
            // dump_node(7);  dump_node(15);
            // dump_node(3);  dump_node(11);
            // // 2/3 add

            // ---- Update nodal velocities ---- //
            update_velocity(rk_alpha,
                            dt,
                            mesh,
                            State.node.vel,
                            State.node.vel_n0,
                            State.node.mass,
                            State.node.force,
                            State.corner.force);

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
                              State.MaterialToMeshMaps.elem,
                              State.MaterialToMeshMaps.num_material_elems.host(mat_id),
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
                             State.MaterialToMeshMaps.elem,
                             time_value,
                             dt,
                             rk_alpha,
                             cycle,
                             State.MaterialToMeshMaps.num_material_elems.host(mat_id),
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
// // THROTTLE OUTPUT COMMENT OUT //
//             printf("Writing outputs to file at %f \n", graphics_time);
// // THROTTLE OUTPUT COMMENT OUT //
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
    const Mesh_t& mesh,
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
double sum_domain_node_mass(const Mesh_t& mesh,
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
void set_corner_force_zero(const Mesh_t& mesh,
    const DCArrayKokkos<double>& corner_force)
{
    // set corner force to zero
    FOR_ALL(corner_gid, 0, mesh.num_corners, {
        for (size_t dim = 0; dim < mesh.num_dims; dim++) {
            corner_force(corner_gid, dim) = 0.0;
        }
    }); // end parallel for corners
} // end function

                    