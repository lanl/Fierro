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
#include <cmath>
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
// REORIENTATION TESTING
// rotation matrix about y axis (x2)
KOKKOS_FUNCTION
void rotation_matrix_y(double angle, double R[3][3]) {
    double c = cos(angle);
    double s = sin(angle);
    R[0][0] = c;   R[0][1] = 0.0; R[0][2] = s;
    R[1][0] = 0.0; R[1][1] = 1.0; R[1][2] = 0.0;
    R[2][0] = -s;  R[2][1] = 0.0; R[2][2] = c;
}

// rotation matrix about z axis (x3)
KOKKOS_FUNCTION
void rotation_matrix_z(double angle, double R[3][3]) {
    double c = cos(angle);
    double s = sin(angle);
    R[0][0] = c;   R[0][1] = -s;  R[0][2] = 0.0;
    R[1][0] = s;   R[1][1] = c;   R[1][2] = 0.0;
    R[2][0] = 0.0; R[2][1] = 0.0; R[2][2] = 1.0;
}

// multiply two 3x3 matrices: C = A * B
KOKKOS_FUNCTION
void mat_mult_3x3(const double A[3][3], const double B[3][3], double C[3][3]) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            C[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// apply rotation matrix to a vector: v_out = R * v_in
KOKKOS_FUNCTION
void rotate_vector(const double R[3][3], const double v_in[3], double v_out[3]) {
    for (int i = 0; i < 3; ++i) {
        v_out[i] = 0.0;
        for (int j = 0; j < 3; ++j) {
            v_out[i] += R[i][j] * v_in[j];
        }
    }
}

// PRESCRIBE NODAL POSITIONS FOR REORIENTATION VALIDATION
// this function overrides the solver's position update with prescribed kinematics
void prescribe_reorientation_kinematics(
    Mesh_t& mesh,
    State_t& State,
    const CArrayKokkos<double>& initial_coords, // (num_nodes,3)
    const CArrayKokkos<int>&    cz_b_side_flag, // (num_nodes) 0/1
    double time_value,
    double dt_stage,
    double omega_y,
    double omega_z,
    double cz_opening_rate
) {
    const double angle_y_t   = omega_y * time_value;
    const double angle_z_t   = omega_z * time_value;
    const double angle_y_tdt = omega_y * (time_value + dt_stage);
    const double angle_z_tdt = omega_z * (time_value + dt_stage);

    const double s_t   = cz_opening_rate * time_value;
    const double s_tdt = cz_opening_rate * (time_value + dt_stage);

    double Ry_t[3][3], Rz_t[3][3], R_t[3][3];
    double Ry_tdt[3][3], Rz_tdt[3][3], R_tdt[3][3];

    rotation_matrix_y(angle_y_t,   Ry_t);
    rotation_matrix_z(angle_z_t,   Rz_t);
    mat_mult_3x3(Rz_t, Ry_t, R_t);

    rotation_matrix_y(angle_y_tdt, Ry_tdt);
    rotation_matrix_z(angle_z_tdt, Rz_tdt);
    mat_mult_3x3(Rz_tdt, Ry_tdt, R_tdt);

    // n(t) = R(t)*[1,0,0] => column 0
    const double n_t[3]   = { R_t[0][0],   R_t[1][0],   R_t[2][0]   };
    const double n_tdt[3] = { R_tdt[0][0], R_tdt[1][0], R_tdt[2][0] };

    FOR_ALL(node_gid, 0, mesh.num_nodes, {

        const double Xx = initial_coords(node_gid,0);
        const double Xy = initial_coords(node_gid,1);
        const double Xz = initial_coords(node_gid,2);

        // active rotation of position vector X0 (rotation about the origin 0,0,0)
        // x(t) = R(t)*X0
        double xt   = R_t[0][0]*Xx + R_t[0][1]*Xy + R_t[0][2]*Xz;
        double yt   = R_t[1][0]*Xx + R_t[1][1]*Xy + R_t[1][2]*Xz;
        double zt   = R_t[2][0]*Xx + R_t[2][1]*Xy + R_t[2][2]*Xz;

        // x(t+dt_stage) = R(t+dt_stage)*X0
        double xtdt = R_tdt[0][0]*Xx + R_tdt[0][1]*Xy + R_tdt[0][2]*Xz;
        double ytdt = R_tdt[1][0]*Xx + R_tdt[1][1]*Xy + R_tdt[1][2]*Xz;
        double ztdt = R_tdt[2][0]*Xx + R_tdt[2][1]*Xy + R_tdt[2][2]*Xz;

        if (cz_b_side_flag(node_gid)) {
            xt   += s_t   * n_t[0];   yt   += s_t   * n_t[1];   zt   += s_t   * n_t[2];
            xtdt += s_tdt * n_tdt[0]; ytdt += s_tdt * n_tdt[1]; ztdt += s_tdt * n_tdt[2];
        }

        // set coords(t) and vel consistent with dt_stage
        State.node.coords(node_gid,0) = xt;
        State.node.coords(node_gid,1) = yt;
        State.node.coords(node_gid,2) = zt;

        State.node.vel(node_gid,0) = (xtdt - xt) / dt_stage;
        State.node.vel(node_gid,1) = (ytdt - yt) / dt_stage;
        State.node.vel(node_gid,2) = (ztdt - zt) / dt_stage;
    });

    Kokkos::fence();
}

// REORIENTATION TESTING
void SGH3D::execute(SimulationParameters_t& SimulationParamaters, 
                    Material_t& Materials, 
                    BoundaryCondition_t& BoundaryConditions, 
                    Mesh_t& mesh, 
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
// =============================================================================
// REORIENTATION VALIDATION MODE FLAG
// set to true to use prescribed kinematics instead of dynamic solver
// =============================================================================
const bool REORIENTATION_VALIDATION_MODE = false;

// Validation parameters (matching Figure 11 setup [1])
const double OMEGA_Y = 0.0792860967; // full 360 deg rotation for 79.247 us //0.02747; //full 360 degree rotation for 228.764 us final time;           // rad/us rotation about y-axis
const double OMEGA_Z = 0.0792860967; // full 360 deg rotation for 79.247 us //0.02747; //full 360 degree rotation for 228.764 us final time;          // rad/us rotation about z-axis  
const double CZ_OPENING_RATE = 1e-5;  // cm/us opening rate

// store initial coordinates for prescribed motion
// always allocate; only fill/use when mode is on
CArrayKokkos<double> initial_coords(mesh.num_nodes, 3, "initial_coords");
CArrayKokkos<int>    cz_b_side_flag(mesh.num_nodes, "cz_b_side_flag");
cz_b_side_flag.set_values(0.0);
if (REORIENTATION_VALIDATION_MODE && doing_fracture) {
    FOR_ALL(n, 0, mesh.num_nodes, {
        initial_coords(n,0) = State.node.coords(n,0);
        initial_coords(n,1) = State.node.coords(n,1);
        initial_coords(n,2) = State.node.coords(n,2);
    });
    Kokkos::fence();
    // mark B-side nodes once (gidB in each pair)
    const double x_interface = 0.5; // your interface plane
    const size_t nne = mesh.num_nodes_in_elem;

    // find the right element by centroid and flag ALL its nodes
    FOR_ALL(e, 0, mesh.num_elems, {
        double xc = 0.0;
        for (size_t a = 0; a < nne; ++a) {
            const size_t gid = mesh.nodes_in_elem(e,a);
            xc += initial_coords(gid,0);
        }
        xc /= (double)nne;

        if (xc > x_interface) {
            for (size_t a = 0; a < nne; ++a) {
                const size_t gid = mesh.nodes_in_elem(e,a);
                cz_b_side_flag(gid) = 1; // OK: all writes set to 1
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

            // 2/2 add
            double dt_stage = rk_alpha * dt;
            // 2/2 add
// REORIENTATION TESTING
const bool reorient_mode = (REORIENTATION_VALIDATION_MODE && doing_fracture);
if (reorient_mode) {
    prescribe_reorientation_kinematics(
        mesh, State,
        initial_coords,
        cz_b_side_flag,
        time_value,      // current time
        dt_stage,        // stage dt
        OMEGA_Y, OMEGA_Z,
        CZ_OPENING_RATE
    );
}
// REORIENTATION TESTING            
            // ---- Calculate velocity gradient for the element ----
// REORIENTATIONT TESTING:
        if (!reorient_mode){ 
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
                
                const int fracture_bdy_set = BoundaryConditions.fracture_bc_id;
                    if (fracture_bdy_set >= 0) {

                    const auto &BC     = BoundaryConditions;
                    const size_t npairs = cohesive_zones_bank.overlapping_node_gids.dims(0);

                    if (npairs > 0) {

                        // sync RaggedRightArrayKokos to host before reading
                        //BoundaryConditions.stress_bc_global_vars.host();
                        
                        // create local copy of stress_bc_global_vars
                        auto stress_bc_vars = BoundaryConditions.stress_bc_global_vars;

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

                        // update device for State.node.coorcs and State.node.vel 
                        State.node.coords.update_device();
                        State.node.vel.update_device(); 

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

                        // extract cohesive zone parameters on HOST before calling device function
                        const double E_inf = BC.stress_bc_global_vars(fracture_bdy_set, fractureStressBC::BCVars::E_inf);
                        const double a1 = BC.stress_bc_global_vars(fracture_bdy_set, fractureStressBC::BCVars::a1);
                        const double n_exp = BC.stress_bc_global_vars(fracture_bdy_set, fractureStressBC::BCVars::n_exp);
                        const double u_n_star = BC.stress_bc_global_vars(fracture_bdy_set, fractureStressBC::BCVars::u_n_star);
                        const double u_t_star = BC.stress_bc_global_vars(fracture_bdy_set, fractureStressBC::BCVars::u_t_star);

                        // extract prony terms into a simple CArrayKokkos for passing to device function
                        DCArrayKokkos<double> prony_params(num_prony_terms, 2, "prony_params");
                        for (int j = 0; j < num_prony_terms; ++j) {
                            const int prony_base = fractureStressBC::BCVars::prony_base + 2*j;
                            prony_params.host(j, 0) = BC.stress_bc_global_vars(fracture_bdy_set, prony_base);     // E_j
                            prony_params.host(j, 1) = BC.stress_bc_global_vars(fracture_bdy_set, prony_base + 1); // tau_j
                        }
                        prony_params.update_device();

                        // calling cohesive_zone_var_update with extracted parameters
                        // 3) cohesive law: update internal_vars + compute increments
                        cohesive_zones_bank.cohesive_zone_var_update(
                            local_opening,
                            dt_stage, 
                            time_value,
                            cohesive_zones_bank.overlapping_node_gids,
                            E_inf, a1, n_exp, u_n_star, u_t_star, // pass scalars directly
                            num_prony_terms,
                            prony_params,
                            cz_internal_vars_view,
                            cz_delta_internal_vars_view);

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
// REORIENTATION TESTING
if (cz_fp &&
    rk_stage == rk_num_stages - 1 &&
    cycle == cz_next_cycle)
{
    const double u_n_star = BC.stress_bc_global_vars(fracture_bdy_set, fractureStressBC::BCVars::u_n_star);
    const double u_t_star = BC.stress_bc_global_vars(fracture_bdy_set, fractureStressBC::BCVars::u_t_star);

    // time,pair,gidA,gidB,un_t,ut_t,un_tdt,ut_tdt,lambda,alpha,Tn,n0x,n0y,n0z,n1x,n1y,n1z

        size_t p = 0; // printing first cohesive zone node pair values
    //for (size_t p = 0; p < npairs; ++p) {
        const size_t gidA = cohesive_zones_bank.overlapping_node_gids(p,0);
        const size_t gidB = cohesive_zones_bank.overlapping_node_gids(p,1);

        const double un_t   = local_opening(p,0);
        const double ut_t   = local_opening(p,1);
        const double un_tdt = local_opening(p,2);
        const double ut_tdt = local_opening(p,3);

        const double lambda_t = std::sqrt(
            (un_t/u_n_star)*(un_t/u_n_star) + (ut_t/u_t_star)*(ut_t/u_t_star));

        const double alpha_t = cz_internal_vars_view(p,1);
        const double Tn_t    = cz_internal_vars_view(p,2);
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