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

#include "tlqs_solver_3D.hpp"

#include "simulation_parameters.hpp"
#include "material.hpp"
#include "boundary_conditions.hpp"
#include "state.hpp"
#include "geometry_new.hpp"
#include "mesh_io.hpp"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn execute
///
/// Evolve the state according to the SGH method
///
/////////////////////////////////////////////////////////////////////////////
void TLQS3D::execute(SimulationParameters_t& SimulationParamaters, 
                    Material_t& Materials, 
                    BoundaryCondition_t& BoundaryConditions, 
                    swage::Mesh_t& mesh, 
                    State_t& State)
{
    // Conveinent local variables
    double fuzz  = SimulationParamaters.DynamicOptions.fuzz;
    double tiny  = SimulationParamaters.DynamicOptions.tiny;
    double small = SimulationParamaters.DynamicOptions.small;

    double graphics_dt_ival  = SimulationParamaters.OutputOptions.graphics_time_step;
    int    graphics_cyc_ival = SimulationParamaters.OutputOptions.graphics_iteration_step;

    // double time_initial = SimulationParamaters.DynamicOptions.time_initial;
    double time_final   = this->time_end;
    double dt_start = SimulationParamaters.DynamicOptions.dt_start;

    int cycle_stop    = SimulationParamaters.DynamicOptions.cycle_stop;

    // initialize time, time_step, and cycles
    double time_value = this->time_start;
    double dt = dt_start;

    // defining displacement at start of tlqs solve
    FOR_ALL(i, 0, static_cast<long long>(mesh.num_nodes),
            j, 0, 3, {
                State.node.displacement(i,j) = State.node.coords(i,j) - State.node.coords_t0(i,j);
            });

    // *******************************
    // local variables for this solver
    // *******************************

    // initializing the reference element
    elements::Quadrature_t Quad;
    Quad.initialize_quadrature(reference_space::GaussLegendre,
                               2*SimulationParamaters.MeshInput.p_order,
                               3);
    elements::ReferenceElement_t ref_elem;
    ref_elem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                  reference_space::LagrangeLegendre,
                                  Quad,
                                  SimulationParamaters.MeshInput.p_order);
    
    // element stiffness and force arrays
    CArrayKokkos <double> K_elem(mesh.num_elems,3*mesh.num_nodes_in_elem,3*mesh.num_nodes_in_elem); /// K1 + K2
    CArrayKokkos <double> F_elem(mesh.num_elems,3*mesh.num_nodes_in_elem); /// F02 - F01

    // additive schwarz preconditioning variables
    /* CArrayKokkos <double> K_elem_inv(mesh.num_elems,3*mesh.num_nodes_in_elem,3*mesh.num_nodes_in_elem);
    CArrayKokkos <double> intermediate_K_elem_inv(mesh.num_elems,3*mesh.num_nodes_in_elem,3*mesh.num_nodes_in_elem);
    CArrayKokkos <double> temporary(mesh.num_elems, 3*mesh.num_nodes_in_elem);
    CArrayKokkos <size_t> perm(mesh.num_elems, 3*mesh.num_nodes_in_elem);
    CArrayKokkos<double> zk(3 * mesh.num_nodes);
    CArrayKokkos<double> zkp1(3 * mesh.num_nodes); */

    // conjugate gradient method vectors
    CArrayKokkos <double> p(3*mesh.num_nodes);
    CArrayKokkos <double> rk(3*mesh.num_nodes);
    CArrayKokkos <double> rkp1(3*mesh.num_nodes);

    // Anderson acceleration variables
    size_t window_size = 1;
    const size_t max_hist = (window_size > 1) ? (window_size - 1) : 1;
    DCArrayKokkos <double> anderson_weights(max_hist);
    CArrayKokkos <double> curr_anderson_residual(3*mesh.num_nodes);
    CArrayKokkos <double> hist_anderson_residual(3*mesh.num_nodes,window_size);
    hist_anderson_residual.set_values(0);
    CArrayKokkos <double> hist_displacement_iter(3*mesh.num_nodes,window_size);
    hist_displacement_iter.set_values(0);
    CArrayKokkos<double> DeltaF(3*mesh.num_nodes, max_hist);
    CArrayKokkos<double> DeltaG(3*mesh.num_nodes, max_hist);
    const int startup_iters = 3;       // Number of initial pure Picard steps
    const double max_weight = 50.0;     // Threshold to catch exploded weights
    const double fine_floor = 1e-10;     // Residual norm below which Anderson is unsafe

    // QR variables
    FArrayKokkos <double> Q(3*mesh.num_nodes, max_hist ,"Q");
    CArrayKokkos <double> R(max_hist, max_hist, "R");
    CArrayKokkos <double> y(max_hist, "y");
    DCArrayKokkos <double> v(max_hist, 3*mesh.num_nodes, "v");

    // Picard iteration vectors
    //
    // displacement_step    : current best estimate of the total displacement increment
    //                        for this load step.  Passed to get_gradients so that K/F
    //                        are evaluated at the current linearisation point.
    //                        Set (not accumulated) to displacement_iter_kp1 after each
    //                        Anderson-accelerated Picard iterate.
    //
    // displacement_iter_k  : Picard iterate coming *into* the current iteration (x_k).
    //                        Initialised to zero at the start of every load step;
    //                        updated to displacement_iter_kp1 at the end of each iter.
    //
    // displacement_iter_kp1: result of the CG solve for this Picard iteration G(x_k),
    //                        then overwritten with the Anderson-accelerated update x_{k+1}.
    //                        Reset to zero before every CG solve.
    CArrayKokkos <double> displacement_step(3*mesh.num_nodes); /// current load-step displacement estimate
    CArrayKokkos <double> displacement_iter_k(3*mesh.num_nodes);   /// x_k  (Picard iterate in)
    CArrayKokkos <double> displacement_iter_kp1(3*mesh.num_nodes); /// G(x_k) then x_{k+1}

    // variables for chebyshev smoothing
    CArrayKokkos<double> D_inv(3 * mesh.num_nodes);
    CArrayKokkos<double> zk(3 * mesh.num_nodes);
    CArrayKokkos<double> zkp1(3 * mesh.num_nodes);
    CArrayKokkos<double> delta_z(3 * mesh.num_nodes);
    CArrayKokkos<double> temporary(3 * mesh.num_nodes);


    // Algebraic Multigrid variables
    /* int max_layers = 10;
    int dof_cutoff = 3000; // NOTE: LOOK INTO THIS VALUE MORE
    double theta = 0.25;
    double omega = 2/3;
    CArrayKokkos <double> A1;
    CArrayKokkos <double> A2;
    CArrayKokkos <double> A3;
    CArrayKokkos <double> A4;
    CArrayKokkos <double> A5;
    CArrayKokkos <double> A6;
    CArrayKokkos <double> A7;
    CArrayKokkos <double> A8;
    CArrayKokkos <double> A9;
    CArrayKokkos <double> A10;
    CArrayKokkos <double> P1;
    CArrayKokkos <double> P2;
    CArrayKokkos <double> P3;
    CArrayKokkos <double> P4;
    CArrayKokkos <double> P5;
    CArrayKokkos <double> P6;
    CArrayKokkos <double> P7;
    CArrayKokkos <double> P8;
    CArrayKokkos <double> P9;
    CArrayKokkos <double> P10; */

    // Create mesh writer
    MeshWriter mesh_writer;

    // --- Graphics vars ----
    CArray<double> graphics_times = CArray<double>(20000);
    graphics_times(0) = this->time_start;
    double graphics_time = this->time_start;
    
    double cached_pregraphics_dt = fuzz;

    // the number of materials specified by the user input
    const size_t num_mats = Materials.num_mats;

    // a flag to exit the calculation
    size_t stop_calc = 0;

    auto time_1 = std::chrono::high_resolution_clock::now();

    // Write initial state at t=0
    printf("Writing outputs to file at %f \n", graphics_time);
    mesh_writer.write_mesh_Pn(mesh,
                              State,
                              SimulationParamaters,
                              dt,
                              time_value,
                              graphics_times,
                              TLQS3D_State::required_node_state,
                              TLQS3D_State::required_gauss_pt_state,
                              TLQS3D_State::required_material_pt_state,
                              this->solver_id,
                              ref_elem);
    
    graphics_time = time_value + graphics_dt_ival;

    // ******************************************************************************************************
    // setting max_iter, need to figure out either a good general number or make it a user set with a default
    // ******************************************************************************************************
    int max_iter = 5000;

    // loop over the max number of load steps
    for (size_t cycle = 0; cycle < cycle_stop; cycle++) {
        // stop calculation if flag
        if (stop_calc == 1) {
            break;
        }

        cached_pregraphics_dt = dt;
        
        // Print the initial time step and time value
        if (cycle == 0) {
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        }
        // print time step every 20 cycles
        else if (cycle % 20 == 0) {
            printf("cycle = %lu, time = %f, time step = %f \n", cycle, time_value, dt);
        } 

        // -----------------------------------------------------------
        // Reset all Picard / Anderson state for this new load step.
        // displacement_step is the running estimate of the load-step
        // displacement, so it also starts from zero.
        // -----------------------------------------------------------
        displacement_step.set_values(0);
        displacement_iter_k.set_values(0);
        displacement_iter_kp1.set_values(0);
        hist_anderson_residual.set_values(0);
        hist_displacement_iter.set_values(0);
        Kokkos::fence();

        // start Picard iteration loop
        for (int iter = 0; iter < max_iter; iter++) {
            //std::cout << "ITER: " << iter << std::endl;
            const auto start_time = std::chrono::steady_clock::now();

            // dirichlet (displacement) type
            boundary_displacement(mesh, BoundaryConditions, K_elem, F_elem, displacement_step, dt, time_value, time_start, time_end);

            auto point_A = std::chrono::steady_clock::now();
            auto elapsed_A = std::chrono::duration_cast<std::chrono::milliseconds>(point_A - start_time).count();
            //std::cout << "Time elapsed for first bc: " << elapsed_A << " ms\n";
            // ***************************************************
            // get element arrays
            // ***************************************************
            K_elem.set_values(0);
            F_elem.set_values(0);
            Kokkos::fence();

            // looping through materials
            for (int mat_id = 0; mat_id < num_mats; mat_id++) {
                
                // parallel loop through elements
                FOR_ALL(elem, 0, State.MaterialToMeshMaps.num_mat_elems.host(mat_id), {

                    // setting up views and temp memory
                    const size_t elem_id = State.MaterialToMeshMaps.elem_in_mat_elem(mat_id, elem);
                    ViewCArrayKokkos<size_t> nodes_in_curr_elem(&mesh.nodes_in_elem(elem_id,0),mesh.num_nodes_in_elem);
                    double grad_u[3][3];
                    double inv_J[3][3];
                    double det_J;
                    double PK2_curr_config[6];
                    double material_matrix[6][6];
                    ViewCArrayKokkos<double> curr_K_elem(&K_elem(elem_id,0,0),3*mesh.num_nodes_in_elem,3*mesh.num_nodes_in_elem);
                    ViewCArrayKokkos<double> curr_F_elem(&F_elem(elem_id,0),3*mesh.num_nodes_in_elem);

                    // looping through material points
                    for (int mat_pt = 0; mat_pt < ref_elem.qpt_grad_basis.dims(0); mat_pt++) {
                        // setting up view and getting material matrix
                        ViewCArrayKokkos<double> curr_grad_basis(&ref_elem.qpt_grad_basis(mat_pt,0,0),ref_elem.qpt_grad_basis.dims(1), mesh.num_dims);
                        Materials.MaterialFunctions(mat_id).fill_C_matrix(Materials.strength_global_vars, material_matrix, mat_id);

                        // tallying to element array
                        TLQS3D::get_gradients(material_matrix, nodes_in_curr_elem, State.node.coords_t0, State.node.displacement, displacement_step, curr_grad_basis, grad_u, inv_J, det_J, PK2_curr_config);
                        double weight = Quad.qpt_weights(mat_pt)*det_J;
                        double local_mat_vol_frac = State.MaterialPoints.mat_volfrac(mat_id, elem);

                        TLQS3D::tally_elem_arrays(material_matrix, grad_u, inv_J, curr_grad_basis, weight, PK2_curr_config, curr_K_elem, curr_F_elem, local_mat_vol_frac);
                    } // end mat_pt

                }); // end elem

                Kokkos::fence();

            } // end mat_id
            auto point_B = std::chrono::steady_clock::now();
            auto elapsed_B = std::chrono::duration_cast<std::chrono::milliseconds>(point_B - point_A).count();
            //std::cout << "Time elapsed to build element arrays: " << elapsed_B << " ms\n";

            // ***************************************************
            // end element arrays
            // ***************************************************

            // ***************************************************
            // apply boundary conditions
            // ***************************************************

            // neumann (traction) type

            // dirichlet (displacement) type
            boundary_displacement(mesh, BoundaryConditions, K_elem, F_elem, displacement_step, dt, time_value, time_start, time_end);
            auto point_C = std::chrono::steady_clock::now();
            auto elapsed_C = std::chrono::duration_cast<std::chrono::milliseconds>(point_C - point_B).count();
            //std::cout << "Time elapsed for second bc: " << elapsed_C << " ms\n";
            // ***************************************************
            // end boundary conditions
            // ***************************************************

            /* // get inverse of K_elem arrays for each element
            // initialize K_elem_inv to K_elem since LU functions are destructive
            FOR_ALL(i, 0, mesh.num_elems,
                    j, 0, 3*mesh.num_nodes_in_elem,
                    k, 0, 3*mesh.num_nodes_in_elem,{
                intermediate_K_elem_inv(i,j,k) = K_elem(i,j,k);
            });

            // Regularize: K_e <- K_e + eps * (trace(K_e)/n_local) * I
            // This lifts the null space without significantly affecting the deformation modes
            FOR_ALL(elem_gid, 0, mesh.num_elems, {
                const size_t n_local = 3 * mesh.num_nodes_in_elem;
                double trace = 0.0;
                for (size_t d = 0; d < n_local; d++) {
                    trace += intermediate_K_elem_inv(elem_gid, d, d);
                }
                const double reg = 1e-6 * trace / (double)n_local;
                for (size_t d = 0; d < n_local; d++) {
                    intermediate_K_elem_inv(elem_gid, d, d) += reg;
                }
            });
            Kokkos::fence();

            FOR_ALL(i, 0, mesh.num_elems,{
                ViewCArrayKokkos<double> curr_intermediate_K_elem_inv(&intermediate_K_elem_inv(i,0,0),3*mesh.num_nodes_in_elem,3*mesh.num_nodes_in_elem);
                ViewCArrayKokkos<double> curr_temporary(&temporary(i, 0), 3*mesh.num_nodes_in_elem);
                ViewCArrayKokkos<size_t> curr_perm(&perm(i, 0), 3*mesh.num_nodes_in_elem);
                ViewCArrayKokkos<double> curr_K_elem_inv(&K_elem_inv(i,0,0),3*mesh.num_nodes_in_elem,3*mesh.num_nodes_in_elem);
                int parity;
                LU_decompose(curr_intermediate_K_elem_inv, curr_perm, curr_temporary, parity);
                LU_invert(curr_intermediate_K_elem_inv, curr_perm, curr_K_elem_inv, curr_temporary);
            }); */

            // ***************************************************
            // begin conjugate gradient solve
            // ***************************************************

            // Reset the CG solution vector each Picard iteration
            displacement_iter_kp1.set_values(0);
            rk.set_values(0);
            Kokkos::fence();

            // getting inverse of the diagonal like diagonal jacobi
            get_diagonal_inverse(D_inv, K_elem, mesh.num_nodes, mesh.elems_in_node, mesh.num_nodes_in_elem, mesh.nodes_in_elem);

            // getting spectral bounds for chebyshev smoothing
            // We pass zk and temporary here safely because they are currently uninitialized 
            // scratch space before the CG loop kicks off.
            double alpha = 0.0;
            double beta = 0.0;
            const int cheb_degree = 3; // Choose your Chebyshev polynomial degree (typically 2 to 5)
            
            get_chebyshev_bounds(alpha, beta, D_inv, K_elem, mesh.num_nodes, 
                                 mesh.elems_in_node, mesh.num_nodes_in_elem, mesh.nodes_in_elem,
                                 zk, temporary, 15); // Running 15 power iterations

            // getting r0 = (02F - 01F) - K * displacement_iter_k
            get_r0(mesh.num_nodes, mesh.elems_in_node, mesh.num_nodes_in_elem, mesh.nodes_in_elem, F_elem, K_elem, displacement_iter_kp1, rk);

            // smoothing with chebyshev polynomial
            apply_chebyshev_preconditioner(rk, zk, D_inv, zk, delta_z, temporary, K_elem, 
                                           mesh.num_nodes, mesh.elems_in_node, mesh.num_nodes_in_elem, mesh.nodes_in_elem, 
                                           alpha, beta, cheb_degree);

            // z0 = M_inv * r0,  p0 = z0
            //get_z0(mesh.num_nodes, mesh.num_nodes_in_elem, mesh.nodes_in_elem, mesh.elems_in_node, K_elem, rk, zk);
            FOR_ALL(i, 0, 3 * mesh.num_nodes, {
                //zk(i) = M_inv(i) * rk(i);
                p(i)  = zk(i);
            });
            Kokkos::fence();

            // start of PCG iteration loop
            for (int cgm_iter = 0; cgm_iter < max_iter; cgm_iter++) {

                // r_k^T * z_k
                double rktzk = 0.0;
                double loc_rktzk = 0.0;
                FOR_REDUCE_SUM(i, 0, 3*mesh.num_nodes, loc_rktzk, {
                    loc_rktzk += rk(i) * zk(i);
                }, rktzk);
                Kokkos::fence();

                // alpha_k = (r_k^T * z_k) / (p_k^T * K * p_k)
                double alpha_k = get_alpha(mesh.num_nodes, mesh.num_nodes_in_elem, mesh.nodes_in_elem, K_elem, rktzk, p);

                // displacement_iter_kp1 = displacement_iter_kp1 + alpha_k * p_k
                FOR_ALL(i, 0, 3*mesh.num_nodes, {
                    displacement_iter_kp1(i) += alpha_k * p(i);
                });
                Kokkos::fence();

                // r_{k+1} = r_k - alpha_k * K * p_k
                get_rkp1(mesh.num_nodes, mesh.elems_in_node, mesh.num_nodes_in_elem, mesh.nodes_in_elem, K_elem, rk, p, alpha_k, rkp1);

                // smoothing with chebyshev polynomial
                apply_chebyshev_preconditioner(rkp1, zkp1, D_inv, zk, delta_z, temporary, K_elem, 
                                               mesh.num_nodes, mesh.elems_in_node, mesh.num_nodes_in_elem, mesh.nodes_in_elem, 
                                               alpha, beta, cheb_degree);

                // z_{k+1} = M_inv * r_{k+1}
                //get_zkp1(mesh.num_nodes, mesh.num_nodes_in_elem, mesh.nodes_in_elem, mesh.elems_in_node, K_elem, rkp1, zkp1);
                /* FOR_ALL(i, 0, 3*mesh.num_nodes, {
                    zkp1(i) = M_inv(i) * rkp1(i);
                });
                Kokkos::fence(); */

                // check convergence on true residual norm
                double rkp1trkp1 = 0.0;
                double loc_rkp1trkp1 = 0.0;
                FOR_REDUCE_SUM(i, 0, 3*mesh.num_nodes, loc_rkp1trkp1, {
                    loc_rkp1trkp1 += rkp1(i) * rkp1(i);
                }, rkp1trkp1);
                Kokkos::fence();
                double norm = sqrt(rkp1trkp1);
                //std::cout << "CGM iter " << cgm_iter << " residual norm: " << norm << "\n";
                if (norm < 1.0/*1E-10*/) {
                    break;
                }

                // r_{k+1}^T * z_{k+1}
                double rkp1tzkp1 = 0.0;
                double loc_rkp1tzkp1 = 0.0;
                FOR_REDUCE_SUM(i, 0, 3*mesh.num_nodes, loc_rkp1tzkp1, {
                    loc_rkp1tzkp1 += rkp1(i) * zkp1(i);
                }, rkp1tzkp1);
                Kokkos::fence();

                // beta_k = (r_{k+1}^T * z_{k+1}) / (r_k^T * z_k)
                double beta_k = rkp1tzkp1 / (rktzk + 1e-16);

                // p_{k+1} = z_{k+1} + beta_k * p_k
                FOR_ALL(i, 0, 3*mesh.num_nodes, {
                    p(i) = zkp1(i) + beta_k * p(i);
                });

                // update rk, zk for next iteration
                FOR_ALL(i, 0, 3*mesh.num_nodes, {
                    rk(i) = rkp1(i);
                    zk(i) = zkp1(i);
                });
                Kokkos::fence();

            } // end PCG iteration loop
            auto point_D = std::chrono::steady_clock::now();
            auto elapsed_D = std::chrono::duration_cast<std::chrono::milliseconds>(point_D - point_C).count();
            //std::cout << "Time elapsed for CGM: " << elapsed_D << " ms\n";
            // ***************************************************
            // end conjugate gradient solve
            // displacement_iter_kp1 now holds G(x_k): the Picard
            // fixed-point map applied to the current iterate x_k.
            // ***************************************************

            // ***************************************************
            // Anderson acceleration of Picard solve
            //
            // Standard AA-m formulation:
            //   f_k   = G(x_k) - x_k          (residual of fixed-point eq.)
            //   x_{k+1} = sum_i( theta_i * G(x_{hist_i}) )
            // where theta minimises ||F * theta|| s.t. sum(theta) = 1,
            // and F = [f_{oldest} ... f_k] (columns, m_anderson of them).
            //
            // The history is stored in a circular buffer of depth window_size.
            // Column overwritten each iteration = iter % window_size.
            // A temporary copy of the residual history is passed to qr_solve()
            // because QR factorisation is destructive.
            // ***************************************************

            // --- Step 1: compute Anderson residual f_k = G(x_k) - x_k ---
            FOR_ALL(i, 0, 3*mesh.num_nodes, {
                curr_anderson_residual(i) = displacement_iter_kp1(i) - displacement_iter_k(i);
            });
            Kokkos::fence();

            // Compute current residual norm for safeguarding checks
            double safe_norm = 0.0;
            double safe_loc_norm = 0.0;
            FOR_REDUCE_SUM(i, 0, 3*mesh.num_nodes, safe_loc_norm, {
                safe_loc_norm += curr_anderson_residual(i) * curr_anderson_residual(i);
            }, safe_norm);

            // --- Step 2: store f_k and G(x_k) in circular buffer ---
            // Column to overwrite is determined by the modulo of the current
            // Picard iteration index against the window size.
            int anderson_col = iter % (int)window_size;

            FOR_ALL(i, 0, 3*mesh.num_nodes, {
                hist_anderson_residual(i, anderson_col) = curr_anderson_residual(i);
                hist_displacement_iter(i, anderson_col) = displacement_iter_kp1(i);
            });
            Kokkos::fence();

            // --- Step 3: determine how many history columns to use ---
            // On early iterations the buffer is only partially filled.
            int m_anderson = std::min(iter + 1, (int)window_size);
            int m_diff = m_anderson - 1;   // number of consecutive differences available

            bool run_acceleration = (iter >= startup_iters) && 
                                    (m_diff > 0) && 
                                    (safe_norm > fine_floor);
            
            if (run_acceleration) {
                std::cout << " [Anderson] Safe to run. Switching from Picard." << std::endl;
                int oldest_col = (iter + 1 <= (int)window_size) ? 0 : (anderson_col + 1) % (int)window_size;

                for (int j = 0; j < m_diff; j++) {
                    int col_j   = (oldest_col + j) % (int)window_size;
                    int col_jp1 = (oldest_col + j + 1) % (int)window_size;
                    FOR_ALL(i, 0, 3*mesh.num_nodes, {
                        DeltaF(i, j) = hist_anderson_residual(i, col_jp1) - hist_anderson_residual(i, col_j);
                        DeltaG(i, j) = hist_displacement_iter(i, col_jp1) - hist_displacement_iter(i, col_j);
                    });
                    Kokkos::fence();
                }

                anderson_weights.set_values(0);
                QR_solver(DeltaF, curr_anderson_residual, anderson_weights, Q, R, y, v, (size_t)m_diff);   // see note below
                anderson_weights.update_host();
                
                // SAFEGUARD: Check for ill-conditioning (exploded weights)
                bool weights_are_valid = true;
                for (int a = 0; a < m_diff; a++) {
                    if (fabs(anderson_weights.host(a)) > max_weight) {
                        weights_are_valid = false;
                        break;
                    }
                }
                /* // --- Diagnostic Verification Block ---
                std::vector<double> alpha(m_diff + 1, 0.0);
                alpha[0] = anderson_weights.host(0);

                for (int i = 1; i < m_diff; i++) {
                    alpha[i] = anderson_weights.host(i) - anderson_weights.host(i - 1);
                }
                alpha[m_diff] = 1.0 - anderson_weights.host(m_diff - 1);

                double alpha_sum = 0.0;
                std::cout << "Reconstructed Alpha Weights: ";
                for (double a : alpha) {
                    std::cout << a << " ";
                    alpha_sum += a;
                }
                std::cout << "\nTotal Alpha Sum (Should be 1.0): " << alpha_sum << std::endl;
                // ------------------------------------- */

                if (weights_are_valid) {
                    // --- Step 4a: Apply Accelerated Update ---
                    FOR_ALL(i, 0, 3*mesh.num_nodes, {
                        double correction = 0.0;
                        for (int j = 0; j < m_diff; j++) {
                            correction += anderson_weights(j) * DeltaG(i, j);
                        }
                        displacement_iter_kp1(i) -= correction;
                    });
                    Kokkos::fence();
                } else {
                    // --- Step 4b: Fallback to Pure Picard (Weights Exploded) ---
                    // By doing nothing here, displacement_iter_kp1 retains its raw Picard value.
                    // We print a diagnostic tracking statement.
                    std::cout << " [Anderson] Ill-conditioned history (Weights Exploded). Falling back to Picard." << std::endl;
                }
            } else {
                // --- Step 4c: Normal Pure Picard Mode ---
                // Active during the startup phase or when the residual drops below fine_floor.
                if (safe_norm <= fine_floor && iter >= startup_iters) {
                    std::cout << " [Anderson] Residual below safety floor (" << safe_norm << "). Finalizing via Picard." << std::endl;
                }
            }

            auto point_E = std::chrono::steady_clock::now();
            auto elapsed_E = std::chrono::duration_cast<std::chrono::milliseconds>(point_E - point_D).count();
            //std::cout << "Time elapsed for anderson: " << elapsed_E << " ms\n";
            // ***************************************************
            // end Anderson acceleration
            // ***************************************************

            // update displacement step vector for convergence check and next iteration
            FOR_ALL(i, 0, 3*mesh.num_nodes, {
                displacement_step(i) += displacement_iter_kp1(i);
            });
            Kokkos::fence();

            // convergence check
            /* double norm_num = 0.0;
            double loc_norm_num = 0.0;
            FOR_REDUCE_SUM(i, 0, 3*mesh.num_nodes, loc_norm_num, {
                loc_norm_num += curr_anderson_residual(i) * curr_anderson_residual(i);
            }, norm_num);

            double norm_den = 0.0;
            double loc_norm_den = 0.0;
            FOR_REDUCE_SUM(i, 0, 3*mesh.num_nodes, loc_norm_den, {
                loc_norm_den += displacement_step(i) * displacement_step(i);
            }, norm_den);

            double norm = sqrt(norm_num / norm_den); */
            double norm = 0.0;
            double loc_norm = 0.0;
            FOR_REDUCE_SUM(i, 0, 3*mesh.num_nodes, loc_norm, {
                loc_norm += curr_anderson_residual(i) * curr_anderson_residual(i);
            }, norm);

            std::cout << "ITER: " << iter << "   ANDERSON RESIDUAL NORM: " << norm << std::endl;
            if (norm < 1E-12 && iter > 1) {
                std::cout << "PICARD CONVERGED AT ITER: " << iter+1 << std::endl;
                break;
            }

            // Update x_k <- x_{k+1} for the next Picard iteration.
            FOR_ALL(i, 0, 3*mesh.num_nodes, {
                displacement_iter_k(i) = displacement_iter_kp1(i);
            });
            Kokkos::fence();
            auto point_F = std::chrono::steady_clock::now();
            auto elapsed_F = std::chrono::duration_cast<std::chrono::milliseconds>(point_F - point_E).count();
            //std::cout << "Time elapsed for picard convergence and update: " << elapsed_F << " ms\n";
            auto elapsed_iter = std::chrono::duration_cast<std::chrono::milliseconds>(point_F - start_time).count();
            //std::cout << "Time elapsed for this full iteration: " << elapsed_iter << " ms\n";
        } // end Picard iteration loop

        // updating total displacement for next load step
        FOR_ALL(i, 0, (int)mesh.num_nodes, 
                j, 0, 3, {
                    State.node.displacement(i,j) += displacement_step(3*i + j);
                    State.node.coords(i,j) = State.node.coords_t0(i,j) + State.node.displacement(i,j);
        });
        Kokkos::fence();

        // filling in stress and strain for output
        for (int mat_id = 0; mat_id < num_mats; mat_id++) {

            FOR_ALL(elem, 0, State.MaterialToMeshMaps.num_mat_elems.host(mat_id), {

                // setting up views and temp memory
                const size_t elem_id = State.MaterialToMeshMaps.elem_in_mat_elem(mat_id, elem);

                ViewCArrayKokkos<size_t> nodes_in_curr_elem(&mesh.nodes_in_elem(elem_id,0),mesh.num_nodes_in_elem);
                double material_matrix[6][6];

                // looping through material points
                for (int mat_pt = 0; mat_pt < ref_elem.qpt_grad_basis.dims(0); mat_pt++) {
                    // setting up view and getting material matrix
                    ViewCArrayKokkos<double> curr_grad_basis(&ref_elem.qpt_grad_basis(mat_pt,0,0),ref_elem.qpt_grad_basis.dims(1), mesh.num_dims);
                    Materials.MaterialFunctions(mat_id).fill_C_matrix(Materials.strength_global_vars, material_matrix, mat_id);

                    // views into stress and strain
                    ViewCArrayKokkos<double> stress_view(&State.MaterialPoints.stress(mat_id, State.points_in_mat_elem(elem,mat_pt),0,0),3,3);
                    ViewCArrayKokkos<double> strain_view(&State.MaterialPoints.strain(mat_id, State.points_in_mat_elem(elem,mat_pt),0,0),3,3);

                    // tallying to element array
                    TLQS3D::post_process(material_matrix, nodes_in_curr_elem, State.node.coords_t0, State.node.displacement, curr_grad_basis, stress_view, strain_view);

                } // end mat_pt

            }); // end elem
            Kokkos::fence();

        } // end mat_id

        // increment the time
        time_value += dt;

        // Manage outputs
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
            mesh_writer.write_mesh_Pn(mesh,
                                   State,
                                   SimulationParamaters,
                                   dt,
                                   time_value,
                                   graphics_times,
                                   TLQS3D_State::required_node_state,
                                   TLQS3D_State::required_gauss_pt_state,
                                   TLQS3D_State::required_material_pt_state,
                                   this->solver_id,
                                   ref_elem);

            graphics_time = time_value + graphics_dt_ival;

            dt = cached_pregraphics_dt;
        } // end if

        // end of calculation
        if (time_value >= time_final) {
            break;
        }
    } // end for cycle loop

    auto time_2    = std::chrono::high_resolution_clock::now();
    auto calc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(time_2 - time_1).count();

    printf("\nCalculation time in seconds: %f \n", calc_time * 1e-9);

} // end of TLQS3D execute