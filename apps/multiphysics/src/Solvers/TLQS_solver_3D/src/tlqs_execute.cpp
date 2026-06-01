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
                    swage::Mesh& mesh, 
                    State_t& State,
                    elements::fe_ref_elem_t& ref_elem)
{
    // Conveinent local variables
    double fuzz  = SimulationParamaters.dynamic_options.fuzz;
    double tiny  = SimulationParamaters.dynamic_options.tiny;
    double small = SimulationParamaters.dynamic_options.small;

    double graphics_dt_ival  = SimulationParamaters.output_options.graphics_time_step;
    int    graphics_cyc_ival = SimulationParamaters.output_options.graphics_iteration_step;

    // double time_initial = SimulationParamaters.dynamic_options.time_initial;
    double time_final   = this->time_end;
    double dt_start = SimulationParamaters.dynamic_options.dt_start;

    int cycle_stop    = SimulationParamaters.dynamic_options.cycle_stop;

    // initialize time, time_step, and cycles
    double time_value = this->time_start;
    double dt = dt_start;

    // *******************************
    // local variables for this solver
    // *******************************
    
    // element stiffness and force arrays
    CArrayKokkos <double> K_elem(mesh.num_elems,3*mesh.num_nodes_in_elem,3*mesh.num_nodes_in_elem); /// K1 + K2
    CArrayKokkos <double> F_elem(mesh.num_elems,3*mesh.num_nodes_in_elem); /// F02 - F01

    // conjugate gradient method vectors
    CArrayKokkos <double> p(3*mesh.num_nodes);
    CArrayKokkos <double> rk(3*mesh.num_nodes);
    CArrayKokkos <double> rkp1(3*mesh.num_nodes);

    // Anderson acceleration variables
    size_t window_size = 5;
    DCArrayKokkos <double> anderson_weights(window_size);
    CArrayKokkos <double> curr_anderson_residual(3*mesh.num_nodes);
    CArrayKokkos <double> hist_anderson_residual(3*mesh.num_nodes,window_size);
    hist_anderson_residual.set_values(0);
    CArrayKokkos <double> hist_displacement_iter(3*mesh.num_nodes,window_size);
    hist_displacement_iter.set_values(0);

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
    int max_iter = 500;

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
            std::cout << "ITER: " << iter << std::endl;
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
                    for (int mat_pt = 0; mat_pt < ref_elem.gauss_point_grad_basis.dims(0); mat_pt++) {
                        // setting up view and getting material matrix
                        ViewCArrayKokkos<double> curr_grad_basis(&ref_elem.gauss_point_grad_basis(mat_pt,0,0),ref_elem.num_basis, mesh.num_dims);
                        Materials.MaterialFunctions(mat_id).fill_C_matrix(Materials.strength_global_vars, material_matrix, mat_id);

                        // tallying to element array
                        get_gradients(material_matrix, nodes_in_curr_elem, State.node.coords_t0, State.node.displacement, displacement_step, curr_grad_basis, grad_u, inv_J, det_J, PK2_curr_config);
                        double weight = ref_elem.gauss_point_weights(mat_pt)*det_J;
                        double local_mat_vol_frac = State.MaterialPoints.volfrac(mat_id, elem);

                        tally_elem_arrays(material_matrix, grad_u, inv_J, curr_grad_basis, weight, PK2_curr_config, curr_K_elem, curr_F_elem, local_mat_vol_frac);
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

            // ***************************************************
            // begin conjugate gradient solve
            // ***************************************************

            // Reset the CG solution vector each Picard iteration
            displacement_iter_kp1.set_values(0);
            rk.set_values(0);
            Kokkos::fence();

            // ***************************************************
            // build Jacobi preconditioner (assembly-free diagonal)
            // ***************************************************
            CArrayKokkos<double> M_inv(3 * mesh.num_nodes);
            FOR_ALL(node_gid, 0, mesh.num_nodes, {
                const size_t num_elems_in_node = mesh.elems_in_node.stride(node_gid);
                for (size_t p_dir = 0; p_dir < 3; p_dir++) {
                    const size_t global_dof = 3 * node_gid + p_dir;
                    double diag = 0.0;
                    for (size_t elem_lid = 0; elem_lid < num_elems_in_node; elem_lid++) {
                        const size_t elem_gid = mesh.elems_in_node(node_gid, elem_lid);
                        size_t local_node_lid = mesh.num_nodes_in_elem;
                        for (size_t a = 0; a < mesh.num_nodes_in_elem; a++) {
                            if (mesh.nodes_in_elem(elem_gid, a) == node_gid) {
                                local_node_lid = a;
                                break;
                            }
                        }
                        const size_t local_dof = 3 * local_node_lid + p_dir;
                        diag += K_elem(elem_gid, local_dof, local_dof);
                    }
                    M_inv(global_dof) = 1.0 / (diag + 1e-16);
                }
            });
            Kokkos::fence();

            // getting r0 = (02F - 01F) - K * displacement_iter_k
            get_r0(mesh.num_nodes, mesh.elems_in_node, mesh.num_nodes_in_elem, mesh.nodes_in_elem, F_elem, K_elem, displacement_iter_kp1, rk);

            // z0 = M_inv * r0,  p0 = z0
            CArrayKokkos<double> zk(3 * mesh.num_nodes);
            CArrayKokkos<double> zkp1(3 * mesh.num_nodes);
            FOR_ALL(i, 0, 3 * mesh.num_nodes, {
                zk(i) = M_inv(i) * rk(i);
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

                            // z_{k+1} = M_inv * r_{k+1}
            FOR_ALL(i, 0, 3*mesh.num_nodes, {
                zkp1(i) = M_inv(i) * rkp1(i);
            });
            Kokkos::fence();

            // check convergence on true residual norm
            double rkp1trkp1 = 0.0;
            double loc_rkp1trkp1 = 0.0;
            FOR_REDUCE_SUM(i, 0, 3*mesh.num_nodes, loc_rkp1trkp1, {
                loc_rkp1trkp1 += rkp1(i) * rkp1(i);
            }, rkp1trkp1);
            Kokkos::fence();
            double norm = sqrt(rkp1trkp1);
            //std::cout << "CGM iter " << cgm_iter << " residual norm: " << norm << "\n";
            if (norm < 1E-10) {
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

            // --- Step 4: build a temporary copy of the residual history ---
            // The copy is ordered oldest -> newest so column 0 of temp is the
            // oldest retained residual and column m_anderson-1 is f_k.
            //
            // Mapping from window slot w (0 = oldest, m_anderson-1 = newest)
            // to circular buffer column:
            //   if buffer not yet full  -> hist_col = w
            //   if buffer full          -> hist_col = (anderson_col - m_anderson + 1 + w
            //                                          + window_size) % window_size
            // Both branches reduce to the same expression when the buffer is not
            // yet full (anderson_col == iter, m_anderson == iter+1), but we keep
            // them explicit for clarity.
            //
            // We allocate a fresh CArrayKokkos each iteration so the size matches
            // m_anderson exactly.  The destructor frees the memory when we leave
            // this scope.
            CArrayKokkos<double> temp_residual_hist(3*mesh.num_nodes, m_anderson);

            for (int w = 0; w < m_anderson; w++) {
                int hist_col;
                if (iter + 1 <= (int)window_size) {
                    // Buffer not yet full: columns 0..iter are valid, in order.
                    hist_col = w;
                } else {
                    // Buffer full: walk forward from the slot after anderson_col
                    // (= oldest) to anderson_col itself (= newest).
                    hist_col = (anderson_col - m_anderson + 1 + w + (int)window_size) % (int)window_size;
                }

                FOR_ALL(i, 0, 3*mesh.num_nodes, {
                    temp_residual_hist(i, w) = hist_anderson_residual(i, hist_col);
                });
                Kokkos::fence();
            }

            // --- Step 5: solve for Anderson mixing weights via QR ---
            // qr_solve() finds theta such that:
            //   min_theta  || temp_residual_hist * theta ||_2
            //   subject to    sum(theta_i) = 1
            // and stores the result in anderson_weights(0..m_anderson-1).
            //
            // The call is destructive with respect to temp_residual_hist; the
            // circular buffer hist_anderson_residual is untouched because we
            // passed a copy.
            anderson_weights.set_values(0);
            QR_solver(temp_residual_hist, curr_anderson_residual, anderson_weights);
            anderson_weights.update_host();

            // --- Step 6: form the Anderson-accelerated iterate ---
            //   x_{k+1} = sum_{w=0}^{m_anderson-1}  theta_w * G(x_{hist_col_w})
            //
            // anderson_weights(w) are assumed to be accessible on the host after
            // qr_solve() returns (sync/copy is the responsibility of qr_solve).
            displacement_iter_kp1.set_values(0);
            Kokkos::fence();

            for (int w = 0; w < m_anderson; w++) {
                int hist_col;
                if (iter + 1 <= (int)window_size) {
                    hist_col = w;
                } else {
                    hist_col = (anderson_col - m_anderson + 1 + w + (int)window_size) % (int)window_size;
                }

                // Read weight on host; capture by value into the kernel.
                // (Requires qr_solve to have made anderson_weights host-accessible.)
                double theta_w = anderson_weights.host(w);

                FOR_ALL(i, 0, 3*mesh.num_nodes, {
                    displacement_iter_kp1(i) += theta_w * hist_displacement_iter(i, hist_col);
                });
                Kokkos::fence();
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
                for (int mat_pt = 0; mat_pt < ref_elem.gauss_point_grad_basis.dims(0); mat_pt++) {
                    // setting up view and getting material matrix
                    ViewCArrayKokkos<double> curr_grad_basis(&ref_elem.gauss_point_grad_basis(mat_pt,0,0),ref_elem.num_basis, mesh.num_dims);
                    Materials.MaterialFunctions(mat_id).fill_C_matrix(Materials.strength_global_vars, material_matrix, mat_id);

                    // views into stress and strain
                    ViewCArrayKokkos<double> stress_view(&State.MaterialPoints.stress(mat_id, State.points_in_mat_elem(elem,mat_pt),0,0),3,3);
                    ViewCArrayKokkos<double> strain_view(&State.MaterialPoints.strain(mat_id, State.points_in_mat_elem(elem,mat_pt),0,0),3,3);

                    // tallying to element array
                    post_process(material_matrix, nodes_in_curr_elem, State.node.coords_t0, State.node.displacement, curr_grad_basis, stress_view, strain_view);

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