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

/**
 * @brief Applies a matrix-free Chebyshev polynomial preconditioner using a thread-safe,
 * node-based gathering approach (no atomic operations required).
 */
void TLQS3D::apply_chebyshev_preconditioner(const MPICArrayKokkos<double>& rk,
                                            const MPICArrayKokkos<double>& zkp1,
                                            const MPICArrayKokkos<double>& D_inv,
                                            const MPICArrayKokkos<double>& zk,
                                            const CArrayKokkos<double>& delta_z,
                                            const MPICArrayKokkos<double>& temporary,
                                            const CArrayKokkos<double>& K_elem,
                                            const size_t num_nodes,
                                            const RaggedRightArrayKokkos<size_t>& elems_in_node,
                                            const size_t num_nodes_in_elem,
                                            const DCArrayKokkos<size_t>& nodes_in_elem,
                                            const double alpha,
                                            const double beta,
                                            const int degree)
{
    const size_t total_dofs = 3 * num_nodes;
    
    // Compute Chebyshev parameters based on spectral bounds
    const double d = (beta + alpha) / 2.0;
    const double c = (beta - alpha) / 2.0;

    // --- Step 1: Initialize the 3-term recurrence (Iteration k = 0) ---
    FOR_ALL(i, 0, (int)num_nodes,
            j, 0, 3, {
        zk(i,j) = 0.0;
        delta_z(i,j) = (1.0 / d) * D_inv(i,j) * rk(i,j);
        zk(i,j) += delta_z(i,j);
    });
    MATAR_FENCE();

    double rho_prev = c / (2.0 * d);
    
    // --- Step 2: Recurrence Loop (Iteration k = 1 to degree-1) ---
    for (int k = 1; k < degree; ++k) {
        double rho_k = 1.0 / (2.0 * d / c - rho_prev);
        double gamma_k = rho_k * 2.0 / c;

        // --- THREAD-SAFE MATRIX-FREE MULTIPLICATION (temporary = K * zk) ---
        // Parallelized by node; each thread computes its exact global DOF value.
        FOR_ALL(node_gid, 0, num_nodes, {
            const size_t num_elems_in_node = elems_in_node.stride(node_gid);

            for (size_t p = 0; p < 3; p++) {
                double val = 0.0;

                // Sum contributions from all elements containing this global node
                for (size_t elem_lid = 0; elem_lid < num_elems_in_node; elem_lid++) {
                    const size_t elem_gid = elems_in_node(node_gid, elem_lid);

                    // Find local index of this node within the current element
                    size_t local_node_lid = num_nodes_in_elem; 
                    for (size_t a = 0; a < num_nodes_in_elem; a++) {
                        if (nodes_in_elem(elem_gid, a) == node_gid) {
                            local_node_lid = a;
                            break;
                        }
                    }

                    const size_t local_dof = 3 * local_node_lid + p;

                    // Contract K_elem row with the incoming Chebyshev vector (zk)
                    for (size_t b = 0; b < num_nodes_in_elem; b++) {
                        const size_t node_gid_b = nodes_in_elem(elem_gid, b);
                        
                        for (size_t q = 0; q < 3; q++) {
                            const size_t local_dof_b = 3 * b + q;
                            
                            val += K_elem(elem_gid, local_dof, local_dof_b) * zk(node_gid_b, q);
                        }
                    }
                }

                // Directly assign to scratch array without atomics or clearing passes
                temporary(node_gid, p) = val;
            }
        });
        MATAR_FENCE();
        // -----------------------------------------------------------------

        // Perform the vector updates using MATAR 1D indexing
        FOR_ALL(i, 0, (int)num_nodes,
                j, 0, 3, {
            delta_z(i,j) = rho_k * delta_z(i,j) + gamma_k * D_inv(i,j) * (rk(i,j) - temporary(i,j));
            zk(i,j) += delta_z(i,j);
        });
        MATAR_FENCE();

        rho_prev = rho_k;
    }

    // --- Step 3: Finalize Output ---
    FOR_ALL(i, 0, (int)num_nodes,
            j, 0, 3, {
        zkp1(i,j) = zk(i,j);
    });
    MATAR_FENCE();
}

/**
 * @brief Computes the inverse of the global stiffness matrix diagonal (D_inv)
 * using a thread-safe, node-based assembly approach.
 * * @param D_inv              Output vector tracking 1.0 / (K_ii + epsilon)
 * @param K_elem             Element stiffness arrays (num_elems, 3*npe, 3*npe)
 * @param num_nodes          Total number of global nodes in the mesh
 * @param elems_in_node      Ragged array mapping global node ID to local elements
 * @param num_nodes_in_elem  Number of nodes per element (e.g., 64 for cubic hex)
 * @param nodes_in_elem      Array mapping element ID to its global node IDs
 */
void TLQS3D::get_diagonal_inverse(MPICArrayKokkos<double>& D_inv,
                                  const CArrayKokkos<double>& K_elem,
                                  const size_t num_nodes,
                                  const RaggedRightArrayKokkos<size_t>& elems_in_node,
                                  const size_t num_nodes_in_elem,
                                  const DCArrayKokkos<size_t>& nodes_in_elem)
{
    FOR_ALL(node_gid, 0, num_nodes, {
        const size_t num_elems_in_node = elems_in_node.stride(node_gid);
        
        for (size_t p_dir = 0; p_dir < 3; p_dir++) {
            double diag = 0.0;
            
            for (size_t elem_lid = 0; elem_lid < num_elems_in_node; elem_lid++) {
                const size_t elem_gid = elems_in_node(node_gid, elem_lid);
                
                // Find local index of this node within the element
                size_t local_node_lid = num_nodes_in_elem;
                for (size_t a = 0; a < num_nodes_in_elem; a++) {
                    if (nodes_in_elem(elem_gid, a) == node_gid) {
                        local_node_lid = a;
                        break;
                    }
                }
                
                const size_t local_dof = 3 * local_node_lid + p_dir;
                
                // Accumulate diagonal entries (local_dof, local_dof)
                diag += K_elem(elem_gid, local_dof, local_dof);
            }
            
            // Invert with safety epsilon protection
            D_inv(node_gid, p_dir) = 1.0 / (diag + 1e-16);
        }
    });
    MATAR_FENCE();
    D_inv.communicate();
}

/**
 * @brief Computes the spectral bounds (alpha and beta) for the Chebyshev preconditioner
 * using a performance-optimized Power Iteration loop with fused kernels and proper MATAR reductions.
 */
void TLQS3D::get_chebyshev_bounds(double& alpha,
                                  double& beta,
                                  const MPICArrayKokkos<double>& D_inv,
                                  const CArrayKokkos<double>& K_elem,
                                  const size_t num_nodes,
                                  const RaggedRightArrayKokkos<size_t>& elems_in_node,
                                  const size_t num_nodes_in_elem,
                                  const DCArrayKokkos<size_t>& nodes_in_elem,
                                  MPICArrayKokkos<double>& v_scratch,
                                  MPICArrayKokkos<double>& w_scratch,
                                  const int max_iters,
                                  const int num_owned_nodes,
                                  const DCArrayKokkos<bool> shared_tally_owned_nodes)
{
    double lambda_max = 0.0;

    // Initialize initial guess vector v to 1.0
    FOR_ALL(i, 0, (int)num_nodes, 
            j, 0, 3, {
        v_scratch(i,j) = 1.0;
    });
    MATAR_FENCE();

    // --- Power Iteration Loop ---
    for (int iter = 0; iter < max_iters; ++iter) {

        // updating stale indices from previous iteration
        v_scratch.communicate();
        
        // 1. FUSED KERNEL: Compute Matrix-Free Product AND Scale by D_inv directly
        FOR_ALL(node_gid, 0, num_nodes, {
            const size_t num_elems_in_node = elems_in_node.stride(node_gid);

            for (size_t p = 0; p < 3; p++) {
                double val = 0.0;

                for (size_t elem_lid = 0; elem_lid < num_elems_in_node; elem_lid++) {
                    const size_t elem_gid = elems_in_node(node_gid, elem_lid);

                    size_t local_node_lid = num_nodes_in_elem; 
                    for (size_t a = 0; a < num_nodes_in_elem; a++) {
                        if (nodes_in_elem(elem_gid, a) == node_gid) {
                            local_node_lid = a;
                            break;
                        }
                    }

                    const size_t local_dof = 3 * local_node_lid + p;

                    for (size_t b = 0; b < num_nodes_in_elem; b++) {
                        const size_t node_gid_b = nodes_in_elem(elem_gid, b);
                        
                        for (size_t q = 0; q < 3; q++) {
                            const size_t local_dof_b = 3 * b + q;
                            
                            val += K_elem(elem_gid, local_dof, local_dof_b) * v_scratch(node_gid_b, q);
                        }
                    }
                }
                
                // Fused operation: Scale the accumulated stiffness action by D_inv 
                w_scratch(node_gid, p) = val * D_inv(node_gid, p);
            }
        });
        MATAR_FENCE();

        // 2. Compute Dot Products for Rayleigh Quotient using 6-arg FOR_REDUCE_SUM
        double v_dot_w = 0.0;
        double v_dot_v = 0.0;

        double local_v_dot_w = 0.0;
        FOR_REDUCE_SUM(i, 0, (int)num_nodes,
                       j, 0, 3, local_v_dot_w, {
            local_v_dot_w += v_scratch(i,j) * w_scratch(i,j);
        }, v_dot_w);

        double local_v_dot_v = 0.0;
        FOR_REDUCE_SUM(i, 0, (int)num_nodes, 
                       j, 0, 3, local_v_dot_v, {
            local_v_dot_v += v_scratch(i,j) * v_scratch(i,j);
        }, v_dot_v);
        MATAR_FENCE();

        // Calculate current estimate of the maximum eigenvalue
        lambda_max = v_dot_w / (v_dot_v + 1e-16);

        // 3. Normalize w vector to update our guess v: v = w / ||w||_2
        double w_norm2 = 0.0;
        double local_w_norm2 = 0.0;
        FOR_REDUCE_SUM(i, 0, (int)num_nodes, 
                       j, 0, 3, local_w_norm2, {
            local_w_norm2 += w_scratch(i,j) * w_scratch(i,j);
        }, w_norm2);
        MATAR_FENCE();
        
        double inv_norm = 1.0 / (sqrt(w_norm2) + 1e-16);

        FOR_ALL(i, 0, (int)num_nodes,
                j, 0, 3, {
            v_scratch(i,j) = w_scratch(i,j) * inv_norm;
        });
        MATAR_FENCE();
    }

    // --- 4. Apply Heuristic Bounding Box ---
    beta  = 1.1 * lambda_max; 
    alpha = beta / 10.0;      
}