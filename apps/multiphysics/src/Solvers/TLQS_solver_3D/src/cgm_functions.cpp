/**********************************************************************************************
� 2020. Triad National Security, LLC. All rights reserved.
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

void TLQS3D::get_r0(
    const size_t num_nodes,
    const RaggedRightArrayKokkos <size_t>& elems_in_node,
    const size_t num_nodes_in_elem,
    const DCArrayKokkos <size_t>& nodes_in_elem,
    const CArrayKokkos <double>& F_elem,
    const CArrayKokkos <double>& K_elem,
    const CArrayKokkos <double>& displacement_iter,
    const CArrayKokkos <double>& r0
)
{
    // getting r0 = (02F - 01F) - K * displacement_iter
    FOR_ALL(node_gid, 0, num_nodes, {
        const size_t num_elems_in_node = elems_in_node.stride(node_gid);

        for (size_t p = 0; p < 3; p++) {
            const size_t global_dof = 3 * node_gid + p;
            double val = 0.0;

            // Sum contributions from all elements containing this node
            for (size_t elem_lid = 0; elem_lid < num_elems_in_node; elem_lid++) {
                const size_t elem_gid = elems_in_node(node_gid, elem_lid);

                // Find local index of this node within the element
                size_t local_node_lid = num_nodes_in_elem; // sentinel
                for (size_t a = 0; a < num_nodes_in_elem; a++) {
                    if (nodes_in_elem(elem_gid, a) == node_gid) {
                        local_node_lid = a;
                        break;
                    }
                }

                const size_t local_dof = 3 * local_node_lid + p;

                // F_elem contribution
                val += F_elem(elem_gid, local_dof);

                // Subtract K_elem * displacement_iter
                for (size_t b = 0; b < num_nodes_in_elem; b++) {
                    const size_t node_gid_b = nodes_in_elem(elem_gid, b);
                    for (size_t q = 0; q < 3; q++) {
                        const size_t local_dof_b = 3 * b + q;
                        const size_t global_dof_b = 3 * node_gid_b + q;
                        val -= K_elem(elem_gid, local_dof, local_dof_b) * displacement_iter(global_dof_b);
                    }
                }
            }

            r0(global_dof) = val;
        }
    });
} // end get_r0

double TLQS3D::get_alpha(
    const size_t num_nodes,
    const size_t num_nodes_in_elem,
    const DCArrayKokkos<size_t>& nodes_in_elem,
    const CArrayKokkos<double>& K_elem,
    const double rktrk,
    const CArrayKokkos<double>& p)
{

    // denominator: p^T * K * p
    // first compute Kp = K * p via assembly-free matvec
    // then dot with p
    double ptkp = 0.0;
    double loc_ptkp = 0.0;
    FOR_REDUCE_SUM(elem_gid, 0, K_elem.dims(0), loc_ptkp, {

        for (size_t a = 0; a < num_nodes_in_elem; a++) {
            const size_t node_gid_a = nodes_in_elem(elem_gid, a);
            for (size_t p_dir = 0; p_dir < 3; p_dir++) {
                const size_t local_dof_a  = 3 * a + p_dir;
                const size_t global_dof_a = 3 * node_gid_a + p_dir;

                double Kp_val = 0.0;
                for (size_t b = 0; b < num_nodes_in_elem; b++) {
                    const size_t node_gid_b = nodes_in_elem(elem_gid, b);
                    for (size_t q = 0; q < 3; q++) {
                        const size_t local_dof_b  = 3 * b + q;
                        const size_t global_dof_b = 3 * node_gid_b + q;
                        Kp_val += K_elem(elem_gid, local_dof_a, local_dof_b) * p(global_dof_b);
                    }
                }
                loc_ptkp += p(global_dof_a) * Kp_val;
            }
        }
    }, ptkp);
    //std::cout << "PTKP: " << ptkp << std::endl;

    return rktrk / (ptkp+1E-16);
} // end get_alpha

void TLQS3D::get_rkp1(
    const size_t num_nodes,
    const RaggedRightArrayKokkos<size_t>& elems_in_node,
    const size_t num_nodes_in_elem,
    const DCArrayKokkos<size_t>& nodes_in_elem,
    const CArrayKokkos<double>& K_elem,
    const CArrayKokkos<double>& rk,
    const CArrayKokkos<double>& p,
    const double alpha,
    const CArrayKokkos<double>& rkp1)
{
    // r_{k+1} = r_k - alpha * K * p
    FOR_ALL(node_gid, 0, num_nodes, {
        const size_t num_elems_in_node = elems_in_node.stride(node_gid);

        for (size_t p_dir = 0; p_dir < 3; p_dir++) {
            const size_t global_dof = 3 * node_gid + p_dir;
            double Kp_val = 0.0;

            for (size_t elem_lid = 0; elem_lid < num_elems_in_node; elem_lid++) {
                const size_t elem_gid = elems_in_node(node_gid, elem_lid);

                // find local index of this node within the element
                size_t local_node_lid = num_nodes_in_elem; // sentinel
                for (size_t a = 0; a < num_nodes_in_elem; a++) {
                    if (nodes_in_elem(elem_gid, a) == node_gid) {
                        local_node_lid = a;
                        break;
                    }
                }

                const size_t local_dof = 3 * local_node_lid + p_dir;

                for (size_t b = 0; b < num_nodes_in_elem; b++) {
                    const size_t node_gid_b = nodes_in_elem(elem_gid, b);
                    for (size_t q = 0; q < 3; q++) {
                        const size_t local_dof_b  = 3 * b + q;
                        const size_t global_dof_b = 3 * node_gid_b + q;
                        Kp_val += K_elem(elem_gid, local_dof, local_dof_b) * p(global_dof_b);
                    }
                }
            }

            rkp1(global_dof) = rk(global_dof) - alpha * Kp_val;
        }
    });
} // end get_rkp1