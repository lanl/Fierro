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

/////////////////////////////////////////////////////////////////////////////
///
/// \fn build_node_lid_in_elem
///
/// \brief Caches the local index of each node inside every element that
///        contains it, so the matrix-free kernels do not have to search
///        nodes_in_elem for it on every call.
///
/////////////////////////////////////////////////////////////////////////////
void TLQS3D::build_node_lid_in_elem(
    const size_t num_nodes,
    const RaggedRightArrayKokkos <size_t>& elems_in_node,
    const size_t num_nodes_in_elem,
    const DCArrayKokkos <size_t>& nodes_in_elem,
    CArrayKokkos <size_t>& num_elems_in_node,
    RaggedRightArrayKokkos <size_t>& node_lid_in_elem
)
{
    num_elems_in_node = CArrayKokkos<size_t>(num_nodes, "num_elems_in_node");
    FOR_ALL(node_gid, 0, num_nodes, {
        num_elems_in_node(node_gid) = elems_in_node.stride(node_gid);
    });
    Kokkos::fence();

    node_lid_in_elem = RaggedRightArrayKokkos<size_t>(num_elems_in_node, "node_lid_in_elem");

    FOR_ALL(node_gid, 0, num_nodes, {
        for (size_t elem_lid = 0; elem_lid < elems_in_node.stride(node_gid); elem_lid++) {
            const size_t elem_gid = elems_in_node(node_gid, elem_lid);

            // Find local index of this node within the element
            size_t local_node_lid = num_nodes_in_elem; // sentinel
            for (size_t a = 0; a < num_nodes_in_elem; a++) {
                if (nodes_in_elem(elem_gid, a) == node_gid) {
                    local_node_lid = a;
                    break;
                }
            }
            node_lid_in_elem(node_gid, elem_lid) = local_node_lid;
        }
    });
    Kokkos::fence();
} // end build_node_lid_in_elem

void TLQS3D::get_r0(
    const size_t num_nodes,
    const RaggedRightArrayKokkos <size_t>& elems_in_node,
    const RaggedRightArrayKokkos <size_t>& node_lid_in_elem,
    const size_t num_nodes_in_elem,
    const DCArrayKokkos <size_t>& nodes_in_elem,
    const CArrayKokkos <double>& F_elem,
    const CArrayKokkos <double>& K_elem,
    const CArrayKokkos <double>& displacement_iter,
    MPICArrayKokkos <double>& r0
)
{
    const size_t num_dof_in_elem = 3 * num_nodes_in_elem;

    // getting r0 = (02F - 01F) - K * displacement_iter
    // the 3 dofs of the node are gathered together (displacement_iter is read once for all 3),
    // each one is still summed in the same order as the original per-dof loop
    FOR_ALL(node_gid, 0, num_nodes, {
        const size_t num_elems_in_node = elems_in_node.stride(node_gid);

        double val0 = 0.0;
        double val1 = 0.0;
        double val2 = 0.0;

        // Sum contributions from all elements containing this node
        for (size_t elem_lid = 0; elem_lid < num_elems_in_node; elem_lid++) {
            const size_t elem_gid  = elems_in_node(node_gid, elem_lid);
            const size_t local_dof = 3 * node_lid_in_elem(node_gid, elem_lid);

            // F_elem contribution
            val0 += F_elem(elem_gid, local_dof);
            val1 += F_elem(elem_gid, local_dof + 1);
            val2 += F_elem(elem_gid, local_dof + 2);

            const double* K_row0 = &K_elem(elem_gid, local_dof, 0);
            const double* K_row1 = K_row0 + num_dof_in_elem;
            const double* K_row2 = K_row1 + num_dof_in_elem;

            // Subtract K_elem * displacement_iter
            for (size_t b = 0; b < num_nodes_in_elem; b++) {
                const size_t node_gid_b = nodes_in_elem(elem_gid, b);
                const double x0 = displacement_iter(node_gid_b, 0);
                const double x1 = displacement_iter(node_gid_b, 1);
                const double x2 = displacement_iter(node_gid_b, 2);

                val0 -= K_row0[3*b] * x0;
                val0 -= K_row0[3*b+1] * x1;
                val0 -= K_row0[3*b+2] * x2;

                val1 -= K_row1[3*b] * x0;
                val1 -= K_row1[3*b+1] * x1;
                val1 -= K_row1[3*b+2] * x2;

                val2 -= K_row2[3*b] * x0;
                val2 -= K_row2[3*b+1] * x1;
                val2 -= K_row2[3*b+2] * x2;
            }
        }

        r0(node_gid, 0) = val0;
        r0(node_gid, 1) = val1;
        r0(node_gid, 2) = val2;
    });
    Kokkos::fence();
    r0.communicate();
} // end get_r0

double TLQS3D::get_alpha(
    const size_t num_nodes,
    const size_t num_nodes_in_elem,
    const size_t num_owned_nodes,
    const RaggedRightArrayKokkos<size_t>& elems_in_node,
    const RaggedRightArrayKokkos<size_t>& node_lid_in_elem,
    const DCArrayKokkos<size_t>& nodes_in_elem,
    const CArrayKokkos<double>& K_elem,
    const double rktrk,
    const CArrayKokkos<double>& p,
    MPICArrayKokkos<double>& temporary,
    const DCArrayKokkos<bool> shared_tally_owned_nodes
    )
{
    // Kernel 1: compute temporary = K * p
    FOR_ALL(node_gid, 0, num_nodes, {
        double val0 = 0.0;
        double val1 = 0.0;
        double val2 = 0.0;
        node_K_times_x(node_gid, elems_in_node, node_lid_in_elem, num_nodes_in_elem, nodes_in_elem, K_elem, p, val0, val1, val2);
        temporary(node_gid, 0) = val0;
        temporary(node_gid, 1) = val1;
        temporary(node_gid, 2) = val2;
    });
    MATAR_FENCE();

    // temporary ghost entries needed for dot product
    temporary.communicate();

    // Kernel 2: p^T * temporary over owned nodes only, then Allreduce
    double ptkp = 0.0;
    double loc_ptkp = 0.0;
    FOR_REDUCE_SUM(node_gid, 0, (int)num_owned_nodes, loc_ptkp, {
        if(shared_tally_owned_nodes(node_gid)){
            for (int j = 0; j < 3; j++) {
                loc_ptkp += p(node_gid, j) * temporary(node_gid, j);
            }
        }
    }, ptkp);
    MATAR_FENCE();

    MPI_Allreduce(MPI_IN_PLACE, &ptkp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return rktrk / (ptkp + 1e-16);
} // end get_alpha

void TLQS3D::get_rkp1(
    const size_t num_nodes,
    const MPICArrayKokkos<double>& rk,
    const MPICArrayKokkos<double>& Kp,
    const double alpha,
    MPICArrayKokkos<double>& rkp1)
{
    // r_{k+1} = r_k - alpha * K * p
    // K * p was already computed (and communicated) by get_alpha and p has not changed since,
    // so this is now a vector update instead of a second matrix-free product. Kp_val is the same
    // value the original recomputed, so rkp1 is unchanged.
    FOR_ALL(node_gid, 0, num_nodes, {
        for (size_t p_dir = 0; p_dir < 3; p_dir++) {
            const double Kp_val = Kp(node_gid, p_dir);
            rkp1(node_gid, p_dir) = rk(node_gid, p_dir) - alpha * Kp_val;
        }
    });
    Kokkos::fence();
    rkp1.communicate();
} // end get_rkp1
