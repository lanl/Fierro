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

#ifndef BOUNDARY_DISP_FIXED_H
#define BOUNDARY_DISP_FIXED_H

#include "boundary_conditions.hpp"

struct BoundaryConditionEnums_t;

/////////////////////////////////////////////////////////////////////////////
///
/// \fn Boundary displacement is set to zero in all directions
///
/// \brief This is a function to set the displacement all directions to a value
///
/// \param Mesh object
/// \param Boundary condition enums to select options
/// \param Boundary condition global variables array
/// \param Boundary condition state variables array
/// \param Node displacement
/// \param Time of the simulation
/// \param Boundary global index for the surface node
/// \param Boundary set local id
///
/////////////////////////////////////////////////////////////////////////////
namespace FixedDisplacementBC
{
// add an enum for boundary statevars and global vars

KOKKOS_FUNCTION
static void displacement(const swage::Mesh& mesh,
        const DCArrayKokkos<BoundaryConditionEnums_t>& BoundaryConditionEnums,
        const RaggedRightArrayKokkos<double>& disp_bc_global_vars,
        const DCArrayKokkos<double>& bc_state_vars,
        const CArrayKokkos<double>& K_elem,
        const CArrayKokkos<double>& F_elem,
        const double time_value,
        const double time_start,
        const double time_end,
        const size_t bdy_node_gid,
        const size_t bdy_set)
{
    const size_t num_elems_in_node = mesh.elems_in_node.stride(bdy_node_gid);
    const size_t num_nodes_in_elem = mesh.num_nodes_in_elem;
    const size_t num_dof_in_elem   = 3 * num_nodes_in_elem;

    for (size_t elem_lid = 0; elem_lid < num_elems_in_node; elem_lid++) {
        const size_t elem_gid = mesh.elems_in_node(bdy_node_gid, elem_lid);

        // Find the local node index for bdy_node_gid within this element
        size_t local_node_lid = num_nodes_in_elem; // sentinel
        for (size_t a = 0; a < num_nodes_in_elem; a++) {
            if (mesh.nodes_in_elem(elem_gid, a) == bdy_node_gid) {
                local_node_lid = a;
                break;
            }
        }

        // Zero rows and columns for all 3 constrained DOFs of this node
        for (size_t p = 0; p < 3; p++) {
            const size_t constrained_dof = 3 * local_node_lid + p;

            for (size_t col = 0; col < num_dof_in_elem; col++)
                K_elem(elem_gid, constrained_dof, col) = 0.0;

            for (size_t row = 0; row < num_dof_in_elem; row++)
                K_elem(elem_gid, row, constrained_dof) = 0.0;
        }
    }

    return;
} // end func
} // end namespace

#endif // end Header Guard