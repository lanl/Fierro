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

#include "sgh_solver_rz.h"
#include "mesh.h"
#include "boundary_conditions.h"


/////////////////////////////////////////////////////////////////////////////
///
/// \fn boundary_velocity
///
/// \brief Evolves the boundary according to a give velocity
///
/// \param The simulation mesh
/// \param Boundary contain arrays of information about BCs
/// \param A view into the nodal velocity array
/// \param The current simulation time
///
/////////////////////////////////////////////////////////////////////////////
void SGHRZ::boundary_velocity_rz(const Mesh_t&      mesh,
                                 const BoundaryCondition_t& BoundaryConditions,
                                 DCArrayKokkos<double>& node_vel,
                                 const double time_value) const
{

    size_t num_vel_bdy_sets = BoundaryConditions.num_vel_bdy_sets_in_solver.host(this->solver_id); 

    // Loop over the velocity boundary sets
    for (size_t bc_lid = 0; bc_lid < num_vel_bdy_sets; bc_lid++) {
        
        size_t bdy_set = BoundaryConditions.vel_bdy_sets_in_solver.host(bc_lid);
        
        // Loop over boundary nodes in a boundary set
        FOR_ALL(bdy_node_lid, 0, mesh.num_bdy_nodes_in_set.host(bdy_set), {

            // get the global index for this node on the boundary
            size_t bdy_node_gid = mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid);

            // evaluate velocity on this boundary node
            BoundaryConditions.BoundaryConditionFunctions(bdy_set).velocity(mesh,
                                                                  BoundaryConditions.BoundaryConditionEnums,
                                                                  BoundaryConditions.velocity_bc_global_vars,
                                                                  BoundaryConditions.bc_state_vars,
                                                                  node_vel,
                                                                  time_value,
                                                                  1, // rk_stage
                                                                  bdy_node_gid,
                                                                  bdy_set);
                
        }); // end for bdy_node_lid

        
    } // end for bdy_set

    return;
} // end boundary_velocity function


/////////////////////////////////////////////////////////////////////////////
///
/// \fn boundary_velocity
///
/// \brief Evolves the boundary according to a give velocity
///
/// \param The simulation mesh
/// \param An array of BoundaryCondition_t that contain information about BCs
/// \param A view into the nodal velocity array
/// \param The current simulation time
///
/////////////////////////////////////////////////////////////////////////////
void SGHRZ::boundary_contact_rz(const Mesh_t& mesh,
                                const BoundaryCondition_t& BoundaryConditions,
                                DCArrayKokkos<double>& node_vel,
                                const double time_value) const
{
    return;
} // end boundary_contact function


