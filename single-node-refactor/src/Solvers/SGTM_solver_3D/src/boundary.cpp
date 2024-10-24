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

#include "sgtm_solver_3D.h"
#include "mesh.h"
#include "boundary_conditions.h"

/////////////////////////////////////////////////////////////////////////////
///
/// \fn boundary_temperature
///
/// \brief Evolves the boundary according to a give temperature
///
/// \param The simulation mesh
/// \param Boundary contain arrays of information about BCs
/// \param A view into the nodal temperature array
/// \param The current simulation time
///
/////////////////////////////////////////////////////////////////////////////
void SGTM3D::boundary_temperature(const Mesh_t& mesh,
                                  const BoundaryCondition_t& BoundaryConditions,
                                  DCArrayKokkos<double>& node_temp,
                                  const double time_value) const
{
    // Loop over boundary sets
    for (size_t bdy_set = 0; bdy_set < mesh.num_bdy_sets; bdy_set++) {
        // Loop over boundary nodes in a boundary set
        FOR_ALL(bdy_node_lid, 0, mesh.num_bdy_nodes_in_set.host(bdy_set), {
            
            // get the global index for this node on the boundary
            size_t bdy_node_gid = mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid);

            // evaluate temperature on this boundary node
            BoundaryConditions.BoundaryConditionFunctions(bdy_set).temperature(mesh,
                                                                  BoundaryConditions.BoundaryConditionEnums,
                                                                  BoundaryConditions.bc_global_vars,
                                                                  BoundaryConditions.bc_state_vars,
                                                                  node_temp,
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
/// \fn boundary_heat_flux
///
/// \brief Evolves the boundary according to a given heat flux (Q)
///
/// \param The simulation mesh
/// \param Boundary contain arrays of information about BCs
/// \param A view into the nodal temperature array
/// \param The current simulation time
///
/////////////////////////////////////////////////////////////////////////////
void SGTM3D::boundary_heat_flux(const Mesh_t& mesh,
                                  const BoundaryCondition_t& BoundaryConditions,
                                  DCArrayKokkos<double>& node_temp,
                                  const double time_value) const
{
    // Loop over boundary sets
    for (size_t bdy_set = 0; bdy_set < mesh.num_bdy_sets; bdy_set++) {
        

        size_t num_bdy_patches_in_set = 2; //mesh.bdy_patches_in_set.stride.host(bdy_set);

        std::cout<<"Num bdy patches in set "<<bdy_set<<" = "<<num_bdy_patches_in_set<<std::endl;

        // Loop over boundary nodes in a boundary set
        FOR_ALL(bdy_patch_lid, 0, num_bdy_patches_in_set, {
            
            // get the global index for this node on the boundary
            size_t bdy_node_gid = mesh.bdy_nodes_in_set(bdy_set, bdy_patch_lid);

            // // evaluate temperature on this boundary node
            // BoundaryConditions.BoundaryConditionFunctions(bdy_set).heat_flux(mesh,
            //                                                       BoundaryConditions.BoundaryConditionEnums,
            //                                                       BoundaryConditions.bc_global_vars,
            //                                                       BoundaryConditions.bc_state_vars,
            //                                                       elem_flux,
            //                                                       time_value,
            //                                                       1, // rk_stage
            //                                                       bdy_node_gid,
            //                                                       bdy_set);
        }); // end for bdy_node_lid
    } // end for bdy_set

    return;
} // end boundary_velocity function
