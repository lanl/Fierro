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

#include "sgtm_solver_3D.h"
#include "state.h"
#include "mesh.h"
#include "simulation_parameters.h"

void SGTM3D::initialize(SimulationParameters_t& SimulationParameters, 
                	   Material_t& Materials, 
                	   Mesh_t& mesh, 
                	   BoundaryCondition_t& Boundary,
                	   State_t& State) const
{
	int num_nodes = mesh.num_nodes;
    int num_gauss_pts = mesh.num_elems;
    int num_corners = mesh.num_corners;
    int num_dim = mesh.num_dims;


    State.node.initialize(num_nodes, num_dim, SGTM3D_State::required_node_state);
    State.GaussPoints.initialize(num_gauss_pts, num_dim, SGTM3D_State::required_gauss_pt_state);
    State.corner.initialize(num_corners, num_dim, SGTM3D_State::required_corner_state);
    
    // NOTE: Material points and material corners are initialize in sgh_setup after calculating the material->mesh maps

        // check that the fills specify the required nodal fields
    bool filled_nodal_state =
        check_fill_node_states(SGTM3D_State::required_fill_node_state,
                               SimulationParameters.region_setups.fill_node_states);
    
    if (filled_nodal_state == false){
        std::cout <<" Missing required nodal state in the fill instructions for the thrmex_FE solver \n";
        std::cout <<" The required nodal state: \n";
        std::cout <<"  - velocity \n";
        std::cout <<"  - temperature \n";
        throw std::runtime_error("**** Provide fill instructions for all required nodal variables ****");
    }
} // end solver initialization


void SGTM3D::initialize_material_state(SimulationParameters_t& SimulationParameters, 
                	                  Material_t& Materials, 
                	                  Mesh_t& mesh, 
                	                  BoundaryCondition_t& Boundary,
                	                  State_t& State) const
{

    const size_t num_nodes = mesh.num_nodes;
    const size_t num_dims = 3;

    const size_t num_mats = Materials.num_mats; // the number of materials on the mesh

    // -----
    //  Allocation of state must include a buffer with ALE
    // -----

    // IMPORTANT, make buffer a parser input variable
    // for ALE, add a buffer to num_elems_for_mat, like 10% of num_elems up to num_elems.
    const size_t buffer = 0; // memory buffer to push back into

    for (int mat_id = 0; mat_id < num_mats; mat_id++) {

        const size_t num_mat_pts_in_elem = mesh.num_leg_gauss_in_elem; 

        size_t num_elems_for_mat = State.MaterialToMeshMaps(mat_id).num_material_elems + buffer; // has a memory buffer for ALE

        size_t num_points_for_mat  = num_elems_for_mat * num_mat_pts_in_elem;
        size_t num_corners_for_mat = num_elems_for_mat * mesh.num_nodes_in_elem;

        State.MaterialToMeshMaps(mat_id).initialize(num_elems_for_mat);
        State.MaterialPoints(mat_id).initialize(num_points_for_mat, 3, SGTM3D_State::required_material_pt_state); // aways 3D, even for 2D-RZ calcs
        State.MaterialCorners(mat_id).initialize(num_corners_for_mat, mesh.num_dims, SGTM3D_State::required_material_corner_state);
        // zones are not used

    } // end for mat_id

    // check that the fills specify the required material point state fields
    bool filled_material_state_A =
        check_fill_mat_states(SGTM3D_State::required_optA_fill_material_pt_state,
                              SimulationParameters.region_setups.fill_gauss_states);
    bool filled_material_state_B =
        check_fill_mat_states(SGTM3D_State::required_optB_fill_material_pt_state,
                              SimulationParameters.region_setups.fill_gauss_states);
    
    // --- full stress tensor is not yet supported in region_fill ---
    //bool filled_material_state_C =
    //    check_fill_mat_states(SGTM3D_State::required_optC_fill_material_pt_state,
    //                                 SimulationParameters.region_setups.fill_gauss_states);

    if (filled_material_state_A == false &&
        filled_material_state_B == false){
        std::cout <<" Missing required material state in the fill instructions for the thrmex_FE solver \n";
        std::cout <<" The required state: \n";
        std::cout <<"  - density \n";
        std::cout <<"  - thermal conductivity \n";
        std::cout <<"  - specific heat \n";
        std::cout <<"  - specific or extensive internal energy \n" << std::endl;
        throw std::runtime_error("**** Provide fill instructions for all required material point variables ****");
    }
    
    // NOTE: Material points are populated in the material_state_setup function

    return;
}