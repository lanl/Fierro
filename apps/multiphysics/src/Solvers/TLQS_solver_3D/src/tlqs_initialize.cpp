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
#include "state.hpp"
#include "simulation_parameters.hpp"


void TLQS3D::initialize(SimulationParameters_t& SimulationParamaters, 
                	   Material_t& Materials, 
                	   swage::Mesh& mesh, 
                	   BoundaryCondition_t& Boundary,
                	   State_t& State) const
{
	const size_t num_nodes = mesh.num_nodes;
    const size_t num_gauss_pts = mesh.num_elems;
    const size_t num_corners = mesh.num_corners;
    const size_t num_dims = mesh.num_dims;

    if (num_dims != 3){
        std::cout << "Wrong dimensions of " << num_dims << "\n";
        throw std::runtime_error("**** Total Lagrangian Quasi-Static solver is for 3D coordinates, wrong dimensions specified in input ****");
    }

    // mesh state
    State.node.initialize(num_nodes, num_dims, TLQS3D_State::required_node_state);
    State.GaussPoints.initialize(num_gauss_pts, num_dims, TLQS3D_State::required_gauss_pt_state);
    State.corner.initialize(num_corners, num_dims, TLQS3D_State::required_corner_state);

    // check that the fills specify the required nodal fields
    bool filled_nodal_state =
        check_fill_node_states(TLQS3D_State::required_fill_node_state,
                               SimulationParamaters.region_setups.fill_node_states);
    
    if (filled_nodal_state == false){
        std::cout <<" Missing required nodal state in the fill instructions for the TLQS solver \n";
        std::cout <<" The nodal velocity must be specified. \n" << std::endl;
        throw std::runtime_error("**** Provide fill instructions for all required nodal variables ****");
    }

    std::cout << "TLQS solver initialized \n";

} // end solver initialization

void TLQS3D::initialize_material_state(SimulationParameters_t& SimulationParamaters, 
    Material_t& Materials, 
    swage::Mesh& mesh, 
    BoundaryCondition_t& Boundary,
    State_t& State) const
{

    // -----
    //  Allocation of state includes the buffer set in region_fill.cpp, it's needed for ALE
    // -----

    State.MaterialToMeshMaps.initialize();
    State.MaterialPoints.initialize(mesh.num_dims, TLQS3D_State::required_material_pt_state); 
    // corners are not used
    // zones are not used

} // end solver initialization


