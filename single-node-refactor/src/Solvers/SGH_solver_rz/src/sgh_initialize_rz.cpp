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

#include "sgh_solver_rz.h"
#include "state.h"
#include "mesh.h"
#include "simulation_parameters.h"

void SGHRZ::initialize(SimulationParameters_t& SimulationParameters, 
                	   Material_t& Materials, 
                	   Mesh_t& mesh, 
                	   BoundaryCondition_t& Boundary,
                	   State_t& State)
{
	size_t num_nodes = mesh.num_nodes;
    size_t num_gauss_pts = mesh.num_elems;
    size_t num_corners = mesh.num_corners;
    size_t num_dim = mesh.num_dims;

    if (num_dim != 2){
        std::cout << "Wrong dimensions of " << num_dim << "\n";
        throw std::runtime_error("**** Solver is for 2D RZ coordinates, wrong dimensions specified in input ****");
    }

    // save the solver_id, which is a pravate class variable
    //this->solver_id = solver_id_inp;

    State.node.initialize(mesh.all_node_map, num_dim, SGHRZ_State::required_node_state, mesh.node_map);
    State.GaussPoints.initialize(num_gauss_pts, 3, SGHRZ_State::required_gauss_pt_state);  // note: dims is always 3 
    State.corner.initialize(num_corners, num_dim, SGHRZ_State::required_corner_state);

    //comms objects
    node_velocity_comms = CommPlan<real_t>(State.node.vel, State.node.local_vel, mesh.node_coords_comms); //copies MPI setup from coordinate comms since the node maps are the same
    node_mass_comms = CommPlan<real_t>(State.node.mass, State.node.local_mass, mesh.node_coords_comms); //copies MPI setup from coordinate comms since the node maps are the same
    
    // NOTE: Material points and material corners are initialize in sgh_setup after calculating the material->mesh maps
}

void SGHRZ::initialize_material_state(SimulationParameters_t& SimulationParameters, 
                	                  Material_t& Materials, 
                	                  Mesh_t& mesh, 
                	                  BoundaryCondition_t& Boundary,
                	                  State_t& State) const
{

    // -----
    //  Allocation of state includes a buffer, set in region_fill.cpp, which is for ALE
    // -----
    
    State.MaterialToMeshMaps.initialize();
    State.MaterialPoints.initialize(3, SGHRZ_State::required_material_pt_state); // aways 3D, even for 2D-RZ calcs
    State.MaterialCorners.initialize(mesh.num_dims, SGHRZ_State::required_material_corner_state); 
    // zones are not used
    
    // NOTE: Material points are populated in the material_state_setup funcion
    
    return;

}