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

#include "sgh_solver_3D.h"
#include "mesh.h"
#include "region_fill.h"
#include "material.h"
#include "boundary_conditions.h"
#include "state.h"
#include "simulation_parameters.h"
#include "geometry_new.h"




/////////////////////////////////////////////////////////////////////////////
///
/// \fn setup the SGH method
///
/// \brief Allocate state, setup models, and fill mesh regions per the YAML input
///
/////////////////////////////////////////////////////////////////////////////
void SGH3D::setup(SimulationParameters_t& SimulationParamaters, 
                Material_t& Materials, 
                Mesh_t& mesh, 
                BoundaryCondition_t& Boundary,
                State_t& State)
{
    // add a flag on whether SGH was set up, if(SGH_setup_already==false)
    
    const size_t num_mats = Materials.num_mats; // the number of materials on the mesh

    // calculate pressure, sound speed, and stress for each material
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {

        // call the initialization function for state vars
        init_state_vars(Materials,
                        mesh,
                        State.MaterialPoints.eos_state_vars,
                        State.MaterialPoints.strength_state_vars,
                        State.MaterialToMeshMaps.elem,
                        State.MaterialPoints.num_material_points.host(mat_id),
                        mat_id);

        // call the init function for pressure, sound speed, and stress
        init_press_sspd_stress(Materials,
                               mesh,
                               State.MaterialPoints.den,
                               State.MaterialPoints.pres,
                               State.MaterialPoints.stress,
                               State.MaterialPoints.sspd,
                               State.MaterialPoints.sie,
                               State.MaterialPoints.eos_state_vars,
                               State.MaterialPoints.strength_state_vars,
                               State.MaterialPoints.shear_modulii,
                               State.MaterialPoints.num_material_points.host(mat_id),
                               mat_id);
    } // for loop over mat_id

    // set corner and node masses to zero
    init_corner_node_masses_zero(mesh, State.node.mass, State.corner.mass);

    // calculate corner and node masses on the mesh
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {

        calc_corner_mass(Materials,
                         mesh,
                         State.node.coords,
                         State.node.mass,
                         State.corner.mass,
                         State.MaterialPoints.mass,
                         State.MaterialToMeshMaps.elem,
                         State.MaterialToMeshMaps.num_material_elems.host(mat_id),
                         mat_id);
    } // end for mat_id

    calc_node_mass(mesh,
                   State.node.coords,
                   State.node.mass,
                   State.corner.mass);

// setting up fracture
    for (size_t i = 0; i < mesh.num_bdy_sets; i++) {
        // if fracture is allowed, then set up the fracture bank
        // note, allow_fracture is set in the parse_bdy_conds_inputs.cpp file and boundary_conditions.h file
        // checking if fracture is allowed... if = 0 then fracture is not enabled; if = 1, then fracture is enabled:
        std::cout << "Boundary.allow_fracture = " << Boundary.allow_fracture << std::endl;
        if (Boundary.allow_fracture) {
            std::cout << "Setting up global fracture (cohesive zones)" << std::endl;
            doing_fracture = true;
        
        // test to see if function is being entered
        std::cout << "Calling initialize()..." << std::endl;
        cohesive_zones_t cohesive_zones_bank;
        cohesive_zones_bank.initialize(mesh, State);
        // to see if function is being entered
        std::cout << "Done calling initialize()..." << std::endl;
        // Example from Gavin's code for running fracture tests:
        // run_fracture_tests(cohesive_zones_bank, mesh, State, sim_param);

        break; 
        }
    }
    // end setting up fracture


// ********************************************************** notes ************** ignore for now ******************* will clean up later *************************************************************
/*     // the following code counts the number of boundary nodes and checks for node overlaps (2 nodes with the same coordinates)
    // this is the beginning step to setting up cohesive zones for fracture
                   
    // counting the number of boundary nodes
    size_t num_bdy_nodes = mesh.num_bdy_nodes;
    std::cout << "Number of boundary nodes: " << num_bdy_nodes << std::endl;
    std::cout << "Total boundary nodes: " << mesh.num_bdy_nodes << std::endl;

    // total number of boundary nodes across all sets
    size_t total_bdy_nodes = 0;
    for (size_t i = 0; i < mesh.num_bdy_sets; ++i) {
        std::cout << "Boundary nodes in set " << i << ": " << mesh.num_bdy_nodes_in_set(i) << std::endl;
    }
    

    size_t overlap_count = 0;

    const double tol = 0.000001; //0.000001; //e-3; // adjust as needed; added just in case coordinate pairs are close but not exactly equal

    for (size_t i = 0; i < mesh.num_bdy_nodes; ++i) { // outer loop goes through each boundary node
        size_t node_i = mesh.bdy_nodes(i);
        for (size_t j = i + 1; j < mesh.num_bdy_nodes; ++j) { // Avoid duplicate pairs and i == j (starts from i +1), for each pair, compares coordinates k = 0, 1, 2
            size_t node_j = mesh.bdy_nodes(j);
            bool overlap = true;

            for (size_t k = 0; k < 3; ++k) {
                double diff = State.node.coords(node_i, k) - State.node.coords(node_j, k);
                if (std::abs(diff) > tol) {
                 overlap = false;
                break;
                }
            }
            if (overlap) {
            std::cout << "Overlap found between node " << node_i << " and node " << node_j << std::endl;
            ++overlap_count;
            }
        }
    }

    std::cout << "Total overlapping node pairs: " << overlap_count << std::endl; //prints overlapping node pairs

    // printing nodal coordinates to check for overlaps 
    for (size_t i = 0; i < mesh.num_bdy_nodes; ++i) {
        size_t node_id = mesh.bdy_nodes(i);
        std::cout << "Node " << node_id << " coordinates: ";
        for (size_t k = 0; k < 3; ++k) { // 3D
            std::cout << State.node.coords(node_id, k) << " ";
        }
        std::cout << std::endl;
    }   
 */
// ********************************************************** notes ************** ignore for now ******************* will clean up later *************************************************************
// ********************************************************** notes ************** ignore for now ******************* will clean up later *************************************************************
    // below is where fracture will be set up
    // if fracture is allowed, then set up the fracture bank
    // this will be a function inside of the fracture code and will have to be called in the SGH setup function
    // fills out UPs
    // UPs are Unique Pairs = unique node pairs with cohesive zone between them
    // this function finds all of the boundary nodes
    // loops again to count how many overlaps

    /*      // setting up fracture
     for (size_t i=0; mesh.num_bdy_nodes; i++) {
        if (Boundary.allow_fracture) {
            std::cout << "Setting up fracture" << std::endl;
            doing_fracture = true;

        //    fracture_bank.initialize(mesh, mesh.num_bdy_nodes, );
        })
            }
        }
     } */
    /*     // this will be a function inside of the fracture code and will have to be called in the SGH setup function
    // fills out UPs
    // UPs are Unique Pairs = unique node pairs with cohesive zone between them
    // this function finds all of the boundary nodes

    // loops again to count how many overlaps
    size_t num_vczs = 0;
    // if (overlap) -> num_vczs += 1;

    CArrayKokkos <size_t> vcz_pairs(num_vczs,2);

    bool overlap = true;
    size_t count = 0;
    for (int i = 0; i < mesh.num_bdy_nodes; i++) {
        for (int j = 0; j < mesh.num_bdy_nodes; j++) {
            if (i != j) {
                for (int k = 0; k < 3; k++) {
                    if (State.node.coords(mesh.bdy_nodes(i),k) != State.node.coords(mesh.bdy_nodes(j),k)) {
                    overlap = false;
                }
            }

            if (overlap) {
                // store mesh.bdy_nodes(i) and mesh.bdy_nodes(j) into vcz_pairs
                vcz_pairs(count,0) = mesh.bdy_nodes(i);
                vcz_pairs(count,1) = mesh.bdy_nodes(j);
                count += 1;
            }
        }
    }
    }  */
    // ********************************************************** notes ************** ignore for now ******************* will clean up later *************************************************************
} // end SGH setup    

 