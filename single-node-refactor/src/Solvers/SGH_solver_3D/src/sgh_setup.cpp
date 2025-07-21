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

    // Setting up contact
    // todo: should this be handled inside of src/boundary_conditions/stress/global_contact ?
    for (size_t i = 0; i < mesh.num_bdy_sets; i++) {
        if (Boundary.allow_contact) {
            std::cout << "Setting up global contact" << std::endl;
            doing_contact = true;

            contact_bank.initialize(mesh, mesh.bdy_patches, State);
            // run_contact_tests(contact_bank, mesh, node, corner, sim_param);
            break;
        }
    }

    
} // end SGH setup
