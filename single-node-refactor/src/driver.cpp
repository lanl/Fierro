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

#include "driver.h"


// Headers for solver classes
#include "sgh_solver_3D.h"
#include "sgh_solver_rz.h"
//#include "sgtm_solver_3D.h"


// Initialize driver data.  Solver type, number of solvers
// Will be parsed from YAML input
void Driver::initialize()
{
	std::cout << "Initializing Driver" << std::endl;
    Yaml::Node root;
    try
    {
        Yaml::Parse(root, yaml_file);
    }
    catch (const Yaml::Exception e)
    {
        std::cout << "Exception " << e.Type() << ": " << e.what() << std::endl;
        exit(0);
    }

    parse_yaml(root, SimulationParamaters, Materials, BoundaryConditions);
    std::cout << "Finished  parsing YAML file" << std::endl;

    if (SimulationParamaters.mesh_input.source == mesh_input::file) {
        // Create and/or read mesh
        std::cout << "Mesh file path: " << SimulationParamaters.mesh_input.file_path << std::endl;
        mesh_reader.set_mesh_file(SimulationParamaters.mesh_input.file_path.data());
        mesh_reader.read_mesh(mesh, 
                              State,
                              SimulationParamaters.mesh_input, 
                              num_dims, 
                              SimulationParamaters.dynamic_options.rk_num_bins);
    }
    else if (SimulationParamaters.mesh_input.source == mesh_input::generate) {
        mesh_builder.build_mesh(mesh, 
                                State.GaussPoints, 
                                State.node, 
                                State.corner, 
                                SimulationParamaters);
    }
    else{
        throw std::runtime_error("**** NO MESH INPUT OPTIONS PROVIDED IN YAML ****");
        return;
    }

    // Build boundary conditions
    int num_bcs = BoundaryConditions.num_bcs;
    printf("Num BC's = %d\n", num_bcs);

    // --- calculate bdy sets ---//
    mesh.init_bdy_sets(num_bcs);
    tag_bdys(BoundaryConditions, mesh, State.node.coords);
    mesh.build_boundry_node_sets(mesh);

    // Setup Solvers
    for (int solver_id = 0; solver_id < SimulationParamaters.solver_inputs.size(); solver_id++) {

        if (SimulationParamaters.solver_inputs[solver_id].method == solver_input::SGH3D) {

            SGH3D* sgh_solver = new SGH3D(); 

            sgh_solver->initialize(SimulationParamaters, 
                                   Materials, 
                                   mesh, 
                                   BoundaryConditions,
                                   State);

            solvers.push_back(sgh_solver);
        } // end if SGH solver
        else if (SimulationParamaters.solver_inputs[solver_id].method == solver_input::SGHRZ) {

            SGHRZ* sgh_solver_rz = new SGHRZ(); 

            sgh_solver_rz->initialize(SimulationParamaters, 
                                   Materials, 
                                   mesh, 
                                   BoundaryConditions,
                                   State);

            solvers.push_back(sgh_solver_rz);
        } // end if SGHRZ solver
        //else if (SimulationParamaters.solver_inputs[solver_id].method == solver_input::SGTM3D) {

        //    SGTM3D* sgtm_solver_3d = new SGTM3D(); 
        //
        //    sgtm_solver_3d->initialize(SimulationParamaters, 
        //                               Materials, 
        //                               mesh, 
        //                               BoundaryConditions,
        //                               State);
        //
        //    solvers.push_back(sgtm_solver_3d);
        //} // end if SGTM solver

    } // end for loop over solvers

}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn setup
///
/// \brief Calls the setup function for each of the created solvers
///
/////////////////////////////////////////////////////////////////////////////
void Driver::setup()
{
    std::cout << "Inside driver setup" << std::endl;

    // allocate state, setup models, and apply fill instructions
    for (auto& solver : solvers) {
        solver->setup(SimulationParamaters, 
                      Materials, 
                      mesh, 
                      BoundaryConditions,
                      State);
    } // end for over solvers

} // end setup function of driver


/////////////////////////////////////////////////////////////////////////////
///
/// \fn execute
///
/// \brief Calls the execute function for each of the created solvers
///
/////////////////////////////////////////////////////////////////////////////
void Driver::execute()
{
    std::cout << "Inside driver execute" << std::endl;
    for (auto& solver : solvers) {
        solver->execute(SimulationParamaters, 
                        Materials, 
                        BoundaryConditions, 
                        mesh, 
                        State);
    } // loop over solvers
}


/////////////////////////////////////////////////////////////////////////////
///
/// \fn finalize
///
/// \brief Calls the finalize function of each of the solvers assuming the 
///        finalize function exists and deletes the solver
///
///
/////////////////////////////////////////////////////////////////////////////
void Driver::finalize()
{
    std::cout << "Inside driver finalize" << std::endl;
    for (auto& solver : solvers) {
        if (solver->finalize_flag) {
            solver->finalize(SimulationParamaters, 
                             Materials, 
                             BoundaryConditions);
        }
    }
    // destroy FEA modules
    for (auto& solver : solvers) {
        std::cout << "Deleting solver" << std::endl;
        delete solver;
    }
}
