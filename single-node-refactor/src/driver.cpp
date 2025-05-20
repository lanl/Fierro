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
#include "sgtm_solver_3D.h"
#include "level_set_solver.h"

#include "region_fill.h"
#include <mpi.h>




// Initialize driver data.  Solver type, number of solvers
// Will be parsed from YAML input
void Driver::initialize()
{
    
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&nranks);
    
    if(myrank == 0){
	    std::cout << "Initializing Driver" << std::endl;
    }
    Yaml::Node root;
    try
    {
        Yaml::Parse(root, yaml_file);
    }
    catch (const Yaml::Exception e)
    {
        if(myrank == 0){
            std::cout << "Exception " << e.Type() << ": " << e.what() << std::endl;
        }
        
        MPI_Finalize();
        MPI_Barrier(MPI_COMM_WORLD);
        exit(0);
    }

    parse_yaml(root, SimulationParameters, Materials, BoundaryConditions);
    
    if(myrank == 0){
        std::cout << "Finished  parsing YAML file" << std::endl;
    }

    if (SimulationParameters.mesh_input.source == mesh_input::file) {
        // Create and/or read mesh
        if(myrank == 0){
            std::cout << "Mesh file path: " << SimulationParameters.mesh_input.file_path << std::endl;
        }
        mesh_reader.set_mesh_file(SimulationParameters.mesh_input.file_path.data());
        mesh_reader.read_mesh(mesh, 
                              State,
                              SimulationParameters.mesh_input, 
                              num_dims);
    }
    else if (SimulationParameters.mesh_input.source == mesh_input::generate) {
        mesh_builder.build_mesh(mesh, 
                                State.GaussPoints, 
                                State.node, 
                                State.corner, 
                                SimulationParameters);
    }
    else{
        throw std::runtime_error("**** NO MESH INPUT OPTIONS PROVIDED IN YAML ****");
        MPI_Finalize();
        MPI_Barrier(MPI_COMM_WORLD);
        exit(0);
    }

    // Build boundary conditions
    const int num_bcs = BoundaryConditions.num_bcs;

    // --- calculate bdy sets ---//
    mesh.init_bdy_sets(num_bcs);
    tag_bdys(BoundaryConditions, mesh, State.node.coords);
    build_boundry_node_sets(mesh);


    // Setup the Solvers
    double time_final = SimulationParameters.dynamic_options.time_final;
    for (size_t solver_id = 0; solver_id < SimulationParameters.solver_inputs.size(); solver_id++) {

        if (SimulationParameters.solver_inputs[solver_id].method == solver_input::SGH3D) {
            
            if(myrank == 0){
                std::cout << "Initializing dynx_FE solver" << std::endl;
            }
            SGH3D* sgh_solver = new SGH3D(); 

            sgh_solver->initialize(SimulationParameters, 
                                   Materials, 
                                   mesh, 
                                   BoundaryConditions,
                                   State);

            // set the variables in the solver class
            setup_solver_vars(sgh_solver, 
                              solver_id);

            solvers.push_back(sgh_solver);


        } // end if SGH solver
        else if (SimulationParameters.solver_inputs[solver_id].method == solver_input::SGHRZ) {
            
            if(myrank == 0){
                std::cout << "Initializing dynx_FE_RZ solver" << std::endl;
            }
            SGHRZ* sgh_solver_rz = new SGHRZ(); 

            sgh_solver_rz->initialize(SimulationParameters, 
                                   Materials, 
                                   mesh, 
                                   BoundaryConditions,
                                   State);

            // set the variables in the solver class
            setup_solver_vars(sgh_solver_rz, 
                              solver_id);

            solvers.push_back(sgh_solver_rz);
        } // end if SGHRZ solver
        else if (SimulationParameters.solver_inputs[solver_id].method == solver_input::SGTM3D) {
            
            if(myrank == 0){
                std::cout << "Initializing thrmex_FE solver" << std::endl;
            }
            SGTM3D* sgtm_solver_3d = new SGTM3D(); 
        
            sgtm_solver_3d->initialize(SimulationParameters, 
                                       Materials, 
                                       mesh, 
                                       BoundaryConditions,
                                       State);

            // set the variables in the solver class
            setup_solver_vars(sgtm_solver_3d, 
                              solver_id);

            solvers.push_back(sgtm_solver_3d);

        } // end if SGTM solver
        else if (SimulationParamaters.solver_inputs[solver_id].method == solver_input::levelSet) {

            std::cout << "Initializing level set solver" << std::endl;
            LevelSet* level_set_solver = new LevelSet(); 
        
            level_set_solver->initialize(SimulationParamaters, 
                                         Materials, 
                                         mesh, 
                                         BoundaryConditions,
                                         State);

            // set the variables in the solver class
            setup_solver_vars(level_set_solver, 
                              solver_id);

            solvers.push_back(level_set_solver);

        } // end if level set solver
        else {
            throw std::runtime_error("**** NO SOLVER INPUT OPTIONS PROVIDED IN YAML, OR OPTION NOT UNDERSTOOD ****");
            MPI_Finalize();
            MPI_Barrier(MPI_COMM_WORLD);
            exit(0);
        }

    } // end for loop over solvers


    // ----
    // setup the simulation by applying all the fills to the mesh

    fillGaussState_t fillGaussState;
    fillElemState_t  fillElemState;

    simulation_setup(SimulationParameters, 
                     Materials, 
                     mesh, 
                     BoundaryConditions,
                     State,
                     fillGaussState,
                     fillElemState);


    // Allocate material state
    for (auto& solver : solvers) {
        solver->initialize_material_state(SimulationParameters, 
                      Materials, 
                      mesh, 
                      BoundaryConditions,
                      State);
    } // end for over solvers


    // populate the material point state
    material_state_setup(SimulationParameters, 
                         Materials, 
                         mesh, 
                         BoundaryConditions,
                         State,
                         fillGaussState,
                         fillElemState);
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
    if(myrank == 0){
        std::cout << "Inside driver setup" << std::endl;
    }

    // allocate state, setup models, and apply fill instructions
    for (auto& solver : solvers) {
        solver->setup(SimulationParameters, 
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
    if(myrank == 0){
        std::cout << "Inside driver execute" << std::endl;
    }
    for (auto& solver : solvers) {
        solver->execute(SimulationParameters, 
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
    if(myrank == 0){
        std::cout << "Inside driver finalize" << std::endl;
    }
    for (auto& solver : solvers) {
        if (solver->finalize_flag) {
            solver->finalize(SimulationParameters, 
                             Materials, 
                             BoundaryConditions);
        }
    }
    // destroy FEA modules
    for (auto& solver : solvers) {
        if(myrank == 0){
            std::cout << "Deleting solver" << std::endl;
        }
        delete solver;
    }
}

/////////////////////////////////////////////////////////////////////////////
///
/// \fn setup_solver_vars
///
/// \brief a function to set the solver variables based on user yaml inputs
///
///
/////////////////////////////////////////////////////////////////////////////
template <typename T>
void Driver::setup_solver_vars(T& a_solver, 
                               const size_t solver_id){


    // the final time of the simulation
    double time_final = this->SimulationParamaters.dynamic_options.time_final;

    // save the solver_id
    a_solver->solver_id = solver_id;

    // add other solver vars here


    // -------
    // setting the ending times are tricky, requiring logic
            
    // set the start and ending times
    double t_end = this->SimulationParamaters.solver_inputs[solver_id].time_end;  // default is t=0
    if(solver_id==0){
        a_solver->time_start = 0.0;

        if(t_end <= 1.0e-14){
            // time wasn't set so end the solver at the final time of the calculation
            a_solver->time_end = fmax(0.0, time_final);
        }
        else {
            // use the specified time in the input file
            a_solver->time_end = t_end;
        } // end if time was set
    }
    else {
        // this is not the first solver run
        double time_prior_solver_ends = this->solvers[solver_id-1]->time_end;

        a_solver->time_start = time_prior_solver_ends;

        if (t_end <= 1.0e-14){
            // time wasn't set so end the solver at the final time of the calculation
            a_solver->time_end = time_final;
        }
        else {
            // use the specified time in the input file
            a_solver->time_end = fmin(t_end, time_final); // ensure t_end is bounded by final time
        } // end if time was set
        
    } // end if solver=0

    std::cout << "Solver " << solver_id << " start time = " << a_solver->time_start << ", ending time = " << a_solver->time_end << "\n";

    return;
} // end setup solver vars function
