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

#include "mesh_io.hpp"
#include "parse_yaml.hpp"
#include "solver.hpp"
#include "simulation_parameters.hpp"

#include "geometry_new.hpp"

// // Headers for solver classes
// #include "sgh_solver_3D.hpp"
// #include "sgh_solver_rz.hpp"
// #include "sgtm_solver_3D.hpp"


// Physical state data
#include "state.hpp"




class Driver
{
public:

    char* mesh_file;
    char* yaml_file;

    // ---------------------------------------------------------------------
    //    input type declarations
    // ---------------------------------------------------------------------

    MeshReader  mesh_reader;
    MeshBuilder mesh_builder;

    SimulationParameters_t SimulationParamaters; ///< the input simulation parameters

    // ---------------------------------------------------------------------
    //    Material and Boundary declarations
    // ---------------------------------------------------------------------

    Material_t Materials;                   ///< Material data for simulation
    BoundaryCondition_t BoundaryConditions; ///< Simulation boundary conditions

    int num_dims = 3;

    // ---------------------------------------------------------------------
    //    mesh data type declarations
    // ---------------------------------------------------------------------
    swage::Mesh mesh;

    // ---------------------------------------------------------------------
    //    state data type declaration
    // ---------------------------------------------------------------------
    State_t  State;

    int num_solvers = 0;

    // set of enabled solvers
    std::vector<Solver*> solvers;

    Driver(char* YAML)
    {
        yaml_file = YAML;
    };
    ~Driver() {};

    // Initialize driver data.  Solver type, number of solvers
    // Will be parsed from YAML input
    void initialize();

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn setup
    ///
    /// \brief Calls the setup function for each of the created solvers
    ///
    /////////////////////////////////////////////////////////////////////////////
    void setup();


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn execute
    ///
    /// \brief Calls the execute function for each of the created solvers
    ///
    /////////////////////////////////////////////////////////////////////////////
    void execute();


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn finalize
    ///
    /// \brief Calls the finalize function of each of the solvers assuming the 
    ///        finalize function exists and deletes the solver
    ///
    ///
    /////////////////////////////////////////////////////////////////////////////
    void finalize();


    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn setup_solver_vars
    ///
    /// \brief a function to set the solver variables based on user yaml inputs
    ///
    ///
    /////////////////////////////////////////////////////////////////////////////
    template <typename T>
    void setup_solver_vars(T& a_solver, 
                           const size_t solver_id);


}; // end driver class





