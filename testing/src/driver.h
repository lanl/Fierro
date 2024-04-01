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

// #include "matar.h"


#include "io_utils.h"

#include "parse_yaml.h"

#include "solver.h"

// Headers for solver classes
#include "sgh_solver.h"


#include "material.h"
#include "region.h"
#include "mesh_inputs.h"
#include "solver_inputs.h"
#include "output_options.h"
#include "boundary_conditions.h"
#include "dynamic_options.h"


class Driver
{
public:

    char* mesh_file;
    char* yaml_file;

    MeshReader mesh_reader;


    mesh_input_t mesh_input;
    output_options_t output_options;
    dynamic_options_t dynamic_options;

    std::vector <solver_input_t> solver_inputs;

    std::vector <boundary_condition_t> boundary_conditions;

    std::vector <reg_fill_t> region_fills;
    std::vector <material_t> materials;
    std::vector <std::vector <double>> eos_global_vars;




    Driver(char* YAML){
        yaml_file = YAML;
    };//Simulation_Parameters& _simparam);
    ~Driver(){};

    // Initialize driver data.  Solver type, number of solvers 
    // Will be parsed from YAML input
    void initialize(){

        std::cout<<"Inside driver initialize"<<std::endl;

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


        parse_yaml(root, solver_inputs, mesh_input, dynamic_options, output_options, region_fills, materials, eos_global_vars, boundary_conditions);

        std::cout << "Done " << std::endl;

        // SGH *sgh_solver = new SGH(mesh_reader);
        // sgh_solver->initialize();
        // solvers.push_back(sgh_solver);
    }
    
    void setup() {
        std::cout<<"Inside driver setup"<<std::endl;
        for (auto & solver : solvers) {
            solver->setup();
        }
    }

    void run(){
        std::cout<<"Inside driver run"<<std::endl;
        for (auto & solver : solvers) {
            solver->execute();
        }
    }

    void finalize(){
        std::cout<<"Inside driver finalize"<<std::endl;
        for (auto & solver : solvers) {
            if (solver->finalize_flag){
                solver->solver_finalize();
            }
        }
        // destroy FEA modules
        for (auto & solver : solvers)
        {
            std::cout<<"Deleting solver"<<std::endl;
            delete solver;
        }
    }


    int num_solvers = 0;


    // set of enabled solvers
    std::vector<Solver*> solvers;


};
