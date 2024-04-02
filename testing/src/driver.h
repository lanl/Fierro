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

#include "simulation_parameters.h"

#include "state.h"



class Driver
{
public:

    char* mesh_file;
    char* yaml_file;

    MeshReader reader;
    simulation_parameters_t sim_param;

    int num_dims = 3;

    // ---------------------------------------------------------------------
    //    mesh data type declarations
    // ---------------------------------------------------------------------
    mesh_t                   mesh;
    CArrayKokkos<reg_fill_t> mat_fill;
    CArrayKokkos<boundary_t> boundary;

    // ---------------------------------------------------------------------
    //    state data type declarations
    // ---------------------------------------------------------------------
    node_t                   node;
    elem_t                   elem;
    corner_t                 corner;
    CArrayKokkos<material_t> material;
    CArrayKokkos<double>     state_vars; // array to hold init model variables


    // ---------------------------------------------------------------------
    //    GPU compatable types
    // ---------------------------------------------------------------------


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

        parse_yaml(root, sim_param);
        std::cout << "Finished  parsing YAML file" << std::endl;


        // Create and/or read mesh
        std::cout<<"Mesh file path: "<<sim_param.mesh_input.file_path<<std::endl;
        reader.set_mesh_file(sim_param.mesh_input.file_path.data());

        reader.read_mesh(mesh, elem, node, corner, num_dims, sim_param.dynamic_options.rk_num_bins);

        std::cout << "Num elements = " << mesh.num_elems << std::endl;
        std::cout << "Num nodes = " << mesh.num_nodes << std::endl;

        // std::cout << "Building corners: " << std::endl;
        mesh.build_corner_connectivity();

        // std::cout << "Building elem elem connectivity: " << std::endl;
        mesh.build_elem_elem_connectivity();

        // std::cout << "Building patches: " << std::endl;
        mesh.build_patch_connectivity();

        // std::cout << "Building node-node connectivity: " << std::endl;
        mesh.build_node_node_connectivity();

        // Setup boundary conditions on mesh

        int num_bcs = sim_param.boundary_conditions.size();

        mesh.init_bdy_sets(num_bcs);
        printf("Num BC's = %lu\n", num_bcs);

        // tag boundary patches in the set
        tag_bdys(boundary, mesh, node.coords);

        build_boundry_node_sets(boundary, mesh);




        // // Create solvers
        // for(int solver_id = 0; solver_id < sim_param.solver_inputs.size(); solver_id++){

        //     if(sim_param.solver_inputs[solver_id].method == solver_input::SGH){
        //         // SGH *sgh_solver = new SGH(sim_param);
        //         // sgh_solver->initialize();
        //         // solvers.push_back(sgh_solver);
        //     }
        // }
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
        // for (auto & solver : solvers) {
        //     if (solver->finalize_flag){
        //         solver->solver_finalize();
        //     }
        // }
        // // destroy FEA modules
        // for (auto & solver : solvers)
        // {
        //     std::cout<<"Deleting solver"<<std::endl;
        //     delete solver;
        // }
    }


    int num_solvers = 0;


    // set of enabled solvers
    std::vector<Solver*> solvers;


};
