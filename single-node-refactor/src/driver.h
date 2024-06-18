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

#include "io_utils.h"
#include "parse_yaml.h"
#include "solver.h"
#include "simulation_parameters.h"

// Headers for solver classes
#include "sgh_solver.h"

// Physical state data
#include "state.h"


void fill_regions(simulation_parameters_t, mesh_t, node_t, elem_t, corner_t);


class Driver
{
public:

    char* mesh_file;
    char* yaml_file;

    MeshReader  mesh_reader;
    MeshBuilder mesh_builder;
    simulation_parameters_t sim_param;

    int num_dims = 3;

    // ---------------------------------------------------------------------
    //    mesh data type declarations
    // ---------------------------------------------------------------------
    mesh_t mesh;

    // ---------------------------------------------------------------------
    //    state data type declarations
    // ---------------------------------------------------------------------
    node_t   node;
    elem_t   elem;
    corner_t corner;

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
    void initialize()
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

        parse_yaml(root, sim_param);
        std::cout << "Finished  parsing YAML file" << std::endl;

        if (sim_param.mesh_input.source == mesh_input::file) {
            // Create and/or read mesh
            std::cout << "Mesh file path: " << sim_param.mesh_input.file_path << std::endl;
            mesh_reader.set_mesh_file(sim_param.mesh_input.file_path.data());
            mesh_reader.read_mesh(mesh, elem, node, corner, num_dims, sim_param.dynamic_options.rk_num_bins);
        }
        else if (sim_param.mesh_input.source == mesh_input::generate) {
            mesh_builder.build_mesh(mesh, elem, node, corner, sim_param);
        }
        else{
            throw std::runtime_error("**** NO MESH INPUT OPTIONS PROVIDED IN YAML ****");
            return;
        }

        // Build boundary conditions
        int num_bcs = sim_param.boundary_conditions.size();
        printf("Num BC's = %d\n", num_bcs);

        // --- calculate bdy sets ---//
        mesh.init_bdy_sets(num_bcs);
        tag_bdys(sim_param.boundary_conditions, mesh, node.coords);
        mesh.build_boundry_node_sets(sim_param.boundary_conditions, mesh);

        // Calculate element volume
        geometry::get_vol(elem.vol, node.coords, mesh);

        // Create memory for state variables
        // look and remove all materials.update_host();
        // CArray<MaterialModelValues_t>& MaterialModelVars
        elem.statev = DCArrayKokkos<double>(mesh.num_elems, sim_param.MaterialModelVars(0).eos_global_vars.size()); // WARNING: HACK

        // --- apply the fill instructions over the Elements---//
        Kokkos::fence();
        
        
        //fill_regions();
        fill_regions(sim_param, mesh, node, elem, corner);


        // --- Move the following sovler setup to yaml parsing routine
        // Create solvers
        for (int solver_id = 0; solver_id < sim_param.solver_inputs.size(); solver_id++) {
            if (sim_param.solver_inputs[solver_id].method == solver_input::SGH) {
                SGH* sgh_solver = new SGH(); // , mesh, node, elem, corner
                sgh_solver->initialize(sim_param);
                solvers.push_back(sgh_solver);
            }
        }

    } // end initialize

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn setup
    ///
    /// \brief Calls the setup function for each of the created solvers
    ///
    /////////////////////////////////////////////////////////////////////////////
    void setup()
    {
        std::cout << "Inside driver setup" << std::endl;
        for (auto& solver : solvers) {
            solver->setup(sim_param, mesh, node, elem, corner);
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn run
    ///
    /// \brief Calls the exectue function for each of the created solvers
    ///
    /////////////////////////////////////////////////////////////////////////////
    void run()
    {
        std::cout << "Inside driver run" << std::endl;
        for (auto& solver : solvers) {
            solver->execute(sim_param, mesh, node, elem, corner);
        }
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
    void finalize()
    {
        std::cout << "Inside driver finalize" << std::endl;
        for (auto& solver : solvers) {
            if (solver->finalize_flag) {
                solver->finalize(sim_param);
            }
        }
        // destroy FEA modules
        for (auto& solver : solvers) {
            std::cout << "Deleting solver" << std::endl;
            delete solver;
        }
    }

    
}; // end driver class

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn fill_regions
    ///
    /// \brief Fills mesh regions based on YAML input
    ///
    /////////////////////////////////////////////////////////////////////////////
    void fill_regions(simulation_parameters_t sim_param, mesh_t mesh, node_t node, elem_t elem, corner_t corner)
    {
        int num_fills = sim_param.region_fills.size();
        printf("Num Fills's = %d\n", num_fills);

        for (int f_id = 0; f_id < num_fills; f_id++) {
            // // voxel mesh setup
            // if (read_voxel_file.host(f_id) == 1)
            // {
            //     // read voxel mesh
            //     user_voxel_init(voxel_elem_values,
            //                     voxel_dx, voxel_dy, voxel_dz,
            //                     orig_x, orig_y, orig_z,
            //                     voxel_num_i, voxel_num_j, voxel_num_k);

            //     // copy values read from file to device
            //     voxel_elem_values.update_device();
            // } // endif

            int num_elems = mesh.num_elems;
            // parallel loop over elements in mesh
            FOR_ALL(elem_gid, 0, num_elems, {
                for (int rk_level = 0; rk_level < 2; rk_level++) {

                    // Set erosion flag to false
                    elem.eroded(elem_gid) = false;
                    
                    // calculate the coordinates and radius of the element
                    double elem_coords[3]; // note:initialization with a list won't work
                    elem_coords[0] = 0.0;
                    elem_coords[1] = 0.0;
                    elem_coords[2] = 0.0;

                    // get the coordinates of the element center
                    for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                        elem_coords[0] += node.coords(rk_level, mesh.nodes_in_elem(elem_gid, node_lid), 0);
                        elem_coords[1] += node.coords(rk_level, mesh.nodes_in_elem(elem_gid, node_lid), 1);
                        if (mesh.num_dims == 3) {
                            elem_coords[2] += node.coords(rk_level, mesh.nodes_in_elem(elem_gid, node_lid), 2);
                        }
                        else{
                            elem_coords[2] = 0.0;
                        }
                    } // end loop over nodes in element (NOTE: Translating using origin so below check works)
                    elem_coords[0] = (elem_coords[0] / mesh.num_nodes_in_elem) - sim_param.region_fills(f_id).origin[0];
                    elem_coords[1] = (elem_coords[1] / mesh.num_nodes_in_elem) - sim_param.region_fills(f_id).origin[1];
                    elem_coords[2] = (elem_coords[2] / mesh.num_nodes_in_elem) - sim_param.region_fills(f_id).origin[2];

                    // spherical radius 
                    double radius = sqrt(elem_coords[0] * elem_coords[0] +
                                          elem_coords[1] * elem_coords[1] +
                                          elem_coords[2] * elem_coords[2]);

                    // cylindrical radius
                    double radius_cyl = sqrt(elem_coords[0] * elem_coords[0] +
                                              elem_coords[1] * elem_coords[1]);

                    // default is not to fill the element
                    size_t fill_this = 0;

                    // check to see if this element should be filled
                    switch (sim_param.region_fills(f_id).volume) {
                        case region::global:
                            {
                                fill_this = 1;
                                break;
                            }
                        case region::box:
                            {

                                double x_lower_bound = sim_param.region_fills(f_id).x1;
                                double x_upper_bound = sim_param.region_fills(f_id).x2;

                                double y_lower_bound = sim_param.region_fills(f_id).y1;
                                double y_upper_bound = sim_param.region_fills(f_id).y2;

                                double z_lower_bound = sim_param.region_fills(f_id).z1;
                                double z_upper_bound = sim_param.region_fills(f_id).z2;


                                if (elem_coords[0] >= x_lower_bound && elem_coords[0] <= x_upper_bound &&
                                    elem_coords[1] >= y_lower_bound && elem_coords[1] <= y_upper_bound &&
                                    elem_coords[2] >= z_lower_bound && elem_coords[2] <= z_upper_bound) {
                                    fill_this = 1;
                                }
                                break;
                            }
                        case region::cylinder:
                            {
                                if (radius_cyl >= sim_param.region_fills(f_id).radius1
                                    && radius_cyl <= sim_param.region_fills(f_id).radius2) {
                                    fill_this = 1;
                                }
                                break;
                            }
                        case region::sphere:
                            {
                                if (radius >= sim_param.region_fills(f_id).radius1
                                    && radius <= sim_param.region_fills(f_id).radius2) {
                                    fill_this = 1;
                                }
                                break;
                            }
                        case region::readVoxelFile:
                            {
                                fill_this = 0; // default is no, don't fill it

                                // // find the closest element in the voxel mesh to this element
                                // double i0_real = (elem_coords[0] - orig_x) / (voxel_dx);
                                // double j0_real = (elem_coords[1] - orig_y) / (voxel_dy);
                                // double k0_real = (elem_coords[2] - orig_z) / (voxel_dz);

                                // int i0 = (int)i0_real;
                                // int j0 = (int)j0_real;
                                // int k0 = (int)k0_real;

                                // // look for the closest element in the voxel mesh
                                // int elem_id0 = get_id(i0, j0, k0, voxel_num_i, voxel_num_j);

                                // // if voxel mesh overlaps this mesh, then fill it if =1
                                // if (elem_id0 < voxel_elem_values.size() && elem_id0 >= 0) {
                                //     // voxel mesh elem values = 0 or 1
                                //     fill_this = voxel_elem_values(elem_id0); // values from file
                                // } // end if
                            } // end case
                    } // end of switch

                    // paint the material state on the element
                    if (fill_this == 1) {
                        // density
                        elem.den(elem_gid) = sim_param.region_fills(f_id).den;

                        // mass
                        elem.mass(elem_gid) = elem.den(elem_gid) * elem.vol(elem_gid);

                        // specific internal energy
                        elem.sie(rk_level, elem_gid) = sim_param.region_fills(f_id).sie;

                        elem.mat_id(elem_gid) = sim_param.region_fills(f_id).material_id;

                        size_t mat_id = elem.mat_id(elem_gid); // short name

                        // get state_vars from the input file or read them in
                        if (false) { // sim_param.materials(mat_id).strength_setup == model_init::user_init) {
                            // use the values read from a file to get elem state vars
                            // for (size_t var = 0; var < sim_param.materials(mat_id).num_eos_state_vars; var++) {
                            //     elem.statev(elem_gid, var) = file_state_vars(mat_id, elem_gid, var);
                            // } // end for
                        }
                        else{
                            // use the values in the input file
                            // set state vars for the region where mat_id resides
                            printf("hello! \n");
                            int num_eos_global_vars = sim_param.MaterialModelVars(mat_id).eos_global_vars.size();
                            printf("ok that worked \n");

                            for (size_t var = 0; var < num_eos_global_vars; var++) {
                                elem.statev(elem_gid, var) = sim_param.MaterialModelVars(mat_id).eos_global_vars(var); // state_vars(mat_id, var); // WARNING: HACK
                            } // end for

                            printf("wow, I made it hear \n");
                        } // end logical on type

                        // --- stress tensor ---
                        // always 3D even for 2D-RZ
                        for (size_t i = 0; i < 3; i++) {
                            for (size_t j = 0; j < 3; j++) {
                                elem.stress(rk_level, elem_gid, i, j) = 0.0;
                            }
                        }  // end for

                        // --- Pressure and stress ---

                        // --- Pressure ---
                        sim_param.materials(mat_id).calc_pressure(elem.pres,
                                                                  elem.stress,
                                                                  elem_gid,
                                                                  elem.mat_id(elem_gid),
                                                                  elem.statev,
                                                                  elem.sspd,
                                                                  elem.den(elem_gid),
                                                                  elem.sie(rk_level, elem_gid));   
                        // --- Sound Speed ---                               
                        sim_param.materials(mat_id).calc_sound_speed(elem.pres,
                                                                     elem.stress,
                                                                     elem_gid,
                                                                     elem.mat_id(elem_gid),
                                                                     elem.statev,
                                                                     elem.sspd,
                                                                     elem.den(elem_gid),
                                                                     elem.sie(rk_level, elem_gid));

                        // loop over the nodes of this element and apply velocity
                        for (size_t node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                            // get the mesh node index
                            size_t node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

                            // --- Velocity ---
                            switch (sim_param.region_fills(f_id).velocity) {
                                case init_conds::cartesian:
                                    {
                                        node.vel(rk_level, node_gid, 0) = sim_param.region_fills(f_id).u;
                                        node.vel(rk_level, node_gid, 1) = sim_param.region_fills(f_id).v;
                                        if (mesh.num_dims == 3) {
                                            node.vel(rk_level, node_gid, 2) = sim_param.region_fills(f_id).w;
                                        }

                                        break;
                                    }
                                // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
                                case init_conds::radial:
                                    {
                                        // Setting up radial
                                        double dir[2];
                                        dir[0] = 0.0;
                                        dir[1] = 0.0;
                                        double radius_val = 0.0;

                                        for (int dim = 0; dim < 2; dim++) {
                                            dir[dim]    = node.coords(rk_level, node_gid, dim);
                                            radius_val += node.coords(rk_level, node_gid, dim) * node.coords(rk_level, node_gid, dim);
                                        } // end for
                                        radius_val = sqrt(radius_val);

                                        for (int dim = 0; dim < 2; dim++) {
                                            if (radius_val > 1.0e-14) {
                                                dir[dim] /= (radius_val);
                                            }
                                            else{
                                                dir[dim] = 0.0;
                                            }
                                        } // end for

                                        node.vel(rk_level, node_gid, 0) = sim_param.region_fills(f_id).speed * dir[0];
                                        node.vel(rk_level, node_gid, 1) = sim_param.region_fills(f_id).speed * dir[1];
                                        if (mesh.num_dims == 3) {
                                            node.vel(rk_level, node_gid, 2) = 0.0;
                                        }

                                        break;
                                    }
                                case init_conds::spherical:
                                    {
                                        // Setting up spherical
                                        double dir[3];
                                        dir[0] = 0.0;
                                        dir[1] = 0.0;
                                        dir[2] = 0.0;
                                        double radius_val = 0.0;

                                        for (int dim = 0; dim < 3; dim++) {
                                            dir[dim]    = node.coords(rk_level, node_gid, dim);
                                            radius_val += node.coords(rk_level, node_gid, dim) * node.coords(rk_level, node_gid, dim);
                                        } // end for
                                        radius_val = sqrt(radius_val);

                                        for (int dim = 0; dim < 3; dim++) {
                                            if (radius_val > 1.0e-14) {
                                                dir[dim] /= (radius_val);
                                            }
                                            else{
                                                dir[dim] = 0.0;
                                            }
                                        } // end for

                                        node.vel(rk_level, node_gid, 0) = sim_param.region_fills(f_id).speed * dir[0];
                                        node.vel(rk_level, node_gid, 1) = sim_param.region_fills(f_id).speed * dir[1];
                                        if (mesh.num_dims == 3) {
                                            node.vel(rk_level, node_gid, 2) = sim_param.region_fills(f_id).speed * dir[2];
                                        }

                                        break;
                                    }
                                case init_conds::radial_linear:
                                    {
                                        printf("**** Radial_linear initial conditions not yet supported ****\n");
                                        break;
                                    }
                                case init_conds::spherical_linear:
                                    {
                                        printf("**** spherical_linear initial conditions not yet supported ****\n");
                                        break;
                                    }
                                case init_conds::tg_vortex:
                                    {
                                        node.vel(rk_level, node_gid, 0) = sin(PI * node.coords(rk_level, node_gid, 0)) * cos(PI * node.coords(rk_level, node_gid, 1));
                                        node.vel(rk_level, node_gid, 1) =  -1.0 * cos(PI * node.coords(rk_level, node_gid, 0)) * sin(PI * node.coords(rk_level, node_gid, 1));
                                        if (mesh.num_dims == 3) {
                                            node.vel(rk_level, node_gid, 2) = 0.0;
                                        }

                                        break;
                                    }
                            } // end of switch
                        } // end loop over nodes of element

                        if (sim_param.region_fills(f_id).velocity == init_conds::tg_vortex) {
                            elem.pres(elem_gid) = 0.25 * (cos(2.0 * PI * elem_coords[0]) + cos(2.0 * PI * elem_coords[1]) ) + 1.0;

                            // p = rho*ie*(gamma - 1)
                            // size_t mat_id = f_id;
                            double gamma  = elem.statev(elem_gid, 4); // gamma value WARNING: BUG HERE
                            elem.sie(rk_level, elem_gid) =
                                elem.pres(elem_gid) / (sim_param.region_fills(f_id).den * (gamma - 1.0));
                        } // end if
                    } // end if fill
                } // end RK loop
            }); // end FOR_ALL element loop
            Kokkos::fence();
        } // end for loop over fills

        // // calculate the corner massess if 2D
        // if (mesh.num_dims == 2) {
        //     FOR_ALL(elem_gid, 0, mesh.num_elems, {
        //         // facial area of the corners
        //         double corner_areas_array[4];

        //         ViewCArrayKokkos<double> corner_areas(&corner_areas_array[0], 4);
        //         ViewCArrayKokkos<size_t> elem_node_gids(&mesh.nodes_in_elem(elem_gid, 0), 4);

        //         geometry::get_area_weights2D(corner_areas, elem_gid, node_coords, elem_node_gids);

        //         // loop over the corners of the element and calculate the mass
        //         for (size_t corner_lid = 0; corner_lid < 4; corner_lid++) {
        //             size_t corner_gid = mesh.corners_in_elem(elem_gid, corner_lid);
        //             corner_mass(corner_gid) = corner_areas(corner_lid) * elem_den(elem_gid); // node radius is added later
        //         } // end for over corners
        //     });
        // } // end of

        // calculate the nodal mass
        FOR_ALL(node_gid, 0, mesh.num_nodes, {
            node.mass(node_gid) = 0.0;

            if (mesh.num_dims == 3) {
                for (size_t elem_lid = 0; elem_lid < mesh.num_corners_in_node(node_gid); elem_lid++) {
                    size_t elem_gid      = mesh.elems_in_node(node_gid, elem_lid);
                    node.mass(node_gid) += 1.0 / 8.0 * elem.mass(elem_gid);
                } // end for elem_lid
            } // end if dims=3
            else{
                // 2D-RZ
                for (size_t corner_lid = 0; corner_lid < mesh.num_corners_in_node(node_gid); corner_lid++) {
                    size_t corner_gid    = mesh.corners_in_node(node_gid, corner_lid);
                    node.mass(node_gid) += corner.mass(corner_gid);  // sans the radius so it is areal node mass

                    corner.mass(corner_gid) *= node.coords(1, node_gid, 1); // true corner mass now
                } // end for elem_lid
            } // end else
        }); // end FOR_ALL
        Kokkos::fence();

    } // end fill regions
