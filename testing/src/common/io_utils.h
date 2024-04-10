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
#ifndef FIERRO_IO_H
#define FIERRO_IO_H

#include "matar.h"
#include "mesh.h"
#include "state.h"
#include "simulation_parameters.h"

#include <map>
#include <memory>

// Reads in mesh connectivity data
class MeshReader
{
public:

    char* mesh_file_ = NULL;

    MeshReader() {} // Simulation_Parameters& _simparam);

    ~MeshReader() = default;

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn set_mesh_file
    ///
    /// \brief Sets the mesh file path for reading in a mesh
    ///
    /// \param Path to mesh file
    ///
    /////////////////////////////////////////////////////////////////////////////
    void set_mesh_file(char* MESH)
    {
        mesh_file_ = MESH;
    }

    // Reads and initializes the mesh and geometric state entities
    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn read_mesh
    ///
    /// \brief Read mesh from file
    ///
    /// \param Simulation mesh
    /// \param Element state struct
    /// \param Node state struct
    /// \param Corner state struct
    /// \param Number of dimensions
    /// \param Number of RK bins
    ///
    ///
    /////////////////////////////////////////////////////////////////////////////
    void read_mesh(mesh_t& mesh, elem_t& elem, node_t& node, corner_t& corner, int num_dims, int rk_num_bins)
    {
        if (mesh_file_ == NULL) {
            printf("No mesh given\n");
            exit(0);
        }

        // Check mesh file extension
        // and read based on extension
        read_ensight_mesh(mesh, elem, node, corner, num_dims, rk_num_bins);
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn read_ensight_mesh
    ///
    /// \brief Read .geo mesh file
    ///
    /// \param Simulation mesh
    /// \param Element state struct
    /// \param Node state struct
    /// \param Corner state struct
    /// \param Number of dimensions
    /// \param Number of RK bins
    ///
    /////////////////////////////////////////////////////////////////////////////
    void read_ensight_mesh(mesh_t& mesh, elem_t& elem, node_t& node, corner_t& corner, int num_dims, int rk_num_bins)
    {
        const size_t rk_level = 0;

        FILE* in;
        char  ch;

        size_t num_nodes_in_elem = 1;
        for (int dim = 0; dim < num_dims; dim++) {
            num_nodes_in_elem *= 2;
        }

        // read the mesh    WARNING: assumes a .geo file
        in = fopen(mesh_file_, "r");

        // skip 8 lines
        for (int j = 1; j <= 8; j++) {
            int i = 0;
            while ((ch = (char)fgetc(in)) != '\n') {
                i++;
                // printf("%c",ch);
            }
            // printf("\n");
        }

        // --- Read in the nodes in the mesh ---

        size_t num_nodes = 0;

        fscanf(in, "%lu", &num_nodes);
        printf("Num nodes read in %lu\n", num_nodes);

        // intialize node variables
        mesh.initialize_nodes(num_nodes);
        node.initialize(rk_num_bins, num_nodes, num_dims);

        // read the initial mesh coordinates
        // x-coords
        for (int node_id = 0; node_id < mesh.num_nodes; node_id++) {
            fscanf(in, "%le", &node.coords(rk_level, node_id, 0));
        }

        // y-coords
        for (int node_id = 0; node_id < mesh.num_nodes; node_id++) {
            fscanf(in, "%le", &node.coords(rk_level, node_id, 1));
        }

        // z-coords
        for (int node_id = 0; node_id < mesh.num_nodes; node_id++) {
            if (num_dims == 3) {
                fscanf(in, "%le", &node.coords(rk_level, node_id, 2));
            }
            else{
                double dummy;
                fscanf(in, "%le", &dummy);

                // printf("dummy = %le\n", dummy);
            }
        } // end for

        ch = (char)fgetc(in);
        // printf("%c",ch);

        // skip 1 line
        for (int j = 1; j <= 1; j++) {
            int i = 0;
            while ((ch = (char)fgetc(in)) != '\n') {
                i++;
                // printf("%c",ch);
            }
            // printf("\n");
        }

        // --- read in the elements in the mesh ---
        size_t num_elem = 0;

        fscanf(in, "%lu", &num_elem);
        printf("Num elements read in %lu\n", num_elem);

        // intialize elem variables
        mesh.initialize_elems(num_elem, num_dims);
        elem.initialize(rk_num_bins, num_nodes, 3); // always 3D here, even for 2D

        // for each cell read the list of associated nodes
        for (int elem_gid = 0; elem_gid < num_elem; elem_gid++) {
            for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++) {
                fscanf(in, "%lu", &mesh.nodes_in_elem.host(elem_gid, node_lid));  // %d vs zu

                // shift to start node index space at 0
                mesh.nodes_in_elem.host(elem_gid, node_lid) -= 1;
            }
        }
        // update device side
        mesh.nodes_in_elem.update_device();

        // intialize corner variables
        int num_corners = num_elem * mesh.num_nodes_in_elem;
        mesh.initialize_corners(num_corners);
        corner.initialize(num_corners, num_dims);

        // save the node coords to the current RK value
        for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
            for (int rk = 1; rk < rk_num_bins; rk++) {
                for (int dim = 0; dim < num_dims; dim++) {
                    node.coords(rk, node_gid, dim) = node.coords(0, node_gid, dim);
                } // end for dim
            } // end for rk
        } // end parallel for

        // Close mesh input file
        fclose(in);

        // Build connectivity
        mesh.build_connectivity();

        return;
    }
};

// Builds a mesh
class MeshBuilder
{
public:

    // MeshBuilder() {}

    ~MeshBuilder()
    {
        // Empty destructor
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_mesh
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_mesh(mesh_t& mesh, elem_t& elem, node_t& node, corner_t& corner, simulation_parameters_t& sim_param)
    {
        if (sim_param.mesh_input.num_dims == 2) {
            if (sim_param.mesh_input.type == mesh_input::Cylinder) {
                build_2d_polar(mesh, elem, node, corner, sim_param);
            }
            else if (sim_param.mesh_input.type == mesh_input::Box) {
                build_2d_box(mesh, elem, node, corner, sim_param);
            }
            else{
                std::cout << "**** 2D MESH TYPE NOT SUPPORTED **** " << std::endl;
                std::cout << "Valid options are: " << std::endl;
                auto map = mesh_input_type_map;
                for (const auto& pair : map) {
                    std::cout << "\t" << pair.first << std::endl;
                }
                throw std::runtime_error("**** 2D MESH TYPE NOT SUPPORTED ****");
                return;
            }
        }
        else if (sim_param.mesh_input.num_dims == 3) {
            build_3d_box(mesh, elem, node, corner, sim_param);
        }
        else{
            throw std::runtime_error("**** ONLY 2D RZ OR 3D MESHES ARE SUPPORTED ****");
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_2d_box
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_2d_box(mesh_t& mesh, elem_t& elem, node_t& node, corner_t& corner, simulation_parameters_t& sim_param)
    {
        printf(" Creating a 2D box mesh \n");

        const int num_dim = 2;

        const double lx = sim_param.mesh_input.length[0];
        const double ly = sim_param.mesh_input.length[1];

        const int num_elems_i = sim_param.mesh_input.num_elems[0];
        const int num_elems_j = sim_param.mesh_input.num_elems[1];

        const int num_points_i = num_elems_i + 1; // num points in x
        const int num_points_j = num_elems_j + 1; // num points in y

        const int num_nodes = num_points_i * num_points_j;

        const double dx = lx / ((double)num_elems_i);  // len/(num_elems_i)
        const double dy = ly / ((double)num_elems_j);  // len/(num_elems_j)

        const int num_elems = num_elems_i * num_elems_j;

        std::vector<double> origin = sim_param.mesh_input.origin;

        // --- 2D parameters ---
        const int num_faces_in_elem  = 4;  // number of faces in elem
        const int num_points_in_elem = 4;  // number of points in elem
        const int num_points_in_face = 2;  // number of points in a face
        const int num_edges_in_elem  = 4;  // number of edges in a elem

        // --- mesh node ordering ---
        // Convert ijk index system to the finite element numbering convention
        // for vertices in elem
        auto convert_point_number_in_quad = CArray<int>(4);
        convert_point_number_in_quad(0) = 0;
        convert_point_number_in_quad(1) = 1;
        convert_point_number_in_quad(2) = 3;
        convert_point_number_in_quad(3) = 2;

        int rk_num_bins = sim_param.dynamic_options.rk_num_bins;

        // intialize node variables
        mesh.initialize_nodes(num_nodes);
        node.initialize(rk_num_bins, num_nodes, num_dim);

        // --- Build nodes ---

        // populate the point data structures
        for (int j = 0; j < num_points_j; j++) {
            for (int i = 0; i < num_points_i; i++) {
                // global id for the point
                int node_gid = get_id(i, j, 0, num_points_i, num_points_j);

                // store the point coordinates
                node.coords(0, node_gid, 0) = origin[0] + (double)i * dx;
                node.coords(0, node_gid, 1) = origin[1] + (double)j * dy;
            } // end for i
        } // end for j

        for (int rk_level = 1; rk_level < rk_num_bins; rk_level++) {
            for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
                node.coords(rk_level, node_gid, 0) = node.coords(0, node_gid, 0);
                node.coords(rk_level, node_gid, 1) = node.coords(0, node_gid, 1);
            }
        }

        // intialize elem variables
        mesh.initialize_elems(num_elems, num_dim);
        elem.initialize(rk_num_bins, num_nodes, 3); // always 3D here, even for 2D

        // populate the elem center data structures
        for (int j = 0; j < num_elems_j; j++) {
            for (int i = 0; i < num_elems_i; i++) {
                // global id for the elem
                int elem_gid = get_id(i, j, 0, num_elems_i, num_elems_j);

                // store the point IDs for this elem where the range is
                // (i:i+1, j:j+1 for a linear quad
                int this_point = 0;

                for (int jcount = j; jcount <= j + 1; jcount++) {
                    for (int icount = i; icount <= i + 1; icount++) {
                        // global id for the points
                        int node_gid = get_id(icount, jcount, 0, num_points_i, num_points_j);

                        // convert this_point index to the FE index convention
                        int this_index = convert_point_number_in_quad(this_point);

                        // store the points in this elem according the the finite
                        // element numbering convention
                        mesh.nodes_in_elem.host(elem_gid, this_index) = node_gid;

                        // increment the point counting index
                        this_point = this_point + 1;
                    } // end for icount
                } // end for jcount
            } // end for i
        } // end for j

        // update device side
        mesh.nodes_in_elem.update_device();

        // intialize corner variables
        int num_corners = num_elems * mesh.num_nodes_in_elem;
        mesh.initialize_corners(num_corners);
        corner.initialize(num_corners, num_dim);

        // Build connectivity
        mesh.build_connectivity();
    } // end build_2d_box

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_2d_polar
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_2d_polar(mesh_t& mesh, elem_t& elem, node_t& node, corner_t& corner, simulation_parameters_t& sim_param)
    {
        printf(" Creating a 2D polar mesh \n");

        int num_dim     = 2;
        int rk_num_bins = sim_param.dynamic_options.rk_num_bins;

        const double inner_radius = sim_param.mesh_input.inner_radius;
        const double outer_radius = sim_param.mesh_input.outer_radius;

        const double start_angle = PI / 180.0 * sim_param.mesh_input.starting_angle;
        const double end_angle   = PI / 180.0 * sim_param.mesh_input.ending_angle;

        const int num_elems_i = sim_param.mesh_input.num_radial_elems;
        const int num_elems_j = sim_param.mesh_input.num_angular_elems;

        const int num_points_i = num_elems_i + 1; // num points in x
        const int num_points_j = num_elems_j + 1; // num points in y

        const int num_nodes = num_points_i * num_points_j;

        const double dx = (outer_radius - inner_radius) / ((double)num_elems_i);  // len/(elems)
        const double dy = (end_angle - start_angle) / ((double)num_elems_j);  // len/(elems)

        const int num_elems = num_elems_i * num_elems_j;

        std::vector<double> origin = sim_param.mesh_input.origin;

        // --- 2D parameters ---
        const int num_faces_in_elem  = 4;  // number of faces in elem
        const int num_points_in_elem = 4;  // number of points in elem
        const int num_points_in_face = 2;  // number of points in a face
        const int num_edges_in_elem  = 4;  // number of edges in a elem

        // --- mesh node ordering ---
        // Convert ijk index system to the finite element numbering convention
        // for vertices in elem
        auto convert_point_number_in_quad = CArray<int>(4);
        convert_point_number_in_quad(0) = 0;
        convert_point_number_in_quad(1) = 1;
        convert_point_number_in_quad(2) = 3;
        convert_point_number_in_quad(3) = 2;

        // intialize node variables
        mesh.initialize_nodes(num_nodes);
        node.initialize(rk_num_bins, num_nodes, num_dim);

        // populate the point data structures
        for (int j = 0; j < num_points_j; j++) {
            for (int i = 0; i < num_points_i; i++) {
                // global id for the point
                int node_gid = get_id(i, j, 0, num_points_i, num_points_j);

                double r_i     = inner_radius + (double)i * dx;
                double theta_j = start_angle + (double)j * dy;

                // store the point coordinates
                node.coords(0, node_gid, 0) = origin[0] + r_i * cos(theta_j);
                node.coords(0, node_gid, 1) = origin[1] + r_i * sin(theta_j);
            } // end for i
        } // end for j

        for (int rk_level = 1; rk_level < rk_num_bins; rk_level++) {
            for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
                node.coords(rk_level, node_gid, 0) = node.coords(0, node_gid, 0);
                node.coords(rk_level, node_gid, 1) = node.coords(0, node_gid, 1);
            }
        }

        // intialize elem variables
        mesh.initialize_elems(num_elems, num_dim);
        elem.initialize(rk_num_bins, num_nodes, 3); // always 3D here, even for 2D

        // populate the elem center data structures
        for (int j = 0; j < num_elems_j; j++) {
            for (int i = 0; i < num_elems_i; i++) {
                // global id for the elem
                int elem_gid = get_id(i, j, 0, num_elems_i, num_elems_j);

                // store the point IDs for this elem where the range is
                // (i:i+1, j:j+1 for a linear quad
                int this_point = 0;

                for (int jcount = j; jcount <= j + 1; jcount++) {
                    for (int icount = i; icount <= i + 1; icount++) {
                        // global id for the points
                        int node_gid = get_id(icount, jcount, 0, num_points_i, num_points_j);

                        // convert this_point index to the FE index convention
                        int this_index = convert_point_number_in_quad(this_point);

                        // store the points in this elem according the the finite
                        // element numbering convention
                        mesh.nodes_in_elem.host(elem_gid, this_index) = node_gid;

                        // increment the point counting index
                        this_point = this_point + 1;
                    } // end for icount
                } // end for jcount
            } // end for i
        } // end for j

        // update device side
        mesh.nodes_in_elem.update_device();

        // intialize corner variables
        int num_corners = num_elems * mesh.num_nodes_in_elem;
        mesh.initialize_corners(num_corners);
        corner.initialize(num_corners, num_dim);

        // Build connectivity
        mesh.build_connectivity();
    } // end build_2d_box

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_3d_box
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_3d_box(mesh_t& mesh, elem_t& elem, node_t& node, corner_t& corner, simulation_parameters_t& sim_param)
    {
        printf(" Creating a 3D box mesh \n");

        const int num_dim = 3;

        const double lx = sim_param.mesh_input.length[0];
        const double ly = sim_param.mesh_input.length[1];
        const double lz = sim_param.mesh_input.length[2];

        const int num_elems_i = sim_param.mesh_input.num_elems[0];
        const int num_elems_j = sim_param.mesh_input.num_elems[1];
        const int num_elems_k = sim_param.mesh_input.num_elems[2];

        const int num_points_i = num_elems_i + 1; // num points in x
        const int num_points_j = num_elems_j + 1; // num points in y
        const int num_points_k = num_elems_k + 1; // num points in y

        const int num_nodes = num_points_i * num_points_j * num_points_k;

        const double dx = lx / ((double)num_elems_i);  // len/(num_elems_i)
        const double dy = ly / ((double)num_elems_j);  // len/(num_elems_j)
        const double dz = lz / ((double)num_elems_k);  // len/(num_elems_k)

        const int num_elems = num_elems_i * num_elems_j * num_elems_k;

        std::vector<double> origin = sim_param.mesh_input.origin;

        // --- 3D parameters ---
        const int num_faces_in_elem  = 6;  // number of faces in elem
        const int num_points_in_elem = 8;  // number of points in elem
        const int num_points_in_face = 4;  // number of points in a face
        const int num_edges_in_elem  = 12; // number of edges in a elem

        // --- mesh node ordering ---
        // Convert ijk index system to the finite element numbering convention
        // for vertices in elem
        auto convert_point_number_in_Hex = CArray<int>(8);
        convert_point_number_in_Hex(0) = 0;
        convert_point_number_in_Hex(1) = 1;
        convert_point_number_in_Hex(2) = 3;
        convert_point_number_in_Hex(3) = 2;
        convert_point_number_in_Hex(4) = 4;
        convert_point_number_in_Hex(5) = 5;
        convert_point_number_in_Hex(6) = 7;
        convert_point_number_in_Hex(7) = 6;

        int rk_num_bins = sim_param.dynamic_options.rk_num_bins;

        // intialize node variables
        mesh.initialize_nodes(num_nodes);
        node.initialize(rk_num_bins, num_nodes, num_dim);

        // --- Build nodes ---

        // populate the point data structures
        for (int k = 0; k < num_points_k; k++) {
            for (int j = 0; j < num_points_j; j++) {
                for (int i = 0; i < num_points_i; i++) {
                    // global id for the point
                    int node_gid = get_id(i, j, k, num_points_i, num_points_j);

                    // store the point coordinates
                    node.coords(0, node_gid, 0) = origin[0] + (double)i * dx;
                    node.coords(0, node_gid, 1) = origin[1] + (double)j * dy;
                    node.coords(0, node_gid, 2) = origin[2] + (double)k * dz;
                } // end for i
            } // end for j
        } // end for k

        for (int rk_level = 1; rk_level < rk_num_bins; rk_level++) {
            for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
                node.coords(rk_level, node_gid, 0) = node.coords(0, node_gid, 0);
                node.coords(rk_level, node_gid, 1) = node.coords(0, node_gid, 1);
                node.coords(rk_level, node_gid, 2) = node.coords(0, node_gid, 2);
            }
        }

        // intialize elem variables
        mesh.initialize_elems(num_elems, num_dim);
        elem.initialize(rk_num_bins, num_nodes, 3); // always 3D here, even for 2D

        // --- Build elems  ---

        // populate the elem center data structures
        for (int k = 0; k < num_elems_k; k++) {
            for (int j = 0; j < num_elems_j; j++) {
                for (int i = 0; i < num_elems_i; i++) {
                    // global id for the elem
                    int elem_gid = get_id(i, j, k, num_elems_i, num_elems_j);

                    // store the point IDs for this elem where the range is
                    // (i:i+1, j:j+1, k:k+1) for a linear hexahedron
                    int this_point = 0;
                    for (int kcount = k; kcount <= k + 1; kcount++) {
                        for (int jcount = j; jcount <= j + 1; jcount++) {
                            for (int icount = i; icount <= i + 1; icount++) {
                                // global id for the points
                                int node_gid = get_id(icount, jcount, kcount,
                                                  num_points_i, num_points_j);

                                // convert this_point index to the FE index convention
                                int this_index = convert_point_number_in_Hex(this_point);

                                // store the points in this elem according the the finite
                                // element numbering convention
                                mesh.nodes_in_elem.host(elem_gid, this_index) = node_gid;

                                // increment the point counting index
                                this_point = this_point + 1;
                            } // end for icount
                        } // end for jcount
                    }  // end for kcount
                } // end for i
            } // end for j
        } // end for k

        // update device side
        mesh.nodes_in_elem.update_device();

        // intialize corner variables
        int num_corners = num_elems * mesh.num_nodes_in_elem;
        mesh.initialize_corners(num_corners);
        corner.initialize(num_corners, num_dim);

        // Build connectivity
        mesh.build_connectivity();
    } // end build_3d_box

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_3d_HexN_box
    ///
    /// \brief <insert brief description>
    ///
    /// <Insert longer more detailed description which
    /// can span multiple lines if needed>
    ///
    /// \param <function parameter description>
    /// \param <function parameter description>
    /// \param <function parameter description>
    ///
    /// \return <return type and definition description if not void>
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_3d_HexN_box(mesh_t& mesh, elem_t& elem, node_t& node, corner_t& corner, simulation_parameters_t& sim_param)
    {
        printf(" ***** WARNING::  build_3d_HexN_box not yet implemented\n");
    }

    // -------------------------------------------------------
    // This gives the index value of the point or the elem
    // the elem = i + (j)*(num_points_i-1) + (k)*(num_points_i-1)*(num_points_j-1)
    // the point = i + (j)*num_points_i + (k)*num_points_i*num_points_j
    // --------------------------------------------------------
    //
    // Returns a global id for a given i,j,k
    int get_id(int i, int j, int k, int num_i, int num_j)
    {
        return i + j * num_i + k * num_i * num_j;
    }

    /**\brief Given (i,j,k) coordinates within the Lagrange hex, return an offset into the local connectivity (PointIds) array.
    *
    * The \a order parameter must point to an array of 3 integers specifying the order
    * along each axis of the hexahedron.
    */
    int PointIndexFromIJK(int i, int j, int k, const int* order)
    {
        bool ibdy = (i == 0 || i == order[0]);
        bool jbdy = (j == 0 || j == order[1]);
        bool kbdy = (k == 0 || k == order[2]);

        // How many boundaries do we lie on at once?
        int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0) + (kbdy ? 1 : 0);

        if (nbdy == 3) { // Vertex DOF
            // ijk is a corner node. Return the proper index (somewhere in [0,7]):
            return (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);
        }

        int offset = 8;
        if (nbdy == 2) { // Edge DOF
            if (!ibdy) { // On i axis
                return (i - 1) +
                       (j ? order[0] - 1 + order[1] - 1 : 0) +
                       (k ? 2 * (order[0] - 1 + order[1] - 1) : 0) +
                       offset;
            }
            if (!jbdy) { // On j axis
                return (j - 1) +
                       (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) +
                       (k ? 2 * (order[0] - 1 + order[1] - 1) : 0) +
                       offset;
            }
            // !kbdy, On k axis
            offset += 4 * (order[0] - 1) + 4 * (order[1] - 1);
            return (k - 1) + (order[2] - 1) * (i ? (j ? 3 : 1) : (j ? 2 : 0)) + offset;
        }

        offset += 4 * (order[0] - 1 + order[1] - 1 + order[2] - 1);
        if (nbdy == 1) { // Face DOF
            if (ibdy) { // On i-normal face
                return (j - 1) + ((order[1] - 1) * (k - 1)) + (i ? (order[1] - 1) * (order[2] - 1) : 0) + offset;
            }
            offset += 2 * (order[1] - 1) * (order[2] - 1);
            if (jbdy) { // On j-normal face
                return (i - 1) + ((order[0] - 1) * (k - 1)) + (j ? (order[2] - 1) * (order[0] - 1) : 0) + offset;
            }
            offset += 2 * (order[2] - 1) * (order[0] - 1);
            // kbdy, On k-normal face
            return (i - 1) + ((order[0] - 1) * (j - 1)) + (k ? (order[0] - 1) * (order[1] - 1) : 0) + offset;
        }

        // nbdy == 0: Body DOF
        offset += 2 * ((order[1] - 1) * (order[2] - 1) + (order[2] - 1) * (order[0] - 1) + (order[0] - 1) * (order[1] - 1));
        return offset + (i - 1) + (order[0] - 1) * ((j - 1) + (order[1] - 1) * ((k - 1)));
    }
};

#endif // end Header Guard