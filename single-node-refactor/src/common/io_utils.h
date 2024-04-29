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
#include <cstring>
#include <sys/stat.h>

/////////////////////////////////////////////////////////////////////////////
///
/// \class MeshReader
///
/// \brief Class for simplifying reading meshes
///
/// This class contains the requisite functions required to read different 
/// mesh formats. The idea is to set the mesh file name, and parse the 
/// extension to decide which reader to use. Currently, only ensight .geo
/// files are supported.
///
/////////////////////////////////////////////////////////////////////////////
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
            }
        }

        // --- Read in the nodes in the mesh ---

        size_t num_nodes = 0;

        fscanf(in, "%lu", &num_nodes);
        printf("Number if nodes read in %lu\n", num_nodes);

        // initialize node variables
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


            }
        } // end for

        ch = (char)fgetc(in);

        // skip 1 line
        for (int j = 1; j <= 1; j++) {
            int i = 0;
            while ((ch = (char)fgetc(in)) != '\n') {
                i++;
            }
        }

        // --- read in the elements in the mesh ---
        size_t num_elem = 0;

        fscanf(in, "%lu", &num_elem);
        printf("Number of elements read in %lu\n", num_elem);

        // initialize elem variables
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

        // initialize corner variables
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

/////////////////////////////////////////////////////////////////////////////
///
/// \class MeshBuilder
///
/// \brief Class for building simple meshes
///
/// This class contains the requisite functions required to build simple
/// 2D and 3D Box meshes as well as 2D polar meshes. It uses the parsed
/// simulation parameters to decide what type of mesh to build.
///
/////////////////////////////////////////////////////////////////////////////
class MeshBuilder
{
public:

    MeshBuilder() {}

    ~MeshBuilder()
    {
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn build_mesh
    ///
    /// \brief Build a mesh for Fierro based on the input instructions
    ///
    /// \param Simulation mesh that is built
    /// \param Element state data
    /// \param Node state data
    /// \param Corner state data
    /// \param Simulation parameters
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
    /// \brief Builds an unstructured 2D rectilinear mesh
    ///
    /// \param Simulation mesh that is built
    /// \param Element state data
    /// \param Node state data
    /// \param Corner state data
    /// \param Simulation parameters
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
    /// \brief Builds an unstructured 2D polar mesh
    ///
    /// \param Simulation mesh that is built
    /// \param Element state data
    /// \param Node state data
    /// \param Corner state data
    /// \param Simulation parameters
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
    /// \brief Builds an unstructured 3D rectilinear mesh
    ///
    /// \param Simulation mesh that is built
    /// \param Element state data
    /// \param Node state data
    /// \param Corner state data
    /// \param Simulation parameters
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
     /// \brief Builds an unstructured high order 3D rectilinear mesh
    ///
    /// \param Simulation mesh that is built
    /// \param Element state data
    /// \param Node state data
    /// \param Corner state data
    /// \param Simulation parameters
    ///
    /////////////////////////////////////////////////////////////////////////////
    void build_3d_HexN_box(mesh_t& mesh, elem_t& elem, node_t& node, corner_t& corner, simulation_parameters_t& sim_param)
    {
        printf(" ***** WARNING::  build_3d_HexN_box not yet implemented\n");
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn get_id
    ///
    /// \brief This gives the index value of the point or the elem
    ///
    /// Assumes that the grid has an i,j,k structure
    /// the elem = i + (j)*(num_points_i-1) + (k)*(num_points_i-1)*(num_points_j-1)
    /// the point = i + (j)*num_points_i + (k)*num_points_i*num_points_j
    ///
    /// \param i index
    /// \param j index
    /// \param k index
    /// \param Number of i indices
    /// \param Number of j indices
    ///
    /////////////////////////////////////////////////////////////////////////////
    int get_id(int i, int j, int k, int num_i, int num_j)
    {
        return i + j * num_i + k * num_i * num_j;
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn PointIndexFromIJK
    ///
    /// \brief Given (i,j,k) coordinates within the Lagrange hex, return an 
    ///        offset into the local connectivity (PointIds) array.
    ///
    /// Assumes that the grid has an i,j,k structure
    /// the elem = i + (j)*(num_points_i-1) + (k)*(num_points_i-1)*(num_points_j-1)
    /// the point = i + (j)*num_points_i + (k)*num_points_i*num_points_j
    ///
    /// \param i index
    /// \param j index
    /// \param k index
    /// \param array of 3 integers specifying the order along each axis of the hexahedron
    ///
    /////////////////////////////////////////////////////////////////////////////
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

/////////////////////////////////////////////////////////////////////////////
///
/// \class MeshWriter
///
/// \brief Class for writing out a mesh with its associated state from Fierro
///
/// This class contains the requisite functions required to write out a mesh
/// with its associated state data from solvers in Fierro. Currently only ensight
/// outputs are supported
///
/////////////////////////////////////////////////////////////////////////////
class MeshWriter
{
private:
    int graphics_id = 0;

public:

    MeshWriter() {}

    ~MeshWriter()
    {
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn writes mesh with the format given in the input.yaml file
    ///
    /// \param Simulation mesh
    /// \param Element related state
    /// \param Node related state
    /// \param Corner related state
    /// \param Simulation input parameters
    ///
    /////////////////////////////////////////////////////////////////////////////
    void write_mesh(
    mesh_t&   mesh,
    elem_t&   elem,
    node_t&   node,
    corner_t& corner,
    simulation_parameters_t& sim_param,
    double time_value,
    CArray<double> graphics_times)
    {
        if (sim_param.output_options.format == output_options::vtk) {
            write_vtk(mesh, elem, node, corner, sim_param, time_value, graphics_times);
        }
        else if (sim_param.output_options.format == output_options::ensight) {
            write_ensight(mesh, elem, node, corner, sim_param, time_value, graphics_times);
        }
        else{
            std::cout << "**** MESH OUTPUT TYPE NOT SUPPORTED **** " << std::endl;
            std::cout << "Valid options are: " << std::endl;
            auto map = output_format_map;
            for (const auto& pair : map) {
                std::cout << "\t" << pair.first << std::endl;
            }
            throw std::runtime_error("**** MESH OUTPUT TYPE NOT SUPPORTED ****");
            return;
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn write_ensight
    ///
    /// \brief Writes an ensight output file
    ///
    /// \param Simulation mesh
    /// \param Element related state
    /// \param Node related state
    /// \param Corner related state
    /// \param Simulation parameters
    /// \param current time value
    /// \param Vector of all graphics output times
    ///
    /////////////////////////////////////////////////////////////////////////////
    void write_ensight(
        mesh_t&   mesh,
        elem_t&   elem,
        node_t&   node,
        corner_t& corner,
        simulation_parameters_t& sim_param,
        double time_value,
        CArray<double> graphics_times)
    {
        // Update host data
        elem.den.update_host();
        elem.pres.update_host();
        elem.stress.update_host();
        elem.sspd.update_host();
        elem.sie.update_host();
        elem.vol.update_host();
        elem.mass.update_host();
        elem.mat_id.update_host();

        node.coords.update_host();
        node.vel.update_host();
        node.mass.update_host();
        Kokkos::fence();

        const int num_scalar_vars = 9;
        const int num_vec_vars    = 2;

        std::string name_tmp;
        name_tmp = "Outputs_SGH";

        char* name = new char [name_tmp.length() + 1];
        std::strcpy(name, name_tmp.c_str());

        const char scalar_var_names[num_scalar_vars][15] = {
            "den", "pres", "sie", "vol", "mass", "sspd", "speed", "mat_id", "elem_switch"
        };

        const char vec_var_names[num_vec_vars][15] = {
            "pos", "vel"
        };

        // short hand
        const size_t num_nodes = mesh.num_nodes;
        const size_t num_elems = mesh.num_elems;
        const size_t num_dims  = mesh.num_dims;

        // save the cell state to an array for exporting to graphics files
        auto elem_fields = CArray<double>(num_elems, num_scalar_vars);
        int  elem_switch = 1;

        DCArrayKokkos<double> speed(num_elems);
        FOR_ALL(elem_gid, 0, num_elems, {
            double elem_vel[3]; // note:initialization with a list won't work
            elem_vel[0] = 0.0;
            elem_vel[1] = 0.0;
            elem_vel[2] = 0.0;
            // get the coordinates of the element center
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                elem_vel[0] += node.vel(1, mesh.nodes_in_elem(elem_gid, node_lid), 0);
                elem_vel[1] += node.vel(1, mesh.nodes_in_elem(elem_gid, node_lid), 1);
                if (mesh.num_dims == 3) {
                    elem_vel[2] += node.vel(1, mesh.nodes_in_elem(elem_gid, node_lid), 2);
                }
                else{
                    elem_vel[2] = 0.0;
                }
            } // end loop over nodes in element
            elem_vel[0] = elem_vel[0] / mesh.num_nodes_in_elem;
            elem_vel[1] = elem_vel[1] / mesh.num_nodes_in_elem;
            elem_vel[2] = elem_vel[2] / mesh.num_nodes_in_elem;

            double speed_sqrd = 0.0;
            for (int dim = 0; dim < num_dims; dim++) {
                speed_sqrd += elem_vel[dim] * elem_vel[dim];
            }
            speed(elem_gid) = sqrt(speed_sqrd);
        }); // end parallel for
        speed.update_host();

        // save the output scale fields to a single 2D array
        double e_switch = 1;
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            // save outputs
            elem_fields(elem_gid, 0) = elem.den.host(elem_gid);
            elem_fields(elem_gid, 1) = elem.pres.host(elem_gid);
            elem_fields(elem_gid, 2) = elem.sie.host(1, elem_gid);
            elem_fields(elem_gid, 3) = elem.vol.host(elem_gid);
            elem_fields(elem_gid, 4) = elem.mass.host(elem_gid);
            elem_fields(elem_gid, 5) = elem.sspd.host(elem_gid);
            elem_fields(elem_gid, 6) = speed.host(elem_gid);
            elem_fields(elem_gid, 7) = (double)elem.mat_id.host(elem_gid);
            elem_fields(elem_gid, 8) = e_switch;
            elem_switch *= -1;
        } // end for elements

        // save the vertex vector fields to an array for exporting to graphics files
        CArray<double> vec_fields(num_nodes, num_vec_vars, 3);

        for (size_t node_gid = 0; node_gid < num_nodes; node_gid++) {
            // position, var 0
            vec_fields(node_gid, 0, 0) = node.coords.host(1, node_gid, 0);
            vec_fields(node_gid, 0, 1) = node.coords.host(1, node_gid, 1);
            if (num_dims == 2) {
                vec_fields(node_gid, 0, 2) = 0.0;
            }
            else{
                vec_fields(node_gid, 0, 2) = node.coords.host(1, node_gid, 2);
            }

            // position, var 1
            vec_fields(node_gid, 1, 0) = node.vel.host(1, node_gid, 0);
            vec_fields(node_gid, 1, 1) = node.vel.host(1, node_gid, 1);
            if (num_dims == 2) {
                vec_fields(node_gid, 1, 2) = 0.0;
            }
            else{
                vec_fields(node_gid, 1, 2) = node.vel.host(1, node_gid, 2);
            }
        } // end for loop over vertices

        //  ---------------------------------------------------------------------------
        //  Setup of file and directoring for exporting
        //  ---------------------------------------------------------------------------
        FILE* out[20];   // the output files that are written to
        char  filename[128];

        struct stat st;

        if (stat("ensight", &st) != 0) {
            system("mkdir ensight");
        }

        if (stat("ensight/data", &st) != 0) {
            system("mkdir ensight/data");
        }

        //  ---------------------------------------------------------------------------
        //  Write the Geometry file
        //  ---------------------------------------------------------------------------
        sprintf(filename, "ensight/data/%s.%05d.geo", name, graphics_id);
        // filename has the full string

        out[0] = fopen(filename, "w");

        fprintf(out[0], "A graphics dump by Fierro \n");
        fprintf(out[0], "%s", "EnSight Gold geometry\n");
        fprintf(out[0], "%s", "node id assign\n");
        fprintf(out[0], "%s", "element id assign\n");

        fprintf(out[0], "part\n");
        fprintf(out[0], "%10d\n", 1);
        fprintf(out[0], "Mesh\n");

        // --- vertices ---
        fprintf(out[0], "coordinates\n");
        fprintf(out[0], "%10lu\n", num_nodes);

        // write all components of the point coordinates
        for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
            fprintf(out[0], "%12.5e\n", node.coords.host(1, node_gid, 0));
        }

        for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
            fprintf(out[0], "%12.5e\n", node.coords.host(1, node_gid, 1));
        }

        for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
            if (num_dims == 3) {
                fprintf(out[0], "%12.5e\n", node.coords.host(1, node_gid, 2));
            }
            else{
                fprintf(out[0], "%12.5e\n", 0.0);
            }
        }

        // --- elements ---
        if (num_dims == 3) {
            fprintf(out[0], "hexa8\n");
        }
        else{
            fprintf(out[0], "quad4\n");
        }
        fprintf(out[0], "%10lu\n", num_elems);

        // write all global point numbers for this cell
        for (int elem_gid = 0; elem_gid < num_elems; elem_gid++) {
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem; node_lid++) {
                fprintf(out[0], "%10lu\t", mesh.nodes_in_elem.host(elem_gid, node_lid) + 1); // note: node_gid starts at 1
            }
            fprintf(out[0], "\n");
        }

        fclose(out[0]);

        // ---------------------------------------------------------------------------
        // Write the Scalar variable files
        // ---------------------------------------------------------------------------

        // ensight_vars = (den, pres,...)
        for (int var = 0; var < num_scalar_vars; var++) {
            // write a scalar value
            sprintf(filename, "ensight/data/%s.%05d.%s", name, graphics_id, scalar_var_names[var]);

            out[0] = fopen(filename, "w");

            fprintf(out[0], "Per_elem scalar values\n");
            fprintf(out[0], "part\n");
            fprintf(out[0], "%10d\n", 1);
            if (num_dims == 3) {
                fprintf(out[0], "hexa8\n");
            }
            else{
                fprintf(out[0], "quad4\n");
            }

            for (int elem_id = 0; elem_id < num_elems; elem_id++) {
                fprintf(out[0], "%12.5e\n", elem_fields(elem_id, var));
            }

            fclose(out[0]);
        } // end for var

        //  ---------------------------------------------------------------------------
        //  Write the Vector variable files
        //  ---------------------------------------------------------------------------

        // ensight vector vars = (position, velocity, force)
        for (int var = 0; var < num_vec_vars; var++) {
            sprintf(filename, "ensight/data/%s.%05d.%s", name, graphics_id, vec_var_names[var]);

            out[0] = fopen(filename, "w");
            // fprintf(out[0],"Per_node vector values\n");
            // fprintf(out[0],"part\n");
            // fprintf(out[0],"%10d \n",1);
            // fprintf(out[0],"hexa8\n"); // WARNING, maybe bug here?

            fprintf(out[0], "Per_node vector values\n");
            fprintf(out[0], "part\n");
            fprintf(out[0], "%10d\n", 1);
            fprintf(out[0], "block\n");

            for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
                fprintf(out[0], "%12.5e\n", vec_fields(node_gid, var, 0));
            }

            for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
                fprintf(out[0], "%12.5e\n", vec_fields(node_gid, var, 1));
            }

            for (int node_gid = 0; node_gid < num_nodes; node_gid++) {
                fprintf(out[0], "%12.5e\n", vec_fields(node_gid, var, 2));
            }

            fclose(out[0]);
        } // end for var

        // ---------------------------------------------------------------------------
        // Write the case file
        // ---------------------------------------------------------------------------

        sprintf(filename, "ensight/%s.case", name);
        out[0] = fopen(filename, "w");

        fprintf(out[0], "FORMAT\n");
        fprintf(out[0], "type: ensight gold\n");
        fprintf(out[0], "GEOMETRY\n");

        sprintf(filename, "model: data/%s.*****.geo\n", name);
        fprintf(out[0], "%s", filename);
        fprintf(out[0], "VARIABLE\n");

        for (int var = 0; var < num_scalar_vars; var++) {
            sprintf(filename,
                    "scalar per element: %s data/%s.*****.%s\n",
                    scalar_var_names[var], name, scalar_var_names[var]);
            fprintf(out[0], "%s", filename);
        }

        for (int var = 0; var < num_vec_vars; var++) {
            sprintf(filename,
                    "vector per node: %s data/%s.*****.%s\n",
                    vec_var_names[var], name, vec_var_names[var]);
            fprintf(out[0], "%s", filename);
        }

        fprintf(out[0], "TIME\n");
        fprintf(out[0], "time set: 1\n");
        fprintf(out[0], "number of steps: %4d\n", graphics_id + 1);
        fprintf(out[0], "filename start number: 0\n");
        fprintf(out[0], "filename increment: 1\n");
        fprintf(out[0], "time values: \n");

        graphics_times(graphics_id) = time_value;

        for (int i = 0; i <= graphics_id; i++) {
            fprintf(out[0], "%12.5e\n", graphics_times(i));
        }
        fclose(out[0]);

        // ---------------------------------------------------------------------------
        // Done writing the graphics dump
        // ---------------------------------------------------------------------------

        // increment graphics id counter
        graphics_id++;

        delete[] name;

        return;
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn write_ensight
    ///
    /// \brief Writes a vtk output file
    ///
    /// \param Simulation mesh
    /// \param Element related state
    /// \param Node related state
    /// \param Corner related state
    /// \param Simulation parameters
    /// \param current time value
    /// \param Vector of all graphics output times
    ///
    /////////////////////////////////////////////////////////////////////////////
    void write_vtk(
        mesh_t&   mesh,
        elem_t&   elem,
        node_t&   node,
        corner_t& corner,
        simulation_parameters_t& sim_param,
        double time_value,
        CArray<double> graphics_times)
    {
        // Not yet supported
        throw std::runtime_error("**** VTK OUTPUT TYPE NOT YET SUPPORTED ****");

    }
};

#endif // end Header Guard