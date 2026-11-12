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
#include "mesh_builder.hpp"

#include "region.hpp"
#include "string_utils.hpp"

#include <map>
#include <memory>
#include <cstring>
#include <sys/stat.h>
#include <iostream>
#include <regex>
#include <fstream>
#include <sstream>
#include <mpi.h>


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
static int get_id(int i, int j, int k, int num_i, int num_j)
{
    return i + j * num_i + k * num_i * num_j;
}


/////////////////////////////////////////////////////////////////////////////
///
/// \fn build_mesh
///
/// \brief Build a mesh for Fierro based on the input instructions
///
/// \param mesh Simulation mesh to be built
/// \param node_coords Node coordinates
/// \param SimulationParamaters Simulation parameters
///
/////////////////////////////////////////////////////////////////////////////
void MeshBuilder::build_mesh(swage::Mesh_t& mesh,
                MPICArrayKokkos<real_t>& node_coords,
                SimulationParameters_t& SimulationParamaters,
                bool HexN)
{
    if (SimulationParamaters.MeshInput.num_dims == 2) {
        if (SimulationParamaters.MeshInput.type == mesh_input::Polar) {
            build_2d_polar(mesh, node_coords, SimulationParamaters);
        }
        else if (SimulationParamaters.MeshInput.type == mesh_input::Box) {
            build_2d_box(mesh, node_coords, SimulationParamaters);
        }
        else{
            std::cout << "**** 2D MESH TYPE NOT SUPPORTED **** " << std::endl;
            std::cout << "Valid options are: " << std::endl;
            auto map = mesh_input_type_map;
            for (const auto& pair : map) {
                std::cout << "\t" << pair.first << std::endl;
            }
            throw std::runtime_error("**** 2D MESH TYPE NOT SUPPORTED ****");
        }
    }
    else if (SimulationParamaters.MeshInput.num_dims == 3) {
        if (!HexN) {
            build_3d_box(mesh, node_coords, SimulationParamaters);
        }
        else {
            build_3d_HexN_box(mesh, node_coords, SimulationParamaters);
        }
    }
    else{
        throw std::runtime_error("**** ONLY 2D RZ OR 3D MESHES ARE SUPPORTED ****");
    }

    int num_corners = mesh.num_elems * mesh.num_nodes_in_elem;
    mesh.build_connectivity();
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
void MeshBuilder::build_2d_box(swage::Mesh_t& mesh,
                  MPICArrayKokkos<real_t>& node_coords,
                  SimulationParameters_t& SimulationParamaters) const
{
    printf("Creating a 2D box mesh \n");

    const int num_dim = 2;

    const real_t lx = SimulationParamaters.MeshInput.length[0];
    const real_t ly = SimulationParamaters.MeshInput.length[1];

    const int num_elems_i = SimulationParamaters.MeshInput.num_elems[0];
    const int num_elems_j = SimulationParamaters.MeshInput.num_elems[1];

    const int num_points_i = num_elems_i + 1; // num points in x
    const int num_points_j = num_elems_j + 1; // num points in y

    const int num_nodes = num_points_i * num_points_j;

    const real_t dx = lx / ((double)num_elems_i);  // len/(num_elems_i)
    const real_t dy = ly / ((double)num_elems_j);  // len/(num_elems_j)

    const int num_elems = num_elems_i * num_elems_j;

    std::vector<double> origin(num_dim);
    // SimulationParamaters.MeshInput.origin.update_host();
    for (int i = 0; i < num_dim; i++) { origin[i] = SimulationParamaters.MeshInput.origin[i]; }

    // --- 2D parameters ---
    // const int num_faces_in_elem  = 4;  // number of faces in elem
    // const int num_points_in_elem = 4;  // number of points in elem
    // const int num_points_in_face = 2;  // number of points in a face
    // const int num_edges_in_elem  = 4;  // number of edges in a elem

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

    node_coords = MPICArrayKokkos<real_t>(num_nodes, num_dim, "node_coordinates_in_mesh_io");

    // --- Build nodes ---

    // populate the point data structures
    for (int j = 0; j < num_points_j; j++) {
        for (int i = 0; i < num_points_i; i++) {
            // global id for the point
            int node_gid = get_id(i, j, 0, num_points_i, num_points_j);

            // store the point coordinates
            node_coords.host(node_gid, 0) = origin[0] + (double)i * dx;
            node_coords.host(node_gid, 1) = origin[1] + (double)j * dy;
        } // end for i
    } // end for j


    node_coords.update_device();

    // initialize elem variables
    mesh.initialize_elems(num_elems);

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
void MeshBuilder::build_2d_polar(swage::Mesh_t& mesh,
                    MPICArrayKokkos<real_t>& node_coords,
                    SimulationParameters_t& SimulationParamaters) const
{
    printf("Creating a 2D polar mesh \n");

    int num_dim     = 2;

    const real_t inner_radius = SimulationParamaters.MeshInput.inner_radius;
    const real_t outer_radius = SimulationParamaters.MeshInput.outer_radius;

    const real_t start_angle = PI / real_t(180.0) * SimulationParamaters.MeshInput.starting_angle;
    const real_t end_angle   = PI / real_t(180.0) * SimulationParamaters.MeshInput.ending_angle;

    const int num_elems_i = SimulationParamaters.MeshInput.num_radial_elems;
    const int num_elems_j = SimulationParamaters.MeshInput.num_angular_elems;

    const int num_points_i = num_elems_i + 1; // num points in x
    const int num_points_j = num_elems_j + 1; // num points in y

    const int num_nodes = num_points_i * num_points_j;

    const real_t dx = (outer_radius - inner_radius) / ((double)num_elems_i);  // len/(elems)
    const real_t dy = (end_angle - start_angle) / ((double)num_elems_j);  // len/(elems)

    const int num_elems = num_elems_i * num_elems_j;

    std::vector<double> origin(num_dim);

    for (int i = 0; i < num_dim; i++) { origin[i] = SimulationParamaters.MeshInput.origin[i]; }

    // --- 2D parameters ---
    // const int num_faces_in_elem  = 4;  // number of faces in elem
    // const int num_points_in_elem = 4;  // number of points in elem
    // const int num_points_in_face = 2;  // number of points in a face
    // const int num_edges_in_elem  = 4;  // number of edges in a elem

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
    node_coords = MPICArrayKokkos<real_t>(num_nodes, num_dim, "node_coordinates_in_mesh_io");

    // populate the point data structures
    for (int j = 0; j < num_points_j; j++) {
        for (int i = 0; i < num_points_i; i++) {
            // global id for the point
            int node_gid = get_id(i, j, 0, num_points_i, num_points_j);

            real_t r_i     = inner_radius + (double)i * dx;
            real_t theta_j = start_angle + (double)j * dy;

            // store the point coordinates
            node_coords.host(node_gid, 0) = origin[0] + r_i * Kokkos::cos(theta_j);
            node_coords.host(node_gid, 1) = origin[1] + r_i * Kokkos::sin(theta_j);

            if(node_coords.host(node_gid, 0) < real_t(0.0)){
                throw std::runtime_error("**** NODE RADIUS FOR RZ MESH MUST BE POSITIVE ****");
            }

        } // end for i
    } // end for j


    node_coords.update_device();

    // initialize elem variables
    mesh.initialize_elems(num_elems);

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
void MeshBuilder::build_3d_box(swage::Mesh_t& mesh,
                  MPICArrayKokkos<real_t>& node_coords,
                  SimulationParameters_t& SimulationParamaters) const
{
    printf("Creating a 3D box mesh \n");

    const int num_dim = 3;

    // SimulationParamaters.MeshInput.length.update_host();
    const real_t lx = SimulationParamaters.MeshInput.length[0];
    const real_t ly = SimulationParamaters.MeshInput.length[1];
    const real_t lz = SimulationParamaters.MeshInput.length[2];

    // SimulationParamaters.MeshInput.num_elems.update_host();
    const int num_elems_i = SimulationParamaters.MeshInput.num_elems[0];
    const int num_elems_j = SimulationParamaters.MeshInput.num_elems[1];
    const int num_elems_k = SimulationParamaters.MeshInput.num_elems[2];

    const int num_points_i = num_elems_i + 1; // num points in x
    const int num_points_j = num_elems_j + 1; // num points in y
    const int num_points_k = num_elems_k + 1; // num points in y

    const int num_nodes = num_points_i * num_points_j * num_points_k;

    const real_t dx = lx / ((double)num_elems_i);  // len/(num_elems_i)
    const real_t dy = ly / ((double)num_elems_j);  // len/(num_elems_j)
    const real_t dz = lz / ((double)num_elems_k);  // len/(num_elems_k)

    const int num_elems = num_elems_i * num_elems_j * num_elems_k;

    std::vector<double> origin(num_dim);
    // SimulationParamaters.MeshInput.origin.update_host();
    for (int i = 0; i < num_dim; i++) { origin[i] = SimulationParamaters.MeshInput.origin[i]; }

    // --- 3D parameters ---
    // const int num_faces_in_elem  = 6;  // number of faces in elem
    // const int num_points_in_elem = 8;  // number of points in elem
    // const int num_points_in_face = 4;  // number of points in a face
    // const int num_edges_in_elem  = 12; // number of edges in a elem


    // initialize mesh node variables
    mesh.initialize_nodes(num_nodes);
    node_coords = MPICArrayKokkos<real_t>(num_nodes, num_dim, "node_coordinates_in_mesh_io");

    // --- Build nodes ---

    // populate the point data structures
    for (int k = 0; k < num_points_k; k++) {
        for (int j = 0; j < num_points_j; j++) {
            for (int i = 0; i < num_points_i; i++) {
                // global id for the point
                int node_gid = get_id(i, j, k, num_points_i, num_points_j);

                // store the point coordinates
                node_coords.host(node_gid, 0) = origin[0] + (double)i * dx;
                node_coords.host(node_gid, 1) = origin[1] + (double)j * dy;
                node_coords.host(node_gid, 2) = origin[2] + (double)k * dz;
            } // end for i
        } // end for j
    } // end for k


    node_coords.update_device();

    // initialize elem variables
    mesh.initialize_elems(num_elems);

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
                            int this_index = this_point; //convert_point_number_in_Hex(this_point);

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
void MeshBuilder::build_3d_HexN_box(swage::Mesh_t& mesh,
                       MPICArrayKokkos<real_t>& node_coords,
                       SimulationParameters_t& SimulationParamaters) const
{
    const int num_dim = 3;

    // SimulationParamaters.MeshInput.length.update_host();
    const real_t lx = SimulationParamaters.MeshInput.length[0];
    const real_t ly = SimulationParamaters.MeshInput.length[1];
    const real_t lz = SimulationParamaters.MeshInput.length[2];

    // SimulationParamaters.MeshInput.num_elems.update_host();
    const int num_elems_i = SimulationParamaters.MeshInput.num_elems[0];
    const int num_elems_j = SimulationParamaters.MeshInput.num_elems[1];
    const int num_elems_k = SimulationParamaters.MeshInput.num_elems[2];

    // creating zones for the Pn order
    const int Pn_order = SimulationParamaters.MeshInput.p_order;
    
    if (Pn_order > 19) {
        printf("Fierro DG and RD solvers are only valid for elements up to Pn = 19 \n");
        return;
    }

    const int num_zones_i = Pn_order*num_elems_i;
    const int num_zones_j = Pn_order*num_elems_j;
    const int num_zones_k = Pn_order*num_elems_k;
    
    const int num_points_i = num_zones_i+1; // num points in x accounting for Pn
    const int num_points_j = num_zones_j+1; // num points in y accounting for Pn
    const int num_points_k = num_zones_k+1; // num points in y accounting for Pn
    
    
    const real_t dx = lx/((double)num_zones_i);  // len/(num_zones_i)
    const real_t dy = ly/((double)num_zones_j);  // len/(num_zones_j)
    const real_t dz = lz/((double)num_zones_k);  // len/(num_zones_k)
    
    const int num_elems = num_elems_i*num_elems_j*num_elems_k;
    // const int num_zones = num_zones_i*num_zones_j*num_zones_k; // accounts for Pn

    std::vector<double> origin(num_dim);
    for (int i = 0; i < num_dim; i++) { origin[i] = SimulationParamaters.MeshInput.origin[i]; }

    // --- 3D parameters ---
    // const int num_faces_in_zone = 6;   // number of faces in zone
    // const int num_points_in_zone = 8;  // number of points in zone
    // const int num_points_in_face = 4;  // number of points in a face
    
    // p_order   = 1, 2, 3, 4, 5
    // num_nodes = 2, 3, 4, 5, 6
    const int num_1D_points = Pn_order+1;
    const int num_points_in_elem = num_1D_points*num_1D_points*num_1D_points;
       
    
    // --- elem ---
    auto elem_coords = CArray<real_t> (num_elems, num_dim);
    auto elem_point_list = CArray <int> (num_elems, num_points_in_elem);
    
    
    // --- point ---
    int num_points = num_points_i * num_points_j * num_points_k;


    // initialize node variables
    mesh.initialize_nodes(num_points);

    node_coords = MPICArrayKokkos<real_t>(num_points, num_dim, "node_coordinates_in_mesh_io");

    // populate the point data structures
    for (int k = 0; k < num_points_k; k++){
        for (int j = 0; j < num_points_j; j++){
            for (int i = 0; i < num_points_i; i++){

            
                // global id for the point
                int node_gid = get_id(i, j, k, num_points_i, num_points_j);

                // store the point coordinates
                node_coords.host(node_gid, 0) = origin[0] + (double)i * dx;
                node_coords.host(node_gid, 1) = origin[1] + (double)j * dy;
                node_coords.host(node_gid, 2) = origin[2] + (double)k * dz;
                
            } // end for k
        } // end for i
    } // end for j


    node_coords.update_device();


    // initialize elem variables
    mesh.initialize_elems_Pn(num_elems, SimulationParamaters.MeshInput.p_order, 2*SimulationParamaters.MeshInput.p_order);

    // --- Build elems  ---
    
    // populate the elem center data structures accounting for Pn
    for (int k=0; k<num_elems_k; k++){
        for (int j=0; j<num_elems_j; j++){
            for (int i=0; i<num_elems_i; i++){
              
                // global id for the elem
                size_t elem_gid = get_id(i, j, k, num_elems_i, num_elems_j);
                
                // store the point IDs for this elem where the range is
                // (i:i+1, j:j+1, k:k+1) for a linear hexahedron
                // (i:(i+1)*Pn_order, j:(j+1)*Pn_order, k:(k+1)*Pn_order) for a Pn hexahedron
                int node_lid = 0;
                
                int k_local = 0;
                for (int kcount=k*Pn_order; kcount<=(k+1)*Pn_order; kcount++){
                    
                    int j_local = 0;
                    for (int jcount=j*Pn_order; jcount<=(j+1)*Pn_order; jcount++){
                        
                        int i_local = 0;
                        for (int icount=i*Pn_order; icount<=(i+1)*Pn_order; icount++){
                            
                            // global id for the points
                            size_t node_gid = get_id(icount, jcount, kcount,
                                              num_points_i, num_points_j);

                            // Saved using i,j,k indexing
                            mesh.nodes_in_elem.host(elem_gid, node_lid) = node_gid;
                            
                            // increment the point counting index
                            node_lid = node_lid + 1;
                            
                            i_local++;
                        } // end for icount
                        
                        j_local++;
                    } // end for jcount
                    
                    k_local ++;
                }  // end for kcount
            } // end for i
        } // end for j
    } // end for k

    // update device side
    mesh.nodes_in_elem.update_device();
}
MeshBuilder::MeshBuilder()
{}

MeshBuilder::~MeshBuilder()
{
}
