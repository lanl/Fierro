
//==============================================================================
/*

 mesh-builder.cpp
 
 This function builds a mesh
 -------------------------------------------------------------------------------
 
 Created by Nathaniel on 12/13/22.
 
 
 Finite Element local point numbering for a Hexahedral
 
 
    7--------6
   /|       /|
  / |      / |
 4--------5  |
 |  |     |  |
 |  |     |  |
 |  3-----|--2
 | /      | /
 |/       |/
 0--------1
 
 
 The number of a Hexahedral element using i,j,k is
 
 
    6--------7
   /|       /|
  / |      / |
 2--------3  |
 |  |     |  |
 |  |     |  |
 |  4-----|--5
 | /      | /
 |/       |/
 0--------1
 
 
 
 The number of a quad element using FE convention is
 
 
        y
       /
      /

    3--------2
   /        /
  /        /
 0--------1  ---> x


The number of a quad element using i,j is

        y
       /
      /
 
    2--------3
   /        /
  /        /
 0--------1  ---> x

 
 */
//==============================================================================
//
//  if building standalone
//
//     g++ --std=c++17 mesh-builder.cpp Yaml.cpp
//
//==============================================================================

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <map>
#include <exception>

#include "matar.h"
#include "yaml-serializable.h"
#include "MeshBuilderInput.h"
#include "MeshBuilder.h"
#include "MeshIO.h"
#include "Mesh.h"

#include <cmath>

using namespace mtr;

//==============================================================================
//   Functions
//==============================================================================

// -------------------------------------------------------
// This gives the index value of the point or the elem
// the elem = i + (j)*(num_points_i-1) + (k)*(num_points_i-1)*(num_points_j-1)
// the point = i + (j)*num_points_i + (k)*num_points_i*num_points_j
//--------------------------------------------------------
//
// Returns a global id for a given i,j,k
int get_id(int i, int j, int k, int num_i, int num_j) {
   return i + j*num_i + k*num_i*num_j;
}

//==============================================================================
//    Main
//==============================================================================

Mesh build_rectilinear(const Input_Rectilinear& input) {
    int num_points_in_elem = std::pow(input.p_order + 1, input.num_dims);
    
    auto mesh = Mesh();
    mesh.num_dim = input.num_dims;
    mesh.points = CArray <double> (input.total_points, input.num_dims);
    mesh.element_point_index = CArray <int> (input.total_elems, num_points_in_elem);
    
    // --- Build nodes ---
    // populate the point data structures
    for (int k = 0; k < input.num_points[2]; k++) {
        for (int j = 0; j < input.num_points[1]; j++){
            for (int i = 0; i < input.num_points[0]; i++){
                // global id for the point
                int point_id = get_id(i, j, k, input.num_points[0], input.num_points[1]);
                
                // store the point coordinates
                mesh.points(point_id, 0) = input.lower_bound[0] + (double)i*input.delta[0];
                mesh.points(point_id, 1) = input.lower_bound[1] + (double)j*input.delta[1];
                if (input.num_dims == 3)
                    mesh.points(point_id, 2) = input.lower_bound[2] + (double)k*input.delta[2];
            }
        }
    }

    // --- Build elems  ---
    int p_order = input.p_order;
    size_t k_local_max = input.num_dims < 3 ? 0 : p_order;
    // populate the elem center data structures
    for (size_t k = 0; k < input.num_elems[2]; k++) {
        for (size_t j = 0; j < input.num_elems[1]; j++){
            for (size_t i = 0; i < input.num_elems[0]; i++){
                // global id for the elem
                size_t elem_id = get_id(i, j, k, input.num_elems[0], input.num_elems[1]);
                
                // store the point IDs for this elem where the range is
                // (i:(i+1)*Pn_order, j:(j+1)*Pn_order, k:(k+1)*Pn_order) for a Pn hexahedron
                size_t point_id_local = 0;
                size_t i_offset = i * p_order,
                       j_offset = j * p_order,
                       k_offset = k * p_order;

                for (size_t k_local = 0; k_local <= k_local_max; k_local++) {
                    for (size_t j_local = 0; j_local <= p_order; j_local++) {
                        for (size_t i_local = 0; i_local <= p_order; i_local++) {
                            // Modulus is included to allow for 
                            // cyclic connectivity structures.
                            int point_id = get_id(
                                (i_local + i_offset) % input.num_points[0], 
                                (j_local + j_offset) % input.num_points[1], 
                                (k_local + k_offset) % input.num_points[2],
                                input.num_points[0], input.num_points[1]
                            );
                            mesh.element_point_index(elem_id, point_id_local) = point_id;
                            point_id_local++;
                        }
                    }
                }
            }
        }
    }

    mesh.p_order = input.p_order;
    mesh.element_types = mtr::CArray<int>(mesh.element_point_index.dims(0));
    int element_type = -1;
    if (input.p_order == 1) {
        if (mesh.num_dim == 2)
            element_type = 9;  // Linear Quad
        if (mesh.num_dim == 3)
            element_type = 12; // Linear Hex
    } else if (input.p_order > 1) {
        if (mesh.num_dim == 2)
            element_type = 70;  // Lagrange Quadrilateral 
        if (mesh.num_dim == 3)
            element_type = 72;  // Lagrange Hexahedron 
    }

    if (element_type < 0)
        throw std::runtime_error(
            "Could not determine element type. mesh.num_dim: "
            + std::to_string(mesh.num_dim) 
            + ", input.p_order: " + std::to_string(input.p_order) + "."
        );

    for (size_t i = 0; i < mesh.element_point_index.dims(0); i++)
        mesh.element_types(i) = element_type;

    return mesh;
}

/**
 * \brief Given a mesh in cylindrical space, convert it to a mesh in rectilinear space.
 * 
 * The coordinates of the points before this function are <r, theta, z>.
 * The coordinates of the points after this function are <x, y, z>.
*/
void cylinder_transform(Mesh& mesh, const Input_Cylinder& input) {
    for (size_t i = 0; i < mesh.points.dims(0); i++) {
        double r_i = mesh.points(i, 0);
        double theta = mesh.points(i, 1);
        mesh.points(i, 0) = input.origin[0] + r_i * cos(theta);
        mesh.points(i, 1) = input.origin[1] + r_i * sin(theta);
    }
}

/**
 * \brief Translates a mesh a distance `tx` along its axes.
*/
void translate(Mesh& mesh, const std::vector<double>& tx) {
    for (size_t i = 0; i < mesh.points.dims(0); i++)
        for (size_t j = 0; j < mesh.points.dims(1); j++)
            mesh.points(i, j) += tx[j];
}

/**
 * \brief Builds a mesh from a mesh builder input file.
*/
Mesh MeshBuilder::build_mesh(std::shared_ptr<MeshBuilderInput> input) {
    Mesh mesh;
    switch (input->type) {
    case MeshType::Box:
        mesh = build_rectilinear(*std::dynamic_pointer_cast<Input_Rectilinear>(input));
        break;
    case MeshType::Cylinder:
        mesh = build_rectilinear(*std::dynamic_pointer_cast<Input_Rectilinear>(input));
        cylinder_transform(mesh, *std::dynamic_pointer_cast<Input_Cylinder>(input));
        break;
    default:
        throw std::runtime_error("Unsupported mesh shape.");
    }

    translate(mesh, input->origin);

    return mesh;
}

void MeshBuilder::build_mesh_from_file(std::string mesh_file){
    MeshBuilderConfig config;
    Yaml::from_file_strict(mesh_file, config);
    
    Mesh mesh = MeshBuilder::build_mesh(config.input);
    
    switch (config.output.file_type) {
        case FileType::Ensight:
            MeshIO::write_ensight(config.output.name, mesh, true);
            break;
        case FileType::VTK:
            MeshIO::write_vtk(config.output.name, config.output.file_location, mesh, true);
            break;
    }
}
