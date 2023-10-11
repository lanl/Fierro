#include <gtest/gtest.h>
#include "MeshBuilder.h"
#include "MeshBuilderInput.h"
#include "MeshIO.h"
#include <string>
#include <sstream>

TEST(MeshBuilderInput, BoxDeserialization) {
    std::string input = R"(
    type: Box
    file_type: VTK
    p_order: 1
    length: [1, 1, 1]
    num_elems: [10, 10, 10]
    origin: [0, 0, 0]
    )";


    std::shared_ptr<MeshBuilderInput> in;
    Yaml::from_string_strict(input, in);


    Mesh mesh = MeshBuilder::build_mesh(in);

    EXPECT_EQ(mesh.p_order, 1);
    EXPECT_EQ(mesh.element_point_index.dims(0), 10*10*10);
    EXPECT_EQ(mesh.element_point_index.dims(1), 8);
    EXPECT_EQ(in->type, MeshType::Box);
}

TEST(MeshBuilderInput, CylinderDeserialization) {
    std::string input = R"(
    type: Cylinder
    file_type: VTK
    p_order: 1
    length: [1, 1, 1]
    num_elems: [10, 10, 10]
    origin: [0, 0, 0]
    inner_radius: 0.1
    start_angle: 1
    )";


    std::shared_ptr<MeshBuilderInput> in;
    Yaml::from_string_strict(input, in);

    Mesh mesh = MeshBuilder::build_mesh(in);

    EXPECT_EQ(mesh.p_order, 1);
    EXPECT_EQ(mesh.element_point_index.dims(0), 10*10*10);
    EXPECT_EQ(mesh.element_point_index.dims(1), 8);
    EXPECT_EQ(in->type, MeshType::Cylinder);
}


TEST(MeshBuilder, WriteRead) {
    std::string input = R"(
    type: Cylinder
    file_type: VTK
    p_order: 1
    length: [1, 1, 1]
    num_elems: [1, 1, 1]
    origin: [0, 0, 0]
    inner_radius: 0.1
    start_angle: 1
    )";

    std::shared_ptr<MeshBuilderInput> in;
    Yaml::from_string_strict(input, in);
    
    Mesh mesh = MeshBuilder::build_mesh(in);
    std::stringstream buffer;

    MeshIO::write_vtk(buffer, mesh);

    Mesh read_back = MeshIO::read_vtk(buffer);

    EXPECT_EQ(mesh.num_dim, read_back.num_dim);
    for (size_t i = 0; i < mesh.points.dims(0); i++)
        for (size_t j = 0; j < mesh.points.dims(1); j++)
            EXPECT_NEAR(mesh.points(i, j), read_back.points(i, j), 1e-5);
    for (size_t i = 0; i < mesh.element_point_index.dims(0); i++)
        for (size_t j = 0; j < mesh.element_point_index.dims(1); j++)
            EXPECT_EQ(mesh.element_point_index(i, j), read_back.element_point_index(i, j));
}