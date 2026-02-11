#include <gtest/gtest.h>
#include "MeshBuilder.h"
#include "MeshBuilderInput.h"
#include "MeshIO.h"
#include <string>
#include <sstream>

TEST(MeshBuilderInput, BoxDeserialization) {
    std::string input = R"(
    output:
        file_type: VTK
    input:
        type: Box
        p_order: 1
        length: [1, 1, 1]
        num_elems: [10, 10, 10]
        origin: [0, 0, 0]
    )";


    MeshBuilderConfig in;
    Yaml::from_string_strict(input, in);


    Mesh mesh = MeshBuilder::build_mesh(in.input);

    EXPECT_EQ(mesh.p_order, 1);
    EXPECT_EQ(mesh.element_point_index.dims(0), 10*10*10);
    EXPECT_EQ(mesh.element_point_index.dims(1), 8);
    EXPECT_EQ(in.input->type, MeshType::Box);
}

TEST(MeshBuilderInput, CylinderDeserialization) {
    std::string input = R"(
    output:
        file_type: VTK
    input:
        type: Cylinder
        p_order: 1
        length: [1, 1, 1]
        num_elems: [10, 10, 10]
        origin: [0, 0, 0]
        inner_radius: 0.1
        start_angle: 1
    )";


    MeshBuilderConfig in;
    Yaml::from_string_strict(input, in);

    Mesh mesh = MeshBuilder::build_mesh(in.input);

    EXPECT_EQ(mesh.p_order, 1);
    EXPECT_EQ(mesh.element_point_index.dims(0), 10*10*10);
    EXPECT_EQ(mesh.element_point_index.dims(1), 8);
    EXPECT_EQ(in.input->type, MeshType::Cylinder);
}


TEST(MeshBuilder, WriteRead) {
    std::string input = R"(
    output:
        file_type: VTK
    input:
        type: Cylinder
        p_order: 1
        length: [1, 1, 1]
        num_elems: [1, 1, 1]
        origin: [0, 0, 0]
        inner_radius: 0.1
        start_angle: 1
    )";

    MeshBuilderConfig in;
    Yaml::from_string_strict(input, in);
    
    
    Mesh mesh = MeshBuilder::build_mesh(in.input);
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


TEST(MeshBuilder, ExampleCylinder) {
    MeshBuilderConfig in;
    Yaml::from_string_strict(MeshBuilderConfig::example_cylinder(), in);

    EXPECT_EQ(in.input->type, MeshType::Cylinder);
}

TEST(MeshBuilder, ExampleBox) {
    
    MeshBuilderConfig in;
    Yaml::from_string_strict(MeshBuilderConfig::example_box(), in);

    EXPECT_EQ(in.input->type, MeshType::Box);
}