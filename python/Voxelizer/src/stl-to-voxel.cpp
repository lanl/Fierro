#include <pybind11/pybind11.h>
#include "stl-to-voxelvtk.h"

int add(int i, int j) {
    return i + j;
}

PYBIND11_MODULE(fierro_voxelizer, m) {
    m.doc() = "backend voxelizer for fierro GUI front-ends.";

    m.def("create_voxel_vtk", &Voxelizer::create_voxel_vtk, "Create a voxelized vtk from an STl");
}