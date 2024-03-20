#include <pybind11/pybind11.h>
#include "stl-to-voxelvtk.h"

PYBIND11_MODULE(fierro_voxelizer, m) {
    m.doc() = "backend voxelizer for fierro GUI front-ends.";

    m.def("create_voxel_vtk", &Voxelizer::create_voxel_vtk_GUI, "Create a voxelized vtk from an STl");
}
