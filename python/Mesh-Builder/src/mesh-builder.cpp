#include <pybind11/pybind11.h>
#include "MeshBuilder.h"


PYBIND11_MODULE(fierro_mesh_builder, m) {
    m.doc() = "backend mesh builder for fierro GUI front-ends.";

    m.def("build_mesh_from_file", &MeshBuilder::build_mesh_from_file, "Build a mesh from .yaml input file");
}
