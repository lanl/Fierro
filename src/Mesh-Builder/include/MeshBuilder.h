#pragma once

#include "Mesh.h"
#include <memory>
#include "MeshBuilderInput.h"

namespace MeshBuilder {
    Mesh build_mesh(std::shared_ptr<MeshBuilderInput> input);
    void build_mesh_from_file(std::string mesh_file);
}
