#pragma once

#include "yaml-serializable.h"
#include <string>

SERIALIZABLE_ENUM(MESH_FORMAT,
    ensight,
    tecplot,
    vtk,
    ansys_dat,
    abaqus_inp
)

SERIALIZABLE_ENUM(ELEMENT_TYPE, 
  quad4, quad8, quad12,
  hex8, hex20, hex32
)

struct Input_Options : Yaml::ValidatedYaml, Yaml::DerivedFields {
  std::string mesh_file_name;
  MESH_FORMAT mesh_file_format;
  int p_order = 2;
  double unit_scaling = 1.0;
  bool topology_optimization_restart = false;
  bool shape_optimization_restart = false;

  ELEMENT_TYPE element_type = ELEMENT_TYPE::hex8;
  bool zero_index_base = false;

  // Non-serialized fields
  int words_per_line;
  int elem_words_per_line;

  /**
   * Determine a couple of file parsing parameters from the specified filetype.
  */
  void derive() {
    if (mesh_file_format == MESH_FORMAT::ansys_dat) {
      words_per_line = 4;
      elem_words_per_line = 11;
    } else {
      switch (mesh_file_format) {
        case MESH_FORMAT::ensight:
          words_per_line = 1;
          break;
        case MESH_FORMAT::abaqus_inp:
          words_per_line = 3;
          break;
        case MESH_FORMAT::vtk:
        case MESH_FORMAT::tecplot:
          words_per_line = 3;
          break;
        default:
          break;
      }

      switch (element_type) {
        case ELEMENT_TYPE::hex8:
          elem_words_per_line = 8;
          break;
        case ELEMENT_TYPE::hex20:
          elem_words_per_line = 20;
          break;
        case ELEMENT_TYPE::hex32:
          elem_words_per_line = 32;
          break;
        case ELEMENT_TYPE::quad4:
          elem_words_per_line = 4;
          break;
        case ELEMENT_TYPE::quad8:
          elem_words_per_line = 8;
          break;
        case ELEMENT_TYPE::quad12:
          elem_words_per_line = 12;
          break;
        // TODO: Implement handling for other element types
        default:
          throw Yaml::ConfigurationException("Unsupported element type `" + to_string(element_type) + "`.");
          break;
      }
    }

    mesh_file_name = std::filesystem::absolute(mesh_file_name).string();
  }
  
  /**
   * Ensures that the provided filepath is valid.
  */
  void validate() {
    Yaml::validate_filepath(mesh_file_name);
  }
};
IMPL_YAML_SERIALIZABLE_FOR(Input_Options, mesh_file_name, mesh_file_format, element_type, zero_index_base, p_order, unit_scaling,
topology_optimization_restart, shape_optimization_restart)
