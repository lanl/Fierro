#pragma once

#include "utilities.h"
#include "error.h"
#include "common.h"

#include "swage.h"

// For VTK input
#include <vtkUnstructuredGridReader.h>     
#include <vtkXMLUnstructuredGridReader.h>

// For VTK output
#include "vtkNew.h"                        
#include "vtkPoints.h"                     
#include "vtkCellType.h"
#include "vtkCell.h"
#include "vtkHexahedron.h"         
#include "vtkLagrangeHexahedron.h"         
#include "vtkCellArray.h"                  
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkUnstructuredGrid.h"           
#include "vtkUnstructuredGridWriter.h"           
#include "vtkXMLUnstructuredGridWriter.h"  

typedef vtkSmartPointer<vtkUnstructuredGrid> VtkGrid;
typedef swage::mesh_t SwageMesh;

namespace swage2vtk {
  // Mapping from VTK local vertex ordering to SWAGE local vertex ordering
  void make_vtk_vert_map(int elem_order, int *map, bool reverse);

  // Reading VTK grid from file and initializing a SWAGE mesh from it
  VtkGrid read_grid(const std::string &file_name);
  void init_swage_mesh(const int &elem_order, const VtkGrid &grid, 
      SwageMesh &mesh);

  // Initializing a VTK grid from a SWAGE mesh and solution data and writing it
  // to a file
  VtkGrid init_vtk_grid(SwageMesh &mesh, const std::string &solution_name, 
      const CArray<double> &solution);
  void write_grid(const VtkGrid &grid, const std::string &file_name);
}
