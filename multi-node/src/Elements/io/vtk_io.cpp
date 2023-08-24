#include "vtk_io.h"

#include <ostream>

/*
 * Generate a map from the local vertex ordering in SWAGE to the local vertex
 * ordering in VTK for high-order elements
 */
void swage2vtk::make_vtk_vert_map(int elem_order, int *map, bool reverse) {
  // Create a VTK N-Hexahedron reference element
  int num_dims = 3;
  int num_verts_per_dim = elem_order + 1;
  int num_verts_per_elem = std::pow(num_verts_per_dim, num_dims);

  vtkNew<vtkLagrangeHexahedron> hex;
  hex->SetOrder(elem_order, elem_order, elem_order);
  hex->GetPointIds()->SetNumberOfIds(num_verts_per_elem);
  hex->GetPoints()->SetNumberOfPoints(num_verts_per_elem);

  SizeType radices[num_dims] = {SizeType(num_verts_per_dim), 
      SizeType(num_verts_per_dim), SizeType(num_verts_per_dim)};

  // Loop over the IJK coordinates to generate the map
  for (int I = 0; I < num_verts_per_dim; I++) {
    for (int J = 0; J < num_verts_per_dim; J++) {
      for (int K = 0; K < num_verts_per_dim; K++) {
        SizeType IJK[num_dims] = {SizeType(I), SizeType(J), SizeType(K)};

        int swage_vert_id = int(common::mixed_radix_to_base_10(num_dims, 
            radices, IJK));

        int vtk_vert_id = hex->PointIndexFromIJK(I, J, K); 

        if (!reverse) {
          map[swage_vert_id] = vtk_vert_id;
        } else {
          map[vtk_vert_id] = swage_vert_id;
        }
      }
    }
  }
}

/*
 * Given the name of a VTK unstructured grid file (extensions .vtu or .vtk),
 * read the file into a VTK unstructured grid object and return it
 */
VtkGrid swage2vtk::read_grid(const std::string &file_name) {
  // Initialize smart pointer to unstructured grid
  VtkGrid grid;

  // Extract the extension from the file name
  std::string extension = "";
  if (file_name.find_last_of(".") != std::string::npos)
    extension = file_name.substr(file_name.find_last_of("."));

  // Drop the case of the extension
  std::transform(extension.begin(), extension.end(), extension.begin(), 
      ::tolower);

  // Initialize the appropriate reader to the type indicated by the extension
  // and read the grid
  if (extension == ".vtu") {
    vtkNew<vtkXMLUnstructuredGridReader> reader;
    reader->SetFileName(file_name.c_str());
    reader->Update();
    grid = reader->GetOutput();
  } else if (extension == ".vtk") {
    vtkNew<vtkUnstructuredGridReader> reader;
    reader->SetFileName(file_name.c_str());
    reader->Update();
    grid = reader->GetOutput();
  }

  return grid;
}

/*
 * Provided a valid VTK unstructured grid object, extract the vertex
 * coordinates and cell-vertex connectivity and store these in a SWAGE mesh
 * object
 *
 * Notes
 * -----
 * This routine currently assumes: 
 *  - all cells in VTK grid are of the same type (VTK_LAGRANGE_HEXAHEDRON or
 *    VTK_TRIQUADRATIC_HEXAHEDRON) and of the same order; and
 *  - all elements in the SWAGE mesh are to be hexahedra of the same order.
 */
void swage2vtk::init_swage_mesh(const int &elem_order, const VtkGrid &grid, 
    SwageMesh &mesh) {
  // Initialize the coordinate data structures in the SWAGE mesh object and
  // store the vertex coordinates from VTK grid in them
  int num_verts = grid->GetPoints()->GetNumberOfPoints();
  mesh.init_nodes(num_verts);
  for (int global_vert_id = 0; global_vert_id < num_verts; global_vert_id++) {
    double x[3]; 
    grid->GetPoints()->GetPoint(global_vert_id, x);
    mesh.node_coords(global_vert_id, 0) = x[0];
    mesh.node_coords(global_vert_id, 1) = x[1];
    mesh.node_coords(global_vert_id, 2) = x[2]; 
  }

  // Loop over the cells in the VTK grid, checking that the cells are all of
  // the correct type. Extract the order of the first cell. If things check
  // out, initialize element data structures in SWAGE mesh object
  int num_cells = grid->GetCells()->GetNumberOfCells();
  for (int cell_id = 0; cell_id < num_cells; cell_id++) {
    int cell_type = grid->GetCellType(cell_id);
    if (cell_type != VTK_LAGRANGE_HEXAHEDRON 
        and cell_type != VTK_TRIQUADRATIC_HEXAHEDRON) {
      std::ostringstream err_oss;
      err_oss << "Error: cell " << cell_id 
          << " is not a VTK_LAGRANGE_HEXAHEDRON" << std::endl;
      throw VtkGridError(err_oss.str());
    }
  }

  int num_dims = 3;
  int num_elems = num_cells;
  mesh.init_element(elem_order, num_dims, num_elems);

  // Map from VTK ordering to SWAGE ordering
  int num_verts_per_dim = elem_order + 1;
  int num_verts_per_elem = std::pow(num_verts_per_dim, num_dims);
  int *vert_map = new int[num_verts_per_elem];
  swage2vtk::make_vtk_vert_map(elem_order, vert_map, true);

  // Loop over the elements in SWAGE and store the cell-vertex connectivity
  // from VTK in the SWAGE mesh object
  for (int elem_id = 0; elem_id < num_elems; elem_id++) {
    // Extract global vert IDs from VTK elem
    vtkNew<vtkIdList> global_vert_ids;
    grid->GetCells()->GetCellAtId(elem_id, global_vert_ids);

    // Loop over VTK vertex IDs 
    for (int vtk_vert_id = 0; vtk_vert_id < num_verts_per_elem; vtk_vert_id++) {
      int global_vert_id = global_vert_ids->GetId(vtk_vert_id);
      int swg_vert_id = vert_map[vtk_vert_id];
      mesh.nodes_in_elem(elem_id, swg_vert_id) = global_vert_id;
    }
  }

  // Deallocate map
  delete[] vert_map;
}

/*
 * Provided a valid VTK unstructured grid object, extract the vertex
 * coordinates and cell-vertex connectivity and store these in a SWAGE mesh
 * object
 *
 * Notes
 * -----
 * This routine currently assumes: 
 *  - all cells in VTK grid are of the same type (VTK_LAGRANGE_HEXAHEDRON) and
 *    of the same order; and
 *  - all elements in the SWAGE mesh are to be hexahedra of the same order.
 */
VtkGrid swage2vtk::init_vtk_grid(SwageMesh &mesh, 
    const std::string &solution_name, const CArray<double> &solution) {
  // Initialize VTK points array
  int num_verts = mesh.num_nodes();
  vtkNew<vtkPoints> points;
  points->SetNumberOfPoints(num_verts);

  // Initialize solution (scalar point data) array
  vtkNew<vtkDoubleArray> solution_array;
  solution_array->SetNumberOfComponents(1);
  solution_array->SetNumberOfTuples(num_verts);
  solution_array->SetName(solution_name.c_str());

  // Initialize elements array
  vtkNew<vtkCellArray> elems;

  // Initialize VTK arbitrary order Lagrange hexahedron object
  int num_dims = mesh.num_dim();
  int elem_order = mesh.elem_order();
  int num_verts_per_dim = elem_order + 1;
  int num_verts_per_elem = std::pow(num_verts_per_dim, num_dims);
  vtkNew<vtkLagrangeHexahedron> hex;
  hex->SetOrder(elem_order, elem_order, elem_order);
  hex->GetPointIds()->SetNumberOfIds(num_verts_per_elem);
  hex->GetPoints()->SetNumberOfPoints(num_verts_per_elem);

  // Map from SWAGE ordering to VTK ordering
  int *vert_map = new int[num_verts_per_elem];
  swage2vtk::make_vtk_vert_map(elem_order, vert_map, false);

  // Loop over the elements in the mesh
  int num_elems = mesh.num_elems();
  for (int elem_id = 0; elem_id < num_elems; elem_id++) {
    // Loop over the vertices in the element
    for (int swage_vert_id = 0; swage_vert_id < num_verts_per_elem; 
        swage_vert_id++) {
      // Obtain the vertex's global ID and its coordinates and add both the ID
      // and the coordinates to the points array
      int global_vert_id = mesh.nodes_in_elem(elem_id, swage_vert_id);
      double x = mesh.node_coords(global_vert_id, 0);
      double y = mesh.node_coords(global_vert_id, 1);
      double z = mesh.node_coords(global_vert_id, 2);
      points->SetPoint(global_vert_id, x, y, z);

      // Get the VTK vertex ID corresponding to the SWAGE vertex ID from the map
      int vtk_vert_id = vert_map[swage_vert_id];
      hex->GetPointIds()->SetId(vtk_vert_id, global_vert_id);

      // Extract the solution at the node from the MATAR array and set its
      // value into the appropriate place in the VTK solution array
      double sol_at_vert = solution(elem_id, swage_vert_id);
      solution_array->SetValue(global_vert_id, sol_at_vert);
    }

    // Add hex to cell array
    elems->InsertNextCell(hex);
  }

  // Add the points, cell array (element connectivity), and solution array to
  // the grid
  VtkGrid grid = VtkGrid::New();
  grid->SetPoints(points);
  grid->SetCells(VTK_LAGRANGE_HEXAHEDRON, elems);
  grid->GetPointData()->SetScalars(solution_array);

  // Deallocate map
  delete[] vert_map;

  return grid;
}

void swage2vtk::write_grid(const VtkGrid &grid, const std::string &file_name) {
  // Give the input data to the writer and write out the file
  // Extract the extension from the file name
  std::string extension = "";
  if (file_name.find_last_of(".") != std::string::npos) {
    extension = file_name.substr(file_name.find_last_of("."));
  }

  // Drop the case of the extension
  std::transform(
      extension.begin(), extension.end(), extension.begin(), ::tolower);

  // Initialize the appropriate writer to the type indicated by the extension
  // and read the grid into a VTK unstructured grid object
  if (extension == ".vtu") {
    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    writer->SetFileName(file_name.c_str());
    writer->SetInputData(grid);
    writer->Write();
  } else if (extension == ".vtk") {
    vtkNew<vtkUnstructuredGridWriter> writer;
    writer->SetFileName(file_name.c_str());
    writer->SetInputData(grid);
    writer->Write();
  }
}
