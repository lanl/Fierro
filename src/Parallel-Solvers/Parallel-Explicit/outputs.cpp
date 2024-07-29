#include "Explicit_Solver.h"
#include "FEA_Module_SGH.h"
#include <map>
#include <fstream>
#include <sys/stat.h>

typedef Kokkos::LayoutRight CArrayLayout;
typedef Kokkos::LayoutLeft FArrayLayout;
typedef Tpetra::MultiVector<>::dual_view_type::t_host::array_layout TpetraHostViewLayout;

MPI_Offset mpi_get_file_position_shared(
    const MPI_Comm comm,
    const MPI_File file_parallel);

void write_string_stream_to_file_mpi_all(
    const std::stringstream& str_stream,
    MPI_Offset               rank_offset_multiplier,
    const MPI_File           file_parallel,
    const MPI_Comm           comm);

template<typename ArrayLayout, typename T>
void write_data_to_string_stream(
    const T*           ptr,
    size_t             dim0,
    size_t             dim1,
    std::stringstream& str_stream);

template<typename ArrayLayout, typename SC, typename LO, typename GO, typename NO>
Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NO>> sort_data(
    const SC*                                 unsorted_data_ptr,
    Teuchos::RCP<Tpetra::Map<LO, GO, NO>>     unsorted_map,
    size_t                                    dim1,
    size_t                                    num_global_unique_elements,
    MPI_Comm                                  comm,
    Teuchos::RCP<Tpetra::Import<LO, GO, NO>>& sorting_importer);

template<typename ArrayLayout, typename SC, typename LO, typename GO, typename NO>
void sort_and_write_data_to_file_mpi_all(
    const SC*                                 unsorted_data_ptr,
    Teuchos::RCP<Tpetra::Map<LO, GO, NO>>     unsorted_map,
    size_t                                    dim1,
    size_t                                    num_global_unique_elements,
    MPI_Comm                                  comm,
    MPI_File                                  file_parallel,
    Teuchos::RCP<Tpetra::Import<LO, GO, NO>>& sorting_importer);

Teuchos::RCP<CArray<int>> get_cell_nodes(
    const CArray<size_t>&            nodes_in_elem,
    size_t                           num_dim,
    Solver::node_ordering_convention active_node_ordering_convention);

Teuchos::RCP<CArray<double>> calculate_elem_speed(
    const Teuchos::RCP<Solver::MV>                             all_node_velocities_distributed,
    const Teuchos::RCP<Solver::MCONN>                          global_nodes_in_elem_distributed,
    const Teuchos::RCP<Solver::MV>                             ghost_node_velocities_distributed,
    const Teuchos::RCP<Tpetra::Import<Solver::LO, Solver::GO>> ghost_importer);

Teuchos::RCP<CArray<int>> calculate_elem_switch(
    Teuchos::RCP<Tpetra::Map<Solver::LO, Solver::GO, Solver::node_type>> all_element_map);

Teuchos::RCP<CArray<int>> get_elem_proc_id(
    Teuchos::RCP<Tpetra::Map<Solver::LO, Solver::GO, Solver::node_type>> all_element_map,
    size_t                                                               myrank);

Teuchos::RCP<CArray<int>> get_elem_gid(
    Teuchos::RCP<Tpetra::Map<Solver::LO, Solver::GO, Solver::node_type>> all_element_map);

Teuchos::RCP<CArray<double>> get_design_density(
    size_t                         rnum_nodes,
    bool                           topology_optimization_on,
    const Teuchos::RCP<Solver::MV> design_node_densities_distributed);

void Explicit_Solver::write_outputs()
{
    // No output for OUTPUT_FORMAT::none
    if (simparam.output_options.output_file_format == OUTPUT_FORMAT::none)
    {
        return;
    }

    const size_t                 rk_level = simparam.dynamic_options.rk_num_bins - 1;
    Teuchos::RCP<CArray<double>> design_density;
    Teuchos::RCP<CArray<int>>    elem_switch;
    Teuchos::RCP<CArray<int>>    elem_proc_id;
    Teuchos::RCP<CArray<int>>    elem_gid;
    Teuchos::RCP<CArray<double>> elem_speed;

    for (const FIELD& field_name : simparam.output_options.output_fields)
    {
        switch (field_name)
        {
        case FIELD::design_density:
            // node "design_density"
            design_density = get_design_density(map->getLocalNumElements(),
          simparam.topology_optimization_on, design_node_densities_distributed);
            point_data_scalars_double["design_density"] = design_density->pointer();
            break;

        case FIELD::speed:
            // element "speed"
            elem_speed = calculate_elem_speed(
                all_node_velocities_distributed,
                global_nodes_in_elem_distributed,
                ghost_node_velocities_distributed,
                ghost_importer);
            cell_data_scalars_double["speed"] = elem_speed->pointer();
            break;

        case FIELD::element_switch:
            // element "element_switch"
            elem_switch = calculate_elem_switch(all_element_map);
            cell_data_scalars_int["element_switch"] = elem_switch->pointer();
            break;

        case FIELD::processor_id:
            // element "processor_id"
            elem_proc_id = get_elem_proc_id(all_element_map, myrank);
            cell_data_scalars_int["processor_id"] = elem_proc_id->pointer();
            break;

        case FIELD::element_id:
            // element "element_id"
            elem_gid = get_elem_gid(all_element_map);
            cell_data_scalars_int["element_id"] = elem_gid->pointer();
            break;

        default:
            break;
        } // end switch
    } // end if

    for (int imodule = 0; imodule < nfea_modules; imodule++)
    {
        fea_modules[imodule]->write_data(
            point_data_scalars_double,
            point_data_vectors_double,
            cell_data_scalars_double,
            cell_data_scalars_int,
            cell_data_fields_double);
    }

    switch (simparam.output_options.output_file_format)
    {
    case OUTPUT_FORMAT::vtk:
        parallel_vtk_writer_new();
        break;

    case OUTPUT_FORMAT::pvtu:
      parallel_vtu_writer_new();
      break;

    default:
      break;
  };
}

void
Explicit_Solver::parallel_vtu_writer_new()
{
  /* to be added... */
  //throw std::runtime_error("parallel_vtu_writer_new() not yet implemented. use parallel_vtk_writer_new()");

  std::string tmp_str;
  std::stringstream mixed_str;

  char name[128];
  char pfilename[128];
  char filename[128];
  char pvtu_dir[128];
  char dirname[128];
  char subdirname[128];
  char tmp[128];

  int num_dims = simparam.num_dims;
  int graphics_idx = simparam.output_options.graphics_id;
  int num_nodes_in_elem = 8;
  int num_points = nnonoverlap_elem_nodes;
  int num_cells = nlocal_elem_non_overlapping;

  all_node_coords_distributed->doImport(*node_coords_distributed, *importer, Tpetra::INSERT);
  host_vec_array node_coords = all_node_coords_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadWrite);

  CArray <size_t> nodes_in_elem (rnum_elem, max_nodes_per_element);
  {
  auto host_view = global_nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
  for (size_t ielem = 0; ielem < nodes_in_elem.dims(0); ielem++) {
    for (size_t inode = 0; inode < nodes_in_elem.dims(1); inode++) {
      nodes_in_elem(ielem, inode) = nonoverlap_element_node_map->getLocalElement(host_view(ielem, inode));
    }
  }
  }

  sprintf(name,"FierroOut");
  sprintf(pvtu_dir, "vtu");
  sprintf(subdirname, "%s_%05lu", name, graphics_idx);

  // mkdir if needed
  struct stat st;
  if (myrank == 0) {
    if (stat(pvtu_dir, &st) != 0) {
      // str_stream.str("");
      // str_stream << "mkdir" << " " << pvtu_dir;
      // system(str_stream.str().c_str());
      //tmp = "";
      sprintf(tmp, "mkdir %s",pvtu_dir);
      system(tmp);
    }
    sprintf(dirname,"%s/%s",pvtu_dir,subdirname);
    if (stat(dirname, &st) != 0) {
      // str_stream.str("");
      // str_stream << "mkdir" << " " << vtu_dir;
      // system(str_stream.str().c_str());
      //tmp = "";
      sprintf(tmp, "mkdir %s", dirname);
      system(tmp);
    }
  }
  MPI_Barrier(world);

  //  ---------------------------------------------------------------------------
  //  Write the PVTU file (only done by a single rank)
  //  ---------------------------------------------------------------------------
  unsigned long int byte_offset = 0;
  int block_header_size = sizeof(unsigned long int);
  unsigned long int data_block_size;

  if (myrank == 0){
    sprintf(pfilename, "%s/%s_%05lu.pvtu",pvtu_dir, name, graphics_idx);
    std::ofstream pout;  // FILE *out;
    pout.open(pfilename,std::ofstream::binary);

    //  Write Header
    tmp_str = "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    tmp_str = "<PUnstructuredGrid GhostLevel=\"1\">\n"; //unsure of exact correct usage of ghostlevel, leaving as 1 for now
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  Write Field Data (only part in pvtu using byte_offset)
    tmp_str = "<FieldData>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    mixed_str.str("");
    mixed_str << "<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(double)+block_header_size;
    tmp_str = mixed_str.str();
    pout.write(tmp_str.c_str(), tmp_str.length());
    tmp_str = "</FieldData>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  Write PPoints Section
    tmp_str = "<PPoints>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    mixed_str.str("");
    mixed_str << "<PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n";
    tmp_str = mixed_str.str();
    pout.write(tmp_str.c_str(), tmp_str.length());
    tmp_str = "</PPoints>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  Write PPoint Data List
    tmp_str = "<PPointData>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    //
    //SCALARS float
    for (auto it = point_data_scalars_double.begin(); it != point_data_scalars_double.end(); it++) {
      mixed_str.str("");
      mixed_str << "<PDataArray type=\"Float64\" Name=\"" << it->first << "\"/>\n";
      tmp_str = mixed_str.str();
      pout.write(tmp_str.c_str(), tmp_str.length());
    }
    //VECTORS float
    for (auto it = point_data_vectors_double.begin(); it != point_data_vectors_double.end(); it++) {
      mixed_str.str("");
      mixed_str << "<PDataArray type=\"Float64\" Name=\"" << it->first << "\" NumberOfComponents=\"" << num_dims << "\"/>\n";
      tmp_str = mixed_str.str();
      pout.write(tmp_str.c_str(), tmp_str.length());
    }
    // //  Velocity
    // mixed_str.str("");
    // mixed_str << "<PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>\n";
    // tmp_str = mixed_str.str();
    // pout.write(tmp_str.c_str(), tmp_str.length());
    //
    tmp_str = "</PPointData>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  Write PCell Data List
    tmp_str = "<PCellData>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    //
    //SCALARS float
    for (auto it = cell_data_scalars_double.begin(); it != cell_data_scalars_double.end(); it++) {
      mixed_str.str("");
      mixed_str << "<PDataArray type=\"Float64\" Name=\"" << it->first << "\"/>\n";
      tmp_str = mixed_str.str();
      pout.write(tmp_str.c_str(), tmp_str.length());
    }
    //SCALARS int
    for (auto it = cell_data_scalars_int.begin(); it != cell_data_scalars_int.end(); it++) {
      mixed_str.str("");
      mixed_str << "<PDataArray type=\"Int32\" Name=\"" << it->first << "\"/>\n";
      tmp_str = mixed_str.str();
      pout.write(tmp_str.c_str(), tmp_str.length());
    }
    // NON-SCALARS float
    for (auto it = cell_data_fields_double.begin(); it != cell_data_fields_double.end(); it++) {
      auto data_name = it->first;
      auto [data_ptr, data_num_comps] = it->second; // Structured binding C++17
      mixed_str.str("");
      mixed_str << "<PDataArray type=\"Float64\" Name=\"" << data_name << "\" NumberOfComponents=\"" << data_num_comps << "\"/>\n";
      tmp_str = mixed_str.str();
      pout.write(tmp_str.c_str(), tmp_str.length());
    }
    //  Rank
    mixed_str.str("");
    mixed_str << "<PDataArray type=\"Int32\" Name=\"rank\"/>\n";
    tmp_str = mixed_str.str();
    pout.write(tmp_str.c_str(), tmp_str.length());
    // //  Stress
    // mixed_str.str("");
    // mixed_str << "<PDataArray type=\"Float64\" Name=\"Stress\" NumberOfComponents=\"9\"/>\n";
    // tmp_str = mixed_str.str();
    // pout.write(tmp_str.c_str(), tmp_str.length());
    // //  Density
    // mixed_str.str("");
    // mixed_str << "<PDataArray type=\"Float64\" Name=\"density\"/>\n";
    // tmp_str = mixed_str.str();
    // pout.write(tmp_str.c_str(), tmp_str.length());
    // //  IPFColor
    // mixed_str.str("");
    // mixed_str << "<PDataArray type=\"unsigned_char\" Name=\"IPFColor\" NumberOfComponents=\"3\"/>\n";
    // tmp_str = mixed_str.str();
    // pout.write(tmp_str.c_str(), tmp_str.length());
    //
    tmp_str = "</PCellData>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  Write Partition Piece List
    char piecename[128];
    for (int rank = 0; rank < nranks; rank++){
      sprintf(piecename, "%s/%s_%05lu_%05lu.vtu",subdirname, name, graphics_idx, rank);
      mixed_str.str("");
      mixed_str << "<Piece Source=\"" << piecename << "\"/>\n";
      tmp_str = mixed_str.str();
      pout.write(tmp_str.c_str(), tmp_str.length());
    }
    //
    tmp_str = "</PUnstructuredGrid>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  Write Appended Data
    tmp_str = "<AppendedData encoding=\"raw\">\n_";
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  Write Time Value
    data_block_size = sizeof(double);
    pout.write((char *) &data_block_size,block_header_size);
    pout.write((char *) &time_value,sizeof(time_value));
    pout.put('\n');
    tmp_str = "</AppendedData>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());  

    tmp_str = "</VTKFile>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());

    pout.close();
  }

  //  ---------------------------------------------------------------------------
  //  Write the VTU file
  //  ---------------------------------------------------------------------------
  sprintf(filename, "%s/%s/%s_%05lu_%05lu.vtu",pvtu_dir ,subdirname, name, graphics_idx, myrank);
  // filename has the full string
  
  std::ofstream out;  // FILE *out;
  out.open(filename,std::ofstream::binary);

  byte_offset = 0;
  double crap = 0.0;
  
  //  Write Header
  //tmp_str = "<?xml version=\"1.0\"?>\n";
  //out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "<UnstructuredGrid>\n";
  out.write(tmp_str.c_str(), tmp_str.length());

  //  Write Time Value Header
  tmp_str = "<FieldData>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
  byte_offset += sizeof(double)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</FieldData>\n";
  out.write(tmp_str.c_str(), tmp_str.length());

  mixed_str.str("");
  mixed_str << "<Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\"" << num_cells << "\">\n";
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());

  //  **Write Points Header**
  tmp_str = "<Points>\n";
  out.write(tmp_str.c_str(), tmp_str.length());   
  //tmp_str = "<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\">\n"; 
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  byte_offset += num_points*num_dims*sizeof(double)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</DataArray>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  // tmp_str = "</DataArray>\n";
  // out.write(tmp_str.c_str(), tmp_str.length());  
  tmp_str = "</Points>\n";
  out.write(tmp_str.c_str(), tmp_str.length());

  //  **Write Cells Header**
  tmp_str = "<Cells>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  //  Connectivity
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
  byte_offset += num_cells*num_nodes_in_elem*sizeof(int)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  //  Offsets
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
  byte_offset += num_cells*sizeof(int)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  //  Types
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n"; 
  byte_offset += num_cells*sizeof(int)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</Cells>\n";
  out.write(tmp_str.c_str(), tmp_str.length());

  //  **Write Point Data Headers**
  tmp_str = "<PointData>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  //
  //SCALARS float
  for (auto it = point_data_scalars_double.begin(); it != point_data_scalars_double.end(); it++) {
    mixed_str.str("");
    mixed_str << "<DataArray type=\"Float64\" Name=\"" << it->first << "\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
    byte_offset += num_points*sizeof(double)+block_header_size;
    tmp_str = mixed_str.str();
    out.write(tmp_str.c_str(), tmp_str.length());
    tmp_str = "</DataArray>\n";
    out.write(tmp_str.c_str(), tmp_str.length());
  }
  //VECTORS float
  for (auto it = point_data_vectors_double.begin(); it != point_data_vectors_double.end(); it++) {
    mixed_str.str("");
    mixed_str << "<DataArray type=\"Float64\" Name=\"" << it->first << "\" NumberOfComponents=\"" << num_dims << "\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
    byte_offset += num_points*num_dims*sizeof(double)+block_header_size;
    tmp_str = mixed_str.str();
    out.write(tmp_str.c_str(), tmp_str.length());
    tmp_str = "</DataArray>\n";
    out.write(tmp_str.c_str(), tmp_str.length());
  }
  // //  Velocity
  // mixed_str.str("");
  // mixed_str << "<DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  // byte_offset += num_points*num_dims*sizeof(double)+block_header_size;
  // tmp_str = mixed_str.str();
  // out.write(tmp_str.c_str(), tmp_str.length());
  // tmp_str = "</DataArray>\n";
  // out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</PointData>\n";
  out.write(tmp_str.c_str(), tmp_str.length());

  //  **Write Cell Data Headers**
  tmp_str = "<CellData>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  //
  //SCALARS float
  for (auto it = cell_data_scalars_double.begin(); it != cell_data_scalars_double.end(); it++) {
    mixed_str.str("");
    mixed_str << "<DataArray type=\"Float64\" Name=\"" << it->first << "\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
    byte_offset += num_cells*sizeof(double)+block_header_size;
    tmp_str = mixed_str.str();
    out.write(tmp_str.c_str(), tmp_str.length());
    tmp_str = "</DataArray>\n";
    out.write(tmp_str.c_str(), tmp_str.length());
  }
  //SCALARS int
  for (auto it = cell_data_scalars_int.begin(); it != cell_data_scalars_int.end(); it++) {
    mixed_str.str("");
    mixed_str << "<DataArray type=\"Int32\" Name=\"" << it->first << "\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
    byte_offset += num_cells*sizeof(int)+block_header_size;
    tmp_str = mixed_str.str();
    out.write(tmp_str.c_str(), tmp_str.length());
    tmp_str = "</DataArray>\n";
    out.write(tmp_str.c_str(), tmp_str.length());
  }
  // NON-SCALARS float
  for (auto it = cell_data_fields_double.begin(); it != cell_data_fields_double.end(); it++) {
    auto data_name = it->first;
    auto [data_ptr, data_num_comps] = it->second; // Structured binding C++17
    mixed_str.str("");
    mixed_str << "<DataArray type=\"Float64\" Name=\"" << data_name << "\" NumberOfComponents=\"" << data_num_comps << "\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
    byte_offset += num_cells*data_num_comps*sizeof(double)+block_header_size;
    tmp_str = mixed_str.str();
    out.write(tmp_str.c_str(), tmp_str.length());
    tmp_str = "</DataArray>\n";
    out.write(tmp_str.c_str(), tmp_str.length());
  }
  //  Rank
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Int32\" Name=\"rank\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  byte_offset += num_cells*sizeof(int)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</DataArray>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  // //  Stress
  // mixed_str.str("");
  // mixed_str << "<DataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"9\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  // byte_offset += num_cells*sizeof(double)+block_header_size;
  // tmp_str = mixed_str.str();
  // out.write(tmp_str.c_str(), tmp_str.length());
  // tmp_str = "</DataArray>\n";
  // out.write(tmp_str.c_str(), tmp_str.length());
  // //  Density
  // mixed_str.str("");
  // mixed_str << "<DataArray type=\"Float64\" Name=\"density\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  // byte_offset += num_cells*sizeof(double)+block_header_size;
  // tmp_str = mixed_str.str();
  // out.write(tmp_str.c_str(), tmp_str.length());
  // tmp_str = "</DataArray>\n";
  // out.write(tmp_str.c_str(), tmp_str.length());
  // //  IPFColor
  // mixed_str.str("");
  // mixed_str << "<DataArray type=\"UInt8\" Name=\"IPFColor\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  // byte_offset += num_cells*num_dims*sizeof(unsigned char)+block_header_size;
  // tmp_str = mixed_str.str();
  // out.write(tmp_str.c_str(), tmp_str.length());
  // tmp_str = "</DataArray>\n";
  // out.write(tmp_str.c_str(), tmp_str.length());
  //
  tmp_str = "</CellData>\n";
  out.write(tmp_str.c_str(), tmp_str.length());

  //  Write Mesh Close
  tmp_str = "</Piece>\n";
  out.write(tmp_str.c_str(), tmp_str.length());  
  tmp_str = "</UnstructuredGrid>\n";
  out.write(tmp_str.c_str(), tmp_str.length());  

  //  Write Appended Data
  tmp_str = "<AppendedData encoding=\"raw\">\n_";
  out.write(tmp_str.c_str(), tmp_str.length());  

  //  Write Time Value
  data_block_size = sizeof(double);
  out.write((char *) &data_block_size,block_header_size);
  out.write((char *) &time_value,sizeof(time_value));

  //  **Write Points Data**
  //double coords[num_points*num_dims];
  data_block_size = num_points*num_dims*sizeof(double);
  out.write((char *) &data_block_size,block_header_size);
  double coord_tmp;
  // for (int node_gid = 0; node_gid < num_points; node_gid++){
  //   for (int dim = 0; dim < num_dims; dim++){
  //     coord_tmp = node_coords(node_gid,dim);
  //     out.write((char *) &coord_tmp,sizeof(coord_tmp));
  //   }
  // }
  for (int node_id = 0; node_id < num_points; node_id++){
    int node_gid = nonoverlap_element_node_map->getGlobalElement(node_id);
    int node_lid = all_node_map->getLocalElement(node_gid);
    for (int dim = 0; dim < num_dims; dim++){
      //printf("rank %d point %d coord %d = %g \n",myrank,node_lid,dim,node_coords(node_lid,dim));
      coord_tmp = node_coords(node_lid,dim);
      out.write((char *) &coord_tmp,sizeof(coord_tmp));
    }
  }

  //  **Write Cells Data**
  //int connect[num_cells*mesh.num_nodes_in_elem];
  data_block_size = num_cells*num_nodes_in_elem*sizeof(int);
  out.write((char *) &data_block_size,block_header_size);
  int connect_tmp;
  for (int elem_gid = 0; elem_gid < num_cells; elem_gid++){
    for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
      connect_tmp = nodes_in_elem(elem_gid, node_lid);
      out.write((char *) &connect_tmp,sizeof(connect_tmp));
    }
  }
  //int offsets[num_cells];
  data_block_size = num_cells*sizeof(int);
  out.write((char *) &data_block_size,block_header_size);
  int offset_tmp;
  for (int elem_gid = 0; elem_gid < num_cells; elem_gid++){
    offset_tmp = (elem_gid+1)*num_nodes_in_elem;
    out.write((char *) &offset_tmp,sizeof(offset_tmp));
  }
  //int types[num_cells];
  data_block_size = num_cells*sizeof(int);
  out.write((char *) &data_block_size,block_header_size);
  int type_tmp;
  for (int elem_gid = 0; elem_gid < num_cells; elem_gid++){
    type_tmp = 12;
    out.write((char *) &type_tmp,sizeof(type_tmp));
  }

  //  **Write Point Data**
  //SCALARS float
  for (auto it = point_data_scalars_double.begin(); it != point_data_scalars_double.end(); it++) {
    auto data_ptr = it->second;
    ViewCArrayKokkos <const double> data_tmp(data_ptr,nall_nodes);

    data_block_size = num_points*sizeof(double);
    out.write((char *) &data_block_size,block_header_size);
    double pscalar_tmp;
    for (int node_id = 0; node_id < num_points; node_id++){
      int node_gid = nonoverlap_element_node_map->getGlobalElement(node_id);
      int node_lid = all_node_map->getLocalElement(node_gid);
      //printf("rank %d point %d coord %d = %g \n",myrank,node_lid,dim,node_coords(node_lid,dim));
      pscalar_tmp = data_tmp(node_lid);
      out.write((char *) &pscalar_tmp,sizeof(pscalar_tmp));
    }
  }
  //VECTORS float
  for (auto it = point_data_vectors_double.begin(); it != point_data_vectors_double.end(); it++) {
    auto data_ptr = it->second;
    ViewCArrayKokkos <const double> data_tmp(data_ptr,nall_nodes,num_dims);

    data_block_size = num_points*num_dims*sizeof(double);
    out.write((char *) &data_block_size,block_header_size);
    double pvector_tmp;
    for (int node_id = 0; node_id < num_points; node_id++){
      int node_gid = nonoverlap_element_node_map->getGlobalElement(node_id);
      int node_lid = all_node_map->getLocalElement(node_gid);
      for (int dim = 0; dim < num_dims; dim++){
        //printf("rank %d point %d coord %d = %g \n",myrank,node_lid,dim,node_coords(node_lid,dim));
        pvector_tmp = data_tmp(node_lid,dim);
        out.write((char *) &pvector_tmp,sizeof(pvector_tmp));
      }
    }
  }
  // //  Velocity
  // data_block_size = num_points*num_dims*sizeof(double);
  // out.write((char *) &data_block_size,block_header_size);
  // for (int node_gid = 0; node_gid < num_points; node_gid++){
  //   for (int dim = 0; dim < num_dims; dim++){
  //     //out.write((char *) &node_vel.host(1, node_gid, dim),sizeof(double));
  //     out.write((char *) &crap,sizeof(double));
  //   }
  // }

  //  **Write Cell Data**
  //SCALARS float
  for (auto it = cell_data_scalars_double.begin(); it != cell_data_scalars_double.end(); it++) {
    auto data_ptr = it->second;
    ViewCArrayKokkos <const double> data_tmp(data_ptr,rnum_elem);

    data_block_size = num_cells*sizeof(double);
    out.write((char *) &data_block_size,block_header_size);
    double cscalar_ftmp;
    for (int elem_gid = 0; elem_gid < num_cells; elem_gid++){
      cscalar_ftmp = data_tmp(elem_gid);
      out.write((char *) &cscalar_ftmp,sizeof(double));
    }
  }
  //SCALARS int
  for (auto it = cell_data_scalars_int.begin(); it != cell_data_scalars_int.end(); it++) {
    auto data_ptr = it->second;
    ViewCArrayKokkos <const int> data_tmp(data_ptr,rnum_elem);

    data_block_size = num_cells*sizeof(int);
    out.write((char *) &data_block_size,block_header_size);
    double cscalar_itmp;
    for (int elem_gid = 0; elem_gid < num_cells; elem_gid++){
      cscalar_itmp = data_tmp(elem_gid);
      out.write((char *) &cscalar_itmp,sizeof(int));
    }
  }
  // NON-SCALARS float
  for (auto it = cell_data_fields_double.begin(); it != cell_data_fields_double.end(); it++) {
    auto data_name = it->first;
    auto [data_ptr, data_num_comps] = it->second; // Structured binding C++17
    ViewCArrayKokkos <const double> data_tmp(data_ptr,rnum_elem,data_num_comps);

    data_block_size = num_cells*data_num_comps*sizeof(double);
    out.write((char *) &data_block_size,block_header_size);
    double cfield_ftmp;
    for (int elem_gid = 0; elem_gid < num_cells; elem_gid++){
      for (int comp = 0; comp < data_num_comps; comp++){
        cfield_ftmp = data_tmp(elem_gid,comp);
        out.write((char *) &cfield_ftmp,sizeof(double));
      }
    }
  }
  //  Rank
  data_block_size = num_cells*sizeof(int);
  out.write((char *) &data_block_size,block_header_size);
  for (int elem_gid = 0; elem_gid < num_cells; elem_gid++){
      out.write((char *) &myrank,sizeof(int));
  }
  // //  Stress
  // data_block_size = num_cells*sizeof(double);
  // out.write((char *) &data_block_size,block_header_size);
  // for (int elem_gid = 0; elem_gid < num_cells; elem_gid++){
  //   for (int ii = 0; ii < 9; ii++){
  //     //out.write((char *) &elem_stress.host(1,elem_gid,ii),sizeof(double));
  //     out.write((char *) &crap,sizeof(double));
  //   }
  // }
  // //  Density
  // data_block_size = num_cells*sizeof(double);
  // out.write((char *) &data_block_size,block_header_size);
  // for (int elem_gid = 0; elem_gid < num_cells; elem_gid++){
  //     out.write((char *) &elem_den.host(elem_gid),sizeof(double));
  // }
  // //  IPFColor
  // data_block_size = num_cells*num_dims*sizeof(unsigned char);
  // out.write((char *) &data_block_size,block_header_size);
  // for (int elem_gid = 0; elem_gid < num_cells; elem_gid++){
  //     for (int dim = 0; dim < num_dims; dim++){
  //         out.write((char *) &elem_IPF_colors.host(elem_gid,dim),sizeof(unsigned char));
  //     }
  // }
  
  out.put('\n');
  tmp_str = "</AppendedData>\n";
  out.write(tmp_str.c_str(), tmp_str.length());  

  //  Write File Close
  tmp_str = "</VTKFile>\n";
  out.write(tmp_str.c_str(), tmp_str.length());  

  out.close();

  simparam.output_options.graphics_id++;

}

std::string construct_file_name(
    size_t file_count,
    int    displacement_module)
{
    bool        displace_geometry = false;
    std::string current_file_name;
    std::string base_file_name = "VTK";
    std::string base_file_name_undeformed = "VTK_undeformed";
    std::string file_extension = ".vtk";

    if (displace_geometry && displacement_module >= 0)
    {
        current_file_name = base_file_name_undeformed + std::to_string(file_count) + file_extension;
    }
    else
    {
        current_file_name = base_file_name + std::to_string(file_count) + file_extension;
    }
    return current_file_name;
}

void Explicit_Solver::parallel_vtk_writer_new()
{
    int               num_dim = simparam.num_dims;
    std::stringstream str_stream;
    MPI_Offset        current_offset;
    MPI_File          myfile_parallel;
    bool              displace_geometry = false;

    std::string vtk_dir      = simparam.output_options.output_file_location;
    std::string vtk_data_dir = vtk_dir + "data/";

    // mkdir if needed
    struct stat st;
    if (myrank == 0)
    {
        if (stat(vtk_dir.c_str(), &st) != 0)
        {
            str_stream.str("");
            str_stream << "mkdir" << " " << vtk_dir;
            system(str_stream.str().c_str());
        }
        if (stat(vtk_data_dir.c_str(), &st) != 0)
        {
            str_stream.str("");
            str_stream << "mkdir" << " " << vtk_data_dir;
            system(str_stream.str().c_str());
        }
    }
    MPI_Barrier(world);

    // construct file name
    std::string current_file_name = construct_file_name(simparam.output_options.graphics_id, displacement_module);
    std::string file_path = vtk_data_dir + current_file_name;
    // open mpi file
    int err = MPI_File_open(MPI_COMM_WORLD, file_path.c_str(),
                          MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &myfile_parallel);
    // allows overwriting the file if it already existed in the directory
    if (err != MPI_SUCCESS)
    {
        if (myrank == 0)
        {
            MPI_File_delete(file_path.c_str(), MPI_INFO_NULL);
        }
        MPI_File_open(MPI_COMM_WORLD, file_path.c_str(),
                  MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &myfile_parallel);
    }

    /*************** write header of the vtk file ***************/
    str_stream.str("");
    str_stream << "# vtk DataFile Version 2.0" << std::endl;
    str_stream << "Mesh for Fierro" << std::endl;
    str_stream << "ASCII" << std::endl;
    str_stream << "DATASET UNSTRUCTURED_GRID" << std::endl;
    current_offset = mpi_get_file_position_shared(world, myfile_parallel);
    if (myrank == 0)
    {
        MPI_File_write_at(
            myfile_parallel,
            current_offset,
            str_stream.str().c_str(),
            str_stream.str().length(),
            MPI_CHAR,
            MPI_STATUS_IGNORE);
    }

    /*************** write POINTS ***************/
    str_stream.str("");
    str_stream << std::endl << "POINTS " << num_nodes << " float" << std::endl;
    current_offset = mpi_get_file_position_shared(world, myfile_parallel);
    if (myrank == 0)
    {
        MPI_File_write_at(
            myfile_parallel,
            current_offset,
            str_stream.str().c_str(),
            str_stream.str().length(),
            MPI_CHAR,
            MPI_STATUS_IGNORE);
    }
    // sgh_module->node_coords.update_host();
    { // view scope
        host_vec_array node_coords = node_coords_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadWrite);
        double*        coord_data  = node_coords.data();
        sort_and_write_data_to_file_mpi_all<array_layout, double, LO, GO, node_type>(
            coord_data,
            map,
            num_dim,
            num_nodes,
            world,
            myfile_parallel,
            node_sorting_importer);
    }

    /*************** write CELLS ***************/
    str_stream.str("");
    str_stream << std::endl << "CELLS " << num_elem << " " << num_elem * (max_nodes_per_element + 1) << std::endl;
    current_offset = mpi_get_file_position_shared(world, myfile_parallel);
    if (myrank == 0)
    {
        MPI_File_write_at(
            myfile_parallel,
            current_offset,
            str_stream.str().c_str(),
            str_stream.str().length(),
            MPI_CHAR,
            MPI_STATUS_IGNORE);
    }
    CArray<size_t> nodes_in_elem(rnum_elem, max_nodes_per_element);
    { // view scope
        auto host_view = global_nodes_in_elem_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
        for (size_t ielem = 0; ielem < nodes_in_elem.dims(0); ielem++)
        {
            for (size_t inode = 0; inode < nodes_in_elem.dims(1); inode++)
            {
                nodes_in_elem(ielem, inode) = host_view(ielem, inode);
            }
        }
    } // end view scope
    auto cell_data = get_cell_nodes(nodes_in_elem, num_dim, active_node_ordering_convention);
    sort_and_write_data_to_file_mpi_all<CArrayLayout, int, LO, GO, node_type>(
        cell_data->pointer(),
        all_element_map,
        cell_data->dims(1),
        num_elem,
        world,
        myfile_parallel,
        element_sorting_importer);

    /*************** write CELL_TYPES ***************/
    str_stream.str("");
    str_stream << std::endl << "CELL_TYPES " << num_elem << std::endl;
    current_offset = mpi_get_file_position_shared(world, myfile_parallel);
    if (myrank == 0)
    {
        MPI_File_write_at(
            myfile_parallel,
            current_offset,
            str_stream.str().c_str(),
            str_stream.str().length(),
            MPI_CHAR,
            MPI_STATUS_IGNORE);
    }
    CArray<int> cell_type(all_element_map->getLocalNumElements(), 1);
    for (int i = 0; i < cell_type.dims(0); i++)
    {
        for (int j = 0; j < cell_type.dims(1); j++)
        {
            cell_type(i, j) = 12;
        }
    }
    sort_and_write_data_to_file_mpi_all<CArrayLayout, int, LO, GO, node_type>(
        cell_type.pointer(),
        all_element_map,
        cell_type.dims(1),
        num_elem,
        world,
        myfile_parallel,
        element_sorting_importer);

    /*************** write POINT_DATA ***************/
    str_stream.str("");
    str_stream << std::endl << "POINT_DATA " << num_nodes;
    current_offset = mpi_get_file_position_shared(world, myfile_parallel);
    if (myrank == 0)
    {
        MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(),
                      str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
    }

    // SCALARS float
    for (auto it = point_data_scalars_double.begin(); it != point_data_scalars_double.end(); it++)
    {
        str_stream.str("");
        str_stream << std::endl << "SCALARS " << it->first << " float 1" << std::endl;
        str_stream << "LOOKUP_TABLE default" << std::endl;
        current_offset = mpi_get_file_position_shared(world, myfile_parallel);
        if (myrank == 0)
        {
            MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(),
                        str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
        }
        sort_and_write_data_to_file_mpi_all<CArrayLayout, double, LO, GO, node_type>(
      it->second, map, 1, num_nodes, world, myfile_parallel, node_sorting_importer);
    }

    // VECTORS float
    for (auto it = point_data_vectors_double.begin(); it != point_data_vectors_double.end(); it++)
    {
        str_stream.str("");
        str_stream << std::endl << "VECTORS " << it->first << " float" << std::endl;
        current_offset = mpi_get_file_position_shared(world, myfile_parallel);
        if (myrank == 0)
        {
            MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(),
                        str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
        }
        sort_and_write_data_to_file_mpi_all<CArrayLayout, double, LO, GO, node_type>(
      it->second, map, num_dim, num_nodes, world, myfile_parallel, node_sorting_importer);
    }

    /*************** write CELL_DATA ***************/
    str_stream.str("");
    str_stream << std::endl << "CELL_DATA " << num_elem;
    current_offset = mpi_get_file_position_shared(world, myfile_parallel);
    if (myrank == 0)
    {
        MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(),
                      str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
    }

    // SCALARS float
    for (auto it = cell_data_scalars_double.begin(); it != cell_data_scalars_double.end(); it++)
    {
        str_stream.str("");
        str_stream << std::endl << "SCALARS " << it->first << " float 1" << std::endl;
        str_stream << "LOOKUP_TABLE default" << std::endl;
        current_offset = mpi_get_file_position_shared(world, myfile_parallel);
        if (myrank == 0)
        {
            MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(),
                        str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
        }
        sort_and_write_data_to_file_mpi_all<CArrayLayout, double, LO, GO, node_type>(
      it->second, all_element_map, 1, num_elem, world, myfile_parallel, element_sorting_importer);
    }

    // SCALARS int
    for (auto it = cell_data_scalars_int.begin(); it != cell_data_scalars_int.end(); it++)
    {
        str_stream.str("");
        str_stream << std::endl << "SCALARS " << it->first << " int 1" << std::endl;
        str_stream << "LOOKUP_TABLE default" << std::endl;
        current_offset = mpi_get_file_position_shared(world, myfile_parallel);
        if (myrank == 0)
        {
            MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(),
                        str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
        }
        sort_and_write_data_to_file_mpi_all<CArrayLayout, int, LO, GO, node_type>(
      it->second, all_element_map, 1, num_elem, world, myfile_parallel, element_sorting_importer);
    }

    // FIELD
    str_stream.str("");
    str_stream << std::endl << "FIELD FieldData " << cell_data_fields_double.size() << std::endl;
    current_offset = mpi_get_file_position_shared(world, myfile_parallel);
    if (myrank == 0)
    {
        MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(),
                      str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    for (auto it = cell_data_fields_double.begin(); it != cell_data_fields_double.end(); it++)
    {
        auto data_name = it->first;
        auto [data_ptr, data_num_comps] = it->second; // Structured binding C++17
        str_stream.str("");
        str_stream << data_name << " " << data_num_comps << " " << num_elem << " float" << std::endl;
        current_offset = mpi_get_file_position_shared(world, myfile_parallel);
        if (myrank == 0)
        {
            MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(),
                        str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
        }
        sort_and_write_data_to_file_mpi_all<CArrayLayout, double, LO, GO, node_type>(
      data_ptr, all_element_map, data_num_comps, num_elem, world, myfile_parallel, element_sorting_importer);
    }

    if(simparam.output_options.optimization_restart_file){
        // Write commented restart data to be used by Fierro
        str_stream.str("");
        str_stream << std::endl << "#RESTART DATA: Objective_Normalization_Constant " <<
                    simparam.optimization_options.objective_normalization_constant << std::endl;
        current_offset = mpi_get_file_position_shared(world, myfile_parallel);
        if (myrank == 0)
        {
            MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(),
                        str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
        }
    }

    MPI_Barrier(world);
    MPI_File_sync(myfile_parallel);
    MPI_File_close(&myfile_parallel);

    /*************** write .vtk.series file ***************/
    simparam.output_options.graphics_times(simparam.output_options.graphics_id) = time_value;
    if (myrank == 0)
    {
        FILE*       myfile;
        std::string filename = vtk_dir + "outputs.vtk.series";
        myfile = fopen(filename.c_str(), "w");

        fprintf(myfile, "{\n");
        fprintf(myfile, "  \"file-series-version\" : \"1.0\",\n");
        fprintf(myfile, "  \"files\" : [\n");
        for (int i = 0; i <= simparam.output_options.graphics_id; i++)
        {
            std::string vtk_filename = construct_file_name(i, displacement_module);
            fprintf(myfile, "    { \"name\" : \"data/%s\", \"time\" : %12.5e },\n",
        vtk_filename.c_str(), simparam.output_options.graphics_times(i));
        }
        fprintf(myfile, "  ]\n");
        fprintf(myfile, "}\n");

        fclose(myfile);
    }

    simparam.output_options.graphics_id++;
}

MPI_Offset mpi_get_file_position_shared(
    const MPI_Comm comm,
    const MPI_File file_parallel)
{
    MPI_Offset position;

    MPI_Barrier(comm);
    MPI_File_seek_shared(file_parallel, 0, MPI_SEEK_END);
    MPI_File_get_position_shared(file_parallel, &position);
    MPI_Barrier(comm);
    return position;
}

void write_string_stream_to_file_mpi_all(
    const std::stringstream& str_stream,
    MPI_Offset               rank_offset_multiplier,
    const MPI_File           file_parallel,
    const MPI_Comm           comm)
{
    /* `str_stream`: must have the same line length
    *  `rank_offset_multiplier`: should be the index of the first node or element in the rank
    *  `file_parallel`: the parallel mpi file to write to
    *  `comm`: MPI_Comm that should do the write
    * */

    MPI_Offset current_offset, file_offset;
    current_offset = mpi_get_file_position_shared(comm, file_parallel);
    std::string line;
    std::getline(const_cast<std::stringstream&>(str_stream), line);
    auto line_length = line.length() + 1; // the +1 is for end of line character
    file_offset = line_length * rank_offset_multiplier + current_offset;
    MPI_Barrier(comm);
    MPI_File_write_at_all(file_parallel, file_offset, str_stream.str().c_str(),
                        str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_Barrier(comm);
    return;
}

template<typename ArrayLayout, typename T>
void write_data_to_string_stream(
    const T*           ptr,
    size_t             dim0,
    size_t             dim1,
    std::stringstream& str_stream)
{
    size_t w;
    size_t p;
    // crude way to ensure that all lines have the same length. required for MPI IO
    if constexpr (std::is_same<T, float>::value or
                  std::is_same<T, double>::value)
    {
        w = 25;
        p = 8;
    }
    else if constexpr (std::is_same<T, int>::value or
                       std::is_same<T, size_t>::value or
                       std::is_same<T, Solver::GO>::value or
                       std::is_same<T, Solver::LO>::value)
    {
        w = 10;
        p = 0;
    }
    else
    {
        throw std::runtime_error(std::string("unknown datatype in function ") +
                             std::string(__FUNCTION__) + std::string("in file ") +
                             std::string(__FILE__));
    }

    // create view of unsorted_data_ptr. not using matar here. array_layout is needed
    Kokkos::View<const T**, ArrayLayout, HostSpace, Kokkos::MemoryUnmanaged>
    data_view(ptr, dim0, dim1);

    str_stream.str("");
    str_stream << std::fixed << std::setprecision(p);
    for (size_t i = 0; i < dim0; i++)
    {
        for (size_t j = 0; j < dim1; j++)
        {
            str_stream << std::left << std::setw(w) << data_view(i, j) << " ";
        }
        str_stream << std::endl;
    }

    return;
}

template<typename ArrayLayout, typename SC, typename LO, typename GO, typename NO>
Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NO>> sort_data(
    const SC*                                 unsorted_data_ptr,
    Teuchos::RCP<Tpetra::Map<LO, GO, NO>>     unsorted_map,
    size_t                                    dim1,
    size_t                                    num_global_unique_elements,
    MPI_Comm                                  comm,
    Teuchos::RCP<Tpetra::Import<LO, GO, NO>>& sorting_importer)
{
    // create view of unsorted_data_ptr. not using matar here. array_layout is needed
    Kokkos::View<const SC**, ArrayLayout, HostSpace, Kokkos::MemoryUnmanaged>
    unsorted_data(unsorted_data_ptr, unsorted_map->getLocalNumElements(), dim1);

    // sorted_map
    Teuchos::RCP<Tpetra::Map<LO, GO, NO>> sorted_map =
        Teuchos::rcp(new Tpetra::Map<LO, GO, NO>(
    num_global_unique_elements, unsorted_map->getIndexBase(), unsorted_map->getComm()));

    // sorted storage
    Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NO>> sorted_storage =
        Teuchos::rcp(new Tpetra::MultiVector<SC, LO, GO, NO>(sorted_map, dim1));

    // unsorted storage
    Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NO>> unsorted_storage =
        Teuchos::rcp(new Tpetra::MultiVector<SC, LO, GO, NO>(unsorted_map, dim1));

    { // view scope
        auto host_view = unsorted_storage->getLocalViewHost(Tpetra::Access::ReadWrite);
        for (int i = 0; i < unsorted_map->getLocalNumElements(); i++)
        {
            for (int j = 0; j < dim1; j++)
            {
                host_view(i, j) = unsorted_data(i, j);
            }
        } // for
    } // end view scope

    sorted_storage->doImport(*unsorted_storage, *sorting_importer, Tpetra::INSERT);

    return sorted_storage;
}

template<typename ArrayLayout, typename SC, typename LO, typename GO, typename NO>
void sort_and_write_data_to_file_mpi_all(
    const SC*                                 unsorted_data_ptr,
    Teuchos::RCP<Tpetra::Map<LO, GO, NO>>     unsorted_map,
    size_t                                    dim1,
    size_t                                    num_global_unique_elements,
    MPI_Comm                                  comm,
    MPI_File                                  file_parallel,
    Teuchos::RCP<Tpetra::Import<LO, GO, NO>>& sorting_importer)
{
    auto sorted_data =
        sort_data<ArrayLayout, SC, LO, GO, NO>(
            unsorted_data_ptr,
            unsorted_map,
            dim1,
            num_global_unique_elements,
            comm,
            sorting_importer);

    { // view scope
        auto              const_host_view = sorted_data->getLocalViewHost(Tpetra::Access::ReadOnly);
        std::stringstream str_stream;
        write_data_to_string_stream<TpetraHostViewLayout, SC>(
            const_host_view.data(),
            sorted_data->getMap()->getLocalNumElements(),
            dim1,
            str_stream);
        write_string_stream_to_file_mpi_all(
            str_stream,
            sorted_data->getMap()->getGlobalElement(0),
            file_parallel,
            comm);
    } // end view scope
}

Teuchos::RCP<CArray<int>> get_cell_nodes(
    const CArray<size_t>&            nodes_in_elem,
    size_t                           num_dim,
    Solver::node_ordering_convention active_node_ordering_convention)
{
    CArray<size_t> convert_ijk_to_ensight(nodes_in_elem.dims(1));
    convert_ijk_to_ensight(0) = 0;
    convert_ijk_to_ensight(1) = 1;
    convert_ijk_to_ensight(2) = 3;
    convert_ijk_to_ensight(3) = 2;
    if (num_dim == 3)
    {
        convert_ijk_to_ensight(4) = 4;
        convert_ijk_to_ensight(5) = 5;
        convert_ijk_to_ensight(6) = 7;
        convert_ijk_to_ensight(7) = 6;
    }

    Teuchos::RCP<CArray<int>> cell_data =
        Teuchos::rcp(new CArray<int>(nodes_in_elem.dims(0), nodes_in_elem.dims(1) + 1) );

    size_t temp_convert;
    for (size_t ielem = 0; ielem < nodes_in_elem.dims(0); ielem++)
    {
        (*cell_data)(ielem, 0) = 8;
        for (int ii = 0; ii < nodes_in_elem.dims(1); ii++)
        {
            if (active_node_ordering_convention == Solver::node_ordering_convention::IJK)
            {
                temp_convert = convert_ijk_to_ensight(ii);
            }
            else
            {
                temp_convert = ii;
            }
            (*cell_data)(ielem, ii + 1) = nodes_in_elem(ielem, temp_convert);
            // (*cell_data)(ielem,ii+1) = nodes_in_elem(ielem,ii);
        }
    }

    return cell_data;
}

Teuchos::RCP<CArray<double>> calculate_elem_speed(
    const Teuchos::RCP<Solver::MV>                             all_node_velocities_distributed,
    const Teuchos::RCP<Solver::MCONN>                          global_nodes_in_elem_distributed,
    const Teuchos::RCP<Solver::MV>                             ghost_node_velocities_distributed,
    const Teuchos::RCP<Tpetra::Import<Solver::LO, Solver::GO>> ghost_importer)
{
    size_t rnum_elems = global_nodes_in_elem_distributed->getLocalLength();
    size_t num_nodes_in_elem = global_nodes_in_elem_distributed->getNumVectors();
    size_t num_dims = all_node_velocities_distributed->getNumVectors();

    Teuchos::RCP<CArray<double>> elem_speed =
        Teuchos::rcp(new CArray<double>(rnum_elems));

    auto nodes_in_elem_hview = global_nodes_in_elem_distributed->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto all_node_vel_hview  = all_node_velocities_distributed->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto vector_map = all_node_velocities_distributed->getMap();
    for (size_t elem_gid = 0; elem_gid < rnum_elems; elem_gid++)
    {
        double elem_vel[3];
        elem_vel[0] = 0.0;
        elem_vel[1] = 0.0;
        elem_vel[2] = 0.0;
        // get the coordinates of the element center
        for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++)
        {
            size_t inode = vector_map->getLocalElement(nodes_in_elem_hview(elem_gid, node_lid));
            elem_vel[0] += all_node_vel_hview(inode, 0);
            elem_vel[1] += all_node_vel_hview(inode, 1);
            if (num_dims == 3)
            {
                elem_vel[2] += all_node_vel_hview(inode, 2);
            }
            else
            {
                elem_vel[2] = 0.0;
            }
        } // end loop over nodes in element

        elem_vel[0] = elem_vel[0] / num_nodes_in_elem;
        elem_vel[1] = elem_vel[1] / num_nodes_in_elem;
        elem_vel[2] = elem_vel[2] / num_nodes_in_elem;

        double speed_sqrd = 0.0;
        for (int dim = 0; dim < num_dims; dim++)
        {
            speed_sqrd += elem_vel[dim] * elem_vel[dim];
        }

        (*elem_speed)(elem_gid) = sqrt(speed_sqrd);
    } // end for

    return elem_speed;
}

Teuchos::RCP<CArray<int>> calculate_elem_switch(
    Teuchos::RCP<Tpetra::Map<Solver::LO, Solver::GO, Solver::node_type>> all_element_map)
{
    Teuchos::RCP<CArray<int>> elem_switch =
        Teuchos::rcp(new CArray<int>(all_element_map->getLocalNumElements()));

    for (size_t ielem = 0; ielem < elem_switch->dims(0); ielem++)
    {
        (*elem_switch)(ielem) = all_element_map->getGlobalElement(ielem) % 2;
    }
    return elem_switch;
}

Teuchos::RCP<CArray<int>> get_elem_proc_id(
    Teuchos::RCP<Tpetra::Map<Solver::LO, Solver::GO, Solver::node_type>> all_element_map,
    size_t                                                               myrank)
{
    Teuchos::RCP<CArray<int>> elem_proc_id =
        Teuchos::rcp(new CArray<int>(all_element_map->getLocalNumElements()));

    for (size_t ielem = 0; ielem < elem_proc_id->dims(0); ielem++)
    {
        (*elem_proc_id)(ielem) = myrank;
    }
    return elem_proc_id;
}

Teuchos::RCP<CArray<int>> get_elem_gid(
    Teuchos::RCP<Tpetra::Map<Solver::LO, Solver::GO, Solver::node_type>> all_element_map)
{
    Teuchos::RCP<CArray<int>> elem_gid =
        Teuchos::rcp(new CArray<int>(all_element_map->getLocalNumElements()));

    for (size_t ielem = 0; ielem < elem_gid->dims(0); ielem++)
    {
        (*elem_gid)(ielem) = all_element_map->getGlobalElement(ielem);
    }
    return elem_gid;
}

Teuchos::RCP<CArray<double>> get_design_density(
    size_t                         rnum_nodes,
    bool                           topology_optimization_on,
    const Teuchos::RCP<Solver::MV> design_node_densities_distributed)
{
    Teuchos::RCP<CArray<double>> design_density =
        Teuchos::rcp(new CArray<double>(rnum_nodes));

    if (topology_optimization_on)
    {
        auto host_view = design_node_densities_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
        for (size_t inode = 0; inode < rnum_nodes; inode++)
        {
            (*design_density)(inode) = host_view(inode, 0);
        }
    }
    else
    {
        for (size_t inode = 0; inode < rnum_nodes; inode++)
        {
            (*design_density)(inode) = 1.0;
        }
    }

    return design_density;
}
