#include "Explicit_Solver_SGH.h"
#include "Simulation_Parameters_SGH.h"
#include "FEA_Module_SGH.h"
#include "Simulation_Parameters_Dynamic_Optimization.h"
#include <map>
#include <fstream>
#include <sys/stat.h>

typedef Kokkos::LayoutRight CArrayLayout;
typedef Kokkos::LayoutLeft FArrayLayout;
typedef Tpetra::MultiVector<>::dual_view_type::t_host::array_layout TpetraHostViewLayout;


MPI_Offset
mpi_get_file_position_shared(
  const MPI_Comm comm,
  const MPI_File file_parallel);

void 
write_string_stream_to_file_mpi_all(
  const std::stringstream& str_stream,
  MPI_Offset rank_offset_multiplier,
  const MPI_File file_parallel,
  const MPI_Comm comm);

template <typename ArrayLayout, typename T>
void
write_data_to_string_stream(
  const T* ptr,
  size_t dim0,
  size_t dim1,
  std::stringstream& str_stream);

template<typename ArrayLayout, typename SC, typename LO, typename GO, typename NO>
Teuchos::RCP<Tpetra::MultiVector<SC,LO,GO,NO>>
sort_data(
  const SC* unsorted_data_ptr,
  Teuchos::RCP<Tpetra::Map<LO,GO,NO>> unsorted_map,
  size_t dim1,
  size_t num_global_unique_elements,
  MPI_Comm comm);

template<typename ArrayLayout, typename SC, typename LO, typename GO, typename NO>
void
sort_and_write_data_to_file_mpi_all(
  const SC* unsorted_data_ptr,
  Teuchos::RCP<Tpetra::Map<LO,GO,NO>> unsorted_map,
  size_t dim1,
  size_t num_global_unique_elements,
  MPI_Comm comm, 
  MPI_File file_parallel);

Teuchos::RCP<CArray<int>>
get_cell_nodes(
  const CArray <size_t>& nodes_in_elem,
  size_t num_dim,
  Solver::node_ordering_convention active_node_ordering_convention);

Teuchos::RCP<CArray<double>>
calculate_elem_speed(
  const Teuchos::RCP<Solver::MV> node_vel,
  const Teuchos::RCP<Solver::MCONN> nodes_in_elem);

Teuchos::RCP<CArray<int>>
calculate_elem_switch(
  Teuchos::RCP<Tpetra::Map<Solver::LO,Solver::GO,Solver::node_type>> all_element_map);

Teuchos::RCP<CArray<int>>
get_elem_proc_id(
  Teuchos::RCP<Tpetra::Map<Solver::LO,Solver::GO,Solver::node_type>> all_element_map,
  size_t myrank);

Teuchos::RCP<CArray<int>>
get_elem_gid(
  Teuchos::RCP<Tpetra::Map<Solver::LO,Solver::GO,Solver::node_type>> all_element_map);

Teuchos::RCP<CArray<double>>
get_design_density(
  size_t rnum_nodes,
  bool topology_optimization_on,
  const Teuchos::RCP<Solver::MV> design_node_densities_distributed);

void
Explicit_Solver_SGH::write_outputs()
{
  // No output for OUTPUT_FORMAT::none
  if (simparam.output_options.output_file_format == OUTPUT_FORMAT::none)
    return;

  const size_t rk_level = simparam.rk_num_bins - 1;

  Teuchos::RCP<CArray<double>> design_density;
  Teuchos::RCP<CArray<int>> elem_switch;
  Teuchos::RCP<CArray<int>> elem_proc_id;
  Teuchos::RCP<CArray<int>> elem_gid;
  Teuchos::RCP<CArray<double>> elem_speed;

  for (const FIELD_OUTPUT_SGH& field_name : simparam.field_output) {
    switch (field_name)
    {
      case FIELD_OUTPUT_SGH::design_density:
        // node "design_density"
        design_density = get_design_density(map->getLocalNumElements(),
          simparam_dynamic_opt.topology_optimization_on, design_node_densities_distributed);
        point_data_scalars_double["design_density"] = design_density->pointer();
        break;

      case FIELD_OUTPUT_SGH::velocity:
        // node "velocity"
        sgh_module->node_vel.update_host();
        point_data_vectors_double["velocity"] = &sgh_module->node_vel.host(rk_level,0,0);
        break;

      case FIELD_OUTPUT_SGH::element_density:
        // element "density"
        sgh_module->elem_den.update_host();
        cell_data_scalars_double["element_density"] = sgh_module->elem_den.host_pointer();
        break;

      case FIELD_OUTPUT_SGH::pressure:  
        // element "pressure"
        sgh_module->elem_pres.update_host();
        cell_data_scalars_double["pressure"] = sgh_module->elem_pres.host_pointer();
        break;

      case FIELD_OUTPUT_SGH::SIE:
        // element "SIE"
        sgh_module->elem_sie.update_host();
        cell_data_scalars_double["SIE"] = &sgh_module->elem_sie.host(rk_level,0);
        break;

      case FIELD_OUTPUT_SGH::volume:
        // element "volume"
        sgh_module->elem_vol.update_host();
        cell_data_scalars_double["volume"] = sgh_module->elem_vol.host_pointer();
        break;

      case FIELD_OUTPUT_SGH::mass:
        // element "mass"
        sgh_module->elem_mass.update_host();
        cell_data_scalars_double["mass"] = sgh_module->elem_mass.host_pointer();
        break;

      case FIELD_OUTPUT_SGH::sound_speed:
        // element "sspd"
        sgh_module->elem_sspd.update_host();
        cell_data_scalars_double["sound_speed"] = sgh_module->elem_sspd.host_pointer();
        break;

      case FIELD_OUTPUT_SGH::speed:
        // element "speed"
        elem_speed = calculate_elem_speed(all_node_velocities_distributed, global_nodes_in_elem_distributed);
        cell_data_scalars_double["speed"] = elem_speed->pointer();
        break;

      case FIELD_OUTPUT_SGH::material_id:
        // element "material_id"
        sgh_module->elem_mat_id.update_host();
        cell_data_scalars_int["material_id"] = reinterpret_cast<int*>(sgh_module->elem_mat_id.host_pointer());
        break;

      case FIELD_OUTPUT_SGH::element_switch:
        // element "element_switch"
        elem_switch = calculate_elem_switch(all_element_map);
        cell_data_scalars_int["element_switch"] = elem_switch->pointer();
        break;

      case FIELD_OUTPUT_SGH::processor_id:
        // element "processor_id"
        elem_proc_id = get_elem_proc_id(all_element_map, myrank);
        cell_data_scalars_int["processor_id"] = elem_proc_id->pointer();
        break;

      case FIELD_OUTPUT_SGH::element_id:
        // element "element_id"
        elem_gid = get_elem_gid(all_element_map);
        cell_data_scalars_int["element_id"] = elem_gid->pointer();
        break;

      case FIELD_OUTPUT_SGH::user_vars:
        // element "user_vars"
        sgh_module->elem_user_output_vars.update_host();
        cell_data_fields_double["user_vars"] = std::make_pair(sgh_module->elem_user_output_vars.host_pointer(), 
                                                                   sgh_module->elem_user_output_vars.dims(1));
      case FIELD_OUTPUT_SGH::stress:
        // element "stress"
        sgh_module->elem_stress.update_host();
        cell_data_fields_double["stress"] = std::make_pair(&sgh_module->elem_stress.host(rk_level,0,0,0), 9);
        break;

      default:
        break;

    } // end switch
  } // end if

  switch (simparam.output_options.output_file_format)
  {
    case OUTPUT_FORMAT::vtk:
      parallel_vtk_writer_new();
      break;

    case OUTPUT_FORMAT::vtu:
      parallel_vtu_writer_new();
      break;

    default:
      break;
  };
}

void
Explicit_Solver_SGH::parallel_vtu_writer_new()
{
  /* to be added... */
  throw std::runtime_error("parallel_vtu_writer_new() not yet implemented. use parallel_vtk_writer_new()");
}

std::string
construct_file_name(
  size_t file_count,
  int displacement_module)
{
  bool displace_geometry = false;
  std::string current_file_name;
  std::string base_file_name= "VTK";
  std::string base_file_name_undeformed= "VTK_undeformed";
	std::string file_extension= ".vtk";

  if(displace_geometry && displacement_module>=0)
    current_file_name = base_file_name_undeformed + std::to_string(file_count) + file_extension;
  else
    current_file_name = base_file_name + std::to_string(file_count) + file_extension;
  return current_file_name;
}

void
Explicit_Solver_SGH::parallel_vtk_writer_new()
{
  const size_t rk_level = simparam.rk_num_bins - 1;

  int num_dim = simparam.num_dims;
  std::stringstream str_stream;
  MPI_Offset current_offset;
  MPI_File myfile_parallel;
  bool displace_geometry = false;

  std::string vtk_dir = "vtk/";
  std::string vtk_data_dir = vtk_dir + "data/";

  // mkdir if needed
  struct stat st;
  if (myrank == 0) {
    if (stat(vtk_dir.c_str(), &st) != 0) {
      str_stream.str("");
      str_stream << "mkdir" << " " << vtk_dir;
      system(str_stream.str().c_str());
    }
    if (stat(vtk_data_dir.c_str(), &st) != 0) {
      str_stream.str("");
      str_stream << "mkdir" << " " << vtk_data_dir;
      system(str_stream.str().c_str());
    }
  }
  MPI_Barrier(world);

  //construct file name
  std::string current_file_name = construct_file_name(simparam.graphics_options.graphics_id, displacement_module);
  std::string file_path = vtk_data_dir + current_file_name;
  //open mpi file
  int err = MPI_File_open(MPI_COMM_WORLD, file_path.c_str(), 
                          MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &myfile_parallel);
  //allows overwriting the file if it already existed in the directory
  if (err != MPI_SUCCESS)  {
    if (myrank == 0){
      MPI_File_delete(file_path.c_str(),MPI_INFO_NULL);
    }
    MPI_File_open(MPI_COMM_WORLD, file_path.c_str(),
                  MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &myfile_parallel);
  }

  /*************** write header of the vtk file ***************/
  str_stream.str("");
  str_stream << "# vtk DataFile Version 2.0" << std::endl;
  str_stream << "Mesh for Fierro" << std::endl;
  str_stream << "ASCII" << std::endl;
  str_stream << "DATASET UNSTRUCTURED_GRID" << std::endl;
  current_offset = mpi_get_file_position_shared(world, myfile_parallel);
  if (myrank == 0) {
    MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(), str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
  }


  /*************** write POINTS ***************/
  str_stream.str("");
  str_stream << std::endl << "POINTS " << num_nodes << " float" << std::endl;
  current_offset = mpi_get_file_position_shared(world, myfile_parallel);
  if(myrank == 0) {
    MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(), str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
  }
  sgh_module->node_coords.update_host();
  sort_and_write_data_to_file_mpi_all <CArrayLayout,double,LO,GO,node_type> (
    &sgh_module->node_coords.host(rk_level,0,0), map, num_dim, num_nodes, world, myfile_parallel);


  /*************** write CELLS ***************/
  str_stream.str("");
  str_stream << std::endl << "CELLS " << num_elem << " " << num_elem*(max_nodes_per_element+1) << std::endl;
  current_offset = mpi_get_file_position_shared(world, myfile_parallel);
  if(myrank == 0) {
    MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(), str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
  }
  CArray <size_t> nodes_in_elem (sgh_module->nodes_in_elem.dims(0), sgh_module->nodes_in_elem.dims(1));
  { //view scope
    auto host_view = global_nodes_in_elem_distributed->getLocalView<HostSpace> (Tpetra::Access::ReadOnly);
    for (size_t ielem = 0; ielem < nodes_in_elem.dims(0); ielem++) {
      for (size_t inode = 0; inode < nodes_in_elem.dims(1); inode++) {
        nodes_in_elem(ielem, inode) = host_view(ielem, inode);
      }
    }
  } //end view scope
  auto cell_data = get_cell_nodes(nodes_in_elem, num_dim, active_node_ordering_convention);
  sort_and_write_data_to_file_mpi_all <CArrayLayout,int,LO,GO,node_type> (
    cell_data->pointer(), all_element_map, cell_data->dims(1), num_elem, world, myfile_parallel);


  /*************** write CELL_TYPES ***************/
  str_stream.str("");
	str_stream << std::endl << "CELL_TYPES " << num_elem << std::endl;
  current_offset = mpi_get_file_position_shared(world, myfile_parallel);
  if(myrank == 0) {
    MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(), str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
  }
  CArray<int> cell_type (all_element_map->getLocalNumElements(), 1);
  for (int i = 0; i < cell_type.dims(0); i++)
    for (int j = 0; j < cell_type.dims(1); j++)
      cell_type(i,j) = 12;
  sort_and_write_data_to_file_mpi_all <CArrayLayout,int,LO,GO,node_type> (
    cell_type.pointer(), all_element_map, cell_type.dims(1), num_elem, world, myfile_parallel);


  /*************** write POINT_DATA ***************/
  str_stream.str("");
	str_stream << std::endl << "POINT_DATA " << num_nodes;
  current_offset = mpi_get_file_position_shared(world, myfile_parallel);
  if(myrank == 0) {
    MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(),
                      str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
  }

  //SCALARS float
  for (auto it = point_data_scalars_double.begin(); it != point_data_scalars_double.end(); it++) {
    str_stream.str("");
    str_stream << std::endl << "SCALARS " << it->first << " float 1" << std::endl;
    str_stream << "LOOKUP_TABLE default" << std::endl;
    current_offset = mpi_get_file_position_shared(world, myfile_parallel);
    if(myrank == 0) {
      MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(),
                        str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    sort_and_write_data_to_file_mpi_all <CArrayLayout,double,LO,GO,node_type> (
      it->second, map, 1, num_nodes, world, myfile_parallel);
  }

  //VECTORS float
  for (auto it = point_data_vectors_double.begin(); it != point_data_vectors_double.end(); it++) {
    str_stream.str("");
    str_stream << std::endl << "VECTORS " << it->first << " float" << std::endl;
    current_offset = mpi_get_file_position_shared(world, myfile_parallel);
    if(myrank == 0) {
      MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(),
                        str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    sort_and_write_data_to_file_mpi_all <CArrayLayout,double,LO,GO,node_type> (
      it->second, map, num_dim, num_nodes, world, myfile_parallel);
  }


  /*************** write CELL_DATA ***************/
  str_stream.str("");
	str_stream << std::endl << "CELL_DATA " << num_elem;
  current_offset = mpi_get_file_position_shared(world, myfile_parallel);
  if(myrank == 0) {
    MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(),
                      str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
  }

  //SCALARS float
  for (auto it = cell_data_scalars_double.begin(); it != cell_data_scalars_double.end(); it++) {
    str_stream.str("");
    str_stream << std::endl << "SCALARS " << it->first << " float 1" << std::endl;
    str_stream << "LOOKUP_TABLE default" << std::endl;
    current_offset = mpi_get_file_position_shared(world, myfile_parallel);
    if(myrank == 0) {
      MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(),
                        str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    sort_and_write_data_to_file_mpi_all <CArrayLayout,double,LO,GO,node_type> (
      it->second, all_element_map, 1, num_elem, world, myfile_parallel);
  }

  //SCALARS int
  for (auto it = cell_data_scalars_int.begin(); it != cell_data_scalars_int.end(); it++) {
    str_stream.str("");
    str_stream << std::endl << "SCALARS " << it->first << " int 1" << std::endl;
    str_stream << "LOOKUP_TABLE default" << std::endl;
    current_offset = mpi_get_file_position_shared(world, myfile_parallel);
    if(myrank == 0) {
      MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(),
                        str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    sort_and_write_data_to_file_mpi_all <CArrayLayout,int,LO,GO,node_type> (
      it->second, all_element_map, 1, num_elem, world, myfile_parallel);
  }

  //FIELD
  str_stream.str("");
  str_stream << std::endl << "FIELD FieldData " << cell_data_fields_double.size() << std::endl;
  current_offset = mpi_get_file_position_shared(world, myfile_parallel);
  if(myrank == 0) {
    MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(),
                      str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
  }
  for (auto it = cell_data_fields_double.begin(); it != cell_data_fields_double.end(); it++) {
    str_stream.str("");
    str_stream << it->first << " " << it->second.second << " " << num_elem << " float" << std::endl;
    current_offset = mpi_get_file_position_shared(world, myfile_parallel);
    if(myrank == 0) {
      MPI_File_write_at(myfile_parallel, current_offset, str_stream.str().c_str(),
                        str_stream.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    sort_and_write_data_to_file_mpi_all <CArrayLayout,double,LO,GO,node_type> (
      it->second.first, all_element_map, it->second.second, num_elem, world, myfile_parallel);
  }
  
  MPI_Barrier(world);
  MPI_File_sync(myfile_parallel);
  MPI_File_close(&myfile_parallel);


  /*************** write .vtk.series file ***************/
  simparam.graphics_options.graphics_times(simparam.graphics_options.graphics_id) = sgh_module->simparam.time_value;
  if (myrank == 0) {
    FILE *myfile;
    std::string filename = vtk_dir + "outputs.vtk.series";
    myfile = fopen(filename.c_str(), "w");
    
    fprintf(myfile, "{\n");
    fprintf(myfile, "  \"file-series-version\" : \"1.0\",\n");
    fprintf(myfile, "  \"files\" : [\n");
    for (int i = 0; i <= simparam.graphics_options.graphics_id; i++) {
      std::string vtk_filename = construct_file_name(i, displacement_module);
      fprintf(myfile, "    { \"name\" : \"data/%s\", \"time\" : %12.5e },\n",
        vtk_filename.c_str(), simparam.graphics_options.graphics_times(i));
    }
    fprintf(myfile, "  ]\n");
    fprintf(myfile, "}\n");

    fclose(myfile);
  }
  
  simparam.graphics_options.graphics_id++;

}


MPI_Offset
mpi_get_file_position_shared(
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

void
write_string_stream_to_file_mpi_all(
  const std::stringstream& str_stream,
  MPI_Offset rank_offset_multiplier,
  const MPI_File file_parallel,
  const MPI_Comm comm)
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

template <typename ArrayLayout, typename T>
void write_data_to_string_stream(
  const T* ptr,
  size_t dim0,
  size_t dim1,
  std::stringstream& str_stream)
{
  size_t w;
  size_t p;
  //crude way to ensure that all lines have the same length. required for MPI IO
  if constexpr (std::is_same<T, float>::value or
                std::is_same<T, double>::value) {
    w = 25;
    p = 8;
  }
  else if constexpr (std::is_same<T, int>::value or 
                     std::is_same<T, size_t>::value or 
                     std::is_same<T, Solver::GO>::value or
                     std::is_same<T, Solver::LO>::value) {
    w = 10;
    p = 0;
  }
  else {
    throw std::runtime_error(std::string("unknown datatype in function ") + 
                             std::string(__FUNCTION__) + std::string("in file ") +
                             std::string(__FILE__));
  }

  // create view of unsorted_data_ptr. not using matar here. array_layout is needed
  Kokkos::View<const T**, ArrayLayout, HostSpace, Kokkos::MemoryUnmanaged>
    data_view(ptr, dim0, dim1);

  str_stream.str("");
  str_stream << std::fixed << std::setprecision(p);
  for (size_t i = 0; i < dim0; i++) {
    for (size_t j = 0; j < dim1; j++) {
      str_stream << std::left << std::setw(w) << data_view(i,j) << " ";
    }
    str_stream << std::endl;
  }

  return;
}

template<typename ArrayLayout, typename SC, typename LO, typename GO, typename NO>
Teuchos::RCP<Tpetra::MultiVector<SC,LO,GO,NO>>
sort_data(
  const SC* unsorted_data_ptr,
  Teuchos::RCP<Tpetra::Map<LO,GO,NO>> unsorted_map,
  size_t dim1,
  size_t num_global_unique_elements,
  MPI_Comm comm)
{

  // create view of unsorted_data_ptr. not using matar here. array_layout is needed
  Kokkos::View<const SC**, ArrayLayout, HostSpace, Kokkos::MemoryUnmanaged>
    unsorted_data(unsorted_data_ptr, unsorted_map->getLocalNumElements(), dim1);

  // sorted_map
  Teuchos::RCP<Tpetra::Map<LO,GO,NO>> sorted_map = 
    Teuchos::rcp(new Tpetra::Map<LO,GO,NO>(
    num_global_unique_elements, unsorted_map->getIndexBase(), unsorted_map->getComm()));

  // importer
  Tpetra::Import<LO,GO,NO> sorting_importer(unsorted_map, sorted_map);

  // sorted storage
  Teuchos::RCP<Tpetra::MultiVector<SC,LO,GO,NO>> sorted_storage = 
    Teuchos::rcp(new Tpetra::MultiVector<SC,LO,GO,NO>(sorted_map, dim1));

  // unsorted storage
  Teuchos::RCP<Tpetra::MultiVector<SC,LO,GO,NO>> unsorted_storage = 
    Teuchos::rcp(new Tpetra::MultiVector<SC,LO,GO,NO>(unsorted_map, dim1));

  { //view scope
    auto host_view = unsorted_storage->getLocalViewHost(Tpetra::Access::ReadWrite);
    for(int i = 0; i < unsorted_map->getLocalNumElements(); i++) {
      for (int j = 0; j < dim1; j++) {
        host_view(i,j) = unsorted_data(i,j);
      }
    } // for
  } //end view scope

  sorted_storage->doImport(*unsorted_storage, sorting_importer, Tpetra::INSERT);

  return sorted_storage;
}

template<typename ArrayLayout, typename SC, typename LO, typename GO, typename NO>
void
sort_and_write_data_to_file_mpi_all(
  const SC* unsorted_data_ptr,
  Teuchos::RCP<Tpetra::Map<LO,GO,NO>> unsorted_map,
  size_t dim1,
  size_t num_global_unique_elements,
  MPI_Comm comm, 
  MPI_File file_parallel)
{

  auto sorted_data = 
    sort_data <ArrayLayout,SC,LO,GO,NO> (unsorted_data_ptr, unsorted_map, dim1, num_global_unique_elements, comm);

  { //view scope
    auto const_host_view = sorted_data->getLocalViewHost(Tpetra::Access::ReadOnly);
    std::stringstream str_stream;
    write_data_to_string_stream <TpetraHostViewLayout,SC> (const_host_view.data(), sorted_data->getMap()->getLocalNumElements(), dim1, str_stream);
    write_string_stream_to_file_mpi_all(str_stream, sorted_data->getMap()->getGlobalElement(0), file_parallel, comm);
  } //end view scope
}

Teuchos::RCP<CArray<int>>
get_cell_nodes(
  const CArray <size_t>& nodes_in_elem,
  size_t num_dim,
  Solver::node_ordering_convention active_node_ordering_convention)
{
  CArray<size_t> convert_ijk_to_ensight(nodes_in_elem.dims(1));
  convert_ijk_to_ensight(0) = 0;
  convert_ijk_to_ensight(1) = 1;
  convert_ijk_to_ensight(2) = 3;
  convert_ijk_to_ensight(3) = 2;
  if(num_dim==3){
    convert_ijk_to_ensight(4) = 4;
    convert_ijk_to_ensight(5) = 5;
    convert_ijk_to_ensight(6) = 7;
    convert_ijk_to_ensight(7) = 6;
  }

  Teuchos::RCP<CArray<int>> cell_data = 
    Teuchos::rcp( new CArray<int>(nodes_in_elem.dims(0), nodes_in_elem.dims(1)+1) );

  size_t temp_convert;
  for (size_t ielem = 0; ielem < nodes_in_elem.dims(0); ielem++) {
    (*cell_data)(ielem,0) = 8;
		for (int ii = 0; ii < nodes_in_elem.dims(1); ii++) {
      if(active_node_ordering_convention == Solver::node_ordering_convention::IJK)
        temp_convert = convert_ijk_to_ensight(ii);
      else
        temp_convert = ii;
      (*cell_data)(ielem,ii+1) = nodes_in_elem(ielem,temp_convert);
      //(*cell_data)(ielem,ii+1) = nodes_in_elem(ielem,ii);
    }
  }

  return cell_data;
}


Teuchos::RCP<CArray<double>>
calculate_elem_speed(
  const Teuchos::RCP<Solver::MV> all_node_vel_distributed,
  const Teuchos::RCP<Solver::MCONN> global_nodes_in_elem_distributed)
{

  size_t rnum_elems = global_nodes_in_elem_distributed->getLocalLength();
  size_t num_nodes_in_elem = global_nodes_in_elem_distributed->getNumVectors();
  size_t num_dims = all_node_vel_distributed->getNumVectors();

  Teuchos::RCP<CArray<double>> elem_speed = 
    Teuchos::rcp(new CArray<double>(rnum_elems));

  auto nodes_in_elem_hview = global_nodes_in_elem_distributed->getLocalViewHost(Tpetra::Access::ReadOnly);
  auto all_node_vel_hview = all_node_vel_distributed->getLocalViewHost(Tpetra::Access::ReadOnly);

  for (size_t elem_gid = 0; elem_gid < rnum_elems; elem_gid++) { 
    double elem_vel[3];
    elem_vel[0] = 0.0;
    elem_vel[1] = 0.0;
    elem_vel[2] = 0.0;
    // get the coordinates of the element center
    for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
      size_t inode = all_node_vel_distributed->getMap()->getLocalElement(nodes_in_elem_hview(elem_gid, node_lid));
      elem_vel[0] += all_node_vel_hview(inode, 0);
      elem_vel[1] += all_node_vel_hview(inode, 1);
      if (num_dims == 3){
        elem_vel[2] += all_node_vel_hview(inode, 2);
      }
      else {
        elem_vel[2] = 0.0;
      }
    } // end loop over nodes in element

    elem_vel[0] = elem_vel[0]/num_nodes_in_elem;
    elem_vel[1] = elem_vel[1]/num_nodes_in_elem;
    elem_vel[2] = elem_vel[2]/num_nodes_in_elem;

    double speed_sqrd = 0.0;
    for (int dim=0; dim<num_dims; dim++){
      speed_sqrd += elem_vel[dim]*elem_vel[dim];
    }

    (*elem_speed)(elem_gid) = sqrt(speed_sqrd);
  } // end for

  return elem_speed;

}


Teuchos::RCP<CArray<int>>
calculate_elem_switch(
  Teuchos::RCP<Tpetra::Map<Solver::LO,Solver::GO,Solver::node_type>> all_element_map)
{
  Teuchos::RCP<CArray<int>> elem_switch = 
    Teuchos::rcp(new CArray<int>(all_element_map->getLocalNumElements()));

  for (size_t ielem = 0; ielem < elem_switch->dims(0); ielem++) {
    (*elem_switch)(ielem) = all_element_map->getGlobalElement(ielem) % 2;
  }
  return elem_switch;
}

Teuchos::RCP<CArray<int>>
get_elem_proc_id(
  Teuchos::RCP<Tpetra::Map<Solver::LO,Solver::GO,Solver::node_type>> all_element_map,
  size_t myrank)
{
  Teuchos::RCP<CArray<int>> elem_proc_id = 
    Teuchos::rcp(new CArray<int>(all_element_map->getLocalNumElements()));

  for (size_t ielem = 0; ielem < elem_proc_id->dims(0); ielem++) {
    (*elem_proc_id)(ielem) = myrank;
  }
  return elem_proc_id;
}

Teuchos::RCP<CArray<int>>
get_elem_gid(
  Teuchos::RCP<Tpetra::Map<Solver::LO,Solver::GO,Solver::node_type>> all_element_map)
{
  Teuchos::RCP<CArray<int>> elem_gid = 
    Teuchos::rcp(new CArray<int>(all_element_map->getLocalNumElements()));

  for (size_t ielem = 0; ielem < elem_gid->dims(0); ielem++) {
    (*elem_gid)(ielem) = all_element_map->getGlobalElement(ielem);
  }
  return elem_gid;
}

Teuchos::RCP<CArray<double>>
get_design_density(
  size_t rnum_nodes,
  bool topology_optimization_on,
  const Teuchos::RCP<Solver::MV> design_node_densities_distributed)
{
  Teuchos::RCP<CArray<double>> design_density =
    Teuchos::rcp(new CArray<double>(rnum_nodes));

  if(topology_optimization_on) {
    auto host_view = design_node_densities_distributed->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
    for (size_t inode = 0; inode < rnum_nodes; inode++) {
      (*design_density)(inode) = host_view(inode,0);
    }
  }
  else {
    for (size_t inode = 0; inode < rnum_nodes; inode++) {
      (*design_density)(inode) = 1.0;
    }
  }

  return design_density;
}

