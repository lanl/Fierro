#include "hdf5_io_functions.h"

hid_t create_hdf5_filespace(int num_dims_, 
    const int *global_array_dims_, const int *local_array_dims_,
    const int *start_coordinates_, const int *stride_)
{
    hsize_t *global_array_dims = new hsize_t[num_dims_];
    hsize_t *local_array_dims = new hsize_t[num_dims_];
    hsize_t *start_coordinates = new hsize_t[num_dims_];
    hsize_t *stride = new hsize_t[num_dims_];

    for (int i = 0; i < num_dims_; i++)
    {
        global_array_dims[i] = global_array_dims_[i];
        local_array_dims[i] = local_array_dims_[i];
        start_coordinates[i] = start_coordinates_[i];
        stride[i] = stride_[i];
    }

    hid_t filespace = H5Screate_simple(num_dims_, global_array_dims, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, 
                        start_coordinates, stride, local_array_dims, NULL);

    delete [] global_array_dims;
    delete [] local_array_dims;
    delete [] start_coordinates;
    delete [] stride;

    return filespace;
}

hid_t create_hdf5_memspace(int num_dims_, const int *local_array_dims_,
    const int *stride_, int num_ghost_nodes_)
{
    /*
        local_array_dims is should not include the ghost_nodes
    */

    hsize_t *global_array_dims = new hsize_t[num_dims_];
    hsize_t *local_array_dims = new hsize_t[num_dims_];
    hsize_t *start_coordinates = new hsize_t[num_dims_];
    hsize_t *stride = new hsize_t[num_dims_];

    for (int i = 0; i < num_dims_; i++)
    {
        global_array_dims[i] = local_array_dims_[i]+2*num_ghost_nodes_;
        local_array_dims[i] = local_array_dims_[i];
        start_coordinates[i] = num_ghost_nodes_;
        stride[i] = stride_[i];
    }

    hid_t memspace = H5Screate_simple(num_dims_, global_array_dims, NULL);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
                        start_coordinates, stride, local_array_dims, NULL);

    delete [] global_array_dims;
    delete [] local_array_dims;
    delete [] start_coordinates;
    delete [] stride;

    return memspace;
}

hid_t create_hdf5_file(const char *filename, MPI_Comm mpi_hdf5_comm)
{
    hid_t file_creation_plist = H5P_DEFAULT;   // File creation property list
    // set the file access template for parallel IO access
    hid_t file_access_plist   = H5P_DEFAULT;   // File access property list
    file_access_plist = H5Pcreate(H5P_FILE_ACCESS);

    // set collective mode for metadata writes
    H5Pset_coll_metadata_write(file_access_plist, true);

    MPI_Info mpi_info = MPI_INFO_NULL; // For MPI IO hints
    MPI_Info_create(&mpi_info);
    MPI_Info_set(mpi_info, "striping_factor", "8");
    MPI_Info_set(mpi_info, "striping_unit", "4194304");

    // tell the HDF5 library that we want to use MPI-IO to do the writing
    H5Pset_fapl_mpio(file_access_plist, mpi_hdf5_comm, mpi_info);

    // Open the file collectively
    // H5F_ACC_TRUNC - overwrite existing file. H5F_ACC_EXCL - no overwrite
    // 3rd argument is file creation property list. Using default here
    // 4th argument is the file access property list identifier
    hid_t file_identifier = H5Fcreate(filename, H5F_ACC_TRUNC, file_creation_plist,
                    file_access_plist);

    // release the file access template
    H5Pclose(file_access_plist);
    MPI_Info_free(&mpi_info);

    return file_identifier;
}

hid_t create_hdf5_dataset(const char *dataset_name, hid_t hdf5_datatype, 
  hid_t file_identifier, hid_t filespace)
{
    // create the dataset
    hid_t link_creation_plist    = H5P_DEFAULT; // Link creation property list
    hid_t dataset_creation_plist = H5P_DEFAULT; // Dataset creation property list
    hid_t dataset_access_plist   = H5P_DEFAULT; // Dataset access property list
    hid_t dataset = H5Dcreate2(
      file_identifier,          // Arg 1: file identifier
      dataset_name,             // Arg 2: dataset name
      hdf5_datatype,            // Arg 3: datatype identifier
      filespace,                // Arg 4: filespace identifier
      link_creation_plist,      // Arg 5: link creation property list
      dataset_creation_plist,   // Arg 6: dataset creation property list
      dataset_access_plist);    // Arg 7: dataset access property list

    return dataset;
}

hid_t open_hdf5_file(const char *filename, MPI_Comm mpi_hdf5_comm)
{
    // set the file access template for parallel IO access
    hid_t file_access_plist = H5P_DEFAULT;   // File access property list
    file_access_plist = H5Pcreate(H5P_FILE_ACCESS);

    // set collective mode for metadata reads (ops)
    H5Pset_all_coll_metadata_ops(file_access_plist, true);

    // tell the HDF5 library that we want to use MPI-IO to do the reading
    H5Pset_fapl_mpio(file_access_plist, mpi_hdf5_comm, MPI_INFO_NULL);

    // Open the file collectively
    // H5F_ACC_RDONLY - sets access to read or write on open of an existing file.
    // 3rd argument is the file access property list identifier
    hid_t file_identifier = H5Fopen(filename, H5F_ACC_RDONLY, file_access_plist);

    // release the file access template
    H5Pclose(file_access_plist);

    return file_identifier;
}

hid_t open_hdf5_dataset(const char *dataset_name, hid_t file_identifier)
{
    
    hid_t dataset_access_plist = H5P_DEFAULT; // Dataset access property list

    hid_t dataset = H5Dopen2(
      file_identifier,        // Arg 1: file identifier
      dataset_name,           // Arg 2: dataset name to match for read
      dataset_access_plist);  // Arg 3: dataset access property list

    return dataset;
}
