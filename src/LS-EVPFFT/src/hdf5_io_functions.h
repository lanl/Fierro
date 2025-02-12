#pragma once

#include "hdf5.h"
#include "mpi.h"
#include <iostream>

hid_t create_hdf5_filespace(int num_dims_, 
    const int *global_array_dims_, const int *local_array_dims_,
    const int *start_coordinates_, const int *stride_);

hid_t create_hdf5_memspace(int num_dims_, const int *local_array_dims_,
const int *stride_, int num_ghost_nodes_);

hid_t create_hdf5_file(const char *filename, MPI_Comm mpi_hdf5_comm);

hid_t create_hdf5_dataset(const char *dataset_name, hid_t hdf5_datatype, 
  hid_t file_identifier, hid_t filespace);

hid_t open_hdf5_file(const char *filename, MPI_Comm mpi_hdf5_comm);

hid_t open_hdf5_dataset(const char *dataset_name, hid_t file_identifier);

template <typename R>
hid_t deduce_hdf5_datatype(R *data_ptr);

template <typename R>
void write_hdf5_dataset(const char *dataset_name, R *data_ptr,
  hid_t file_identifier, hid_t memspace, hid_t filespace, MPI_Comm mpi_hdf5_comm);

template <typename R>
void read_hdf5_dataset(const char *dataset_name, R *data_ptr,
  hid_t file_identifier, hid_t memspace, hid_t filespace, MPI_Comm mpi_hdf5_comm);


template <typename R>
hid_t deduce_hdf5_datatype(R *data_ptr)
{
    hid_t hdf5_datatype;

    if constexpr (std::is_same<R, int>::value) {
        hdf5_datatype = H5T_NATIVE_INT;
    }
    else if constexpr (std::is_same<R, float>::value) {
        hdf5_datatype = H5T_NATIVE_FLOAT;
    }
    else if constexpr (std::is_same<R, double>::value) {
        hdf5_datatype = H5T_NATIVE_DOUBLE;
    }
    else if constexpr (std::is_same<R, long double>::value) {
        hdf5_datatype = H5T_NATIVE_LDOUBLE;
    }
    else {
        throw std::runtime_error("HDF5 datatype cannot be deduced.\n");
    }

    return hdf5_datatype;
}

template <typename R>
void write_hdf5_dataset(const char *dataset_name, R *data_ptr,
  hid_t file_identifier, hid_t memspace, hid_t filespace, MPI_Comm mpi_hdf5_comm)
{
    hid_t hdf5_datatype = deduce_hdf5_datatype <R> (data_ptr);

    hid_t dataset = create_hdf5_dataset(dataset_name, hdf5_datatype, file_identifier, filespace);

    // Create property list for collective dataset write.
    hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);

    // write the data to disk using both the memory space and the data space.
    H5Dwrite(dataset, hdf5_datatype, memspace, filespace, xfer_plist,
             data_ptr);

    H5Dclose(dataset);

    H5Pclose(xfer_plist);
}

template <typename R>
void read_hdf5_dataset(const char *dataset_name, R *data_ptr,
  hid_t file_identifier, hid_t memspace, hid_t filespace, MPI_Comm mpi_hdf5_comm)
{
    hid_t hdf5_datatype = deduce_hdf5_datatype <R> (data_ptr);

    hid_t dataset = open_hdf5_dataset(dataset_name, file_identifier);

    // Create property list for collective dataset write.
    hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);

    // read the data from disk using both the memory space and the data space.
    H5Dread(dataset, hdf5_datatype, memspace, filespace, H5P_DEFAULT, data_ptr);

    H5Dclose(dataset);

    H5Pclose(xfer_plist);
}
