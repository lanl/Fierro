#include "hdf5_reader.h"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

//mpicc -O3 -Wall -shared -std=c++17 -undefined -fPIC $(python3 -m pybind11 --includes) hdf5_reader.cpp -o open_hdf5_file$(python3-config --extension-suffix)

//CC=/usr/bin/mpicc ./configure --enable-parallel --prefix=/usr/local/hdf5

namespace py = pybind11;

hid_t open_hdf5_file(const char *filename, const char *dataset_name, MPI_Comm mpi_hdf5_comm)
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

    hid_t dataset_access_plist = H5P_DEFAULT; // Dataset access property list

    hid_t dataset = H5Dopen2(
      file_identifier,        // Arg 1: file identifier
      dataset_name,           // Arg 2: dataset name to match for read
      dataset_access_plist);  // Arg 3: dataset access property list

    return dataset;
}

PYBIND11_MODULE(open_hdf5_file, m) {
    m.doc() = "";

    m.def("open_file", &PyInit_open_hdf5_file, "Open hdf5 file");
}