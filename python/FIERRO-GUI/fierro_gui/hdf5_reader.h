#pragma once

#include "hdf5.h"
#include "mpi.h"
#include <iostream>

hid_t open_hdf5_file(const char *filename, MPI_Comm mpi_hdf5_comm);

hid_t open_hdf5_dataset(const char *dataset_name, hid_t file_identifier);
