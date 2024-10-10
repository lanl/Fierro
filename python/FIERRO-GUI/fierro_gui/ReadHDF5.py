def ReadHDF5(filename, dataset_name):
    import open_hdf5_file # type: ignore
    #hid_t open_hdf5_file(const char *filename, const char *dataset_name, MPI_Comm mpi_hdf5_comm)
    print ("Opening HDF5 file: " + filename)
    open_hdf5_file.open_file(filename, dataset_name, mpi_hdf5_comm)

    
