#pragma once

#include "hdf5_io_functions.h"

template <int ARRAY_LAYOUT>
struct MicroOutputWriter
{
    MicroOutputWriter(MPI_Comm mpi_hdf5_comm, const char *filename,
        const int num_dims, const int *global_array_dims,
        const int *local_array_dims, const int *start_coordinates,
        const int *stride=nullptr);
    MicroOutputWriter(const MicroOutputWriter &) = delete;
    ~MicroOutputWriter();

    template <typename T>
    void reverse_array(T *arr, int start, int end);

    template<typename T>
    void write(const char *dataset_name, T *data_ptr);

    template<typename T>
    void write(const char *dataset_name, T *data_ptr, hid_t group_id);

//private:
    MPI_Comm mpi_hdf5_comm_;
    int num_dims_;
    int *global_array_dims_;
    int *local_array_dims_;
    int *start_coordinates_;
    int *stride_;
    hid_t filespace_;
    hid_t memspace_;
    hid_t file_id_;
    const char *filename_;
};


template <int ARRAY_LAYOUT>
MicroOutputWriter<ARRAY_LAYOUT>::MicroOutputWriter(MPI_Comm mpi_hdf5_comm, const char *filename,
        const int num_dims, const int *global_array_dims, const int *local_array_dims,
        const int *start_coordinates, const int *stride)
    : mpi_hdf5_comm_ (mpi_hdf5_comm)
    , filename_ (filename)
    , num_dims_ (num_dims)
    , filespace_ (H5S_NULL)
    , memspace_ (H5S_NULL)
{
#ifndef ABSOLUTE_NO_OUTPUT
    global_array_dims_ = new int[num_dims_];
    local_array_dims_ = new int[num_dims_];
    start_coordinates_ = new int[num_dims_];
    stride_ = new int[num_dims_];

    for (int i = 0; i < num_dims_; i++)
    {
        global_array_dims_[i] = global_array_dims[i];
        local_array_dims_[i] = local_array_dims[i];
        start_coordinates_[i] = start_coordinates[i];
        if (stride == nullptr)
        {
          stride_[i] = 1;
        }
        else
        {
          stride_[i] = stride[i];
        }
    }

    if constexpr (ARRAY_LAYOUT == MPI_ORDER_FORTRAN)
    {
        // reverse arrays because hdf5 called from c code writes in c order
        reverse_array(global_array_dims_, 0, num_dims_-1);
        reverse_array(local_array_dims_, 0, num_dims_-1);
        reverse_array(start_coordinates_, 0, num_dims_-1);
        reverse_array(stride_, 0, num_dims_-1);
    }

    filespace_ = create_hdf5_filespace(3, global_array_dims_, local_array_dims_,
                                       start_coordinates_, stride_);
    memspace_ = create_hdf5_memspace(3, local_array_dims_, stride_, 0);
    file_id_ = create_hdf5_file(filename, mpi_hdf5_comm_);
#endif
}


template <int ARRAY_LAYOUT>
MicroOutputWriter<ARRAY_LAYOUT>::~MicroOutputWriter()
{
#ifndef ABSOLUTE_NO_OUTPUT
    delete [] global_array_dims_;
    delete [] local_array_dims_;
    delete [] start_coordinates_;
    delete [] stride_;

    H5Sclose(filespace_);
    filespace_ = H5S_NULL;

    H5Sclose(memspace_);
    memspace_ = H5S_NULL;

    H5Fclose(file_id_);
#endif
};


template <int ARRAY_LAYOUT>
template<typename T>
void MicroOutputWriter<ARRAY_LAYOUT>::reverse_array(T *arr, int start, int end)
{
    /* Function to reverse T *arr from start to end*/

    while (start < end)
    {
        T temp = arr[start];
        arr[start] = arr[end];
        arr[end] = temp;
        start++;
        end--;
    }
}


template <int ARRAY_LAYOUT>
template<typename T>
void MicroOutputWriter<ARRAY_LAYOUT>::write(const char *dataset_name, T *data_ptr)
{
    /* dataset_name: can be absolute link path or a relative link path */

    write_hdf5_dataset(dataset_name, data_ptr, file_id_,
                       memspace_, filespace_, mpi_hdf5_comm_);
}

template <int ARRAY_LAYOUT>
template<typename T>
void MicroOutputWriter<ARRAY_LAYOUT>::write(const char *dataset_name, T *data_ptr, hid_t group_id)
{
    /* dataset_name: relative link path to group_id */

    write_hdf5_dataset(dataset_name, data_ptr, group_id,
                       memspace_, filespace_, mpi_hdf5_comm_);
}
