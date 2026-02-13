#pragma once

#include <mpi.h>
#include <array>
#include <string>
#include <cstring>
#include <iostream>

template <typename R, int ARRAY_LAYOUT>
class Manager_MPI_IO
{
private:
    MPI_Comm mpi_io_comm_;
    int root;
    int my_rank_;
    int num_ranks_;
    std::array<int,3> dimensions_full_array_;
    std::array<int,3> dimensions_subarray_;
    std::array<int,3> start_coordinates_;
    std::array<int,3> local_dimensions_full_array_;  // with ghost nodes
    std::array<int,3> local_dimensions_subarray_;    // without ghost nodes
    std::array<int,3> local_start_coordinates_;      // without ghost nodes
    const int num_dims_;
    const int num_ghost_nodes_;
    MPI_Datatype file_space_type_;
    MPI_Datatype memory_space_type_;
    MPI_Datatype mpi_native_type_;

    // for writing ASCII
    std::string format_;
    int chars_per_num_;
    MPI_Datatype chars_per_num_type_;
    MPI_Datatype ascii_file_space_type_;
    MPI_Datatype ascii_memory_space_type_;

    MPI_File create_mpi_io_file_(const char *filename) const;
    MPI_File open_mpi_io_file_(const char *filename) const;

public:
    Manager_MPI_IO(MPI_Comm mpi_io_comm,
                   const std::array<int,3> & dimensions_full_array,
                   const std::array<int,3> & dimensions_subarray,
                   const std::array<int,3> & start_coordinates,
                   int num_dims,
                   int num_ghost_nodes);
    Manager_MPI_IO(const Manager_MPI_IO &) = delete;
    ~Manager_MPI_IO();
    void write_binary(const char *filename, R *data, int file_offset) const;
    void read_binary(const char *filename, R *data, int file_offset) const;
    void write_ascii(const char *filename, R *data, const char* header_text) const;
    void write_ascii(const char *filename, const char *data_as_chars, const char* header_text, int chars_per_line) const;
    
};

template <typename R, int ARRAY_LAYOUT>
Manager_MPI_IO<R,ARRAY_LAYOUT>::Manager_MPI_IO(MPI_Comm mpi_io_comm,
                   const std::array<int,3> & dimensions_full_array,
                   const std::array<int,3> & dimensions_subarray,
                   const std::array<int,3> & start_coordinates,
                   int num_dims,
                   int num_ghost_nodes)
    : mpi_io_comm_ (mpi_io_comm)
    , root (0)
    , dimensions_full_array_ (dimensions_full_array)
    , dimensions_subarray_ (dimensions_subarray)
    , start_coordinates_ (start_coordinates)
    , num_dims_(num_dims)
    , num_ghost_nodes_ (num_ghost_nodes)
    , mpi_native_type_ (MPI_DATATYPE_NULL)
    , format_ ("%013.6E\n")
{
    /*
    Note that dimensions_full_array, etc. should be in the order the original data was
    allocated. For example if C layout, {nz,ny,nx} and if FORTRAN layout, {nx,ny,nz}.
    More example: if the array to be used with mpi_io was allocated as {nz,ny,nx}
                  then dimensions_full_array, etc. should be in the form {z,y,x}.
    */

    MPI_Comm_rank(mpi_io_comm_, &my_rank_);
    MPI_Comm_size(mpi_io_comm_, &num_ranks_);

    local_dimensions_full_array_ = {dimensions_subarray_[0]+2*num_ghost_nodes_,
                                    dimensions_subarray_[1]+2*num_ghost_nodes_,
                                    dimensions_subarray_[2]+2*num_ghost_nodes_};
    local_dimensions_subarray_ = dimensions_subarray_;
    local_start_coordinates_ = {num_ghost_nodes_, num_ghost_nodes_, num_ghost_nodes_};

    // set correct mpiType dependent on typename R
    if constexpr (std::is_same<R, int>::value)
    {
        mpi_native_type_ = MPI_INT;
    }
    else if constexpr (std::is_same<R, float>::value)
    {
        mpi_native_type_ = MPI_FLOAT;
    }
    else if constexpr (std::is_same<R, double>::value)
    {
        mpi_native_type_ = MPI_DOUBLE;
    }
    else if constexpr (std::is_same<R, long double>::value)
    {
        mpi_native_type_ = MPI_LONG_DOUBLE;
    }
    else
    {
        throw std::runtime_error("mpi_native_type_ can not be deduced in Manager_MPI_IO. Please add the mpi type.\n");
    }

    // create file_space_type_
    MPI_Type_create_subarray(num_dims_, dimensions_full_array_.data(), dimensions_subarray_.data(),
                             start_coordinates_.data(), ARRAY_LAYOUT, mpi_native_type_,
                             &file_space_type_);
    MPI_Type_commit(&file_space_type_);

    // create memory_space_type_
    MPI_Type_create_subarray(num_dims_, local_dimensions_full_array_.data(), local_dimensions_subarray_.data(),
                             local_start_coordinates_.data(), ARRAY_LAYOUT, mpi_native_type_, &memory_space_type_);
    MPI_Type_commit(&memory_space_type_);


    // for writing ASCII
    // calculating chars_per_num based on format specified
    char s[100];
    sprintf(s, format_.c_str(), R(0));
    chars_per_num_ = strlen(s);

    // create chars_per_num_type_
    MPI_Type_contiguous(chars_per_num_, MPI_CHAR, &chars_per_num_type_);
    MPI_Type_commit(&chars_per_num_type_);

    // create ascii_file_space_type_
    MPI_Type_create_subarray(num_dims_, dimensions_full_array_.data(), dimensions_subarray_.data(),
                             start_coordinates_.data(), ARRAY_LAYOUT, chars_per_num_type_,
                             &ascii_file_space_type_);
    MPI_Type_commit(&ascii_file_space_type_);

    // create ascii_memory_space_type_
    MPI_Type_create_subarray(num_dims_, local_dimensions_full_array_.data(), local_dimensions_subarray_.data(),
                             local_start_coordinates_.data(), ARRAY_LAYOUT, chars_per_num_type_, &ascii_memory_space_type_);
    MPI_Type_commit(&ascii_memory_space_type_);
}

template <typename R, int ARRAY_LAYOUT>
Manager_MPI_IO<R,ARRAY_LAYOUT>::~Manager_MPI_IO()
{
    MPI_Type_free(&file_space_type_);
    MPI_Type_free(&memory_space_type_);
    MPI_Type_free(&chars_per_num_type_);
    MPI_Type_free(&ascii_file_space_type_);
    MPI_Type_free(&ascii_memory_space_type_);
}

template <typename R, int ARRAY_LAYOUT>
MPI_File Manager_MPI_IO<R,ARRAY_LAYOUT>::create_mpi_io_file_(const char *filename) const
{
    int err;
    int file_mode = MPI_MODE_UNIQUE_OPEN | MPI_MODE_WRONLY | MPI_MODE_CREATE;

    MPI_Info mpi_info = MPI_INFO_NULL; // For MPI IO hints
//#if 0
    MPI_Info_create(&mpi_info);
    MPI_Info_set(mpi_info, "collective_buffering", "true");
    MPI_Info_set(mpi_info, "striping_factor", "8");
    MPI_Info_set(mpi_info, "striping_unit", "4194304");
//#endif

    MPI_File file_handle = NULL;
    err = MPI_File_open(mpi_io_comm_, filename, file_mode, mpi_info, &file_handle);
    if (err)
    {
        std::string error_msg = std::string("MPI_File_open failed to open ") + std::string(filename);
        throw std::runtime_error(error_msg);
    }
    return file_handle;
}

template <typename R, int ARRAY_LAYOUT>
MPI_File Manager_MPI_IO<R,ARRAY_LAYOUT>::open_mpi_io_file_(const char *filename) const
{
    int err;
    int file_mode = MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN;

    MPI_Info mpi_info = MPI_INFO_NULL; // For MPI IO hints
    MPI_Info_create(&mpi_info);
    MPI_Info_set(mpi_info, "collective_buffering", "true");

    MPI_File file_handle = NULL;
    err = MPI_File_open(mpi_io_comm_, filename, file_mode, mpi_info, &file_handle);
    if (err)
    {
        std::string error_msg = std::string("MPI_File_open failed to open ") + std::string(filename);
        throw std::runtime_error(error_msg);
    }
    return file_handle;
}

template <typename R, int ARRAY_LAYOUT>
void Manager_MPI_IO<R,ARRAY_LAYOUT>::write_binary(const char *filename, R *data, int file_offset) const
{
    MPI_File file_handle = create_mpi_io_file_(filename);

    MPI_File_set_view(file_handle, file_offset, mpi_native_type_, file_space_type_, "native", MPI_INFO_NULL);
    MPI_File_write_all(file_handle, data, 1, memory_space_type_, MPI_STATUS_IGNORE);
  
    MPI_File_close(&file_handle);
}

template <typename R, int ARRAY_LAYOUT>
void Manager_MPI_IO<R,ARRAY_LAYOUT>::read_binary(const char *filename, R *data, int file_offset) const
{
    // open file
    MPI_File file_handle = open_mpi_io_file_(filename);

    // set view and read file
    MPI_File_set_view(file_handle, file_offset, mpi_native_type_, file_space_type_, "native", MPI_INFO_NULL);
    MPI_File_read_all(file_handle, data, 1, memory_space_type_, MPI_STATUS_IGNORE);
    MPI_File_close(&file_handle);
}

template <typename R, int ARRAY_LAYOUT>
void Manager_MPI_IO<R,ARRAY_LAYOUT>::write_ascii(const char *filename, R *data, const char* header_text) const
{
    // for holding data converted to chars
    const int size = local_dimensions_full_array_[0]*local_dimensions_full_array_[1]*local_dimensions_full_array_[2];
    char *data_as_chars = new char[size*chars_per_num_];

    // write data into data_as_chars
    for (int i = 0; i < size; i++)
    {
        sprintf(&data_as_chars[i*chars_per_num_], format_.c_str(), (double)data[i]);
    }

    int header_text_size = strlen(header_text);
    int data_as_chars_size = strlen(data_as_chars);

    // open file
    MPI_File file_handle = create_mpi_io_file_(filename);

    // my_rank_ == 0 writes header of file
    if (my_rank_ == root && header_text_size > 0)
    {
        MPI_File_write(file_handle, header_text, header_text_size, MPI_CHAR, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(mpi_io_comm_);

    // set view and write data
    MPI_File_set_view(file_handle, header_text_size, MPI_CHAR, ascii_file_space_type_, "native", MPI_INFO_NULL);
    MPI_File_write_all(file_handle, data_as_chars, 1, ascii_memory_space_type_, MPI_STATUS_IGNORE);

    MPI_File_close(&file_handle);

    delete[] data_as_chars;
}

template <typename R, int ARRAY_LAYOUT>
void Manager_MPI_IO<R,ARRAY_LAYOUT>::write_ascii(const char *filename, const char *data_as_chars, const char* header_text, int chars_per_line) const
{
    /*
      - header_text can be set to "" if no header_text is needed.
      - All lines in data_as_chars must have the same chars_per_line.
      - Do not forget to count the \n char when counting chars_per_line.
    */

    // create chars_per_line_type
    MPI_Datatype chars_per_line_type;
    MPI_Type_contiguous(chars_per_line, MPI_CHAR, &chars_per_line_type);
    MPI_Type_commit(&chars_per_line_type);

    // create custom_file_space_type
    MPI_Datatype custom_file_space_type;
    MPI_Type_create_subarray(num_dims_, dimensions_full_array_.data(), dimensions_subarray_.data(),
                             start_coordinates_.data(), ARRAY_LAYOUT, chars_per_line_type,
                             &custom_file_space_type);
    MPI_Type_commit(&custom_file_space_type);

    int header_text_size = strlen(header_text);
    int data_as_chars_size = strlen(data_as_chars);

    // open file
    MPI_File file_handle = create_mpi_io_file_(filename);

    // my_rank_ == 0 writes header of file
    if (my_rank_ == root && header_text_size > 0)
    {
        MPI_File_write(file_handle, header_text, header_text_size, MPI_CHAR, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(mpi_io_comm_);

    // set view and write data
    MPI_File_set_view(file_handle, header_text_size, MPI_CHAR, custom_file_space_type, "native", MPI_INFO_NULL);
    MPI_File_write_all(file_handle, data_as_chars, data_as_chars_size, MPI_CHAR, MPI_STATUS_IGNORE);

    MPI_File_close(&file_handle);

    // free created MPI_Datatype
    MPI_Type_free(&chars_per_line_type);
    MPI_Type_free(&custom_file_space_type);
}
