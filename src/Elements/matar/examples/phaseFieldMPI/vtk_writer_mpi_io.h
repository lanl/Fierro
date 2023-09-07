#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <array>
#include <string>
#include <iostream>


class VTK_Writer_MPI_IO
{
private:
    MPI_Comm mpi_io_comm_;
    const std::array<int,3> dimensions_full_array_;
    const std::array<int,3> dimensions_subarray_;
    const std::array<int,3> start_coordinates_;
    std::string format_;
    int chars_per_num_;
    MPI_Datatype chars_per_num_type_;
    MPI_Datatype file_space_type_;

    MPI_File create_mpi_io_file(const char *filename);
    void write_mpi_io_file(const char *filename, const char* header_text, 
                           const char *data_as_chars);

public:
    VTK_Writer_MPI_IO(MPI_Comm mpi_io_comm, const std::array<int,3> & dimensions_full_array,
                      const std::array<int,3> & dimensions_subarray,
                      const std::array<int,3> & start_coordinates,
                      const char *format);
    ~VTK_Writer_MPI_IO();

    void write(int iter, const double *data);
};
