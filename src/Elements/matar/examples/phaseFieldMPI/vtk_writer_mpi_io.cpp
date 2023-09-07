#include <vtk_writer_mpi_io.h>

VTK_Writer_MPI_IO::VTK_Writer_MPI_IO(MPI_Comm mpi_io_comm, const std::array<int,3> & dimensions_full_array,
    const std::array<int,3> & dimensions_subarray,
    const std::array<int,3> & start_coordinates,
    const char *format) :
mpi_io_comm_(mpi_io_comm),
dimensions_full_array_(dimensions_full_array),
dimensions_subarray_(dimensions_subarray),
start_coordinates_(start_coordinates),
format_(format),
chars_per_num_type_(MPI_DATATYPE_NULL),
file_space_type_(MPI_DATATYPE_NULL)
{
    // calculating chars_per_num based on format specified
    char s[100];
    sprintf(s, format_.c_str(), double(0));
    chars_per_num_ = strlen(s);

    // create chars_per_num_type_
    MPI_Type_contiguous(chars_per_num_, MPI_CHAR, &chars_per_num_type_);
    MPI_Type_commit(&chars_per_num_type_);

    // create file_space_type_
    int dimensions_full_array_reordered[3] = { dimensions_full_array_[2], dimensions_full_array_[1], dimensions_full_array_[0] };
    int dimensions_subarray_reordered[3] = { dimensions_subarray_[2], dimensions_subarray_[1], dimensions_subarray_[0] };
    int start_coordinates_reordered[3] = { start_coordinates_[2], start_coordinates_[1], start_coordinates_[0] };
    MPI_Type_create_subarray(3, dimensions_full_array_reordered, dimensions_subarray_reordered,
                             start_coordinates_reordered, MPI_ORDER_C, chars_per_num_type_,
                             &file_space_type_);
    MPI_Type_commit(&file_space_type_);
}

VTK_Writer_MPI_IO::~VTK_Writer_MPI_IO()
{
    MPI_Type_free(&chars_per_num_type_);
    MPI_Type_free(&file_space_type_);
}

void VTK_Writer_MPI_IO::write(int iter, const double *data)
{
    // global array dimensions
    int nx = dimensions_full_array_[0];
    int ny = dimensions_full_array_[1];
    int nz = dimensions_full_array_[2];

    // create name of output vtk file
    char filename[50];
    sprintf(filename, "outputComp_%d.vtk", iter);
    
    // for storing header_text
    std::string header_text;
    
    // write vtk file heading
    char buff[300];
    sprintf(buff, "%s\n", "# vtk DataFile Version 3.0");
    header_text += buff;
    sprintf(buff, "%s\n", filename);
    header_text += buff;
    sprintf(buff, "%s\n", "ASCII");
    header_text += buff;
    sprintf(buff, "%s\n", "DATASET STRUCTURED_POINTS");
    header_text += buff;
    sprintf(buff, "%s %d  %d  %d\n", "DIMENSIONS", nx, ny, nz);
    header_text += buff;
    sprintf(buff, "%s %d %d %d\n", "ORIGIN", 0, 0, 0);
    header_text += buff;
    sprintf(buff, "%s %d %d %d\n", "SPACING", 1, 1, 1);
    header_text += buff;
    sprintf(buff, "%s %d\n", "POINT_DATA", nx*ny*nz);
    header_text += buff;
    sprintf(buff, "%s\n", "SCALARS data double");
    header_text += buff;
    sprintf(buff, "%s\n", "LOOKUP_TABLE default");
    header_text += buff;

    // for holding data converted to chars
    const int subarray_size = dimensions_subarray_[0]*dimensions_subarray_[1]*dimensions_subarray_[2];
    char *data_as_chars = new char[subarray_size*chars_per_num_];

    // write data into data_as_chars
    for (int i = 0; i < subarray_size; i++)
    {
        sprintf(&data_as_chars[i*chars_per_num_], format_.c_str(), data[i]);
    }

    write_mpi_io_file(filename, header_text.c_str(), data_as_chars);

    delete[] data_as_chars;
}

MPI_File VTK_Writer_MPI_IO::create_mpi_io_file(const char *filename)
{
    int file_mode = MPI_MODE_UNIQUE_OPEN | MPI_MODE_WRONLY | MPI_MODE_CREATE;

    MPI_Info mpi_info = MPI_INFO_NULL; // For MPI IO hints
#if 0
    MPI_Info_create(&mpi_info);
    MPI_Info_set(mpi_info, "collective_buffering", "true");
    MPI_Info_set(mpi_info, "striping_factor", "8");
    MPI_Info_set(mpi_info, "striping_unit", "4194304");
#endif

    MPI_File file_handle = NULL;
    MPI_File_open(mpi_io_comm_, filename, file_mode, mpi_info, &file_handle);

    return file_handle;
}

void VTK_Writer_MPI_IO::write_mpi_io_file(const char *filename, const char* header_text, const char *data_as_chars)
{
    int my_rank;
    MPI_Comm_rank(mpi_io_comm_, &my_rank);

    int header_text_size = strlen(header_text);
    int data_as_chars_size = strlen(data_as_chars);

    // open file
    MPI_File file_handle = create_mpi_io_file(filename);

    // my_rank == 0 writes header of file
    if (my_rank == 0 && header_text_size > 0)
    {
        MPI_File_write(file_handle, header_text, header_text_size, MPI_CHAR, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(mpi_io_comm_);

    // set view and write data
    MPI_File_set_view(file_handle, header_text_size, MPI_CHAR, file_space_type_, "native", MPI_INFO_NULL);
    MPI_File_write_all(file_handle, data_as_chars, data_as_chars_size, MPI_CHAR, MPI_STATUS_IGNORE);

    MPI_File_close(&file_handle);
}
