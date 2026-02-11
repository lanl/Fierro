#pragma once

#include <cstdio>
#include <array>
#include <string>
#include <stdexcept>


class VTKUniformGridWriter
{
public:
  size_t N1;
  size_t N2;
  size_t N3;
  FILE *vtk_file;

  VTKUniformGridWriter(size_t N1_, size_t N2_, size_t N3_);
  void open_file(const char* vtk_filename);
  void write_header();
  template <typename T>
  void write_scalar_dataset(T* data_ptr, const char* dataName, const char* dataType);
  void close_file();
};


VTKUniformGridWriter::VTKUniformGridWriter(size_t N1_, size_t N2_, size_t N3_)
  : N1 (N1_)
  , N2 (N2_)
  , N3 (N3_)
  , vtk_file (nullptr)
{
}

void VTKUniformGridWriter::open_file(const char* vtk_filename)
{
  vtk_file = fopen(vtk_filename, "w");
}

void VTKUniformGridWriter::write_header()
{
  // write vtk file heading
  fprintf(vtk_file, "%s\n", "# vtk DataFile Version 3.0");
  fprintf(vtk_file, "%s\n", "Simulation outputs");
  fprintf(vtk_file, "%s\n", "ASCII");
  fprintf(vtk_file, "%s\n", "DATASET STRUCTURED_POINTS");

  fprintf(vtk_file, "%s %d  %d  %d\n", "DIMENSIONS", N1, N2, N3);
  fprintf(vtk_file, "%s %d %d %d\n", "ORIGIN", 0, 0, 0);
  fprintf(vtk_file, "%s %d %d %d\n", "SPACING", 1, 1, 1);
  fprintf(vtk_file, "%s %d\n", "POINT_DATA", N1*N2*N3);
}

template <typename T>
void VTKUniformGridWriter::write_scalar_dataset(T* data_ptr, const char* dataName, const char* dataType)
{
  fprintf(vtk_file, "%s %s %s\n", "SCALARS", dataName, dataType);
  fprintf(vtk_file, "%s\n", "LOOKUP_TABLE default");

  for (int i = 0; i < N1*N2*N3; i++) {
    fprintf(vtk_file, " %12.6E\n", data_ptr[i]);
  }
}

void VTKUniformGridWriter::close_file()
{
  if (vtk_file != nullptr) {
    fclose(vtk_file);
  }
}
