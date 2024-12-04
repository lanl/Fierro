#include <stdio.h>
#include <iostream>
#include "utilities.h"

int get_mpi_comm_rank(const MPI_Comm mpi_comm)
{
    int rank;
    MPI_Comm_rank(mpi_comm, &rank);
    return rank;
}

int get_mpi_comm_size(const MPI_Comm mpi_comm)
{
    int size;
    MPI_Comm_size(mpi_comm, &size);
    return size;
}

// FOR DEBUGING **************************
void check_that_file_is_open(const std::ifstream & filestream, const char* filename)
{
  if (!filestream)
  {
    throw std::runtime_error(std::string("Can't open ") + filename);
  }
}

void check_that_file_is_open(FILE *file_ptr, const char* filename)
{
  if (filename == NULL)
  {
    throw std::runtime_error(std::string("Can't open ") + filename);
  }
}


void print_array_to_file(real_t *ptr, int N, int my_rank, const char* filename)
{
  std::string filename_per_rank = "rank";
  filename_per_rank += std::to_string(my_rank);
  filename_per_rank += "_";
  filename_per_rank += filename;

  FILE *my_file = NULL;
  my_file = fopen (filename_per_rank.c_str(), "w");

  for (int i = 0; i < N; i++) {
    fprintf(my_file, "%24.14E\n", ptr[i]);
  }

  fclose(my_file);
}



void print_array_to_file(int *ptr, int N, int my_rank, const char* filename)
{
  std::string filename_per_rank = "rank";
  filename_per_rank += std::to_string(my_rank);
  filename_per_rank += "_";
  filename_per_rank += filename;

  FILE *my_file = NULL;
  my_file = fopen (filename_per_rank.c_str(), "w");

  for (int i = 0; i < N; i++) {
    fprintf(my_file, "%13d\n", ptr[i]);
  }

  fclose(my_file);
}
