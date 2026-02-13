#pragma once
#include "definitions.h"
#include <fstream>
#include <iomanip>
#include "mpi.h"

using namespace utils;

// MPI utility functions
int get_mpi_comm_rank(const MPI_Comm mpi_comm);
int get_mpi_comm_size(const MPI_Comm mpi_comm);

// FOR DEBUGING **************************
void check_that_file_is_open(const std::ifstream & filestream, const char* filename);
void check_that_file_is_open(FILE *file_ptr, const char* filename);

// minval
template <typename T>
T minval(T *ptr, int N);

// maxval
template <typename T>
T maxval(T *ptr, int N);

// sum_array
template <typename T>
T sum_array(T *ptr, int N);

template <typename T>
void print_array(T* ptr, size_t dim0, size_t dim1);

template <typename T>
void print_array_to_file(T *ptr, size_t dim0, size_t dim1, size_t my_rank, const char* filename);

// minval
template <typename T>
T minval(T *ptr, int N)
{
  T minval_ = ptr[0];
  for (int i = 0; i < N; i++) {
    if (ptr[i] < minval_) {
      minval_ = ptr[i];
    }
  }
  return minval_;
}

// maxval
template <typename T>
T maxval(T *ptr, int N)
{
  T maxval_ = ptr[0];
  for (int i = 0; i < N; i++) {
    if (ptr[i] > maxval_) {
      maxval_ = ptr[i];
    }
  }
  return maxval_;
}

// sum_array
template <typename T>
T sum_array(T *ptr, int N)
{
  T sum_ = 0.0;
  for (int i = 0; i < N; i++) {
    sum_ += ptr[i];
  }
  return sum_;
}


template <typename T>
void print_array(T* ptr, size_t dim0, size_t dim1)
{
  int p = 4;
  std::cout << std::setprecision(p) << std::scientific;
  for (size_t i = 0; i < dim0; i++) {
    for (size_t j = 0; j < dim1; j++) {
      std::cout << "  " << ptr[j+(i*dim1)];
    }
    std::cout << std::endl;
  }
}

template <typename T>
void print_array_to_file(T *ptr, size_t dim0, size_t dim1, size_t my_rank, const char* filename)
{
  std::string filename_per_rank = "rank";
  filename_per_rank += std::to_string(my_rank);
  filename_per_rank += "_";
  filename_per_rank += filename;

  std::ofstream myfile;
  myfile.open(filename_per_rank);

  int p = 4;
  myfile << std::setprecision(p) << std::scientific;
  for (size_t i = 0; i < dim0; i++) {
    for (size_t j = 0; j < dim1; j++) {
      myfile << "  " << ptr[j+(i*dim1)];
    }
    myfile << std::endl;
  }

  myfile.close();
}

