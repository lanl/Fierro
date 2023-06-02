#pragma once
#include "definitions.h"
#include <fstream>

using namespace utils;

// FOR DEBUGING **************************
void check_that_file_is_open(const std::ifstream & filestream, const char* filename);
void check_that_file_is_open(FILE *file_ptr, const char* filename);
void print_array(real_t *ptr, int N);
void print_array(int *ptr, int N);
void print_array_to_file(real_t *ptr, int N, int my_rank, const char* filename);
void print_array_to_file(int *ptr, int N, int my_rank, const char* filename);


// minval
template <typename T>
T minval(T *ptr, int N);

// maxval
template <typename T>
T maxval(T *ptr, int N);

// sum_array
template <typename T>
T sum_array(T *ptr, int N);

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
