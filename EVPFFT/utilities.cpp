#include <stdio.h>
#include <iostream>
#include "utilities.h"


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

void print_array(real_t *ptr, int N)
{
  for (int i = 0; i < N; i++) {
    printf("%24.14E", ptr[i]);
    if ((i+1)%3 == 0) printf("\n");
  }
  printf("\n");

}

void print_array(int *ptr, int N)
{
  for (int i = 0; i < N; i++) {
    printf("%13d", ptr[i]);
    if ((i+1)%6 == 0) printf("\n");
  }
  printf("\n");

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
