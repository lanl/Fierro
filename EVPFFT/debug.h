#pragma once

#include <stdio.h>

#if 0
#include <iostream>
#include <iomanip>


#define PAUSE std::cout << "\nC++ PAUSE: enter <return> or <ctrl>d to continue>\n"; std::cin.get();


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
#endif


void print_array(double* ptr, size_t dim0, size_t dim1)
{
  for (size_t i = 0; i < dim0; i++) {
    for (size_t j = 0; j < dim1; j++) {
      printf( "%24.14E", ptr[j+(i*dim1)] );
    }
    printf("\n");
  }

}
