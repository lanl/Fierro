#include <stdio.h>
#include <iostream>

#include "matar.h"
#include "Kokkos_DualView.hpp"

using namespace mtr; // matar namespace

void DViewCArrayKokkosTwoDimensionExample();
void DCArrayKokkosTwoDimensionExample();

int main() {

    Kokkos::initialize();
    { 

    // Run DViewCArrayKokkos 2D example
    DViewCArrayKokkosTwoDimensionExample();

    // Run DCArrayKokkos 2D example
    DCArrayKokkosTwoDimensionExample();

    } // end of kokkos scope
    Kokkos::finalize();
}


void DViewCArrayKokkosTwoDimensionExample()
{
    printf("\n====================Running 2D DViewCArrayKokkos example====================\n");

    int nx = 2;
    int ny = 2;

    // CPU arr
    int arr[nx*ny];
    
    for (int i = 0; i < nx*ny; i++){
        arr[i] = 1;
    }

    // Create A_2D
    auto A_2D = DViewCArrayKokkos <int> (&arr[0], nx, ny);
    
    // Print host copy of data
    printf("Printing host copy of data (should be all 1s):\n");
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            printf("%d\n", A_2D.host(i,j));
        }
    }
    
    // Print device copy of data
    printf("Printing device copy of data (should be all 1s):\n");
    FOR_ALL(i, 0, nx,
            j, 0, ny, {
        printf("%d\n", A_2D(i,j));
    });
    Kokkos::fence();

    // Manupulate data on device and update host
    FOR_ALL(i, 0, nx,
            j, 0, ny, {
        A_2D(i,j) = 2;
    });
    A_2D.update_host();
    Kokkos::fence();
    printf("---Data updated to 2 on device---\n");

    // Print host copy of data
    printf("Printing host copy of data (should be all 2s):\n");
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            printf("%d\n", A_2D.host(i,j));
        }
    }

    // Print device copy of data
    printf("Printing device copy of data (should be all 2s):\n");
    FOR_ALL(i, 0, nx,
            j, 0, ny, {
        printf("%d\n", A_2D(i,j));
    });
    Kokkos::fence();

    // Print pointer to data on host and device
    printf("\nPrinting pointer to data on host and device.\n");
    printf("Should be same address if using OpenMP backend.\n");
    printf("Should be different addresses if using GPU backend.\n");
    printf("Host data pointer: %p\n", A_2D.host_pointer());
    printf("Device data pointer: %p\n", A_2D.device_pointer());
}



void DCArrayKokkosTwoDimensionExample()
{
    printf("\n====================Running 2D DCArrayKokkos example====================\n");

    int nx = 2;
    int ny = 2;

    // Create A_2D
    auto A_2D = DCArrayKokkos <int> (nx, ny);

    // Set data to one on host and updata device
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            A_2D.host(i,j) = 1;
        }
    }
    A_2D.update_device();
    Kokkos::fence();

    // Print host copy of data
    printf("Printing host copy of data (should be all 1s):\n");
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            printf("%d\n", A_2D.host(i,j));
        }
    }
    
    // Print device copy of data
    printf("Printing device copy of data (should be all 1s):\n");
    FOR_ALL(i, 0, nx,
            j, 0, ny, {
        printf("%d\n", A_2D(i,j));
    });
    Kokkos::fence();

    // Manupulate data on device and update host
    FOR_ALL(i, 0, nx,
            j, 0, ny, {
        A_2D(i,j) = 2;
    });
    A_2D.update_host();
    Kokkos::fence();
    printf("---Data updated to 2 on device---\n");

    // Print host copy of data
    printf("Printing host copy of data (should be all 2s):\n");
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            printf("%d\n", A_2D.host(i,j));
        }
    }

    // Print device copy of data
    printf("Printing device copy of data (should be all 2s):\n");
    FOR_ALL(i, 0, nx,
            j, 0, ny, {
        printf("%d\n", A_2D(i,j));
    });
    Kokkos::fence();

    // Print pointer to data on host and device
    printf("\nPrinting pointer to data on host and device.\n");
    printf("Should be same address if using OpenMP backend.\n");
    printf("Should be different addresses if using GPU backend.\n");
    printf("Host data pointer: %p\n", A_2D.host_pointer());
    printf("Device data pointer: %p\n", A_2D.device_pointer());

}
