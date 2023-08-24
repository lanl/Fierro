#include <stdio.h>
#include <math.h>
#include <chrono>
#include <matar.h>

using namespace mtr; // matar namespace

const int width = 1000;
const int height = 1000;
const double temp_tolerance = 0.01;

void initialize(CArray<double> &temperature_previous);
void track_progress(int iteration, CArray<double> &temperature);

int main() {
    int i, j;
    int iteration = 1;
    double worst_dt = 100;

    auto temperature = CArray <double> (height+2, width+2);
    auto temperature_previous = CArray <double> (height+2, width+2);

    // Start measuring time
    auto begin = std::chrono::high_resolution_clock::now();

    // initialize temperature profile
    initialize(temperature_previous);

    while (worst_dt > temp_tolerance) {
        // finite difference
        for (i = 1; i <= height; i++) {
            for (j = 1; j <= width; j++) {
                temperature(i,j) = 0.25 * (temperature_previous(i+1,j)
                                        + temperature_previous(i-1,j)
                                        + temperature_previous(i,j+1)
                                        + temperature_previous(i,j-1));         
            }
        }
        
        // calculate max difference between temperature and temperature_previous
        worst_dt = 0.0;
        for (i = 1; i <= height; i++) {
            for (j = 1; j <= width; j++) {
                worst_dt = fmax(fabs(temperature(i,j) - 
                                temperature_previous(i,j)), 
                                worst_dt);
            }
        }

        // update temperature_previous
        for (i = 1; i <= height; i++) {
            for (j = 1; j <= width; j++) {
                temperature_previous(i,j) = temperature(i,j);
            }
        }

        // track progress
        if (iteration % 100 == 0) {
            track_progress(iteration, temperature);
        }

        iteration++;
    }

    // Stop measuring time and calculate the elapsed time
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

    printf("Total time was %f seconds.\n", elapsed.count() * 1e-9);
    printf("\nMax error at iteration %d was %f\n", iteration-1, worst_dt);      

    return 0;
}

void initialize(CArray<double> &temperature_previous) {
    int i, j;

    // initialize temperature_previous to 0.0
    for (i = 0; i <= height+1; i++) {
        for (j = 0; j <= width+1; j++) {
            temperature_previous(i,j) = 0.0;            
        }
    }

    // setting the left and right boundary conditions
    for (i = 0; i <= height+1; i++) {
        temperature_previous(i,0) = 0.0;
        temperature_previous(i,width+1) = (100.0/height)*i;
    }

    // setting the top and bottom boundary condition
    for (j = 0; j <= width+1; j++) {
        temperature_previous(0,j) = 0.0;
        temperature_previous(height+1,j) = (100.0/width)*j; 
    }
}

void track_progress(int iteration, CArray<double> &temperature) {
    int i;

    printf("---------- Iteration number: %d ----------\n", iteration);
    for (i = height-5; i <= height; i++) {
        printf("[%d,%d]: %5.2f  ", i,i, temperature(i,i));
    }
    printf("\n");
}
