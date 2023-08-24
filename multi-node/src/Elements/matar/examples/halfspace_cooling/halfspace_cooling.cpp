// Populate a dynamic ragged down array with the 
// temperatures from the halfspace cooling model as a function of
// age and depth

// T = Tm erf(z/sqrt(4kt))

// The mantle temperature is 1350, temp will be 1350 at the ridge axis
// the x direction is controlled by increasing age
// the y direction is depth
// k = thermal diffusivity (1 x 10^-6)
// calculate every 1 million years for 1000 Ma years
// have depth increase by 1 km

#include <stdio.h>
#include <math.h>
#include <cmath>
#include <chrono>
#include <iostream>
#include <matar.h>

using namespace mtr; // matar namespace

// set up constant parameters
const int max_age = 1000;
const double mantle_temp = 1350.0;
const double thermal_diff = 0.000001;

int main() {

    Kokkos::initialize();
{
    // depth will need to be adjusted for larger max ages
    // age 2000 Ma, depth 250
    // age 3000 Ma, depth 280
    // age 4000 Ma, depth 320
    int depth = 200; 
    auto begin = std::chrono::high_resolution_clock::now(); // start clock

    DynamicRaggedDownArrayKokkos <double> dyn_ragged_down(max_age+1, depth+1); // create array

    DO_ALL(i, 0, max_age, {
            for (int j = 0; j <= depth; j++) {
                if (i == 0 && j == 0)
                { // when depth and age are 0, give mantle_temp
                    dyn_ragged_down.stride(j)++;
                    dyn_ragged_down(i, j) = mantle_temp;
                }
                double temp = mantle_temp * erf(j / (2.0 * sqrt(thermal_diff * (i * 1e6))));
                dyn_ragged_down.stride(j)++;
                dyn_ragged_down(i, j) = temp;

                // check if we have reached the mantle, if yes, move on to next age
                if (round(dyn_ragged_down(i, j)) == 1350)
                {
                    printf("Depth to mantle %d km, age of lithosphere %d Ma \n", j, i);
                    break;
                }
            }
        });

    // Stop counting time and calculate elapsed time
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin);
    
    printf("Total time was %f seconds.\n", elapsed.count() * 1e-9);
}
    Kokkos::finalize();

    return 0;
}
