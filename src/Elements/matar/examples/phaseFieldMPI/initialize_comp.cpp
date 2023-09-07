#include <iostream>
#include "initialize_comp.h"


void initialize_comp(const SimParameters &sp, DCArrayKokkos<double> &comp, CArray<double> &comp_all)
{
if (0 == rank)
{
    // seed random number generator
    srand(sp.iseed);

    // to hold random number
    double r;

    for (int i = 0; i < sp.nn[0]; ++i) {
        for (int j = 0; j < sp.nn[1]; ++j) {
            for (int k = 0; k < sp.nn[2]; ++k) {
                // random number between 0.0 and 1.0
                r = (double) rand()/RAND_MAX;

                // initialize "comp" with stochastic thermal fluctuations
                comp_all(i,j,k) = sp.c0 + (2.0*r - 1.0)*sp.noise;
            }
        }
    }
}


}
