#include <iostream>
#include "initialize_comp.h"


void initialize_comp(SimParameters& sp, DCArrayKokkos<double> &comp)
{
    // unpack simimulation parameters needed 
    // for calculations in this function
    int nx       = sp.nn[0];
    int ny       = sp.nn[1];
    int nz       = sp.nn[2];
    int iseed    = sp.iseed;
    double c0    = sp.c0;
    double noise = sp.noise;

    // seed random number generator
    srand(iseed);

    // to hold random number
    double r;

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                // random number between 0.0 and 1.0
                r = (double) rand()/RAND_MAX;

                // initialize "comp" with stochastic thermal fluctuations
                comp.host(i,j,k) = c0 + (2.0*r - 1.0)*noise;
            }
        }
    }

    // update host copy of device
    comp.update_device();
}
