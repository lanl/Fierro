#include <iostream>
#include <stdio.h>
#include "complex_arrays.h"
#include "mpi.h"

ComplexArrays::ComplexArrays(const SimParameters & sp, const std::array<int,3> & loc_nn_img, const std::array<int,3> & loc_start_index) :
comp_img(loc_nn_img[2], loc_nn_img[1], loc_nn_img[0], 2),
dfdc_img(loc_nn_img[2], loc_nn_img[1], loc_nn_img[0], 2),
kpow2(loc_nn_img[2], loc_nn_img[1], loc_nn_img[0]),
denominator(loc_nn_img[2], loc_nn_img[1], loc_nn_img[0]),
fs(sp.nn, loc_nn_img, loc_start_index, sp.delta)
{
    // set values of kpow2
    set_kpow2();

    // set values of denominator
    set_denominator(sp);
}


void ComplexArrays::set_kpow2()
{   
    // calculate kpow2
    FOR_ALL_CLASS(k, 0, kpow2.dims(0),
                  j, 0, kpow2.dims(1),
                  i, 0, kpow2.dims(2), {
        kpow2(k,j,i) =   fs.kx(i) * fs.kx(i)
                       + fs.ky(j) * fs.ky(j)
                       + fs.kz(k) * fs.kz(k);
    });
}


void ComplexArrays::set_denominator(const SimParameters & sp)
{
    double dt = sp.dt;
    double M = sp.M;
    double kappa = sp.kappa;

    // calculate denominator_
    FOR_ALL_CLASS(k, 0, denominator.dims(0),
                  j, 0, denominator.dims(1),
                  i, 0, denominator.dims(2), {
        denominator(k,j,i) = 1.0 + (dt * M * kappa * kpow2(k,j,i) * kpow2(k,j,i));
    });
}
