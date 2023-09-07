#pragma once

struct SimParameters
{
    int     nn[3];
    int     ndim;
    int     num_steps; 
    int     print_rate;
    int     iseed;
    double  dx;
    double  delta[3];
    double  dt;
    double  kappa;
    double  M;
    double  c0; 
    double  noise;

    SimParameters();
    void set_ndim();
    void print();
};