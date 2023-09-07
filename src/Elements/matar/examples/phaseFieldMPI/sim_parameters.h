#pragma once
#include <array>


struct SimParameters
{
    std::array<int,3> nn;
    int ndim;
    int num_steps; 
    int print_rate;
    int iseed;
    double dx;
    std::array<double,3> delta;
    double dt;
    double kappa;
    double M;
    double c0; 
    double noise;

    SimParameters();
    void print() const;

private:
    void set_ndim();
};
