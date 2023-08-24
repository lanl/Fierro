#include <iostream>
#include <fstream>
#include <limits>

#include "sim_parameters.h"


SimParameters::SimParameters()
{
    // set default simulation parameters
    this->nn[0]      = 32;        // nx
    this->nn[1]      = 32;        // ny
    this->nn[2]      = 32;        // nz
    this->delta[0]   = 1.0;       // dx
    this->delta[1]   = 1.0;       // dy
    this->delta[2]   = 1.0;       // dz
    this->dt         = 5.0E-2;    // dt
    this->num_steps  = 1000;      // total number of time steps
    this->print_rate = 100;       // time step interval for output file
    this->iseed      = 456;       // random number seed
    this->kappa      = 1.0;       // gradient energy coefficient
    this->M          = 1.0;       // mobility
    this->c0         = 5.0E-1;    // critical composition
    this->noise      = 5.0E-3;    // noise term for thermal fluctuations

    // set number of dimensions
    set_ndim();
}


void SimParameters::set_ndim()
{
    ndim = 0;
    for (int i = 0; i < 3; i++) {
        if (nn[i] > 1) ++ndim;
    }
}


void SimParameters::print() const
{
    std::cout << " nx = "         << nn[0]      << std::endl;
    std::cout << " ny = "         << nn[1]      << std::endl;
    std::cout << " nz = "         << nn[2]      << std::endl;
    std::cout << " dx = "         << delta[0]   << std::endl; 
    std::cout << " dy = "         << delta[1]   << std::endl; 
    std::cout << " dz = "         << delta[2]   << std::endl; 
    std::cout << " dt = "         << dt         << std::endl;
    std::cout << " num_steps = "  << num_steps  << std::endl; 
    std::cout << " print_rate = " << print_rate << std::endl; 
    std::cout << " iseed = "      << iseed      << std::endl; 
    std::cout << " kappa = "      << kappa      << std::endl; 
    std::cout << " M = "          << M          << std::endl; 
    std::cout << " c0 = "         << c0         << std::endl; 
    std::cout << " noise = "      << noise      << std::endl; 
}
