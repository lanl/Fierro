#include <iostream>
#include <stdio.h>
#include <chrono>

#include "sim_parameters.h"
#include "global_arrays.h"
#include "initialize_comp.h"
#include "CH_fourier_spectral_solver.h"
#include "local_free_energy.h"
#include "outputs.h"


int main(int argc, char* argv[])
{
    Kokkos::initialize(argc, argv);
    {

    // simulation parameters
    SimParameters sp;
    sp.print();

    // global arrays needed for simulation
    GlobalArrays ga = GlobalArrays(sp.nn);

    // setup initial composition profile
    initialize_comp(sp, ga.comp);

    // initialize solver
    CHFourierSpectralSolver CH_fss(sp);

    // Start measuring time
    auto begin = std::chrono::high_resolution_clock::now();

    // time stepping loop
    for (int iter = 1; iter <= sp.num_steps; iter++) {
        // calculate df/dc
        calculate_dfdc(sp.nn, ga.comp, ga.dfdc);

        // Cahn Hilliard equation solver
        CH_fss.time_march(ga.comp, ga.dfdc);

        // report simulation progress and output vtk files
        if (iter % sp.print_rate == 0) {

            track_progress(iter, sp.nn, ga.comp);

            write_vtk(iter, sp.nn, sp.delta, ga.comp);

            output_total_free_energy(iter, sp.print_rate, sp.num_steps, 
                                     sp.nn, sp.delta, sp.kappa, 
                                     ga.comp);
        }
    }

    // Stop measuring time and calculate the elapsed time
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf("Total time was %f seconds.\n", elapsed.count() * 1e-9); 


    }
    Kokkos::finalize();

    return 0;
}
