
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <Kokkos_Core.hpp>

#include "matar.h"
#include "driver.h"


int main(int argc, char* argv[])
{

    // check to see of a mesh was supplied when running the code
    if (argc == 1)
    {
        std::cout << "\n\n**********************************\n\n";
        std::cout << " ERROR:\n";
        std::cout << " Please supply a mesh \n";
        std::cout << "   ./fierro my_mesh.geo \n\n";
        std::cout << "**********************************\n\n" << std::endl;
        std::exit(EXIT_FAILURE);
    } // end if

    int num_solvers = 1;


    Kokkos::initialize();

    Driver driver(argv[1]);

    driver.initialize(num_solvers);
    driver.setup();
    driver.run();
    driver.finalize();

    Kokkos::finalize();

    std::cout << "**** End of main **** " << std::endl;
    return 0;
}

