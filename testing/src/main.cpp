
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <Kokkos_Core.hpp>
#include <sys/stat.h>


#include "matar.h"
#include "driver.h"


int main(int argc, char* argv[])
{
    // check to see of an input file was supplied when running the code
    if (argc == 1) {
        std::cout << "\n\n**********************************\n\n";
        std::cout << " ERROR:\n";
        std::cout << " Please supply a YAML input, \n";
        std::cout << "   ./mesh-builder input.yaml \n\n";
        std::cout << "**********************************\n\n" << std::endl;
        return 0;
    } // end if

    Kokkos::initialize();

    // Create driver on heap
    Driver* driver;
    driver = new Driver(argv[1]);
    

    driver->initialize();
    driver->setup();
    driver->run();
    driver->finalize();

    // Delete driver
    delete driver;

    Kokkos::finalize();

    std::cout << "**** End of main **** " << std::endl;
    return 0;
}

