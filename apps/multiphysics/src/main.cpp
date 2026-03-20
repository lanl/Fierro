/**********************************************************************************************
 � 2020. Triad National Security, LLC. All rights reserved.
 This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
 National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
 Department of Energy/National Nuclear Security Administration. All rights in the program are
 reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
 Security Administration. The Government is granted for itself and others acting on its behalf a
 nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
 derivative works, distribute copies to the public, perform publicly and display publicly, and
 to permit others to do so.
 This program is open source under the BSD-3 License.
 Redistribution and use in source and binary forms, with or without modification, are permitted
 provided that the following conditions are met:
 1.  Redistributions of source code must retain the above copyright notice, this list of
 conditions and the following disclaimer.
 2.  Redistributions in binary form must reproduce the above copyright notice, this list of
 conditions and the following disclaimer in the documentation and/or other materials
 provided with the distribution.
 3.  Neither the name of the copyright holder nor the names of its contributors may be used
 to endorse or promote products derived from this software without specific prior
 written permission.
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **********************************************************************************************/
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
// #include <Kokkos_Core.hpp>
#include <sys/stat.h>

#include "ELEMENTS.h"
#include "driver.hpp"
#include "parse_tools.hpp"


/////////////////////////////////////////////////////////////////////////////
///
/// \fn main
///
/// \brief Takes in a YAML input, creates a driver, and build simulation based
///        on the YAML input
///
/////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
    // check to see of an input file was supplied when running the code
    if (argc == 1) {
        std::cout << "\n\n**********************************\n\n";
        std::cout << " ERROR:\n";
        std::cout << " Please supply a YAML input, \n";
        std::cout << "   ./Fierro input.yaml \n\n";
        std::cout << "**********************************\n\n" << std::endl;
        return 0;
    } // end if

    if (std::string(argv[1]) == "--help"){

        print_inputs();

        return 0;
    }

    Kokkos::initialize();
    {

        // Create driver
        Driver* driver = new Driver(argv[1]);


        // Timing data for each step
        auto time_start = std::chrono::high_resolution_clock::now();
        auto time_init = std::chrono::high_resolution_clock::now();
        
        driver->initialize();
        
        auto time_now = std::chrono::high_resolution_clock::now();
        auto calc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(time_now - time_init).count();
        printf("\n**** Total time to initialize driver in seconds  %f ****\n\n", calc_time * 1e-9);

        auto time_setup = std::chrono::high_resolution_clock::now();
        driver->setup();
        time_now = std::chrono::high_resolution_clock::now();
        calc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(time_now - time_setup).count();
        printf("\n**** Total time to setup driver in seconds  %f ****\n\n", calc_time * 1e-9);
        

        auto time_run = std::chrono::high_resolution_clock::now();
        driver->execute();
        time_now = std::chrono::high_resolution_clock::now();
        calc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(time_now - time_setup).count();
        printf("\n**** Total time to execute driver in seconds  %f ****\n\n", calc_time * 1e-9);


        driver->finalize();

        time_now = std::chrono::high_resolution_clock::now();
        calc_time = std::chrono::duration_cast<std::chrono::nanoseconds>(time_now - time_start).count();

        printf("\n**** Total time to run simulation in seconds  %f ****\n\n", calc_time * 1e-9);

        // Delete driver
        delete driver;
    }

    Kokkos::finalize();

    std::cout << "**** End of main **** " << std::endl;
    return 0;
}
