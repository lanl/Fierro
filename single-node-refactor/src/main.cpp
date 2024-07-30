/**********************************************************************************************
 © 2020. Triad National Security, LLC. All rights reserved.
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
#include <Kokkos_Core.hpp>
#include <sys/stat.h>

#include "matar.h"
#include "driver.h"

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
        std::cout << "   ./mesh-builder input.yaml \n\n";
        std::cout << "**********************************\n\n" << std::endl;
        return 0;
    } // end if

    //Kokkos::initialize(argc, argv);
    MATAR_MPI_INIT
    MATAR_KOKKOS_INIT
    { // kokkos scope

        // Create driver
        Driver* driver = new Driver(argv[1]);

        driver->initialize();
        driver->setup();
        driver->run();
        driver->finalize();

        // Delete driver
        delete driver;
    } // end kokkos scope
    MATAR_KOKKOS_FINALIZE
    MATAR_MPI_FINALIZE
    //Kokkos::finalize();

    std::cout << "**** End of main **** " << std::endl;
    return 0;
}
