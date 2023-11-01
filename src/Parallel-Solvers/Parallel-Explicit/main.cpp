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
#include <sys/stat.h>
#include <mpi.h>
#include <Kokkos_Core.hpp>
#include "Explicit_Solver.h"
#include "Simulation_Parameters/Simulation_Parameters.h"
#include "yaml-serializable.h"
#include <memory>
#include <string>

void solver_setup(const char* filename){
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

  std::shared_ptr<Simulation_Parameters> simparam;
  Yaml::from_file_strict(filename, simparam);

  if (simparam->type != SOLVER_TYPE::Explicit) {
      if (myrank == 0)
        std::cerr << "Invalid solver type" << std::endl;
      exit(1);
  }

  Explicit_Solver solver(*std::dynamic_pointer_cast<Simulation_Parameters_Explicit>(simparam));
  //checks for optional solver routines
  if(solver.setup_flag) solver.solver_setup();
  // invoke solver's run function (should perform most of the computation)
  solver.run();
  //invoke optional finalize function
  if(solver.finalize_flag) solver.solver_finalize();
}

int main(int argc, char *argv[]){
  //initialize MPI
  MPI_Init(&argc, &argv);

  Kokkos::initialize();
  if(argc<2){
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    if(myrank==0)
      std::cout << "Fierro requires a yaml setup file as a command line argument" << std::endl;
  }
  else{
    solver_setup(argv[1]);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  Kokkos::finalize();
  MPI_Finalize();

  return 0;
}
