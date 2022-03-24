
#include "Solver.h"
#include <stdlib.h> 
#include <mpi.h>

Solver::Solver(){
  //default flags assume optional routines are off
  setup_flag = finalize_flag = 0;
}

void Solver::exit_solver(int status){
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  exit(status);
}

Solver::~Solver(){}
