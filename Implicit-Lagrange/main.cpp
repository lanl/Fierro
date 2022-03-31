#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <mpi.h>
#include "Solver.h"
#include "Implicit_Solver.h"

void solver_setup(int argc, char *argv[]);

//==============================================================================
//    Main
//==============================================================================

int main(int argc, char *argv[]){
  
  //initialize MPI
  MPI_Init(&argc,&argv);

  solver_setup(argc, argv);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}

void solver_setup(int argc, char *argv[]){
  /*General strategy: process initial input here to determine which solver
    object to construct
  */
  
  //base solver class pointer
  //Solver *solver;
  
  //set base pointer to the chosen solver
  //solver = new Static_Solver_Parallel();
  //Static_Solver_Parallel solver;
  Implicit_Solver solver;
  //solver = new Static_Solver();
  //solver = new Pseudo_Laplacian();

  //checks for optional solver routines
  if(solver.setup_flag) solver.setup();

  // invoke solver's run function (should perform most of the computation)//
  solver.run(argc,argv);
  
  //invoke optional finalize function
  if(solver.finalize_flag) solver.finalize();
 
  //delete solver;
}
