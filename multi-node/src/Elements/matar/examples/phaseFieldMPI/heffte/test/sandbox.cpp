
#include "heffte.h"
#include "test_fft3d.h"

/*
 * This method should be EMPTY in every pull request.
 * The goal is to have a file that is part of the build system
 * but it is not en encumbered by template parameters or complex text logic.
 * A single test for just one backend or one case of options/inputs
 * can be easily isolated here and tested.
 * This can also be used for profiling and running benchmarks
 * on sub-modules, e.g., a single reshape operation
 * or direct call to a backend.
 */
void test_sandbox(MPI_Comm const){
    // add code here to be tested by the system
}

void test_sandbox(){
    // same as above, but no MPI will be used
}

int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);

    test_sandbox();
    test_sandbox(MPI_COMM_WORLD);

    MPI_Finalize();
}
