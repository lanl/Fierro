#include "heffte.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct{
    double r;
    double i;
} example_complex;

/*!
 * \brief HeFFTe example C, simple DFT using two MPI ranks, FFTW backend and the C API.
 *
 * Performing DFT on three dimensional data in a box of 4 by 4 by 4 split
 * across the third dimension between two MPI ranks.
 * See the C++ FFTW example for reference.
 */
void compute_dft(MPI_Comm comm){

    int me; // this process rank within the comm
    MPI_Comm_rank(comm, &me);

    int num_ranks; // total number of ranks in the comm
    MPI_Comm_size(comm, &num_ranks);

    if (num_ranks != 2){
        if (me == 0) printf(" heffte_example_c is set to 2 mpi ranks, exiting \n");
        return;
    }

    int i = 0;
    int box_low[3] = {0, 0, 0};
    int box_high[3] = {3, 3, 3};
    if (me == 0) // split across the last dimension
        box_high[2] = 1;
    else
        box_low[2] = 2;

    heffte_plan plan;
    int status = heffte_plan_create(Heffte_BACKEND_FFTW, box_low, box_high, NULL, box_low, box_high, NULL, comm, NULL, &plan);
    if (status != Heffte_SUCCESS){
        printf("Failed at heffte_plan_create() with error code: %d\n", status);
        MPI_Abort(comm, 1);
    }

    int size_inbox  = heffte_size_inbox(plan);
    int size_outbox = heffte_size_outbox(plan);

    double *input           = malloc(size_inbox * sizeof(double));
    for(i=0; i<size_inbox; i++) input[i] = (double) i;

    example_complex *output = calloc(size_outbox, sizeof(example_complex));

    heffte_forward_d2z(plan, input, output, Heffte_SCALE_NONE);

    if (me == 1){
        int pass = 0;
        if (fabs(output[0].r + 512.0) > 1.E-11 || fabs(output[0].i) > 1.E-11)
            pass = 1;
        for(i=1; i<size_outbox; i++)
            if (fabs(output[i].r) > 1.E-11 || fabs(output[i].i) > 1.E-11)
                pass = 1;
        if (pass){
            printf("The computed transform deviates by more than the tolerance.\n");
            MPI_Abort(comm, 1);
        }
    }

    for(i=0; i<size_inbox; i++) input[i] = 0.0;

    heffte_backward_z2d(plan, output, input, Heffte_SCALE_FULL);

    status = heffte_plan_destroy(plan);
    if (status != Heffte_SUCCESS){
        printf("Failed at heffte_plan_destroy() with error code: %d\n", status);
        MPI_Abort(comm, 1);
    }

    double err = 0.0;
    for(i=0; i<size_inbox; i++)
        if (fabs(input[i] - (double) i) > err)
            err = fabs(input[i] - (double) i);

    if (me == 0) printf("rank 0 computed error: %1.6le\n", err);
    MPI_Barrier(comm);
    if (me == 1) printf("rank 1 computed error: %1.6le\n", err);

    free(input);
    free(output);
}

int main(int argc, char** argv){

    MPI_Init(&argc, &argv);

    compute_dft(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
