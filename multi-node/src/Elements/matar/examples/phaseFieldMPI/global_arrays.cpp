#include "global_arrays.h"
#include "mpi.h"

GlobalArrays::GlobalArrays(const std::array<int,3> & nn_all, const std::array<int,3> & nn) :
comp(nn[2], nn[1], nn[0]),
dfdc(nn[2], nn[1], nn[0])
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (0 == rank) {
        this->comp_all = CArray<double>(nn_all[2], nn_all[1], nn_all[0]);
    }
}
