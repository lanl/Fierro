#include "global_arrays.h"

GlobalArrays::GlobalArrays(int* nn)
{   
    this->comp = DCArrayKokkos<double>(nn[0], nn[1], nn[2]);
    this->dfdc = CArrayKokkos<double>(nn[0], nn[1], nn[2]);
}
