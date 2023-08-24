#pragma once
#include "matar.h"

using namespace mtr; // matar namespace

struct GlobalArrays
{
    DCArrayKokkos<double> comp;
    CArrayKokkos<double> dfdc;

    GlobalArrays(int* nn);
};
