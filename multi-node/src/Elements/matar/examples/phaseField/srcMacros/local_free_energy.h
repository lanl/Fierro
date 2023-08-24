#pragma once
#include "matar.h"

using namespace mtr; // matar namespace

double calculate_total_free_energy(int* nn, double* delta, double kappa, DCArrayKokkos<double> &comp);

void calculate_dfdc(int* nn, DCArrayKokkos<double> &comp, CArrayKokkos<double> &dfdc);
