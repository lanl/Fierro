#pragma once

#include "definitions.h"

using namespace utils;

KOKKOS_FUNCTION
double optimizedPow(double base, int exponent);

KOKKOS_FUNCTION
void inverse_gj(real_t *a_, int n);

KOKKOS_FUNCTION
void lu_inverse(real_t *a_, int n);
KOKKOS_FUNCTION
void ludcmp(real_t *a_, int n, int np, int *indx_, real_t d, int &isingular);
KOKKOS_FUNCTION
void lubksb(real_t *a_, int n, int np, int *indx_, real_t *b_);

void inverse(real_t *a, int n, int np, real_t *y);
void ludcmpc(real_t *a_, int n, int np, int *indx_, real_t d);
void lubksbc(real_t *a_, int n, int np, int *indx_, real_t *b_);
