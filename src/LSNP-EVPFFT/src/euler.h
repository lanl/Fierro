#pragma once

#include "definitions.h"

using namespace utils;

KOKKOS_FUNCTION
void euler(int iopt, real_t &ph, real_t &th, real_t &tm, real_t *a_);
