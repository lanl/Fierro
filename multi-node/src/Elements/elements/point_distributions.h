#pragma once

#include "common.h"

template <typename NumType>
void equispaced_points(SizeType N, NumType &zl, NumType &zr, NumType *z);

template <typename NumType>
void chebyshev_points(SizeType N, NumType &zl, NumType &zr, NumType *z);
