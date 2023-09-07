/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_BACKENDS_H
#define HEFFTE_BACKENDS_H

// vector for RAII memory management (GPU only)
#include "heffte_backend_vector.h"

// the individual backends
#include "heffte_backend_stock.h"
#include "heffte_backend_fftw.h"
#include "heffte_backend_mkl.h"

#include "heffte_backend_cuda.h"
#include "heffte_backend_rocm.h"
#include "heffte_backend_oneapi.h"

// helpers to move data between the GPU and CPU
#include "heffte_backend_data_transfer.h"

#endif   /* HEFFTE_BACKENDS_H */
