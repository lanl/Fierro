#pragma once

#include "heffte.h"

#ifdef USE_CUFFT
    using heffte_backend = heffte::backend::cufft;
#elif USE_ROCFFT
    using heffte_backend = heffte::backend::rocfft;
#elif USE_FFTW
    using heffte_backend = heffte::backend::fftw;
#elif USE_MKL
    using heffte_backend = heffte::backend::mkl;
#endif
