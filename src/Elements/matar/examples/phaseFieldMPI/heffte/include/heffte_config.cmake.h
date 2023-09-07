/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_CONFIG_H
#define HEFFTE_CONFIG_H

#define Heffte_VERSION_MAJOR @Heffte_VERSION_MAJOR@
#define Heffte_VERSION_MINOR @Heffte_VERSION_MINOR@
#define Heffte_VERSION_PATCH @Heffte_VERSION_PATCH@

#define Heffte_GIT_HASH "@Heffte_GIT_HASH@"

#cmakedefine Heffte_ENABLE_AVX
#cmakedefine Heffte_ENABLE_AVX512

#cmakedefine Heffte_ENABLE_FFTW
#cmakedefine Heffte_ENABLE_MKL
#cmakedefine Heffte_ENABLE_CUDA
#cmakedefine Heffte_ENABLE_ROCM
#cmakedefine Heffte_ENABLE_ONEAPI

#cmakedefine Heffte_ENABLE_MAGMA

#cmakedefine Heffte_ENABLE_TRACING

#cmakedefine Heffte_DISABLE_GPU_AWARE_MPI

#if defined(Heffte_ENABLE_CUDA) || defined(Heffte_ENABLE_ROCM) || defined(Heffte_ENABLE_ONEAPI)
#define Heffte_ENABLE_GPU
#endif

#endif  /* HEFFTE_CONFIG_H */
