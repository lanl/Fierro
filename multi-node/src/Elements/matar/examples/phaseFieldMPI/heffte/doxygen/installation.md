# Installation

[TOC]

### Requirements

At the minimum, heFFTe requires a C++11 capable compiler,
an implementation of the Message Passing Library (MPI),
and at least one backend FFT library.
The heFFTe library can be build with either CMake 3.10 or newer,
or a simple GNU Make build engine.
CMake is the recommended way to use heFFTe since dependencies and options
are much easier to export to user projects and not all options
could be cleanly implemented in the rigid Makefile.

| Compiler | Tested versions |
|----|----|
| gcc      | 7 - 11          |
| clang    | 5 - 14          |
| icc      | 18              |
| dpcpp    | 2021.2          |
| OpenMPI  | 4.0.3           |

Tested backend libraries:

| Backend        | Tested versions |
|----|----|
| stock          | all             |
| fftw3          | 3.3.7 - 3.3.10  |
| mkl            | 2016            |
| oneapi/onemkl  | 2021.4          |
| cuda/cufft     | 10 - 11         |
| rocm/rocfft    | 4.0 - 5.6       |

The listed tested versions are part of the continuous integration and nightly build systems,
but heFFTe may yet work with other compilers and backend versions.

### CMake Installation

Typical CMake build follows the steps:
```
mkdir build
cd build
cmake <cmake-build-command-options> <path-to-heffte-source>
make
make install
make test_install
```

Typical CMake build command:
```
cmake \
    -D CMAKE_BUILD_TYPE=Release \
    -D BUILD_SHARED_LIBS=ON     \
    -D CMAKE_INSTALL_PREFIX=<path-for-installation> \
    -D Heffte_ENABLE_AVX=ON \
    -D Heffte_ENABLE_FFTW=ON \
    -D FFTW_ROOT=<path-to-fftw3-installation> \
    -D Heffte_ENABLE_CUDA=ON \
    -D CUDA_TOOLKIT_ROOT_DIR=<path-to-cuda-installation> \
    <path-to-heffte-source-code>
```

The standard CMake options are also accepted:
```
    CMAKE_CXX_COMPILER=<path-to-suitable-cxx-compiler>        (sets the C++ compiler)
    CMAKE_CXX_FLAGS="<additional cxx flags>"                  (adds flags to the build process)
    MPI_CXX_COMPILER=<path-to-suitable-mpi-compiler-wrapper>  (specifies the MPI compiler wrapper)
```

Additional heFFTe options:
```
    Heffte_ENABLE_ROCM=<ON/OFF>      (enable the rocFFT backend)
    Heffte_ENABLE_ONEAPI=<ON/OFF>    (enable the oneMKL backend)
    Heffte_ENABLE_MKL=<ON/OFF>       (enable the MKL backend)
    Heffte_ENABLE_AVX512=<ON/OFF>    (enable AVX512 instructions in the stock backend)
    MKL_ROOT=<path>                  (path to the MKL folder)
    Heffte_ENABLE_DOXYGEN=<ON/OFF>   (build the documentation)
    Heffte_ENABLE_TRACING=<ON/OFF>   (enable the even logging engine)
```

Additional language interfaces and helper methods:
```
    -D Heffte_ENABLE_PYTHON=<ON/OFF>   (configure the Python module)
    -D Heffte_ENABLE_FORTRAN=<ON/OFF>  (build the Fortran modules)
    -D Heffte_ENABLE_SWIG=<ON/OFF>     (generate new Fortrans source files)
    -D Heffte_ENABLE_MAGMA=<ON/OFF>    (link to MAGMA for helper methods)
```
See the Fortran and Python sections for details.

### List of Available Backend Libraries

* **Stock:** the stock backend library is a CPU-based 1-D FFT library built *into* HeFFTe, and requires no additional external libraries or hardware. This is enabled by default, and can optionally perform vectorized complex arithmetic using AVX and AVX512 instruction sets (this accelerates the backend by up to an order of magnitude, but includes no vectorization by default). This vectorization is enabled with
```
    -D Heffte_ENABLE_AVX=<ON/OFF>
    -D Heffte_ENABLE_AVX512=<ON/OFF>
```
Just using `Heffte_ENABLE_AVX512` enables all the vectorization flags. Using vectorization is optional and not required for using the stock backend.

* **FFTW3:** the [fftw3](http://www.fftw.org/) library is the de-facto standard for open source FFT library and is distributed under the GNU General Public License, the fftw3 backend is enabled with:
```
    -D Heffte_ENABLE_FFTW=ON
    -D FFTW_ROOT=<path-to-fftw3-installation>
```
Note that fftw3 uses two different libraries for single and double precision, while HeFFTe handles all precisions in a single template library; thus both the fftw3 (double-precision) and fftw3f (single-precision) variants are needed and those can be installed in the same path.

* **MKL:** the [Intel Math Kernel Library](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html) provides optimized FFT implementation targeting Intel processors and can be enabled within heFFTe with:
```
    -D Heffte_ENABLE_MKL=ON
    -D MKL_ROOT=<path-to-mkl-installation>
```
The `MKL_ROOT` default to the environment variable `MKLROOT` (chosen by Intel). The additional variable `Heffte_MKL_THREAD_LIBS` allows to choose the MKL threaded backend, tested with `mkl_gnu_thread` and `mkl_intel_thread`, the default is to use GNU-threads on GCC compiler and Intel otherwise. Note that `mkl_intel_thread` also requires `libiomp5.so` and heFFTe will search for it in the system paths and `LD_LIBRARY_PATH`, unless the variable `Heffte_MKL_IOMP5` is defined and pointing to `libiomp5.so`. GNU-threads do not use `libiomp5.so` but the GNU `libgomp.so` which CMake finds automatically.

* **oneMKL** the [Intel oneMKL Library](https://spec.oneapi.com/versions/latest/elements/oneMKL/source/index.html) provides optimized FFT implementation targeting Intel GPUs and can be enabled within heFFTe with:
```
    -D CMAKE_CXX_COMPILER=dpcpp
    -D Heffte_ENABLE_ONEAPI=ON
    -D Heffte_ONEMKL_ROOT=<path-to-onemkl-installation>
```
The `Heffte_ONEMKL_ROOT` defaults to the environment variable `ONEAPI_ROOT`. The Intel oneAPI framework provides support for both CPU and GPU devices from the same code-base, thus the oneMKL libraries require the CPU MKL libraries and the MKL options will be enabled automatically if oneMKL is selected.

* **CUFFT:** the [Nvidia CUDA framework](https://developer.nvidia.com/cuda-zone) provides a GPU accelerated FFT library [cuFFT](https://docs.nvidia.com/cuda/cufft/index.html), which can be enabled in heFFTe with:
```
    -D Heffte_ENABLE_CUDA=ON
    -D CUDA_TOOLKIT_ROOT_DIR=<path-to-cuda-installation>
```

* **ROCFFT:**  the [AMD ROCm framework](https://github.com/RadeonOpenCompute/ROCm) provides a GPU accelerated FFT library [rocFFT](https://github.com/ROCmSoftwarePlatform/rocFFT), which can be enabled in heFFTe with:
```
    -D CMAKE_CXX_COMPILER=hipcc
    -D Heffte_ENABLE_ROCM=ON
```

**Note:** Only one of the GPU backends can be enabled (CUDA, ROCM, or ONEAPI) since the three backends operate with arrays allocated in GPU device memory (or alternatively shared/managed memory). By default when using either GPU backend, heFFTe assumes that the MPI implementation is GPU-Aware, see the next section.


### GPU-Aware MPI

Different implementations of MPI can provide GPU-Aware capabilities, where data can be send/received directly in GPU memory. OpenMPI provided CUDA aware capabilities if compiled with the corresponding options, e.g., see [CUDA-Aware OpenMPI](https://www.open-mpi.org/faq/?category=buildcuda). Both CUDA and ROCm support such API; however, the specific implementation available to the user may not be available for various reasons, e.g., insufficient hardware support. HeFFTe can be compiled without GPU-Aware capabilities with the CMake option:
```
    -D Heffte_DISABLE_GPU_AWARE_MPI=ON
```
**Note:** The GPU-Aware capabilities can also be disabled at runtime by setting the corresponding option in `heffte::plan_options`. On some platforms, the GPU-Aware MPI calls have significantly larger latency and moving the buffers to the CPU before ini


### Linking to HeFFTe

HeFFTe installs a CMake package-config file in
```
    <install-prefix>/lib/cmake/
```
Typical project linking to HeFFTe will look like this:
```
    project(heffte_user VERSION 1.0 LANGUAGES CXX)

    find_package(Heffte PATHS <install-prefix>)

    add_executable(foo ...)
    target_link_libraries(foo Heffte::Heffte)
```
An example is installed in `<install-prefix>/share/heffte/examples/`.

The package-config also provides a set of components corresponding to the different compile options, specifically:
```
    FFTW MKL CUDA ROCM ONEAPI PYTHON Fortran GPUAWARE
```


### GNU Make Installation
HeFFTe supports a GNU Make build engine, where dependencies and compilers
are set manually in the included Makefile.
Selecting the backends is done with:
```
    make backends=fftw,cufft
```
The `backends` should be separated by commas and must have correctly selected
compilers, includes, and libraries. Additional options are available, see
```
    make help
```
and see also the comments inside the Makefile.

Testing is invoked with:
```
    make ctest
```
The library will be build in `./lib/`

The GNU Make build engine does not support all options, e.g., MAGMA or disabling the GPU-Aware calls,
but is provided for testing and debugging purposes.


### Known Issues

* the current testing suite requires about 3GB of free GPU RAM
    * CUDA seem to reserve 100-200MB of RAM per MPI rank and some tests use 12 ranks
    * the GPU handles used by MAGMA-CUDA and MAGMA-HIP are not destroyed on time
      (both backends seem to implement some garbage collection mechanism),
      thus even an 8GB GPU may run out of memory when using the 12 rank test
