# Basic Usage

HeFFTe provides a modern C++ interface for distributed Fast Fourier Transforms. Algorithms for data-transformations are implemented using MPI calls and third-patry backend libraries are used for the actual FFT operations. Templated front-end interface allows the use of multiple backends through a single API. Currently the supported backends include FFTW, MKL, oneMKL, cuFFT and rocFFT.

The interface and all associated backends are included in the `libheffte` library which can be linked against using the CMake target `Heffte::Heffte`. The `heffte.h` header includes all needed definitions and the HeFFTe methods are wrapped in the `heffte::` namespace. The user describes a distributed FFT problem with either a `heffte::fft3d` or a `heffte::fft3d_r2c` objects. The global set of indexes is distributed across the ranks of an MPI communicators where each rank has an input and output sub-box of the indexes, those are described with `box3d` structures that list the lowest and highest index in each direction. In addition to the distribution of the indexes described by the boxes and the associated MPI communicator, the the constructors of the fft classes accept a `plan_options` struct set of options about using algorithm details, e.g., point-to-point or all-to-all communication algorithm, strided vs contiguous 1d transforms, and others.

Transforms are supported in both 2D and 3D where a 2D problem is simply a global set of indexes with size 1 in one direction. The `heffte::fft2d` and `heffte::box2d` aliases are supported for more expressive user code in the 2D case.

In addition to the standard Fourier transform, heFFTe also supports the Sine and Cosine transforms that use real-to-real inputs and outputs, see `heffte::rtransform` for details.

The actual transforms are computed with the `forward()` and `backward()` methods that provide overloads that work with either vector or vector-like containers or raw-arrays, and optional external work-buffers. The result can also be scaled, per user request. A single class can handle transforms with different precision and data-types, so long as the distribution of the indexes is the same. The C++11 API also provides type safety and compatibility with `std::complex` types as well as the types provided by the backends, e.g., `fftw_complex`.

See also the examples in `<install-prefix>/share/heffte/examples/` or `<source-prefix>/examples/`.
