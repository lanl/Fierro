# Python Wrappers

HeFFTe provides a set of Python wrappers that use the C interface but provide object-oriented API mimicking C++. Thus, the C++ documentation remains relevant.

#### Requirements

The Python wrappers require the following modules:
```
    ctypes   # standard module included in most python distributions
    numpy    # used at the core of most Python numerical packages
    mpi3py   # provides MPI bindings for Python
```
In addition, using the CUDA backend requires:
```
    from numba import cuda
```
The ROCm bindings using `numba.hsa` are still under development.


#### Enabling Python in CMake

The heFFTe Python bindings are enabled with:
```
    -D Heffte_ENABLE_PYTHON=ON
    -D PYTHON_EXECUTABLE=<path-to-python>
```
The `PYTHON_EXECUTABLE` is used for testing purposes only, the `ctypes` interface allows for multiple Python versions to use the same heFFTe installation.

The module files will be installed in two locations:
```
    <install-prefix>/lib/pythonX.Y/site-packages
    <install-prefix>/share/heffte/python
```
The `site-packages` is the standard location for Python modules where `X` and `Y` represent the Python version. The location in `share` allows for a more version-independent access to the modules.

Using the heFFTe Python module requires that at least one location is accessible to the interpreter by either setting environment variable:
```
    export PYTHONPATH=<install-prefix>/share/heffte/python:$PYTHONPATH
```
or by adding the path to your `sys.path`
```
    import sys
    sys.path.append("<install-prefix>/share/heffte/python")
```


#### Basic Usage

HeFFTe requires at least two modules:
```
    import mpi4py  # gives access to the MPI capabilities
    import heffte  # the heFFTe wrappers
```
In addition, either `numpy` or `numba.cuda` modules are required to work with the CPU or GPU backends.

The Python equivalanets to the `heffte::fft3d` and `heffte::fft3d_r2c` are created with:
```
    fft = heffte.fft3d(backend_tag, inbox, outbox, comm)
    fft = heffte.fft3d_r2c(backend_tag, inbox, oubox, r2c_direction, comm)
```
While the methods have signature similar to Python constructors, the `fft3d` and `fft3d_r2c` are factory methods that create objects of type `heffte_fft_plan` which in turn wraps around a heFFTe plan from the C interface.

The `backend_tag` is a runtime variable and takes one of the values:
```
    heffte.backend.fftw    heffte.backend.cufft
    heffte.backend.mkl     heffte.backend.rocfft
```
The corresponding backend has to be also enabled in CMake and C++.

The `inbox` and `outbox` are instances of:
```
    heffte.box3d
```
The Python `box3d` class offers only the constructor that uses as inputs objects convertible to `numpy.ndarray` with three entries. The rest of the C++ `box3d` API is not available under Python.

The `comm` variable is an MPI communicator from the `mpi4py` module.

The `fft3d` objects provide the same methods for `size_inbox()`, `size_outbox()`, and `size_workspace()`.

The FFT transforms are performed with the `forward()` and `backward()` methods:
```
    fft3d.forward(inarray, outarray, scaling)
    fft3d.backward(outarray, inarray, scaling)
```
The `inarray` and `outarray` are instances of either `numpy.ndarray` or `DeviceNDArray` from the `numba.cuda` module. The sizes and types must match those in the calls to the C++ interface, see the documentation of `heffte::fft3d`.

The `scaling` variable is one of the constants:
```
    heffte.scale.none
    heffte.scale.symmetric
    heffte.scale.full
```


#### More to come soon

* Python examples for CPU and GPU
* Support for ROCm backend and the corresponding array types
* User allocated workspace buffers
