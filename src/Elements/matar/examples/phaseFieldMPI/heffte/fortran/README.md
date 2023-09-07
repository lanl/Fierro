# Fortran Wrappers

HeFFTe provides Fortran 2003 wrappers for the C++ classes using [Swig-Fortran](https://github.com/swig-fortran). However, users don't need to use Swig as the generated files are included in heFFTe and Swig is needed only if changes to the C++ code is made.


#### Building the pre-generated wrappers

The Fortran wrappers for all enabled backends can be compiled with:
```
    -D Heffte_ENABLE_FORTRAN=ON
```
Each backend will generate a separate Fortrain module, e.g.,
```
    use heffte_fftw
    use heffte_cufft
    use heffte_mkl
```
The modules will be installed in `<prefix>\include` and multiple modules can be used simultaneously without causing a conflicts. Since the modules are just wrappers around C++, the corresponding C++ backend must also be enabled.


#### Linking with CMake

The heFFTe package-config imports component `Fortrain` and `Heffte_Fortran_FOUND` variable. All modules can be linked using the `Heffte::Fortrain` target:
```
    target_link_libraries(<user-target> Heffte::Fortran)
```
Other build systems needs the mod-files in `<prefix>\include` and the C++ and Fortrain libraries in `<prefix>\lib`.


#### Fortran 2003 API

The wrappers are designed to mimic the C++ API using Fortran 2003 classes, see the Fortran example installed in:
```
    <prefix>\share\heffte\examples\
```
Each module introduces two classes:
```
    heffte_fft3d_<backend>       e.g., heffte_fft3d_fftw
    heffte_fft3d_r2c_<backend>   e.g., heffte_fft3d_r2c_fftw
```
The constructors of the classes do not use heffte::box3d instead replace each box with 6 or 9 integers indicating the low and high indexes and the optional order. In Swig-Fortran, objects cannot be constructed inline without causing memory leaks and fully copying the C++ interface would require a potentially clumsy construction/destruction of objects. The integer constructors are a better alternative. Two constructors use the 6 and 9 integers, and a third one accepts three booleans corresponding to the plan options for reorder, all-to-all, and pencils. See the constructors for heffte::fft3d marked for internal use.

The Swig-Fortran objects must be manually destroyed with a call to:
```
    call fft%release()
```
This will invoke the C++ destructor.

The transform operations are called with:
```
    call fft%forward(input, output)
    call fft%backward(output, input)
```
The methods use Fortran generic overloads and can automatically handle all valid real and complex inputs with single and double precision. User defined workspace vector and a scaling parameter are optional inputs for the overloads.

The enumerated class used by the scaling parameter is replaced by constants:
```
    scale_<backend>_<none,symmetric,full>  e.g., call fft%forward(input, output, scale_fftw_full)
```
The backend name is added to the constant to avoid conflicts when multiple modules are using in the same Fortran file; however, the constants are just integers and are interchangeable.

The methods to query various sizes work the same as in C++:
```
    in_size   = fft%size_inbox()
    out_size  = fft%size_output()
    work_size = fft%size_workspace()
```

#### Rebuilding the Swig wrappers

Swig-Fortran can be easily installed using [spack](https://spack.io/) with the command:
```
    spack install swig@fortran
```
The swig instance is installed in `<spack-root>/opt/spack/<arch>/<compiler>/swig-fortran-<hash>/`.

Regenerating the heFFTe wrappers can be enabled withing the heFFTe CMake system with:
```
    -D Heffte_ENABLE_FORTRAN=ON
    -D Heffte_ENABLE_SWIG=ON
    -D SWIG_EXECUTABLE=<path-to-spack-install-of-swig>/bin/swigfortran
```
Using the above commands will update the wrappers for the currently selected heFFTe backends.

Developers of heFFTe need to rebuild the wrappers for all backends even if those that are not currently selected and that can be achieved with adding the option:
```
    -D Heffte_regenerate_all_swig=ON
```
Since all heFFTe backends cannot be enabled simultaneously, e.g., CUDA and ROCM, the build with the above option will fail; however, the wrappers are updated before the failure.
