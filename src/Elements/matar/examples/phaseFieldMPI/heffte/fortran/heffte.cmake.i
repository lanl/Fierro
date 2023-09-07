/*!
 * \file heffte.cmake.i
 *
 * Copyright (c) 2020, University of Tennessee.
 */

%module "heffte_@heffte_backend@"

/* -------------------------------------------------------------------------
 * Header definition macros
 * ------------------------------------------------------------------------- */

%define %heffte_add_header
%insert("fbegin") %{
! heFFTe project, https://bitbucket.org/icl/heffte/
! Copyright (c) 2020, University of Tennessee.
! Distributed under a BSD 3-Clause license: see LICENSE for details.
%}
%insert("begin") %{
/*
 * heFFTe project, https://bitbucket.org/icl/heffte/
 * Copyright (c) 2020, University of Tennessee.
 * Distributed under a BSD 3-Clause license: see LICENSE for details.
 */
%}
%enddef

%heffte_add_header

/* -------------------------------------------------------------------------
 * Exception handling
 * ------------------------------------------------------------------------- */

// Rename the error variables' internal C symbols
#define SWIG_FORTRAN_ERROR_INT heffte_@heffte_backend@_ierr
#define SWIG_FORTRAN_ERROR_STR heffte_@heffte_backend@_get_serr

// Restore names in the wrapper code
%rename(ierr) heffte_@heffte_backend@_ierr;
%rename(get_serr) heffte_@heffte_backend@_get_serr;

%include <exception.i>
%include "mpi.i"

/* -------------------------------------------------------------------------
 * Data types and instantiation
 * ------------------------------------------------------------------------- */

// Note: stdint.i inserts #include <stdint.h>
%include <stdint.i>

// Unless otherwise specified, ignore functions that use unknown types
%fortranonlywrapped;

// Backend name appended to each enum, e.g., scale_fftw_full
%rename("@heffte_backend@_%s",%$isenumitem) "";

/* -------------------------------------------------------------------------
 * Wrapped files
 * ------------------------------------------------------------------------- */

%{
#include "heffte.h"
%}

// Allow native fortran arrays to be passed to pointer/arrays
%include <typemaps.i>
%include <complex.i>
%apply SWIGTYPE ARRAY[] {
    int*,
    float*,
    double*,
    std::complex<float>*,
    std::complex<double>*,
    int[],
    float[],
    double[],
    std::complex<float>[],
    std::complex<double>[]
};

// Wrap clases/methods found in this file
%include "heffte_fft3d.h"
%include "heffte_fft3d_r2c.h"

namespace heffte {
    namespace backend {
        %ignore @heffte_backend@;
        struct @heffte_backend@{};
    }

    %extend fft3d {
        %template(forward) forward<float, std::complex<float>>;
        %template(forward) forward<std::complex<float>, std::complex<float>>;
        %template(forward) forward<double, std::complex<double>>;
        %template(forward) forward<std::complex<double>, std::complex<double>>;
        %template(backward) backward<std::complex<float>, float>;
        %template(backward) backward<std::complex<float>, std::complex<float>>;
        %template(backward) backward<std::complex<double>, double>;
        %template(backward) backward<std::complex<double>, std::complex<double>>;
    }
    %extend fft3d_r2c {
        %template(forward) forward<float, std::complex<float>>;
        %template(forward) forward<double, std::complex<double>>;
        %template(backward) backward<std::complex<float>, float>;
        %template(backward) backward<std::complex<double>, double>;
    }

    %template(heffte_fft3d_@heffte_backend@) fft3d<backend::@heffte_backend@>;
    %template(heffte_fft3d_r2c_@heffte_backend@) fft3d_r2c<backend::@heffte_backend@>;
}
