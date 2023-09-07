/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#include "heffte.h"

/*
 * This template will compile if-and-only-if the class T is
 * movable and move constructible, but not copyable and not copy-constructible.
 */
template<typename T>
void assert_move_and_not_copy(){
    static_assert(std::is_move_constructible<T>::value, "class T is not move constructible");
    static_assert(std::is_move_assignable<T>::value, "class T is not move assignable");

    static_assert(not std::is_copy_constructible<T>::value, "class T is copy constructible (will leak memory)");
    static_assert(not std::is_copy_assignable<T>::value, "class T is copy assignable (will leak memory)");
}

int main(int, char**){

    /*
     * Tests if the single header can be included and compiled without
     * any other includes. Guards against bugs due to headers included
     * in the test common header but missing from the main heffte.h.
     * It also does several compile-time checks for class properties.
     */

    assert_move_and_not_copy<heffte::fft3d<heffte::backend::stock>>();
    assert_move_and_not_copy<heffte::fft3d_r2c<heffte::backend::stock>>();
    #ifdef Heffte_ENABLE_FFTW
    assert_move_and_not_copy<heffte::fft3d<heffte::backend::fftw>>();
    assert_move_and_not_copy<heffte::fft3d_r2c<heffte::backend::fftw>>();
    #endif
    #ifdef Heffte_ENABLE_CUDA
    assert_move_and_not_copy<heffte::fft3d<heffte::backend::cufft>>();
    assert_move_and_not_copy<heffte::fft3d_r2c<heffte::backend::cufft>>();
    #endif
    #ifdef Heffte_ENABLE_MKL
    assert_move_and_not_copy<heffte::fft3d<heffte::backend::mkl>>();
    assert_move_and_not_copy<heffte::fft3d_r2c<heffte::backend::mkl>>();
    #endif

    return 0;
}
