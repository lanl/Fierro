/** @class */
/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#include "test_fft3d.h"

template<typename backend_tag>
void test_fft3d_cases(MPI_Comm const comm){
    int const num_ranks = mpi::comm_size(comm);

    if (num_ranks != 12)
        throw std::runtime_error("No test for the given number of ranks!");

    test_fft3d_arrays<backend_tag, float, 19, 20, 21>(comm);
    test_fft3d_arrays<backend_tag, double, 19, 20, 21>(comm);
    test_fft3d_arrays<backend_tag, std::complex<float>, 19, 15, 25>(comm);
    test_fft3d_arrays<backend_tag, std::complex<double>, 19, 19, 17>(comm);
    test_batch_cases<backend_tag>(comm);
}

void perform_tests(MPI_Comm const comm){
    all_tests<> name("heffte::fft class");

    test_fft3d_cases<backend::stock>(comm);
    #ifdef Heffte_ENABLE_FFTW
    test_fft3d_cases<backend::fftw>(comm);
    #endif
    #ifdef Heffte_ENABLE_MKL
    test_fft3d_cases<backend::mkl>(comm);
    #endif
    #ifdef Heffte_ENABLE_CUDA
    test_fft3d_cases<backend::cufft>(comm);
    #endif
    #ifdef Heffte_ENABLE_ROCM
    test_fft3d_cases<backend::rocfft>(comm);
    #endif
    #ifdef Heffte_ENABLE_ONEAPI
    test_fft3d_cases<backend::onemkl>(comm);
    #endif
}

int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);

    perform_tests(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
