/** @class */
/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#include "test_fft3d.h"

template<typename backend_tag>
void test_fft3d_queues(MPI_Comm const comm){
    int const num_ranks = mpi::comm_size(comm);

    switch(num_ranks){
        case 6:
            test_fft3d_queues<backend_tag, std::complex<float>, 11, 11, 20>(comm);
            test_fft3d_queues<backend_tag, std::complex<double>, 11, 11, 20>(comm);
            test_fft3d_r2c_queues<backend_tag, float, 16, 11, 13>(comm);
            test_fft3d_r2c_queues<backend_tag, double, 16, 11, 9>(comm);
            break;
        default:
            throw std::runtime_error("No test for the given number of ranks!");
    }
}

void perform_tests(MPI_Comm const comm){
    all_tests<> name("heffte::fft streams");

    test_fft3d_queues<backend::stock>(comm);
    #ifdef Heffte_ENABLE_FFTW
    test_fft3d_queues<backend::fftw>(comm);
    #endif
    #ifdef Heffte_ENABLE_MKL
    test_fft3d_queues<backend::mkl>(comm);
    #endif
    #ifdef Heffte_ENABLE_CUDA
    test_fft3d_queues<backend::cufft>(comm);
    #endif
    #ifdef Heffte_ENABLE_ROCM
    test_fft3d_queues<backend::rocfft>(comm);
    #endif
    #ifdef Heffte_ENABLE_ONEAPI
    test_fft3d_queues<backend::onemkl>(comm);
    #endif
}

int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);

    perform_tests(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
