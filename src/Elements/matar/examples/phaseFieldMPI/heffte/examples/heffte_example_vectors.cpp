#include "heffte.h"

/*!
 * \brief HeFFTe example 3, simple DFT using arbitrary MPI ranks and vector containers.
 *
 * Performing DFT on three dimensional data in a box of 10 by 20 by 30
 * split across several processors (less than 10).
 *
 * HeFFTe commands can be issues with either arrays or vector containers.
 * The containers allow for more expressive calls and simplify types and memory management,
 * although the vector overloads do not accept scratch buffers.
 * Backends with CPU implementations (e.g., fftw and mkl) use std::vector containers,
 * while the GPU implementation (e.g., cufft) uses heffte::cuda::vector which is similar
 * to the STL vector but wraps around an array allocated on the GPU device.
 */
void compute_dft(MPI_Comm comm){

    // select the default CPU backend, first available in the following order MKL, FFTW, Stock
    using backend_tag = heffte::backend::default_backend<heffte::tag::cpu>::type;

    // wrapper around MPI_Comm_rank() and MPI_Comm_size(), using this is optional
    int const me        = heffte::mpi::comm_rank(comm);
    int const num_ranks = heffte::mpi::comm_size(comm);

    if (num_ranks > 9){
        if (me == 0) std::cout << " heffte_example_vectors should use less than 10 ranks, exiting \n";
        return;
    }

    if (me == 0) std::cout << "using backend: " << heffte::backend::name<backend_tag>() << "\n";

    // using problem with size 10x20x30 problem
    heffte::box3d<> all_indexes({0, 0, 0}, {9, 19, 29});

    // see the heffte_example_options for comments on the proc_grid and boxes
    std::array<int, 3> proc_grid = heffte::proc_setup_min_surface(all_indexes, num_ranks);

    std::vector<heffte::box3d<>> all_boxes = heffte::split_world(all_indexes, proc_grid);

    heffte::box3d<> const inbox  = all_boxes[me];
    heffte::box3d<> const outbox = all_boxes[me]; // same inbox and outbox

    // at this stage we can manually adjust some HeFFTe options
    heffte::plan_options options = heffte::default_options<backend_tag>();

    // using slab decomposition
    options.use_pencils = false;

    // define the heffte class and the input and output geometry
    heffte::fft3d<backend_tag> fft(inbox, outbox, comm, options);

    // vectors with the correct sizes to store the input and output data
    std::vector<double> input(fft.size_inbox());
    std::iota(input.begin(), input.end(), 0); // put some data in the input

    // output has std::vector<std::complex<double>>
    auto output = fft.forward(input, heffte::scale::symmetric);

    // in a backward transform, the result could be either real or complex
    // since the input was real, we know that the result will be real
    // but this is not guaranteed in general
    // backward() returns std::vector<std::complex<double>>
    // backward_real() returns std::vector<double> using truncation of the complex part
    auto complex_inverse = fft.backward(output, heffte::scale::symmetric);
    auto real_inverse    = fft.backward_real(output, heffte::scale::symmetric);

    // double-check the types
    static_assert(std::is_same<decltype(complex_inverse), std::vector<std::complex<double>>>::value,
                  "the complex inverse should be a vector of std::complex<double>");
    static_assert(std::is_same<decltype(real_inverse), std::vector<double>>::value,
                  "the complex inverse should be a vector of double");

    // compare the computed entries to the original input data
    // the error is the max difference between the input and the real parts of the inverse
    // or the absolute value of the complex numbers
    double err = 0.0;
    for(size_t i=0; i<input.size(); i++){
        err = std::max(err, std::abs(std::real(complex_inverse[i]) - input[i]));
        err = std::max(err, std::abs(std::imag(complex_inverse[i])));
    }
    for(size_t i=0; i<input.size(); i++)
        err = std::max(err, std::abs(real_inverse[i] - input[i]));

    // print the error for each MPI rank
    std::cout << std::scientific;
    for(int i=0; i<num_ranks; i++){
        MPI_Barrier(comm);
        if (me == i) std::cout << "rank " << i << " computed error: " << err << std::endl;
    }
}

int main(int argc, char** argv){

    MPI_Init(&argc, &argv);

    compute_dft(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
