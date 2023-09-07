#include "heffte.h"

/*!
 * \brief HeFFTe example 4, simple DFT from real to complex computing only the non-conjugate coefficients
 *
 * Performing DFT on three dimensional data in a box of 20 by 20 by 20
 * split across several processors (less than 10).
 *
 * If the input of an FFT transform consists of all real numbers,
 * the output comes in conjugate pairs which can be exploited to reduce
 * both the floating point operations and MPI communications.
 * Given a global set of indexes, HeFFTe can compute the corresponding DFT
 * and exploit the real-to-complex symmetry by selecting a dimension
 * and reducing the indexes by roughly half (the exact formula is floor(n / 2) + 1).
 */
void compute_dft(MPI_Comm comm){

    // select the default CPU backend, first available in the following order MKL, FFTW, Stock
    using backend_tag = heffte::backend::default_backend<heffte::tag::cpu>::type;

    // wrapper around MPI_Comm_rank() and MPI_Comm_size(), using this is optional
    int const me        = heffte::mpi::comm_rank(comm);
    int const num_ranks = heffte::mpi::comm_size(comm);

    if (num_ranks > 9){
        if (me == 0) std::cout << " heffte_example_r2c should use less than 10 ranks, exiting \n";
        return;
    }

    if (me == 0) std::cout << "using backend: " << heffte::backend::name<backend_tag>() << "\n";

    // the dimension where the data will shrink
    int r2c_direction = 0;

    // using problem with size 20x20x20 problem, the computed indexes are 11x20x20
    // direction 0 is chosen to reduce the number of indexes
    heffte::box3d<> real_indexes({0, 0, 0}, {19, 19, 19});
    heffte::box3d<> complex_indexes({0, 0, 0}, {10, 19, 19});

    // check if the complex indexes have correct dimension
    assert(real_indexes.r2c(r2c_direction) == complex_indexes);

    // report the indexes
    if (me == 0){
        std::cout << "The global input contains " << real_indexes.count() << " real indexes.\n";
        std::cout << "The global output contains " << complex_indexes.count() << " complex indexes.\n";
    }

    // see the heffte_example_options for comments on the proc_grid and boxes
    // the proc_grid is chosen to minimize the real data, but use for both real and complex cases
    std::array<int, 3> proc_grid = heffte::proc_setup_min_surface(real_indexes, num_ranks);

    std::vector<heffte::box3d<>> real_boxes    = heffte::split_world(real_indexes,    proc_grid);
    std::vector<heffte::box3d<>> complex_boxes = heffte::split_world(complex_indexes, proc_grid);

    heffte::box3d<> const inbox  = real_boxes[me];
    heffte::box3d<> const outbox = complex_boxes[me];

    // define the heffte class and the input and output geometry
    heffte::fft3d_r2c<backend_tag> fft(inbox, outbox, r2c_direction, comm);

    // vectors with the correct sizes to store the input and output data
    std::vector<float> input(fft.size_inbox());
    std::iota(input.begin(), input.end(), 0); // put some data in the input

    // output has std::vector<std::complex<float>>
    auto output = fft.forward(input);

    // verify that the output has the correct size
    assert(output.size() == static_cast<size_t>(fft.size_outbox()));

    // in the r2c case, the result from a backward transform is always real
    // thus inverse is std::vector<float>
    auto inverse = fft.backward(output, heffte::scale::full);

    // double-check the types
    static_assert(std::is_same<decltype(output), std::vector<std::complex<float>>>::value,
                  "the output should be a vector of std::complex<float>");
    static_assert(std::is_same<decltype(inverse), std::vector<float>>::value,
                  "the inverse should be a vector of float");

    // compare the computed entries to the original input data
    // the error is the max difference between the input and the real parts of the inverse
    // or the absolute value of the complex numbers
    float err = 0.0;
    for(size_t i=0; i<input.size(); i++)
        err = std::max(err, std::abs(inverse[i] - input[i]));

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
