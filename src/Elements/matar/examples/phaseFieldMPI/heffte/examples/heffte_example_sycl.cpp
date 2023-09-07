#include "heffte.h"

/*!
 * \brief HeFFTe example 6, using the OneMKL backend.
 *
 * The OneAPI backend of heFFTe uses unified memory buffers (UMC) and it is
 * recommended that the user provides their own sycl::queue.
 */
void compute_dft(MPI_Comm comm){

    // wrapper around MPI_Comm_rank() and MPI_Comm_size(), used here for convenience
    int const my_rank   = heffte::mpi::comm_rank(comm);
    int const num_ranks = heffte::mpi::comm_size(comm);

    // it is highly recommended to pass a user created sycl::queue to the constructor of fft3d
    // otherwise heFFTe will use an internal queue which may create a second device context
    // which may in turn cause all sorts of problems
    sycl::queue q;

    heffte::box3d<> const world = {{0, 0, 0}, {63, 63, 63}};

    // find the three dimensional grid of boxes with minimal surface area
    std::array<int,3> proc_grid = heffte::proc_setup_min_surface(world, num_ranks);

    // using the grid dimensions, generate the boxes
    std::vector<heffte::box3d<>> boxes = heffte::split_world(world, proc_grid);

    // create the heFFTe plan and associate it with the given queue
    heffte::fft3d<heffte::backend::onemkl> fft(q, boxes[my_rank], boxes[my_rank], comm);

    // allocate memory on the CPU
    std::vector<std::complex<double>> input(fft.size_inbox());
    std::vector<std::complex<double>> output(fft.size_outbox());
    std::vector<std::complex<double>> inverse(fft.size_inbox());

    // fill in some data, the data is not important for the example
    std::iota(input.begin(), input.end(), 0);

    // allocate memory on the device
    // it is preferable to wrap these in some container for RAII resource management
    // heFFTe already does that with the heffte::gpu::vector template
    // this example uses raw UMC arrays for illustrative purposes
    std::complex<double> *sycl_input = sycl::malloc_device<std::complex<double>>(input.size(), q);
    std::complex<double> *sycl_output = sycl::malloc_device<std::complex<double>>(output.size(), q);
    std::complex<double> *sycl_inverse = sycl::malloc_device<std::complex<double>>(inverse.size(), q);

    // copy the buffer to the device
    q.memcpy(sycl_input, input.data(), input.size() * sizeof(std::complex<double>)).wait();

    // perform a forward and backward FFT
    // note that the calls to heFFTe are actually blocking, which is primarily dictated by MPI
    fft.forward(sycl_input, sycl_output, heffte::scale::full);
    fft.backward(sycl_output, sycl_inverse);

    // get back the result into the std::vector buffers
    q.memcpy(output.data(), sycl_inverse, output.size() * sizeof(std::complex<double>)).wait();
    q.memcpy(inverse.data(), sycl_inverse, inverse.size() * sizeof(std::complex<double>)).wait();

    // without a RAII container, the memory has to be manually deallocated
    sycl::free(sycl_input, q);
    sycl::free(sycl_output, q);
    sycl::free(sycl_inverse, q);

    double err = 0.0;
    for(size_t i=0; i<input.size(); i++)
        err = std::max(err, std::abs(input[i] - inverse[i]));

    std::cout << std::scientific;
    for(int i=0; i<num_ranks; i++){
        MPI_Barrier(comm);
        if (my_rank == i) std::cout << "rank " << i << " computed error: " << err << std::endl;
    }
}

int main(int argc, char** argv){

    MPI_Init(&argc, &argv);

    compute_dft(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
