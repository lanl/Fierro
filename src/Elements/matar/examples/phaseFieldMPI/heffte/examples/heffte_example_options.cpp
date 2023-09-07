#include "heffte.h"

/*!
 * \brief HeFFTe example 2, simple DFT using four MPI ranks and multiple options.
 *
 * Performing DFT on three dimensional data in a box of 10 by 20 by 30
 * split across four processors.
 *
 * This is a more advanced example demonstrating how to adjust algorithm options
 * and how to use externally allocated scratch buffers.
 */
void compute_dft(MPI_Comm comm){

    // select the default CPU backend, first available in the following order MKL, FFTW, Stock
    using backend_tag = heffte::backend::default_backend<heffte::tag::cpu>::type;

    int me; // this process rank within the comm
    MPI_Comm_rank(comm, &me);

    if (me == 0) std::cout << "using backend: " << heffte::backend::name<backend_tag>() << "\n";

    int num_ranks; // total number of ranks in the comm
    MPI_Comm_size(comm, &num_ranks);

    if (num_ranks != 4){
        if (me == 0) std::cout << " heffte_example_options is set to 4 mpi ranks, exiting \n";
        return;
    }

    // using problem with size 10x20x30 problem
    heffte::box3d<> all_indexes({0, 0, 0}, {9, 19, 29});

    // create a processor grid with minimum surface (measured in number of indexes)
    std::array<int,3> proc_grid = heffte::proc_setup_min_surface(all_indexes, num_ranks);

    // split all indexes across the processor grid, defines a set of boxes
    std::vector<heffte::box3d<>> all_boxes = heffte::split_world(all_indexes, proc_grid);

    // pick the box corresponding to this rank
    heffte::box3d<> const inbox  = all_boxes[me];
    // pick a random different box
    // this demonstrates that the input and output boxes can be different
    // rank "me" will use box "(me + 1)%num_ranks" and will reorder
    // the fast, mid, and slow indexes
    // {2, 1, 0} means that the first index is not the slowest while the last is the fastest
    heffte::box3d<> const outbox(all_boxes[(me+1)%num_ranks].low,
                                 all_boxes[(me+1)%num_ranks].high,
                                 {2, 1, 0}
                                );

    // at this stage we can manually adjust some HeFFTe options
    heffte::plan_options options = heffte::default_options<backend_tag>();

    // use strided 1-D FFT operations
    // some backends work just as well when the entries of the data are not contiguous
    // then there is no need to reorder the data in the intermediate stages which saves time
    options.use_reorder = false;

    // use point-to-point communications
    // collaborative all-to-all and individual point-to-point communications are two alternatives
    // one may be better than the other depending on
    // the version of MPI, the hardware interconnect, and the problem size
    options.algorithm = heffte::reshape_algorithm::alltoallv;

    // in the intermediate steps, the data can be shapes as either 2-D slabs or 1-D pencils
    // for sufficiently large problem, it is expected that the pencil decomposition is better
    // but for smaller problems, the slabs may perform better (depending on hardware and backend)
    options.use_pencils = true;

    // define the heffte class and the input and output geometry
    heffte::fft3d<backend_tag> fft(inbox, outbox, comm, options);

    // vectors with the correct sizes to store the input and output data
    // taking the size of the input and output boxes
    // note that the input is real, the output is complex
    //      and we are not taking advantage of symmetry in the indexes
    std::vector<float> input(fft.size_inbox());
    std::vector<std::complex<float>> output(fft.size_outbox());

    // fill the input vector with data that looks like 0, 1, 2, ...
    std::iota(input.begin(), input.end(), 0); // put some data in the input

    // each HeFFTe forward() and backward() operations requires additional memory buffers
    // normally HeFFTe will allocate and free those buffers internally
    // but this can degrade performance and HeFFTe can also accept external scratch buffers
    // the buffer needs to have fft.size_workspace() complex entries of the correct precision
    std::vector<std::complex<float>> workspace(fft.size_workspace());

    // perform a forward DFT, scale the result
    fft.forward(input.data(), output.data(), workspace.data(), heffte::scale::full);

    // store the result from an inverse operation
    std::vector<float> inverse(input.size());

    // perform a backward DFT
    fft.backward(output.data(), inverse.data(), workspace.data());

    // compare the computed entries to the original input data
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
