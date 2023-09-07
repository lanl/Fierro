#include "heffte.h"

/*!
 * \brief HeFFTe example 5, computing the Cosine Transform DCT using an arbitrary number of MPI ranks.
 *
 * Performing the discrete Cosine Transform (DCT) on three dimensional data in a box of 64 by 64 by 64
 * split across an arbitrary number of MPI ranks.
 */
void compute_dct(MPI_Comm comm){

    int me; // this process rank within the comm
    MPI_Comm_rank(comm, &me);

    int num_ranks; // total number of ranks in the comm
    MPI_Comm_size(comm, &num_ranks);

    // Using input configuration with pencil data format in X direction
    // and output configuration with pencil data in the Z direction.
    // This format uses only two internal reshape operation.
    std::array<int, 2> proc_grid = heffte::make_procgrid(num_ranks);
    std::array<int, 3> input_grid = {1, proc_grid[0], proc_grid[1]};
    std::array<int, 3> output_grid = {proc_grid[0], proc_grid[1], 1};

    // Describe all the indexes across all ranks
    heffte::box3d<> const world = {{0, 0, 0}, {63, 63, 63}};

    // Split the world box into a 2D grid of boxes
    std::vector<heffte::box3d<>> inboxes  = heffte::split_world(world, input_grid);
    std::vector<heffte::box3d<>> outboxes = heffte::split_world(world, output_grid);

    // Select the backend to use, prefer FFTW and fallback to the stock backend
    // The real-to-real transforms have _cos and _sin appended
    #ifdef Heffte_ENABLE_FFTW
    using backend_tag = heffte::backend::fftw_cos;
    #else
    using backend_tag = heffte::backend::stock_cos;
    #endif

    // define the heffte class and the input and output geometry
    // note that rtransform is just an alias to fft3d
    heffte::rtransform<backend_tag> tcos(inboxes[me], outboxes[me], comm);

    // vectors with the correct sizes to store the input and output data
    // taking the size of the input and output boxes
    std::vector<double> input(tcos.size_inbox());
    std::vector<double> output(tcos.size_outbox());

    // the workspace vector is of a real type too
    std::vector<double> workspace(tcos.size_workspace());

    // fill the input with random data
    // see the FFTW example 1 on how to properly distribute the data to the boxes
    for(size_t i=0; i<input.size(); i++)
        input[i] = static_cast<double>(i+1);

    // perform a forward DCT
    tcos.forward(input.data(), output.data());

    // compute the inverse or backward transform
    std::vector<double> inverse(tcos.size_inbox());
    tcos.backward(output.data(), inverse.data(), workspace.data(), heffte::scale::full);

    double err = 0.0;
    for(size_t i=0; i<inverse.size(); i++)
        err = std::max(err, std::abs(inverse[i] - input[i]));

    // print the error for each MPI rank
    std::cout << std::scientific;
    for(int i=0; i<num_ranks; i++){
        if (me == i) std::cout << "rank " << i << " error: " << err << std::endl;
        MPI_Barrier(comm);
    }
}

int main(int argc, char** argv){

    MPI_Init(&argc, &argv);

    compute_dct(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
