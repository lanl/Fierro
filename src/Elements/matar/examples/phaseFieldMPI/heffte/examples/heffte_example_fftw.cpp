#include "heffte.h"

/*!
 * \brief HeFFTe example 1, simple DFT using two MPI ranks and FFTW backend.
 *
 * Performing DFT on three dimensional data in a box of 4 by 4 by 4 split
 * across the third dimension between two MPI ranks.
 */
void compute_dft(MPI_Comm comm){

    int me; // this process rank within the comm
    MPI_Comm_rank(comm, &me);

    int num_ranks; // total number of ranks in the comm
    MPI_Comm_size(comm, &num_ranks);

    if (num_ranks != 2){
        if (me == 0) std::cout << " heffte_example_fftw is set to 2 mpi ranks, exiting \n";
        return;
    }

    // define the domain split between two boxes
    heffte::box3d<> const left_box  = {{0, 0, 0}, {3, 3, 1}};
    heffte::box3d<> const right_box = {{0, 0, 2}, {3, 3, 3}};

    // the box associated with this MPI rank
    heffte::box3d<> const my_box = (me == 0) ? left_box : right_box;

    // define the heffte class and the input and output geometry
    heffte::fft3d<heffte::backend::fftw> fft(my_box, my_box, comm);

    // vectors with the correct sizes to store the input and output data
    // taking the size of the input and output boxes
    std::vector<std::complex<double>> input(fft.size_inbox());
    std::vector<std::complex<double>> output(fft.size_outbox());

    // form a global input vector and copy the data to the local input
    // this serves as an example of the relation between global and local indexes
    std::vector<std::complex<double>> world_input(4 * 4 * 4);
    std::iota(world_input.begin(), world_input.end(), 0); // fills with 0, 1, 2, 3, ...

    // set the strides for the triple indexes
    int world_plane = 4 * 4;
    int world_stride = 4;
    int local_plane = my_box.size[0] * my_box.size[1];
    int local_stride = my_box.size[0];
    // note the order of the loops corresponding to the default order (0, 1, 2)
    // order (0, 1, 2) means that the data in dimension 0 is contiguous
    for(int i=my_box.low[2]; i <= my_box.high[2]; i++)
        for(int j=my_box.low[1]; j <= my_box.high[1]; j++)
            for(int k=my_box.low[0]; k <= my_box.high[0]; k++)
                input[(i - my_box.low[2]) * local_plane
                      + (j - my_box.low[1]) * local_stride + k - my_box.low[0]]
                    = world_input[i * world_plane + j * world_stride + k];

    // perform a forward DFT
    fft.forward(input.data(), output.data());

    // check the accuracy
    // Compute an FFT using only a single rank
    MPI_Comm single_rank_comm;
    MPI_Comm_split(comm, me, me, &single_rank_comm);
    std::vector<std::complex<double>> world_output(world_input.size());
    heffte::box3d<> const world_box  = {{0, 0, 0}, {3, 3, 3}};
    // create a heFFTe object associated with the single rank and perform forward transform
    heffte::fft3d<heffte::backend::fftw>(world_box, world_box, single_rank_comm)
        .forward(world_input.data(), world_output.data());

    // check the difference between the local and global outputs
    // using the same indexing scheme as when assigning the inputs
    double err = 0.0;
    for(int i=my_box.low[2]; i <= my_box.high[2]; i++)
        for(int j=my_box.low[1]; j <= my_box.high[1]; j++)
            for(int k=my_box.low[0]; k <= my_box.high[0]; k++)
                err = std::max(err,
                               std::abs(output[(i - my_box.low[2]) * local_plane
                                                + (j - my_box.low[1]) * local_stride + k - my_box.low[0]]
                                        - world_output[i * world_plane + j * world_stride + k]));

    // print the error for each MPI rank
    std::cout << std::scientific;
    if (me == 0) std::cout << "rank 0 forward error: " << err << std::endl;
    MPI_Barrier(comm);
    if (me == 1) std::cout << "rank 1 forward error: " << err << std::endl;

    // reset the input to zero
    std::fill(input.begin(), input.end(), std::complex<double>(0.0, 0.0));

    // perform a backward DFT
    fft.backward(output.data(), input.data());

    // rescale the result
    // alternatively: fft.backward(output.data(), input.data(), heffte::scale::full);
    for(auto &i : input) i /= 64.0;

    // compare the computed entries to the original input data
    err = 0.0; // reset the error
    for(int i=my_box.low[2]; i <= my_box.high[2]; i++)
        for(int j=my_box.low[1]; j <= my_box.high[1]; j++)
            for(int k=my_box.low[0]; k <= my_box.high[0]; k++)
                err = std::max(err,
                               std::abs(input[(i - my_box.low[2]) * local_plane
                                              + (j - my_box.low[1]) * local_stride + k - my_box.low[0]]
                                        - world_input[i * world_plane + j * world_stride + k]));

    // print the error for each MPI rank
    if (me == 0) std::cout << "rank 0 backward error: " << err << std::endl;
    MPI_Barrier(comm);
    if (me == 1) std::cout << "rank 1 backward error: " << err << std::endl;

    MPI_Comm_free(&single_rank_comm);
}

int main(int argc, char** argv){

    MPI_Init(&argc, &argv);

    compute_dft(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
