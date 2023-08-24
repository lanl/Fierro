#pragma once

#include <mpi.h>
#include "heffte.h"
#include <array>

/**************************************************
    FFTBase
***************************************************
*/

template <typename HEFFTE_BACKEND, typename R>
class FFTBase
{
public:
    MPI_Comm comm;
    const int root = 0;
    int my_rank;
    int num_ranks;
    std::array<int,3> globalRealBoxSize;
    std::array<int,3> globalComplexBoxSize;
    heffte::box3d<> globalRealBox;
    heffte::box3d<> globalComplexBox;
    std::array<int, 3> procGrid;
    std::vector<heffte::box3d<>> localRealBoxes;
    std::vector<heffte::box3d<>> localComplexBoxes;
    heffte::box3d<> myRealBox;
    heffte::box3d<> myComplexBox;
    heffte::plan_options options;
    std::vector<std::array<int,3>> localRealBoxSizes;
    std::vector<std::array<int,3>> localComplexBoxSizes;
    MPI_Datatype mpiType;
    std::vector<MPI_Datatype> mpiSubarrayTypes;

    FFTBase(MPI_Comm comm_, const std::array<int,3> & globalRealBoxSize_, const std::array<int,3> & globalComplexBoxSize_);
    virtual ~FFTBase();
    virtual void forward(const R *input, std::complex<R> *output) = 0;
    virtual void forward(const R *input, R *output) = 0;
    virtual void backward(const std::complex<R> *input, R *output) = 0;
    virtual void backward(const R *input, R *output) = 0;
};

template <typename HEFFTE_BACKEND, typename R>
FFTBase<HEFFTE_BACKEND,R>::FFTBase(MPI_Comm comm_, const std::array<int,3> & globalRealBoxSize_, const std::array<int,3> & globalComplexBoxSize_) :
comm(comm_),
my_rank(heffte::mpi::comm_rank(comm)),
num_ranks(heffte::mpi::comm_size(comm)),
globalRealBoxSize(globalRealBoxSize_),
globalComplexBoxSize(globalComplexBoxSize_),
globalRealBox({0, 0, 0}, {globalRealBoxSize[0]-1, globalRealBoxSize[1]-1, globalRealBoxSize[2]-1}),
globalComplexBox({0, 0, 0}, {globalComplexBoxSize[0]-1, globalComplexBoxSize[1]-1, globalComplexBoxSize[2]-1}),
procGrid(heffte::proc_setup_min_surface(globalRealBox, num_ranks)),
localRealBoxes(heffte::split_world(globalRealBox, procGrid)),
localComplexBoxes(heffte::split_world(globalComplexBox, procGrid)),
myRealBox(localRealBoxes[my_rank]),
myComplexBox(localComplexBoxes[my_rank]),
options(heffte::default_options<HEFFTE_BACKEND>()),
localRealBoxSizes(num_ranks),
localComplexBoxSizes(num_ranks),
mpiType(MPI_DATATYPE_NULL),
mpiSubarrayTypes(num_ranks, MPI_DATATYPE_NULL)
{
    // use strided 1-D FFT operations
    // some backends work just as well when the entries of the data are not contiguous
    // then there is no need to reorder the data in the intermediate stages which saves time
    options.use_reorder = false;

    // use point-to-point communications
    // collaborative all-to-all and individual point-to-point communications are two alternatives
    // one may be better than the other depending on
    // the version of MPI, the hardware interconnect, and the problem size
    //options.use_alltoall = true;

    // in the intermediate steps, the data can be shapes as either 2-D slabs or 1-D pencils
    // for sufficiently large problem, it is expected that the pencil decomposition is better
    // but for smaller problems, the slabs may perform better (depending on hardware and backend)
    options.use_pencils = true;

    // calculate sizes of real and complex domains(boxes) for each rank
    for (size_t i = 0; i < num_ranks; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            localRealBoxSizes[i][j]  = localRealBoxes[i].high[j] - localRealBoxes[i].low[j] + 1;
            localComplexBoxSizes[i][j]  = localComplexBoxes[i].high[j] - localComplexBoxes[i].low[j] + 1;
        }
    }

    // set correct mpiType dependent on typename R
    if constexpr (std::is_same<R, float>::value)
    {
        mpiType = MPI_FLOAT;
    }
    else if constexpr (std::is_same<R, double>::value)
    {
        mpiType = MPI_DOUBLE;
    }
    else if constexpr (std::is_same<R, long double>::value)
    {
        mpiType = MPI_LONG_DOUBLE;
    }
    else
    {
        throw std::runtime_error("mpiType can not be deduced in FFTBase. Please add the mpi type.\n");
    }

    // create MPI subarray type for each rank
    for (size_t i = 0; i < num_ranks; i++)
    {
        int dimensions_full_array[3] = { globalRealBoxSize[2], globalRealBoxSize[1], globalRealBoxSize[0] };
        int dimensions_subarray[3] = { localRealBoxSizes[i][2], localRealBoxSizes[i][1], localRealBoxSizes[i][0] };
        int start_coordinates[3] = { localRealBoxes[i].low[2], localRealBoxes[i].low[1], localRealBoxes[i].low[0] };
        MPI_Type_create_subarray(3, dimensions_full_array, dimensions_subarray, start_coordinates, MPI_ORDER_C, mpiType, &mpiSubarrayTypes[i]);
        MPI_Type_commit(&mpiSubarrayTypes[i]);
    }

#if defined(CUDA) || defined(HIP)
    // on a multi-gpu system, distribute the devices across the mpi ranks
    if (heffte::gpu::device_count() > 1) {
        heffte::gpu::device_set(heffte::mpi::comm_rank(comm) % heffte::gpu::device_count());
    }
#endif
}

template <typename HEFFTE_BACKEND, typename R>
FFTBase<HEFFTE_BACKEND,R>::~FFTBase()
{
    for (size_t i = 0; i < num_ranks; i++)
    {
        MPI_Type_free(&mpiSubarrayTypes[i]);
    }
}

/**************************************************
    FFT3D_R2C
***************************************************
*/
template <typename HEFFTE_BACKEND, typename R>
class FFT3D_R2C : public FFTBase<HEFFTE_BACKEND,R>
{
public:
    int r2c_direction;
    heffte::fft3d_r2c<HEFFTE_BACKEND> fft; // heffte class for performing the fft
    typename heffte::fft3d<HEFFTE_BACKEND>::template buffer_container<std::complex<R>> workspace;

    FFT3D_R2C(MPI_Comm comm, const std::array<int,3> & globalRealBoxSize);
    ~FFT3D_R2C();
    void forward(const R *input, std::complex<R> *output) override;
    void forward(const R *input, R *output) override;
    void backward(const std::complex<R> *input, R *output) override;
    void backward(const R *input, R *output) override;
};

template <typename HEFFTE_BACKEND, typename R>
FFT3D_R2C<HEFFTE_BACKEND,R>::FFT3D_R2C(MPI_Comm comm, const std::array<int,3> & globalRealBoxSize) :
FFTBase<HEFFTE_BACKEND,R>(comm, globalRealBoxSize, {globalRealBoxSize[0]/2+1, globalRealBoxSize[1], globalRealBoxSize[2]}),
r2c_direction(0),
fft(this->myRealBox, this->myComplexBox, r2c_direction, this->comm, this->options),
workspace(fft.size_workspace())
{
    // check if the complex indexes have correct dimension
    assert(this->globalRealBox.r2c(r2c_direction) == this->globalComplexBox);
}

template <typename HEFFTE_BACKEND, typename R>
void FFT3D_R2C<HEFFTE_BACKEND,R>::forward(const R *input, std::complex<R> *output)
{
    fft.forward(input, output, workspace.data());
}

template <typename HEFFTE_BACKEND, typename R>
void FFT3D_R2C<HEFFTE_BACKEND,R>::forward(const R *input, R *output)
{
    fft.forward(input, (std::complex<R>*)output, workspace.data());
}

template <typename HEFFTE_BACKEND, typename R>
void FFT3D_R2C<HEFFTE_BACKEND,R>::backward(const std::complex<R> *input, R *output)
{
    fft.backward(input, output, workspace.data(), heffte::scale::full);
}

template <typename HEFFTE_BACKEND, typename R>
void FFT3D_R2C<HEFFTE_BACKEND,R>::backward(const R *input, R *output)
{
    fft.backward((std::complex<R>*)input, output, workspace.data(), heffte::scale::full);
}

template <typename HEFFTE_BACKEND, typename R>
FFT3D_R2C<HEFFTE_BACKEND,R>::~FFT3D_R2C()
{
}

/**************************************************
    FFT3D
***************************************************
*/
template <typename HEFFTE_BACKEND, typename R>
class FFT3D : public FFTBase<HEFFTE_BACKEND,R>
{
public:
    heffte::fft3d<HEFFTE_BACKEND> fft; // heffte class for performing the fft
    typename heffte::fft3d<HEFFTE_BACKEND>::template buffer_container<std::complex<R>> workspace;

    FFT3D(MPI_Comm comm, const std::array<int,3> & globalRealBoxSize);
    ~FFT3D();
    void forward(const R *input, std::complex<R> *output) override;
    void forward(const R *input, R *output) override;
    void backward(const std::complex<R> *input, R *output) override;
    void backward(const R *input, R *output) override;
};

template <typename HEFFTE_BACKEND, typename R>
FFT3D<HEFFTE_BACKEND,R>::FFT3D(MPI_Comm comm, const std::array<int,3> & globalRealBoxSize) :
FFTBase<HEFFTE_BACKEND,R>(comm, globalRealBoxSize, globalRealBoxSize),
fft(this->myRealBox, this->myComplexBox, this->comm, this->options),
workspace(fft.size_workspace())
{
}

template <typename HEFFTE_BACKEND, typename R>
void FFT3D<HEFFTE_BACKEND,R>::forward(const R *input, std::complex<R> *output)
{
    fft.forward(input, output, workspace.data());
}

template <typename HEFFTE_BACKEND, typename R>
void FFT3D<HEFFTE_BACKEND,R>::forward(const R *input, R *output)
{
    fft.forward(input, (std::complex<R>*)output, workspace.data());
}

template <typename HEFFTE_BACKEND, typename R>
void FFT3D<HEFFTE_BACKEND,R>::backward(const std::complex<R> *input, R *output)
{
    fft.backward(input, output, workspace.data(), heffte::scale::full);
}

template <typename HEFFTE_BACKEND, typename R>
void FFT3D<HEFFTE_BACKEND,R>::backward(const R *input, R *output)
{
    fft.backward((std::complex<R>*)input, output, workspace.data(), heffte::scale::full);
}

template <typename HEFFTE_BACKEND, typename R>
FFT3D<HEFFTE_BACKEND,R>::~FFT3D()
{
}
