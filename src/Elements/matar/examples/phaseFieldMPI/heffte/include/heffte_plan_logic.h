/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_PLAN_LOGIC_H
#define HEFFTE_PLAN_LOGIC_H

#include "heffte_common.h"

/*!
 * \ingroup fft3d
 * \addtogroup fft3dplan Plan transformation logic
 *
 * Implements the analysis of the input and output distribution of the boxes
 * and creates the corresponding plan for reshape and 1-D FFT operations.
 */

namespace heffte {

/*!
 * \ingroup fft3d
 * \brief Defines list of potential communication algorithms.
 *
 * Depending on the size of the data and the number of MPI ranks used in the FFT transform,
 * the problems can be classified as either bandwidth-bound or latency-bound.
 * The bandwidth-bound case hits pretty close to the maximum throughput of the MPI interconnect
 * while the latency-bound case is more affected by the latency of the large number of small communications.
 * As a short-hand we can call these small-problems (latency-bound) or large-problems (bandwidth-bound),
 * although the specific cutoff point is dependent on the backend (and the version of the backend),
 * the version of MPI, the machine interconnect, and the specific optimizations that have been implemented in MPI.
 *
 * There is a plan of adding an auto-tuning framework in heFFTe to help users select the best
 * possible set of options; however, currently the users have to manually find the best option for their hardware.
 * The expected "best" algorithm is:
 * \code
 *      reshape_algorithm::alltoallv          : for larger FFT, many MPI ranks
 *      reshape_algorithm::alltoall           : for smaller FFT, many MPI ranks
 *      reshape_algorithm::p2p_plined         : for larger FFT, fewer MPI ranks
 *      reshape_algorithm::p2p                : for smaller FFT, fewer MPI ranks
 * \endcode
 *
 * Note that in the GPU case, the above algorithms are also affected by the GPU latency
 * if MPI calls are made directly from the GPU. This can be controlled with the use_gpu_aware
 * variable of the heffte::plan_options.
 */
enum class reshape_algorithm{
    //! \brief Using the MPI_Alltoallv options, no padding on the data (default option).
    alltoallv = 0,
    //! \brief Using the MPI_Alltoall options, with padding on the data.
    alltoall = 3,
    //! \brief Using MPI_Isend and MPI_Irecv, all sending receiving packing and unpacking are pipelined.
    p2p_plined = 1,
    //! \brief Using MPI_Send and MPI_Irecv, receive is pipelined with packing and sending.
    p2p = 2
};

/*!
 * \ingroup fft3d
 * \brief Defines a set of tweaks and options to use in the plan generation.
 *
 * Example usage:
 * \code
 *  heffte::plan_options options = heffte::default_options<heffte::backend::fftw>();
 *  options.algorithm = reshape_algorithm::p2p; // forces the use of point-to-point communication
 *  heffte::fft3d<heffte::backend::fftw> fft3d(inbox, outbox, comm, options);
 * \endcode
 *
 * \par Option use_reorder
 * Controls whether the backends should be called with strided or contiguous data.
 * If the option is enabled then during the reshape operations heFFTe will reorder
 * the data so that the backend is called for contiguous batch of 1D FFTs.
 * Otherwise the strided call will be performed. Depending on the size and
 * the specific backend (or version of the backend), one or the other may improve performance.
 * The reorder is applied during the unpacking stage of an MPI communication and
 * will be applied even if no MPI communication is used.
 * Note that some backends don't currently support strided transforms, e.g.,
 * the Sine and Cosine transforms, in which case this option will have no effect.
 *
 * \par Option algorithm
 * Specifies the combination of MPI calls to use in the communication.
 * See `heffte::reshape_algorithm` for details.
 *
 * \par Option use_pencils
 * Indicates whether the intermediate steps of the computation should be done
 * either in pencil or slab format. Slabs work better for problems with fewer
 * MPI ranks, while pencils work better when the number ranks increases.
 * The specific cutoff depends on the hardware and MPI implementation.
 * Note that is the input or output shape of the data is in slab format,
 * then this option will be ignored.
 *
 * \par Option use_gpu_aware
 * Applied only when using one of the GPU backends, indicates whether MPI communication
 * should be initiated from the GPU device or if the data has to be moved the CPU first.
 * MPI calls from the GPU have faster throughput but larger latency, thus initiating
 * the calls from the CPU (e.g., setting use_gpu_aware to false) can be faster
 * when using smaller problems compared to the number of MPI ranks.
 *
 * \par Option use_subcomm or use_num_subranks
 * Restricts the intermediate reshape and FFT operations to a subset of the ranks
 * specified by the communicator given in the construction of heffte::fft3d and heffte::fft3d_r2c.
 * By default, heFFTe will use all of the available MPI ranks but this is not always optimal
 * (see the two examples below).
 * The other options are defined as member variables, but the subcomm option is specified
 * with member functions that accept either an integer or an MPI communicator.
 * Using an integer will specify ranks \b 0 to \b num_subranks -1, while using a communicator
 * can define an arbitrary subset. MPI ranks that don't belong to the subcomm should pass MPI_COMM_NULL.
 * The plan_options class will hold a non-owning reference to the MPI subcomm
 * but heffte::fft3d and heffte::fft3d_r2c will use the subcomm only in the constructors,
 * i.e., the subcomm can be safely discarded/freed after the fft3d classes are constructed.
 *
 * \par
 * For example, if the input and output shapes of the data do not form pencils in any direction
 * (i.e., using a brick decomposition),
 * then heFFTe has to perform 4 reshape operations (3 if using slabs) and if the problem size
 * is small relative to the number of ranks this results in 4 sets of small messages which increases
 * the latency and reduces performance.
 * However, if the problem can fit on a single node (e.g., single GPU), then gathering all the data
 * to a single rank then performing a single 3D FFT and scattering the data back will result
 * in only 2 communications involving the small messages. Thus, so long as the two operations
 * are less expensive than the 4, using a subcomm will result in an overall performance boost.
 *
 * \par
 * Similar to the previous example, if we are using a CPU backend with multiple MPI ranks per node,
 * then reducing the MPI ranks to one-per-node can effectively coalesce smaller messages from multiple
 * ranks to larger messages and thus reduce latency. If the CPU backend supports multi-threading,
 * then all CPU cores can still be used by calls from the single rank without reduction of performance.
 *
 */
struct plan_options{
    //! \brief Constructor, initializes all options with the default values for the given backend tag.
    template<typename backend_tag> plan_options(backend_tag const)
        : use_reorder(default_plan_options<backend_tag>::use_reorder),
          algorithm(reshape_algorithm::alltoallv),
          use_pencils(true),
          use_gpu_aware(true),
          num_sub(-1),
          subcomm(MPI_COMM_NULL)
    {}
    //! \brief Constructor, initializes each variable, primarily for internal use.
    plan_options(bool reorder, reshape_algorithm alg, bool pencils)
        : use_reorder(reorder), algorithm(alg), use_pencils(pencils), use_gpu_aware(true), num_sub(-1), subcomm(MPI_COMM_NULL)
    {}
    //! \brief Defines whether to transpose the data on reshape or to use strided 1-D ffts.
    bool use_reorder;
    //! \brief Defines the communication algorithm.
    reshape_algorithm algorithm;
    //! \brief Defines whether to use pencil or slab data distribution in the reshape steps.
    bool use_pencils;
    //! \brief Defines whether to use MPI calls directly from the GPU or to move to the CPU first.
    bool use_gpu_aware;
    //! \brief Defines the number of ranks to use for the internal reshapes, set to -1 to use all ranks.
    void use_num_subranks(int num_subranks){ num_sub = num_subranks; }
    /*!
     * \brief Set sub-communicator to use in the intermediate reshape operations.
     *
     * The ranks defined by \b comm must be a subset of the communicator that will be used
     * in the future call to heffte::fft3d or heffte::fft3d_r2c.
     * The ranks that are not associated with the comm should pass in MPI_COMM_NULL.
     * The plan_options object will take a non-owning reference to \b comm
     * but the reference will not be passed into heffte::fft3d or heffte::fft3d_r2c.
     *
     * This method takes precedence over use_num_subranks() if both methods are called.
     * Avoid calling both methods.
     */
    void use_subcomm(MPI_Comm comm){
        num_sub = 1;
        subcomm = comm;
    }
    //! \brief Return the set number of sub-ranks.
    int get_subranks() const{ return num_sub; }
private:
    int num_sub;
    MPI_Comm subcomm;
};

/*!
 * \ingroup fft3d
 * \brief Simple I/O for the plan options struct.
 */
inline std::ostream & operator << (std::ostream &os, plan_options const options){
    std::string algorithm = "";
    switch (options.algorithm){
        case reshape_algorithm::alltoallv  : algorithm = "mpi:alltoallv"; break;
        case reshape_algorithm::alltoall   : algorithm = "mpi:alltoall"; break;
        case reshape_algorithm::p2p_plined : algorithm = "mpi:point-to-point-pipelined"; break;
        case reshape_algorithm::p2p        : algorithm = "mpi:point-to-point"; break;
    };
    os << "options = ("
       << ((options.use_reorder) ? "fft1d:contiguous" : "fft1d:strided") << ", "
       << algorithm << ", "
       << ((options.use_pencils) ? "decomposition:pencil" : "decomposition:slab") << ", "
       << ((options.use_gpu_aware) ? "mpi:from-gpu" : "mpi:from-cpu") << ")";
    return os;
}

/*!
 * \ingroup fft3d
 * \brief Adjusts the user provided options to what can be handled by the backend.
 *
 * Some backends do not support all available options, e.g., they require the use_reorder
 * option to be set on. This template makes the necessary adjustments so that the correct
 * answer is always computed even if the user provides unsupported options.
 */
template<typename backend_tag, bool use_r2c = false>
plan_options set_options(plan_options opts){
    if (std::is_same<backend_tag, backend::stock_cos>::value
        or std::is_same<backend_tag, backend::mkl_cos>::value
        or std::is_same<backend_tag, backend::cufft_cos>::value
        or std::is_same<backend_tag, backend::rocfft_cos>::value
        or std::is_same<backend_tag, backend::onemkl_cos>::value
        or std::is_same<backend_tag, backend::stock_sin>::value
        or std::is_same<backend_tag, backend::mkl_sin>::value
        or std::is_same<backend_tag, backend::cufft_sin>::value
        or std::is_same<backend_tag, backend::rocfft_sin>::value
        or std::is_same<backend_tag, backend::onemkl_sin>::value
    ){
        // currently the cosine options work only with reorder.
        opts.use_reorder = true;
        return opts;
    }else if (use_r2c and std::is_same<backend_tag, backend::rocfft>::value){
        // the rocfft backend with r2c requires the reorder (problem with the strides)
        opts.use_reorder = true;
        return opts;
    }else{
        return opts; // all options are supported for this backend
    }
}

/*!
 * \ingroup heffterocm
 * \brief Forces the reorder logic for the ROCM r2c variant.
 */
inline plan_options force_reorder(plan_options opts){
    opts.use_reorder = true;
    return opts;
}

/*!
 * \ingroup fft3d
 * \brief Returns the default backend options associated with the given backend.
 */
template<typename backend_tag>
plan_options default_options(){
    return plan_options(backend_tag());
}

/*!
 * \ingroup fft3dplan
 * \brief The logic plan incorporates the order and types of operations in a transform.
 *
 * The logic_plan is used to separate the logic of the order of basic operations (reshape or fft execute)
 * from the constructor of the fft3d and fft3d_r2c classes.
 * In this manner, detection of pencils vs. brick distribution of the data and/or making decisions regarding
 * the transposition of indexing can be done in sufficiently robust and complex logic without
 * clutter of the main classes or unnecessary repetition of code.
 *
 * Node that a reshape operation \b i will be performed only if in_shape[i] and out_shape[i] are different.
 *
 * Specifically:
 * - in_shape[0] is the set of input boxes
 * - out_shape[0] is the geometry to be used for the fist 1-D FFT operation which will happen in fft_direction[0]
 * - in_shape[1] is either out_shape[0] or, in the r2c case, the boxes with reduced dimension
 * - in_shape[i] -> out_shape[i] for i = 1, 2, will hold the intermediate shapes
 * - in_shape[i] may be equal to out_shape[i] for any i = 0, 1, 2, 3
 * - 1-D FFT transforms will be applied to out_shape[1] and out_shape[2] in directions fft_direction[0] and fft_direction[1]
 * - out_shape[3] is the set of output boxes
 * - index_count is the produce of all indexes, used in scaling
 * - options will hold a copy of the set of options used in the construction
 */
template<typename index>
struct logic_plan3d{
    //! \brief Holds the input shapes for the 4 forward reshapes (backwards reverses in and out).
    std::vector<box3d<index>> in_shape[4];
    //! \brief Holds the output shapes for the 4 forward reshapes (backwards reverses in and out).
    std::vector<box3d<index>> out_shape[4];
    //! \brief Sizes for the 1-D transforms.
    std::array<index, 3> fft_sizes;
    //! \brief Direction of the 1-D FFT transforms.
    std::array<int, 3> fft_direction;
    //! \brief The total number of indexes in all directions.
    long long index_count;
    //! \brief Extra options used in the plan creation.
    plan_options const options;
    //! \brief MPI rank used in the plan creation.
    int const mpi_rank;
};

/*!
 * \ingroup fft3dplan
 * \brief Returns true for each direction where the boxes form pencils (i.e., where the size matches the world size).
 */
template<typename index>
inline std::array<bool, 3> pencil_directions(box3d<index> const world, std::vector<box3d<index>> const &boxes){
    std::array<bool, 3> is_pencil = {true, true, true};
    for(auto const &b : boxes){
        for(int i=0; i<3; i++)
            is_pencil[i] = is_pencil[i] and (world.size[i] == b.size[i]);
    }
    return is_pencil;
}

/*!
 * \ingroup fft3dplan
 * \brief Creates the logic plan with the provided user input.
 *
 * \param boxes is the current distribution of the data across the MPI comm
 * \param r2c_direction is the direction is the direction of shrinking of the data for an r2c transform
 *              the c2c case should use -1
 * \param options is a set of heffte::plan_options to use
 *
 * \returns the plan for reshape and 1-D fft transformations
 */
template<typename index>
logic_plan3d<index> plan_operations(ioboxes<index> const &boxes, int r2c_direction, plan_options const options, int const mpi_rank);

/*!
 * \ingroup fft3dplan
 * \brief Assuming the shapes in the plan form grids, reverse engineer the grid dimensions (used in the benchmark).
 */
template<typename index>
std::vector<std::array<int, 3>> compute_grids(logic_plan3d<index> const &plan);

}

#endif
