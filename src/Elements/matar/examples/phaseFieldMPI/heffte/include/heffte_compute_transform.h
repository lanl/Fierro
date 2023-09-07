/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_COMPUTE_TRANSFORM_H
#define HEFFTE_COMPUTE_TRANSFORM_H

#include "heffte_reshape3d.h"

namespace heffte {

    /*!
     * \internal
     * \ingroup fft3d
     * \brief Implements a series of reshape and executor operations.
     *
     * The generic template is explicitly instantiated to work with the standard types,
     * float, double, std::complex<float> and std::complex<double>, e.g.,
     * no custom complex numbers are allowed here.
     * The algorithm implements a series of reshape and executor operations following
     * the logic reshape-executor-reshape-executor-reshape-executor.
     * Any of the reshape/executor operations can be skipped if the corresponding
     * pointers are not set. The algorithm is required to avoid any extraneous
     * data movement or operations.
     *
     * This variant works with the cases when the input and output data-types are identical,
     * e.g., complex-to-complex or real-to-complex transforms.
     *
     * \tparam location_tag must be either heffte::tag::cpu or heffte::tag::gpu
     *                      the tag and the corresponding stream are needed for
     *                      the rare occasion when a copy cannot be avoided,
     *                      e.g., out-of-place transform with only one executor and no reshapes.
     * \tparam index is either int or long long, indicates the index of the box used to
     *               make the reshape operators.
     * \tparam scalar_type is either float, double, std::complex<float> or std::complex<double>,
     *
     * \param stream is the device stream, e.g., cudaStream_t, sycl::queue& or void* (cpu case)
     * \param input is the input for the forward or backward transform
     * \param output is the output for the forward or backward transform
     * \param workspace is the pre-allocated buffer with sufficient size to accommodate all operations,
     * \param shaper are the four stages of the reshape operations
     * \param executor holds the three stages of the one dimensional FFT algorithm
     * \param dir indicates whether to use the forward or backward method of the executor
     * \endinternal
     */
    template<typename location_tag, typename index, typename scalar_type>
    void compute_transform(typename backend::data_manipulator<location_tag>::stream_type stream,
                        int const batch_size,
                        scalar_type const input[], scalar_type output[], scalar_type workspace[],
                        size_t executor_buffer_offset, size_t size_comm_buffers,
                        std::array<std::unique_ptr<reshape3d_base<index>>, 4> const &shaper,
                        std::array<executor_base*, 3> const &executor,
                        direction dir);
    /*!
     * \internal
     * \ingroup fft3d
     * \brief Implements a series of reshape and executor operations, forward real-to-complex case.
     *
     * Identical to compute_transform() but applies to the real-to-complex forward transform.
     * The difference is that the \b scalar_type is float or double and the output is std::complex<scalar_type>.
     * The first executor must be set and it will be called with the method:
     * \code
     *      executor_base::forward(scalar_type input[], std::complex<scalar_type> output[],
     *                             std::complex<scalar_type> *workspace);
     * \endcode
     * \endinternal
     *
     * The \b direction parameter is ignored.
     */
    template<typename location_tag, typename index, typename scalar_type>
    void compute_transform(typename backend::data_manipulator<location_tag>::stream_type stream,
                        int const batch_size,
                        scalar_type const input[], std::complex<scalar_type> output[],
                        std::complex<scalar_type> workspace[],
                        size_t executor_buffer_offset, size_t size_comm_buffers,
                        std::array<std::unique_ptr<reshape3d_base<index>>, 4> const &shaper,
                        std::array<executor_base*, 3> const &executor, direction);
    /*!
     * \internal
     * \ingroup fft3d
     * \brief Implements a series of reshape and executor operations, forward complex-to-real case.
     *
     * Identical to compute_transform() but applies to the complex-to-real backward transform.
     * The difference is that the \b scalar_type is float or double and the output is std::complex<scalar_type>.
     * The first executor must be set and it will be called with the method:
     * \code
     *      executor_base::forward(std::complex<scalar_type> input[], scalar_type output[],
     *                             std::complex<scalar_type> *workspace);
     * \endcode
     *
     * The \b direction parameter is ignored.
     * \endinternal
     */
    template<typename location_tag, typename index, typename scalar_type>
    void compute_transform(typename backend::data_manipulator<location_tag>::stream_type stream,
                        int const batch_size,
                        std::complex<scalar_type> const input[], scalar_type output[],
                        std::complex<scalar_type> workspace[],
                        size_t executor_buffer_offset, size_t size_comm_buffers,
                        std::array<std::unique_ptr<reshape3d_base<index>>, 4> const &shaper,
                        std::array<executor_base*, 3> const &executor, direction);

}

#endif // HEFFTE_COMPUTE_TRANSFORM_H
