/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_FFT3D_R2C_H
#define HEFFTE_FFT3D_R2C_H

#include "heffte_fft3d.h"

namespace heffte {

/*!
 * \ingroup fft3d
 * \brief Similar to heffte::fft3d, but computed fewer redundant coefficients when the input is real.
 *
 * \par Overview
 * Given real input data, there is no unambiguous way to distinguish between the positive and negative
 * direction in the complex plane; therefore, by an argument of symmetry, all complex output must come
 * in conjugate pairs. The heffte::fft3d computes both numbers for each conjugate pair,
 * this class aims at computing fewer redundant coefficients and thus reducing both flops and data movement.
 * This is achieved by selecting one of the three dimensions and the data is shortened in that dimensions
 * to contain only the unique (non-conjugate) coefficients.
 *
 * \par Boxes and Data Distribution
 * Similar to heffte::fft3d the data is organized in boxes using the heffte::box3d structs;
 * however, in the real-to-complex case the global input and output domains do not match.
 * If the original data sits in a box {0, 0, 0}, {x, y, z}, then depending on the dimensions
 * chosen for the shortening, the output data will form the box:
 * \code
 *  {{0, 0, 0}, {x/2 + 1, y,       z}}        // if chosen dimension 0
 *  {{0, 0, 0}, {x,       y/2 + 1, z}}        // if chosen dimension 1
 *  {{0, 0, 0}, {x,       y,       z/2 + 1}}  // if chosen dimension 2
 * // note that x/2 indicates the standard C++ integer division
 * \endcode
 * Thus, the union of the inboxes across all MPI ranks must add up to the global input box,
 * and the union of the outboxes must add up to the shortened global box.
 *
 * \par Compatible Types
 * The real-to-complex variant does not support the cases when the input is complex,
 * the supported types are the ones with real input in
 * \ref HeffteFFT3DCompatibleTypes "the table of compatible types".
 */
template<typename backend_tag, typename index = int>
class fft3d_r2c : public backend::device_instance<typename backend::buffer_traits<backend_tag>::location>{
public:
    //! \brief FFT executor for the complex-to-complex dimensions.
    using backend_executor_c2c = typename one_dim_backend<backend_tag>::executor;
    //! \brief FFT executor for the real-to-complex dimension.
    using backend_executor_r2c = typename one_dim_backend<backend_tag>::executor_r2c;
    /*!
     * \brief Type-tag that is either tag::cpu or tag::gpu to indicate the location of the data.
     */
    using location_tag = typename backend::buffer_traits<backend_tag>::location;

    /*!
     * \brief Alias to the container template associated with the backend (allows for RAII memory management).
     */
    template<typename T> using buffer_container = typename backend::buffer_traits<backend_tag>::template container<T>;
    //! \brief Container of real values corresponding to the complex type T.
    template<typename T> using real_buffer_container = buffer_container<typename define_standard_type<T>::type::value_type>;
    //! \brief Container of the output type corresponding to T, see \ref HeffteFFT3DCompatibleTypes "the table of compatible input and output types".
    template<typename T> using output_buffer_container = buffer_container<typename fft_output<T>::type>;

    /*!
     * \brief Constructor creating a plan for FFT transform across the given communicator and using the box geometry.
     *
     * \param inbox is the box for the non-transformed data, i.e., the input for the forward() transform and the output of the backward() transform.
     * \param outbox is the box for the transformed data, i.e., the output for the forward() transform and the input of the backward() transform.
     * \param r2c_direction indicates the direction where the total set of coefficients will be reduced to hold only the non-conjugate pairs;
     *        selecting a dimension with odd number of indexes will result in (slightly) smaller final data set.
     * \param comm is the MPI communicator with all ranks that will participate in the FFT.
     * \param options is a set of options that define the FFT plan, see heffte::plan_options for details.
     */
    fft3d_r2c(box3d<index> const inbox, box3d<index> const outbox, int r2c_direction, MPI_Comm const comm,
              plan_options const options = default_options<backend_tag>()) :
        fft3d_r2c(plan_operations(mpi::gather_boxes(inbox, outbox, comm), r2c_direction, set_options<backend_tag, true>(options), mpi::comm_rank(comm)), comm){
        assert(r2c_direction == 0 or r2c_direction == 1 or r2c_direction == 2);
        static_assert(backend::is_enabled<backend_tag>::value, "The requested backend is invalid or has not been enabled.");
    }
    /*!
     * \brief See the documentation for fft3d::fft3d()
     */
    fft3d_r2c(typename backend::device_instance<location_tag>::stream_type gpu_stream,
              box3d<index> const inbox, box3d<index> const outbox, int r2c_direction, MPI_Comm const comm,
              plan_options const options = default_options<backend_tag>()) :
        fft3d_r2c(gpu_stream,
                  plan_operations(mpi::gather_boxes(inbox, outbox, comm), r2c_direction, set_options<backend_tag, true>(options), mpi::comm_rank(comm)),
                  comm){
        assert(r2c_direction == 0 or r2c_direction == 1 or r2c_direction == 2);
        static_assert(backend::is_enabled<backend_tag>::value, "The requested backend is invalid or has not been enabled.");
    }

    //! \brief Internal use only, used by the Fortran interface
    fft3d_r2c(int il0, int il1, int il2, int ih0, int ih1, int ih2, int io0, int io1, int io2,
              int ol0, int ol1, int ol2, int oh0, int oh1, int oh2, int oo0, int oo1, int oo2,
              int r2c_direction, MPI_Comm const comm,
              bool use_reorder, int algorithm, bool use_pencils)
        : fft3d_r2c(box3d<index>({il0, il1, il2}, {ih0, ih1, ih2}, {io0, io1, io2}),
                    box3d<index>({ol0, ol1, ol2}, {oh0, oh1, oh2}, {oo0, oo1, oo2}),
                    r2c_direction, comm,
                plan_options(use_reorder, static_cast<reshape_algorithm>(algorithm), use_pencils))
    {}
    //! \brief Internal use only, used by the Fortran interface
    fft3d_r2c(int il0, int il1, int il2, int ih0, int ih1, int ih2, int io0, int io1, int io2,
              int ol0, int ol1, int ol2, int oh0, int oh1, int oh2, int oo0, int oo1, int oo2,
              int r2c_direction, MPI_Comm const comm)
        : fft3d_r2c(box3d<index>({il0, il1, il2}, {ih0, ih1, ih2}, {io0, io1, io2}),
                    box3d<index>({ol0, ol1, ol2}, {oh0, oh1, oh2}, {oo0, oo1, oo2}),
                    r2c_direction, comm)
    {}
    //! \brief Internal use only, used by the Fortran interface
    fft3d_r2c(int il0, int il1, int il2, int ih0, int ih1, int ih2,
              int ol0, int ol1, int ol2, int oh0, int oh1, int oh2,
              int r2c_direction, MPI_Comm const comm)
        : fft3d_r2c(box3d<index>({il0, il1, il2}, {ih0, ih1, ih2}), box3d<index>({ol0, ol1, ol2}, {oh0, oh1, oh2}), r2c_direction, comm)
    {}

    //! \brief Returns the size of the inbox defined in the constructor.
    long long size_inbox() const{ return pinbox->count(); }
    //! \brief Returns the size of the outbox defined in the constructor.
    long long size_outbox() const{ return poutbox->count(); }
    //! \brief Returns the inbox.
    box3d<index> inbox() const{ return *pinbox; }
    //! \brief Returns the outbox.
    box3d<index> outbox() const{ return *poutbox; }
    //! \brief Returns the workspace size that will be used, size is measured in complex numbers.
    size_t size_workspace() const{ return size_buffer_work; }
    //! \brief Returns the size used by the communication workspace buffers (internal use).
    size_t size_comm_buffers() const{ return comm_buffer_offset; }

    /*!
     * \brief Performs a forward Fourier transform using two arrays.
     *
     * \tparam input_type is either float or double type.
     * \tparam output_type is a type compatible with the output of a forward FFT,
     *         see \ref HeffteFFT3DCompatibleTypes "the table of compatible types".
     *
     * \param input is an array of size at least size_inbox() holding the input data corresponding
     *          to the inbox
     * \param output is an array of size at least size_outbox() and will be overwritten with
     *          the result from the transform corresponding to the outbox
     * \param scaling defines the type of scaling to apply (default no-scaling).
     */
    template<typename input_type, typename output_type>
    void forward(input_type const input[], output_type output[], scale scaling = scale::none) const{
        static_assert((std::is_same<input_type, float>::value and is_ccomplex<output_type>::value)
                   or (std::is_same<input_type, double>::value and is_zcomplex<output_type>::value),
                "Using either an unknown complex type or an incompatible pair of types!");

        auto workspace = make_buffer_container<output_type>(this->stream(), size_workspace());
        forward(input, output, workspace.data(), scaling);
    }

    //! \brief Overload utilizing a user provided buffer.
    template<typename input_type, typename output_type>
    void forward(input_type const input[], output_type output[], output_type workspace[], scale scaling = scale::none) const{
        static_assert((std::is_same<input_type, float>::value and is_ccomplex<output_type>::value)
                   or (std::is_same<input_type, double>::value and is_zcomplex<output_type>::value),
                "Using either an unknown complex type or an incompatible pair of types!");

        compute_transform<location_tag, index>(this->stream(), 1, convert_to_standard(input), convert_to_standard(output),
                                               convert_to_standard(workspace),
                                               executor_buffer_offset, size_comm_buffers(), forward_shaper,
                                               forward_executors(), direction::forward);
        apply_scale(1, direction::forward, scaling, output);
    }
    //! \brief Overload utilizing a batch transform.
    template<typename input_type, typename output_type>
    void forward(int batch_size, input_type const input[], output_type output[],
                 output_type workspace[], scale scaling = scale::none) const{
        static_assert((std::is_same<input_type, float>::value and is_ccomplex<output_type>::value)
                   or (std::is_same<input_type, double>::value and is_zcomplex<output_type>::value),
                "Using either an unknown complex type or an incompatible pair of types!");

        compute_transform<location_tag, index>(this->stream(), batch_size, convert_to_standard(input), convert_to_standard(output),
                                               convert_to_standard(workspace),
                                               executor_buffer_offset, size_comm_buffers(), forward_shaper,
                                               forward_executors(), direction::forward);
        apply_scale(batch_size, direction::forward, scaling, output);
    }
    //! \brief Overload utilizing a batch transform using internally allocated workspace.
    template<typename input_type, typename output_type>
    void forward(int batch_size, input_type const input[], output_type output[], scale scaling = scale::none) const{
        static_assert((std::is_same<input_type, float>::value and is_ccomplex<output_type>::value)
                   or (std::is_same<input_type, double>::value and is_zcomplex<output_type>::value),
                "Using either an unknown complex type or an incompatible pair of types!");

        auto workspace = make_buffer_container<output_type>(this->stream(), batch_size * size_workspace());
        forward(batch_size, input, output, workspace.data(), scaling);
    }

    /*!
     * \brief Vector variant of forward() using input and output buffer_container classes.
     *
     * Returns either std::vector or heffte::cuda:vector using only the C++ standard types.
     * Allows for more C++-ish calls and RAII memory management, see the heffte::fft3d equivalent.
     *
     * \tparam input_type is either float or double.
     *
     * \param input is a std::vector or heffte::cuda::vector with size at least size_inbox() corresponding to the input of forward().
     * \param scaling defines the type of scaling to apply (default no-scaling).
     *
     * \returns std::vector or heffte::cuda::vector with entries corresponding to the output type and with size equal to size_outbox()
     *          corresponding to the output of forward().
     *
     * \throws std::invalid_argument is the size of the \b input is less than size_inbox().
     */
    template<typename input_type>
    output_buffer_container<input_type> forward(buffer_container<input_type> const &input, scale scaling = scale::none){
        if (input.size() < static_cast<size_t>(size_inbox()))
            throw std::invalid_argument("The input vector is smaller than size_inbox(), i.e., not enough entries provided to fill the inbox.");
        static_assert(std::is_same<input_type, float>::value or std::is_same<input_type, double>::value,
                      "The input to forward() must be real, i.e., either float or double.");
        auto output = make_buffer_container<typename fft_output<input_type>::type>(this->stream(), size_outbox());
        forward(input.data(), output.data(), scaling);
        return output;
    }

    /*!
     * \brief Performs a backward Fourier transform using two arrays.
     *
     * \tparam input_type is either float or double.
     * \tparam output_type is a type compatible with the output of a backward FFT,
     *                     see \ref HeffteFFT3DCompatibleTypes "the table of compatible types".
     *
     * \param input is an array of size at least size_outbox() holding the input data corresponding
     *          to the outbox
     * \param output is an array of size at least size_inbox() and will be overwritten with
     *          the result from the transform corresponding to the inbox
     * \param scaling defines the type of scaling to apply (default no-scaling)
     */
    template<typename input_type, typename output_type>
    void backward(input_type const input[], output_type output[], scale scaling = scale::none) const{
        static_assert((std::is_same<output_type, float>::value and is_ccomplex<input_type>::value)
                   or (std::is_same<output_type, double>::value and is_zcomplex<input_type>::value),
                "Using either an unknown complex type or an incompatible pair of types!");

        auto workspace = make_buffer_container<input_type>(this->stream(), size_workspace());
        backward(input, output, workspace.data(), scaling);
    }

    //! \brief Overload utilizing a user provided buffer.
    template<typename input_type, typename output_type>
    void backward(input_type const input[], output_type output[], input_type workspace[], scale scaling = scale::none) const{
        static_assert((std::is_same<output_type, float>::value and is_ccomplex<input_type>::value)
                   or (std::is_same<output_type, double>::value and is_zcomplex<input_type>::value),
                "Using either an unknown complex type or an incompatible pair of types!");

        compute_transform<location_tag, index>(this->stream(), 1, convert_to_standard(input), convert_to_standard(output),
                                               convert_to_standard(workspace),
                                               executor_buffer_offset, size_comm_buffers(), backward_shaper,
                                               backward_executors(), direction::backward);
        apply_scale(1, direction::backward, scaling, output);
    }
    //! \brief Overload that performs a batch transform.
    template<typename input_type, typename output_type>
    void backward(int batch_size, input_type const input[], output_type output[],
                  input_type workspace[], scale scaling = scale::none) const{
        static_assert((std::is_same<output_type, float>::value and is_ccomplex<input_type>::value)
                   or (std::is_same<output_type, double>::value and is_zcomplex<input_type>::value),
                "Using either an unknown complex type or an incompatible pair of types!");

        compute_transform<location_tag, index>(this->stream(), batch_size, convert_to_standard(input), convert_to_standard(output),
                                               convert_to_standard(workspace),
                                               executor_buffer_offset, size_comm_buffers(), backward_shaper,
                                               backward_executors(), direction::backward);
        apply_scale(batch_size, direction::backward, scaling, output);
    }
    //! \brief Overload that performs a batch transform using internally allocated workspace.
    template<typename input_type, typename output_type>
    void backward(int batch_size, input_type const input[], output_type output[], scale scaling = scale::none) const{
        static_assert((std::is_same<output_type, float>::value and is_ccomplex<input_type>::value)
                   or (std::is_same<output_type, double>::value and is_zcomplex<input_type>::value),
                "Using either an unknown complex type or an incompatible pair of types!");

        auto workspace = make_buffer_container<input_type>(this->stream(), batch_size * size_workspace());
        backward(batch_size, input, output, workspace.data(), scaling);
    }

    /*!
     * \brief Variant of backward() that uses buffer_container for RAII style of resource management.
     */
    template<typename scalar_type>
    real_buffer_container<scalar_type> backward(buffer_container<scalar_type> const &input, scale scaling = scale::none){
        static_assert(is_ccomplex<scalar_type>::value or is_zcomplex<scalar_type>::value,
                      "Either calling backward() with non-complex input or using an unknown complex type.");
        auto result = make_buffer_container<typename define_standard_type<scalar_type>::type::value_type>(this->stream(), size_inbox());
        backward(input.data(), result.data(), scaling);
        return result;
    }

    /*!
     * \brief Returns the scale factor for the given scaling.
     */
    double get_scale_factor(scale scaling) const{ return (scaling == scale::symmetric) ? std::sqrt(scale_factor) : scale_factor; }

private:
    //! \brief Same as in the fft3d case.
    fft3d_r2c(logic_plan3d<index> const &plan, MPI_Comm const comm) :
        pinbox(new box3d<index>(plan.in_shape[0][plan.mpi_rank])), poutbox(new box3d<index>(plan.out_shape[3][plan.mpi_rank])),
        scale_factor(1.0 / static_cast<double>(plan.index_count))
        #ifdef Heffte_ENABLE_MAGMA
        , hmagma(this->stream())
        #endif
    {
        setup(plan, comm);
    }

    //! \brief Same as in the fft3d case.
    fft3d_r2c(typename backend::device_instance<location_tag>::stream_type gpu_stream,
              logic_plan3d<index> const &plan, MPI_Comm const comm) :
        backend::device_instance<location_tag>(gpu_stream),
        pinbox(new box3d<index>(plan.in_shape[0][plan.mpi_rank])), poutbox(new box3d<index>(plan.out_shape[3][plan.mpi_rank])),
        scale_factor(1.0 / static_cast<double>(plan.index_count))
        #ifdef Heffte_ENABLE_MAGMA
        , hmagma(this->stream())
        #endif
    {
        setup(plan, comm);
    }

    //! \brief Setup the executors and the reshapes.
    void setup(logic_plan3d<index> const &plan, MPI_Comm const comm){
        for(int i=0; i<4; i++){
            forward_shaper[i]    = make_reshape3d<backend_tag>(this->stream(), plan.in_shape[i], plan.out_shape[i], comm, plan.options);
            backward_shaper[3-i] = make_reshape3d<backend_tag>(this->stream(), plan.out_shape[i], plan.in_shape[i], comm, plan.options);
        }

        executors[0] = make_executor_r2c<backend_tag>(this->stream(), plan.out_shape[0][mpi::comm_rank(comm)], plan.fft_direction[0]);
        executors[1] = make_executor<backend_tag>(this->stream(), plan.out_shape[1][mpi::comm_rank(comm)], plan.fft_direction[1]);
        executors[2] = make_executor<backend_tag>(this->stream(), plan.out_shape[2][mpi::comm_rank(comm)], plan.fft_direction[2]);

        size_t executor_workspace_size = get_max_work_size(executors);
        comm_buffer_offset = std::max(get_workspace_size(forward_shaper), get_workspace_size(backward_shaper));
        size_buffer_work = comm_buffer_offset
                        + get_max_box_size_r2c(executors) + executor_workspace_size;
        executor_buffer_offset = (executor_workspace_size == 0) ? 0 : size_buffer_work - executor_workspace_size;
    }
    //! \brief Return references to the executors in forward order.
    std::array<executor_base*, 3> forward_executors() const{
        return std::array<executor_base*, 3>{executors[0].get(), executors[1].get(), executors[2].get()};
    }
    //! \brief Return references to the executors in forward order.
    std::array<executor_base*, 3> backward_executors() const{
        return std::array<executor_base*, 3>{executors[2].get(), executors[1].get(), executors[0].get()};
    }

    //! \brief Applies the scaling factor to the data.
    template<typename scalar_type>
    void apply_scale(int const batch_size, direction dir, scale scaling, scalar_type data[]) const{
        if (scaling != scale::none){
            add_trace name("scale");
            #ifdef Heffte_ENABLE_MAGMA
            if (std::is_same<location_tag, tag::gpu>::value){
                hmagma.scal(batch_size * (dir == direction::forward) ? size_outbox() : size_inbox(), get_scale_factor(scaling), data);
                return;
            }
            #endif
            data_scaling::apply(
                this->stream(),
                batch_size * ((dir == direction::forward) ? size_outbox() : size_inbox()),
                data, get_scale_factor(scaling));
        }
    }

    std::unique_ptr<box3d<index>> pinbox, poutbox;
    double scale_factor;
    std::array<std::unique_ptr<reshape3d_base<index>>, 4> forward_shaper;
    std::array<std::unique_ptr<reshape3d_base<index>>, 4> backward_shaper;

    std::array<std::unique_ptr<executor_base>, 3> executors;
    #ifdef Heffte_ENABLE_MAGMA
    gpu::magma_handle<location_tag> hmagma;
    #endif

    // cache some values for faster read
    size_t size_buffer_work, comm_buffer_offset, executor_buffer_offset;
};

/*!
 * \ingroup fft3d
 * \brief Alias of heffte::fft2d to be used for a two dimensional problem.
 */
template<typename backend_tag, typename index = int>
using fft2d_r2c = fft3d_r2c<backend_tag, index>;

/*!
 * \ingroup fft3d
 * \brief Factory method that auto-detects the index type based on the box.
 */
template<typename backend_tag, typename index>
fft3d_r2c<backend_tag, index> make_fft3d_r2c(box3d<index> const inbox, box3d<index> const outbox,
                                             int r2c_direction, MPI_Comm const comm,
                                             plan_options const options = default_options<backend_tag>()){
    static_assert(std::is_same<index, int>::value or std::is_same<index, long long>::value,
                  "heFFTe works with 'int' and 'long long' indexing only");
    static_assert(backend::is_enabled<backend_tag>::value,
                  "The backend_tag is not valid, perhaps it needs to be enabled in the build system");
    return fft3d_r2c<backend_tag, index>(inbox, outbox, r2c_direction, comm, options);
}

}

#endif
