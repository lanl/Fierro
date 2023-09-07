/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_BACKEND_ROCM_H
#define HEFFTE_BACKEND_ROCM_H

#include "heffte_r2r_executor.h"

#ifdef Heffte_ENABLE_ROCM

#ifndef __HIP_PLATFORM_HCC__
#define __HIP_PLATFORM_HCC__
#endif
#include <hip/hip_runtime.h>
#include <rocfft.h>
#include "heffte_backend_vector.h"

#ifdef Heffte_ENABLE_MAGMA
#include "heffte_magma_helpers.h"
#endif

/*!
 * \ingroup fft3d
 * \addtogroup heffterocm Backend rocfft
 *
 * Wrappers and template specializations related to the rocFFT backend.
 * Requires CMake option:
 * \code
 *  -D Heffte_ENABLE_ROCM=ON
 * \endcode
 *
 * In addition to the rocFFT wrappers, this also includes a series of kernels
 * for packing/unpacking/scaling the data, as well as a simple container
 * that wraps around ROCM arrays for RAII style of resource management.
 */

namespace heffte{

/*!
 * \ingroup heffterocm
 * \brief ROCM specific methods, vector-like container, error checking, etc.
 */
namespace rocm {
    /*!
     * \ingroup heffterocm
     * \brief Checks the status of a ROCm command and in case of a failure, converts it to a C++ exception.
     */
    inline void check_error(hipError_t status, const char *function_name){
        if (status != hipSuccess)
            throw std::runtime_error(std::string(function_name) + " failed with message: " + std::string(hipGetErrorString(status)));
    }
    /*!
     * \ingroup heffterocm
     * \brief Checks the status of a cufft command and in case of a failure, converts it to a C++ exception.
     */
    inline void check_error(rocfft_status status, const char *function_name){
        if (status != rocfft_status_success)
            throw std::runtime_error(std::string(function_name) + " failed with error code: " + std::to_string(status));
    }

    /*!
     * \ingroup heffterocm
     * \brief Convert real numbers to complex when both are located on the GPU device.
     *
     * Launches a ROCM kernel.
     */
    template<typename precision_type, typename index>
    void convert(hipStream_t stream, index num_entries, precision_type const source[], std::complex<precision_type> destination[]);
    /*!
     * \ingroup heffterocm
     * \brief Convert complex numbers to real when both are located on the GPU device.
     *
     * Launches a ROCM kernel.
     */
    template<typename precision_type, typename index>
    void convert(hipStream_t stream, index num_entries, std::complex<precision_type> const source[], precision_type destination[]);

    /*!
     * \ingroup heffterocm
     * \brief Scales real data (double or float) by the scaling factor.
     */
    template<typename scalar_type, typename index>
    void scale_data(hipStream_t stream, index num_entries, scalar_type *data, double scale_factor);

    /*!
     * \ingroup heffterocm
     * \brief Performs a direct-pack operation for data sitting on the GPU device.
     *
     * Launches a HIP kernel.
     */
    template<typename scalar_type, typename index>
    void direct_pack(hipStream_t stream, index nfast, index nmid, index nslow, index line_stride, index plane_stide,
                    scalar_type const source[], scalar_type destination[]);
    /*!
     * \ingroup heffterocm
     * \brief Performs a direct-unpack operation for data sitting on the GPU device.
     *
     * Launches a HIP kernel.
     */
    template<typename scalar_type, typename index>
    void direct_unpack(hipStream_t stream, index nfast, index nmid, index nslow, index line_stride, index plane_stide,
                    scalar_type const source[], scalar_type destination[]);
    /*!
     * \ingroup heffterocm
     * \brief Performs a transpose-unpack operation for data sitting on the GPU device.
     *
     * Launches a HIP kernel.
     */
    template<typename scalar_type, typename index>
    void transpose_unpack(hipStream_t stream, index nfast, index nmid, index nslow, index line_stride, index plane_stide,
                        index buff_line_stride, index buff_plane_stride, int map0, int map1, int map2,
                        scalar_type const source[], scalar_type destination[]);

    /*!
     * \ingroup hefftecuda
     * \brief Implementation of Cosine Transform pre-post processing methods using CUDA.
     */
    struct cos_pre_pos_processor{
        //! \brief Pre-process in the forward transform.
        template<typename precision>
        static void pre_forward(hipStream_t, int length, precision const input[], precision fft_signal[]);
        //! \brief Post-process in the forward transform.
        template<typename precision>
        static void post_forward(hipStream_t, int length, std::complex<precision> const fft_result[], precision result[]);
        //! \brief Pre-process in the inverse transform.
        template<typename precision>
        static void pre_backward(hipStream_t, int length, precision const input[], std::complex<precision> fft_signal[]);
        //! \brief Post-process in the inverse transform.
        template<typename precision>
        static void post_backward(hipStream_t, int length, precision const fft_result[], precision result[]);
    };
    /*!
     * \ingroup hefftecuda
     * \brief Implementation of Sine Transform pre-post processing methods using CUDA.
     */
    struct sin_pre_pos_processor{
        //! \brief Pre-process in the forward transform.
        template<typename precision>
        static void pre_forward(hipStream_t, int length, precision const input[], precision fft_signal[]);
        //! \brief Post-process in the forward transform.
        template<typename precision>
        static void post_forward(hipStream_t, int length, std::complex<precision> const fft_result[], precision result[]);
        //! \brief Pre-process in the inverse transform.
        template<typename precision>
        static void pre_backward(hipStream_t, int length, precision const input[], std::complex<precision> fft_signal[]);
        //! \brief Post-process in the inverse transform.
        template<typename precision>
        static void post_backward(hipStream_t, int length, precision const fft_result[], precision result[]);
    };
}

namespace backend{
    /*!
     * \ingroup heffterocm
     * \brief Indicate that the rocFFT backend has been enabled.
     */
    template<> struct is_enabled<rocfft> : std::true_type{};
    /*!
     * \ingroup heffterocm
     * \brief Indicate that the rocFFT backend has been enabled.
     */
    template<> struct is_enabled<rocfft_cos> : std::true_type{};
    /*!
     * \ingroup heffterocm
     * \brief Indicate that the rocFFT backend has been enabled.
     */
    template<> struct is_enabled<rocfft_sin> : std::true_type{};

    /*!
     * \ingroup heffterocm
     * \brief The ROCm backend uses a HIP stream.
     */
    template<>
    struct device_instance<tag::gpu>{
        //! \brief Constructor, sets up the stream.
        device_instance(hipStream_t new_stream = nullptr) : _stream(new_stream){}
        //! \brief Returns the nullptr.
        hipStream_t stream(){ return _stream; }
        //! \brief Returns the nullptr (const case).
        hipStream_t stream() const{ return _stream; }
        //! \brief Syncs the execution with the queue.
        void synchronize_device() const{ rocm::check_error(hipStreamSynchronize(_stream), "device sync"); }
        //! \brief The CUDA stream to be used in all operations.
        mutable hipStream_t _stream;
        //! \brief The type for the internal stream.
        using stream_type = hipStream_t;
    };

    /*!
     * \ingroup heffterocm
     * \brief In ROCm mode, the default GPU backend is rocfft.
     */
    template<> struct default_backend<tag::gpu>{
        //! \brief Set the rocfft tag.
        using type = rocfft;
    };

    /*!
     * \ingroup heffterocm
     * \brief Specialization for the data operations in ROCm mode.
     */
    template<> struct data_manipulator<tag::gpu> {
        //! \brief The stream type for the device.
        using stream_type = hipStream_t;
        //! \brief Defines the backend_device.
        using backend_device = backend::device_instance<tag::gpu>;
        //! \brief Allocate memory.
        template<typename scalar_type>
        static scalar_type* allocate(hipStream_t stream, size_t num_entries){
            void *new_data;
            if (stream != nullptr) rocm::check_error( hipStreamSynchronize(stream), "hipStreamSynchronize()");
            rocm::check_error(hipMalloc(&new_data, num_entries * sizeof(scalar_type)), "hipMalloc()");
            return reinterpret_cast<scalar_type*>(new_data);
        }
        //! \brief Free memory.
        template<typename scalar_type>
        static void free(hipStream_t stream, scalar_type *pntr){
            if (pntr == nullptr) return;
            if (stream != nullptr) rocm::check_error( hipStreamSynchronize(stream), "hipStreamSynchronize()");
            rocm::check_error(hipFree(pntr), "hipFree()");
        }
        //! \brief Equivalent to std::copy_n() but using CUDA arrays.
        template<typename scalar_type>
        static void copy_n(hipStream_t stream, scalar_type const source[], size_t num_entries, scalar_type destination[]){
            if (stream == nullptr)
                rocm::check_error(hipMemcpy(destination, source, num_entries * sizeof(scalar_type), hipMemcpyDeviceToDevice), "data_manipulator::copy_n()");
            else
                rocm::check_error(hipMemcpyAsync(destination, source, num_entries * sizeof(scalar_type), hipMemcpyDeviceToDevice, stream), "data_manipulator::copy_n()");
        }
        //! \brief Copy-convert complex-to-real.
        template<typename scalar_type>
        static void copy_n(hipStream_t stream, std::complex<scalar_type> const source[], size_t num_entries, scalar_type destination[]){
            rocm::convert(stream, static_cast<long long>(num_entries), source, destination);
        }
        //! \brief Copy-convert real-to-complex.
        template<typename scalar_type>
        static void copy_n(hipStream_t stream, scalar_type const source[], size_t num_entries, std::complex<scalar_type> destination[]){
            rocm::convert(stream, static_cast<long long>(num_entries), source, destination);
        }
        //! \brief Copy the date from the device to the host.
        template<typename scalar_type>
        static void copy_device_to_host(hipStream_t stream, scalar_type const source[], size_t num_entries, scalar_type destination[]){
            rocm::check_error(hipMemcpyAsync(destination, source, num_entries * sizeof(scalar_type), hipMemcpyDeviceToHost, stream),
                            "device_to_host (rocm)");
        }
        //! \brief Copy the date from the device to the device.
        template<typename scalar_type>
        static void copy_device_to_device(hipStream_t stream, scalar_type const source[], size_t num_entries, scalar_type destination[]){
                rocm::check_error(hipMemcpyAsync(destination, source, num_entries * sizeof(scalar_type), hipMemcpyDeviceToDevice, stream),
                                "device_to_device (rocm)");
        }
        //! \brief Copy the date from the host to the device.
        template<typename scalar_type>
        static void copy_host_to_device(hipStream_t stream, scalar_type const source[], size_t num_entries, scalar_type destination[]){
            rocm::check_error(hipMemcpyAsync(destination, source, num_entries * sizeof(scalar_type), hipMemcpyHostToDevice, stream),
                            "host_to_device (rocm)");
        }
    };

    /*!
     * \ingroup heffterocm
     * \brief Defines the location type-tag and the cuda container.
     */
    template<>
    struct buffer_traits<rocfft>{
        //! \brief The rocfft library uses data on the gpu device.
        using location = tag::gpu;
        //! \brief The data is managed by the ROCm vector container.
        template<typename T> using container = heffte::gpu::device_vector<T, data_manipulator<tag::gpu>>;
    };
    /*!
     * \ingroup heffterocm
     * \brief Defines the location type-tag and the cuda container.
     */
    template<>
    struct buffer_traits<rocfft_cos>{
        //! \brief The rocfft library uses data on the gpu device.
        using location = tag::gpu;
        //! \brief The data is managed by the ROCm vector container.
        template<typename T> using container = heffte::gpu::device_vector<T, data_manipulator<tag::gpu>>;
    };
    /*!
     * \ingroup heffterocm
     * \brief Defines the location type-tag and the cuda container.
     */
    template<>
    struct buffer_traits<rocfft_sin>{
        //! \brief The rocfft library uses data on the gpu device.
        using location = tag::gpu;
        //! \brief The data is managed by the ROCm vector container.
        template<typename T> using container = heffte::gpu::device_vector<T, data_manipulator<tag::gpu>>;
    };
}

/*!
 * \ingroup heffterocm
 * \brief Plan for the r2c single precision transform.
 *
 * Note, this is a base template and does not specialize,
 * the complex case is handled in a specialization.
 */
template<typename precision_type, direction dir>
struct plan_rocfft{
    /*!
     * \brief Constructor and initializer of the plan.
     *
     * \param size is the number of entries in a 1-D transform
     * \param batch is the number of transforms in the batch
     * \param stride is the distance between entries of the same transform
     * \param rdist is the distance between the first entries of consecutive real sequences
     * \param cdist is the distance between the first entries of consecutive complex sequences
     */
    plan_rocfft(size_t size, size_t batch, size_t stride, size_t rdist, size_t cdist){

        rocfft_plan_description desc = nullptr;
        rocm::check_error( rocfft_plan_description_create(&desc), "rocm plan create");

        rocm::check_error(
            rocfft_plan_description_set_data_layout(
                desc,
                (dir == direction::forward) ? rocfft_array_type_real : rocfft_array_type_hermitian_interleaved,
                (dir == direction::forward) ? rocfft_array_type_hermitian_interleaved : rocfft_array_type_real,
                nullptr, nullptr,
                1, &stride, (dir == direction::forward) ? rdist : cdist,
                1, &stride, (dir == direction::forward) ? cdist : rdist
            ),
            "plan layout"
        );

        rocm::check_error(
        rocfft_plan_create(&plan, rocfft_placement_notinplace,
                           (dir == direction::forward) ? rocfft_transform_type_real_forward : rocfft_transform_type_real_inverse,
                           (std::is_same<precision_type, float>::value)? rocfft_precision_single : rocfft_precision_double,
                           1, &size, batch, desc),
        "plan create");

        rocm::check_error( rocfft_plan_get_work_buffer_size(plan, &worksize), "get_worksize");

        rocm::check_error( rocfft_plan_description_destroy(desc), "rocm plan destroy");
    }
    //! \brief Destructor, deletes the plan.
    ~plan_rocfft(){ rocfft_plan_destroy(plan); }
    //! \brief Custom conversion to the rocfft_plan.
    operator rocfft_plan() const{ return plan; }
    //! \brief Return the worksize.
    size_t size_work() const{ return worksize; }

private:
    //! \brief The rocfft opaque structure (pointer to struct).
    rocfft_plan plan;
    //! \brief The size of the scratch workspace.
    size_t worksize;
};

/*!
 * \ingroup heffterocm
 * \brief Plan for the single precision complex transform.
 */
template<typename precision_type, direction dir>
struct plan_rocfft<std::complex<precision_type>, dir>{
    /*!
     * \brief Constructor, takes inputs identical to cufftMakePlanMany().
     *
     * \param size is the number of entries in a 1-D transform
     * \param batch is the number of transforms in the batch
     * \param stride is the distance between entries of the same transform
     * \param dist is the distance between the first entries of consecutive sequences
     */
    plan_rocfft(size_t size, size_t batch, size_t stride, size_t dist) : plan(nullptr), worksize(0){
        rocfft_plan_description desc = nullptr;
        rocm::check_error( rocfft_plan_description_create(&desc), "rocm plan create");

        rocm::check_error(
            rocfft_plan_description_set_data_layout(
                desc,
                rocfft_array_type_complex_interleaved,
                rocfft_array_type_complex_interleaved,
                nullptr, nullptr,
                1, &stride, dist, 1, &stride, dist
            ),
            "plan layout"
        );

        rocm::check_error(
        rocfft_plan_create(&plan, rocfft_placement_inplace,
                           (dir == direction::forward) ? rocfft_transform_type_complex_forward : rocfft_transform_type_complex_inverse,
                           (std::is_same<precision_type, float>::value)? rocfft_precision_single : rocfft_precision_double,
                           1, &size, batch, desc),
        "plan create");

        rocm::check_error( rocfft_plan_get_work_buffer_size(plan, &worksize), "get_worksize");

        rocm::check_error( rocfft_plan_description_destroy(desc), "rocm plan destroy");
    }
    /*!
     * \brief Constructor, takes inputs identical to cufftMakePlanMany().
     *
     * \param size1 is the number of entries in a 2-D transform, direction 1
     * \param size2 is the number of entries in a 2-D transform, direction 2
     * \param embed is the stride between entries in each dimension
     * \param batch is the number of transforms in the batch
     * \param dist is the distance between the first entries of consecutive sequences
     */
    plan_rocfft(size_t size1, size_t size2, std::array<size_t, 2> const &embed, size_t batch, size_t dist) : plan(nullptr), worksize(0){
        size_t size[2] = {size1, size2};

        rocfft_plan_description desc = nullptr;
        rocm::check_error( rocfft_plan_description_create(&desc), "rocm plan create");

        rocm::check_error(
            rocfft_plan_description_set_data_layout(
                desc,
                rocfft_array_type_complex_interleaved,
                rocfft_array_type_complex_interleaved,
                nullptr, nullptr,
                2, embed.data(), dist, 2, embed.data(), dist
            ),
            "plan layout"
        );

        rocm::check_error(
        rocfft_plan_create(&plan, rocfft_placement_inplace,
                           (dir == direction::forward) ? rocfft_transform_type_complex_forward : rocfft_transform_type_complex_inverse,
                           (std::is_same<precision_type, float>::value) ? rocfft_precision_single : rocfft_precision_double,
                           2, size, batch, desc),
        "plan create");

        rocm::check_error( rocfft_plan_get_work_buffer_size(plan, &worksize), "get_worksize");

        rocm::check_error( rocfft_plan_description_destroy(desc), "rocm plan destroy");
    }
    //! \brief Constructor, takes inputs identical to cufftPlan3d()
    plan_rocfft(size_t size1, size_t size2, size_t size3){
        std::array<size_t, 3> size = {size1, size2, size3};
        rocfft_plan_description desc = nullptr;
        rocm::check_error( rocfft_plan_description_create(&desc), "rocm plan create");

        rocm::check_error(
            rocfft_plan_description_set_data_layout(
                desc,
                rocfft_array_type_complex_interleaved,
                rocfft_array_type_complex_interleaved,
                nullptr, nullptr, 3, nullptr, 1, 3, nullptr, 1
            ),
            "plan layout"
        );

        rocm::check_error(
        rocfft_plan_create(&plan, rocfft_placement_inplace,
                           (dir == direction::forward) ? rocfft_transform_type_complex_forward : rocfft_transform_type_complex_inverse,
                           (std::is_same<precision_type, float>::value) ? rocfft_precision_single : rocfft_precision_double,
                           3, size.data(), 1, desc),
        "plan create 3d");

        rocm::check_error( rocfft_plan_get_work_buffer_size(plan, &worksize), "get_worksize");

        rocm::check_error( rocfft_plan_description_destroy(desc), "rocm plan destroy");
    }
    //! \brief Destructor, deletes the plan.
    ~plan_rocfft(){ rocm::check_error( rocfft_plan_destroy(plan), "plan destory"); }
    //! \brief Custom conversion to the rocfft_plan.
    operator rocfft_plan() const{ return plan; }
    //! \brief Return the worksize.
    size_t size_work() const{ return worksize; }

private:
    //! \brief The rocfft opaque structure (pointer to struct).
    rocfft_plan plan;
    size_t worksize;
};

/*!
 * \ingroup heffterocm
 * \brief Wrapper around the rocFFT API.
 *
 * A single class that manages the plans and executions of rocFFT
 * so that a single API is provided for all backends.
 * The executor operates on a box and performs 1-D FFTs
 * for the given dimension.
 * The class silently manages the plans and buffers needed
 * for the different types.
 * All input and output arrays must have size equal to the box.
 */
class rocfft_executor : public executor_base{
public:
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::forward;
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::backward;
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::complex_size;
    //! \brief Constructor, specifies the box and dimension.
    template<typename index>
    rocfft_executor(hipStream_t active_stream, box3d<index> const box, int dimension) :
        stream(active_stream),
        size(box.size[dimension]), size2(0),
        howmanyffts(fft1d_get_howmany(box, dimension)),
        stride(fft1d_get_stride(box, dimension)),
        dist((dimension == box.order[0]) ? size : 1),
        blocks((dimension == box.order[1]) ? box.osize(2) : 1),
        block_stride(box.osize(0) * box.osize(1)),
        total_size(box.count()),
        embed({0, 0}),
        worksize(compute_workspace_size())
    {}
    //! \brief Merges two FFTs into one.
    template<typename index>
    rocfft_executor(hipStream_t active_stream, box3d<index> const box, int dir1, int dir2) :
        stream(active_stream),
        size(box.size[std::min(dir1, dir2)]), size2(box.size[std::max(dir1, dir2)]),
        blocks(1), block_stride(0), total_size(box.count()),
        worksize(0)
    {
        int odir1 = box.find_order(dir1);
        int odir2 = box.find_order(dir2);

        if (std::min(odir1, odir2) == 0 and std::max(odir1, odir2) == 1){
            stride = 1;
            dist = size * size2;
            embed = {static_cast<size_t>(stride), static_cast<size_t>(size)};
            howmanyffts = box.size[2];
        }else if (std::min(odir1, odir2) == 1 and std::max(odir1, odir2) == 2){
            stride = box.size[0];
            dist = 1;
            embed = {static_cast<size_t>(stride), static_cast<size_t>(size) * static_cast<size_t>(stride)};
            howmanyffts = box.size[0];
        }else{ // case of directions (0, 2)
            stride = 1;
            dist = size;
            embed = {static_cast<size_t>(stride), static_cast<size_t>(box.size[1]) * static_cast<size_t>(box.size[0])};
            howmanyffts = box.size[1];
        }
        worksize = compute_workspace_size();
    }
    //! \brief Merges three FFTs into one.
    template<typename index>
    rocfft_executor(hipStream_t active_stream, box3d<index> const box) :
        stream(active_stream),
        size(box.size[0]), size2(box.size[1]), howmanyffts(box.size[2]),
        stride(0), dist(0),
        blocks(1), block_stride(0),
        total_size(box.count()),
        embed({0, 0}),
        worksize(compute_workspace_size())
    {}

    //! \brief Perform an in-place FFT on the data in the given direction.
    template<typename precision_type, direction dir>
    void execute(std::complex<precision_type> data[], std::complex<precision_type> *workspace) const{
        if (std::is_same<precision_type, float>::value){
            if (dir == direction::forward)
                make_plan(ccomplex_forward);
            else
                make_plan(ccomplex_backward);
        }else{
            if (dir == direction::forward)
                make_plan(zcomplex_forward);
            else
                make_plan(zcomplex_backward);
        }
        rocfft_execution_info info;
        rocfft_execution_info_create(&info);
        rocfft_execution_info_set_stream(info, stream);

        size_t wsize = (std::is_same<precision_type, float>::value) ?
                            ((dir == direction::forward) ? ccomplex_forward->size_work() : ccomplex_backward->size_work()) :
                            ((dir == direction::forward) ? zcomplex_forward->size_work() : zcomplex_backward->size_work());

        if (wsize > 0)
            rocfft_execution_info_set_work_buffer(info, reinterpret_cast<void*>(workspace), wsize);

        for(int i=0; i<blocks; i++){
            void* block_data = reinterpret_cast<void*>(data + i * block_stride);
            rocm::check_error( rocfft_execute(
                (std::is_same<precision_type, float>::value) ?
                    ((dir == direction::forward) ? *ccomplex_forward : *ccomplex_backward) :
                    ((dir == direction::forward) ? *zcomplex_forward : *zcomplex_backward),
                &block_data, nullptr, info), "rocfft execute");
        }
        rocfft_execution_info_destroy(info);
    }

    //! \brief Forward fft, float-complex case.
    void forward(std::complex<float> data[], std::complex<float> *workspace) const override{
        execute<float, direction::forward>(data, workspace);
    }
    //! \brief Forward fft, double-complex case.
    void forward(std::complex<double> data[], std::complex<double> *workspace) const override{
        execute<double, direction::forward>(data, workspace);
    }
    //! \brief Backward fft, float-complex case.
    void backward(std::complex<float> data[], std::complex<float> *workspace) const override{
        execute<float, direction::backward>(data, workspace);
    }
    //! \brief Backward fft, double-complex case.
    void backward(std::complex<double> data[], std::complex<double> *workspace) const override{
        execute<double, direction::backward>(data, workspace);
    }

    //! \brief Converts the deal data to complex and performs float-complex forward transform.
    void forward(float const indata[], std::complex<float> outdata[], std::complex<float> *workspace) const override{
        rocm::convert(stream, total_size, indata, outdata);
        forward(outdata, workspace);
    }
    //! \brief Converts the deal data to complex and performs double-complex forward transform.
    void forward(double const indata[], std::complex<double> outdata[], std::complex<double> *workspace) const override{
        rocm::convert(stream, total_size, indata, outdata);
        forward(outdata, workspace);
    }
    //! \brief Performs backward float-complex transform and truncates the complex part of the result.
    void backward(std::complex<float> indata[], float outdata[], std::complex<float> *workspace) const override{
        backward(indata, workspace);
        rocm::convert(stream, total_size, indata, outdata);
    }
    //! \brief Performs backward double-complex transform and truncates the complex part of the result.
    void backward(std::complex<double> indata[], double outdata[], std::complex<double> *workspace) const override{
        backward(indata, workspace);
        rocm::convert(stream, total_size, indata, outdata);
    }

    //! \brief Returns the size of the box.
    int box_size() const override{ return total_size; }
    //! \brief Return the size of the needed workspace.
    size_t workspace_size() const override{ return worksize; }
    //! \brief Computes the size of the needed workspace.
    size_t compute_workspace_size() const{
        make_plan(ccomplex_forward);
        make_plan(ccomplex_backward);
        make_plan(zcomplex_forward);
        make_plan(zcomplex_backward);
        return
        std::max( std::max(ccomplex_forward->size_work(), ccomplex_backward->size_work()) / sizeof(std::complex<float>),
                  std::max(zcomplex_forward->size_work(), zcomplex_backward->size_work()) / sizeof(std::complex<double>) ) + 1;
        return 0;
    }

private:
    //! \brief Helper template to create the plan.
    template<typename scalar_type, direction dir>
    void make_plan(std::unique_ptr<plan_rocfft<scalar_type, dir>> &plan) const{
        if (not plan){
            if (dist == 0)
                plan = std::unique_ptr<plan_rocfft<scalar_type, dir>>(new plan_rocfft<scalar_type, dir>(size, size2, howmanyffts));
            else if (size2 == 0)
                plan = std::unique_ptr<plan_rocfft<scalar_type, dir>>(new plan_rocfft<scalar_type, dir>(size, howmanyffts, stride, dist));
            else
                plan = std::unique_ptr<plan_rocfft<scalar_type, dir>>(new plan_rocfft<scalar_type, dir>(size, size2, embed, howmanyffts, dist));
        }
    }

    mutable hipStream_t stream;

    int size, size2, howmanyffts, stride, dist, blocks, block_stride, total_size;
    std::array<size_t, 2> embed;
    mutable std::unique_ptr<plan_rocfft<std::complex<float>, direction::forward>> ccomplex_forward;
    mutable std::unique_ptr<plan_rocfft<std::complex<float>, direction::backward>> ccomplex_backward;
    mutable std::unique_ptr<plan_rocfft<std::complex<double>, direction::forward>> zcomplex_forward;
    mutable std::unique_ptr<plan_rocfft<std::complex<double>, direction::backward>> zcomplex_backward;

    size_t worksize;
};

/*!
 * \ingroup heffterocm
 * \brief Wrapper to rocFFT API for real-to-complex transform with shortening of the data.
 *
 * Serves the same purpose of heffte::rocfft_executor but only real input is accepted
 * and only the unique (non-conjugate) coefficients are computed.
 * All real arrays must have size of real_size() and all complex arrays must have size complex_size().
 */
class rocfft_executor_r2c : public executor_base{
public:
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::forward;
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::backward;
    /*!
     * \brief Constructor defines the box and the dimension of reduction.
     *
     * Note that the result sits in the box returned by box.r2c(dimension).
     */
    template<typename index>
    rocfft_executor_r2c(hipStream_t active_stream, box3d<index> const box, int dimension) :
        stream(active_stream),
        size(box.size[dimension]),
        howmanyffts(fft1d_get_howmany(box, dimension)),
        stride(fft1d_get_stride(box, dimension)),
        blocks((dimension == box.order[1]) ? box.osize(2) : 1),
        rdist((dimension == box.order[0]) ? size : 1),
        cdist((dimension == box.order[0]) ? size/2 + 1 : 1),
        rblock_stride(box.osize(0) * box.osize(1)),
        cblock_stride(box.osize(0) * (box.osize(1)/2 + 1)),
        rsize(box.count()),
        csize(box.r2c(dimension).count()),
        worksize(compute_workspace_size())
    {}

    //! \brief Forward transform, single precision.
    template<typename precision_type>
    void forward(precision_type const indata[], std::complex<precision_type> outdata[], std::complex<precision_type> *workspace) const{
        if (std::is_same<precision_type, float>::value){
            make_plan(sforward);
        }else{
            make_plan(dforward);
        }

        rocfft_execution_info info;
        rocfft_execution_info_create(&info);
        rocfft_execution_info_set_stream(info, stream);

        size_t wsize = (std::is_same<precision_type, float>::value) ? sforward->size_work() : dforward->size_work();
        if (wsize > 0)
            rocfft_execution_info_set_work_buffer(info, reinterpret_cast<void*>(workspace), wsize);

        precision_type *copy_indata = reinterpret_cast<precision_type*>(
            reinterpret_cast<unsigned char *>(workspace) + wsize);
        backend::data_manipulator<tag::gpu>::copy_n(stream, indata, box_size(), copy_indata);

        for(int i=0; i<blocks; i++){
            void *rdata = const_cast<void*>(reinterpret_cast<void const*>(copy_indata + i * rblock_stride));
            void *cdata = reinterpret_cast<void*>(outdata + i * cblock_stride);
            rocm::check_error( rocfft_execute(
                (std::is_same<precision_type, float>::value) ? *sforward : *dforward,
                &rdata, &cdata, info), "rocfft execute");
        }
        rocfft_execution_info_destroy(info);
    }
    //! \brief Backward transform, single precision.
    template<typename precision_type>
    void backward(std::complex<precision_type> indata[], precision_type outdata[], std::complex<precision_type> *workspace) const{
        if (std::is_same<precision_type, float>::value){
            make_plan(sbackward);
        }else{
            make_plan(dbackward);
        }

        rocfft_execution_info info;
        rocfft_execution_info_create(&info);
        rocfft_execution_info_set_stream(info, stream);

        size_t wsize = (std::is_same<precision_type, float>::value) ? sbackward->size_work() : dbackward->size_work();
        if (wsize > 0)
            rocfft_execution_info_set_work_buffer(info, reinterpret_cast<void*>(workspace), wsize);

        std::complex<precision_type> *copy_indata = reinterpret_cast<std::complex<precision_type>*>(
            reinterpret_cast<unsigned char *>(workspace) + wsize);
        backend::data_manipulator<tag::gpu>::copy_n(stream, indata, complex_size(), copy_indata);

        for(int i=0; i<blocks; i++){
            void *cdata = const_cast<void*>(reinterpret_cast<void const*>(copy_indata + i * cblock_stride));
            void *rdata = reinterpret_cast<void*>(outdata + i * rblock_stride);
            rocm::check_error( rocfft_execute(
                (std::is_same<precision_type, float>::value) ? *sbackward : *dbackward,
                &cdata, &rdata, info), "rocfft execute");
        }
        rocfft_execution_info_destroy(info);
    }
    //! \brief Forward transform, single precision.
    void forward(float const indata[], std::complex<float> outdata[], std::complex<float> *workspace) const override{
        forward<float>(indata, outdata, workspace);
    }
    //! \brief Backward transform, single precision.
    void backward(std::complex<float> indata[], float outdata[], std::complex<float> *workspace) const override{
        backward<float>(indata, outdata, workspace);
    }
    //! \brief Forward transform, double precision.
    void forward(double const indata[], std::complex<double> outdata[], std::complex<double> *workspace) const override{
        forward<double>(indata, outdata, workspace);
    }
    //! \brief Backward transform, double precision.
    void backward(std::complex<double> indata[], double outdata[], std::complex<double> *workspace) const override{
        backward<double>(indata, outdata, workspace);
    }

    //! \brief Returns the size of the box with real data.
    int box_size() const override{ return rsize; }
    //! \brief Returns the size of the box with complex coefficients.
    int complex_size() const override{ return csize; }
    //! \brief Return the size of the needed workspace.
    size_t workspace_size() const override{ return worksize; }
    //! \brief Computes the size of the needed workspace.
    size_t compute_workspace_size() const{
        make_plan(sforward);
        make_plan(dforward);
        make_plan(sbackward);
        make_plan(dbackward);
        // Temporary copies have to be made, request that from user in addition, to what rocFFT requires.
        return
        std::max( std::max(sforward->size_work() + box_size() * sizeof(float),  sbackward->size_work() + complex_size() * sizeof(std::complex<float>))  / sizeof(std::complex<float>),
                  std::max(dforward->size_work() + box_size() * sizeof(double), dbackward->size_work() + complex_size() * sizeof(std::complex<double>)) / sizeof(std::complex<double>) ) + 1;
    }

private:
    //! \brief Helper template to initialize the plan.
    template<typename scalar_type, direction dir>
    void make_plan(std::unique_ptr<plan_rocfft<scalar_type, dir>> &plan) const{
        if (!plan) plan = std::unique_ptr<plan_rocfft<scalar_type, dir>>(new plan_rocfft<scalar_type, dir>(size, howmanyffts, stride, rdist, cdist));
    }

    mutable hipStream_t stream;

    int size, howmanyffts, stride, blocks;
    int rdist, cdist, rblock_stride, cblock_stride, rsize, csize;
    mutable std::unique_ptr<plan_rocfft<float, direction::forward>> sforward;
    mutable std::unique_ptr<plan_rocfft<double, direction::forward>> dforward;
    mutable std::unique_ptr<plan_rocfft<float, direction::backward>> sbackward;
    mutable std::unique_ptr<plan_rocfft<double, direction::backward>> dbackward;

    size_t worksize;
};

/*!
 * \ingroup heffterocm
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 *
 * The struct is specialized for each backend.
 */
template<> struct one_dim_backend<backend::rocfft>{
    //! \brief Defines the complex-to-complex executor.
    using executor = rocfft_executor;
    //! \brief Defines the real-to-complex executor.
    using executor_r2c = rocfft_executor_r2c;
};

/*!
 * \ingroup heffterocm
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 *
 * The struct is specialized for each backend.
 */
template<> struct one_dim_backend<backend::rocfft_cos>{
    //! \brief Defines the complex-to-complex executor.
    using executor = real2real_executor<backend::rocfft, rocm::cos_pre_pos_processor>;
    //! \brief Defines the real-to-complex executor.
    using executor_r2c = void;
};
/*!
 * \ingroup heffterocm
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 *
 * The struct is specialized for each backend.
 */
template<> struct one_dim_backend<backend::rocfft_sin>{
    //! \brief Defines the complex-to-complex executor.
    using executor = real2real_executor<backend::rocfft, rocm::sin_pre_pos_processor>;
    //! \brief Defines the real-to-complex executor.
    using executor_r2c = void;
};

/*!
 * \ingroup hefftepacking
 * \brief Simple packer that copies sub-boxes without transposing the order of the indexes.
 */
template<> struct direct_packer<tag::gpu>{
    //! \brief Execute the planned pack operation.
    template<typename scalar_type, typename index>
    void pack(hipStream_t stream, pack_plan_3d<index> const &plan, scalar_type const data[], scalar_type buffer[]) const{
        rocm::direct_pack(stream, plan.size[0], plan.size[1], plan.size[2], plan.line_stride, plan.plane_stride, data, buffer);
    }
    //! \brief Execute the planned unpack operation.
    template<typename scalar_type, typename index>
    void unpack(hipStream_t stream, pack_plan_3d<index> const &plan, scalar_type const buffer[], scalar_type data[]) const{
        rocm::direct_unpack(stream, plan.size[0], plan.size[1], plan.size[2], plan.line_stride, plan.plane_stride, buffer, data);
    }
};

/*!
 * \ingroup hefftepacking
 * \brief GPU version of the transpose packer.
 */
template<> struct transpose_packer<tag::gpu>{
    //! \brief Execute the planned pack operation.
    template<typename scalar_type, typename index>
    void pack(hipStream_t stream, pack_plan_3d<index> const &plan, scalar_type const data[], scalar_type buffer[]) const{
        direct_packer<tag::gpu>().pack(stream, plan, data, buffer); // packing is done the same way as the direct_packer
    }
    //! \brief Execute the planned transpose-unpack operation.
    template<typename scalar_type, typename index>
    void unpack(hipStream_t stream, pack_plan_3d<index> const &plan, scalar_type const buffer[], scalar_type data[]) const{
        rocm::transpose_unpack<scalar_type>(stream, plan.size[0], plan.size[1], plan.size[2], plan.line_stride, plan.plane_stride,
                                            plan.buff_line_stride, plan.buff_plane_stride, plan.map[0], plan.map[1], plan.map[2], buffer, data);
    }
};

namespace data_scaling {
    /*!
     * \ingroup heffterocm
     * \brief Simply multiply the \b num_entries in the \b data by the \b scale_factor.
     */
    template<typename scalar_type, typename index>
    static void apply(hipStream_t stream, index num_entries, scalar_type *data, double scale_factor){
        rocm::scale_data<scalar_type, long long>(stream, static_cast<long long>(num_entries), data, scale_factor);
    }
    /*!
     * \ingroup heffterocm
     * \brief Complex by real scaling.
     */
    template<typename precision_type, typename index>
    static void apply(hipStream_t stream, index num_entries, std::complex<precision_type> *data, double scale_factor){
        apply<precision_type>(stream, 2*num_entries, reinterpret_cast<precision_type*>(data), scale_factor);
    }
};

/*!
 * \ingroup heffterocm
 * \brief Sets the default options for the cufft backend.
 */
template<> struct default_plan_options<backend::rocfft>{
    //! \brief The reshape operations will not transpose the data.
    static const bool use_reorder = true;
};
/*!
 * \ingroup heffterocm
 * \brief Sets the default options for the cufft backend.
 */
template<> struct default_plan_options<backend::rocfft_cos>{
    //! \brief The reshape operations will not transpose the data.
    static const bool use_reorder = true;
};
/*!
 * \ingroup heffterocm
 * \brief Sets the default options for the cufft backend.
 */
template<> struct default_plan_options<backend::rocfft_sin>{
    //! \brief The reshape operations will not transpose the data.
    static const bool use_reorder = true;
};

}

#endif

#endif   /* HEFFTE_BACKEND_FFTW_H */
