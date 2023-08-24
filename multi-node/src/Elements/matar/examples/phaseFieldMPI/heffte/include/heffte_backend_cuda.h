/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_BACKEND_CUDA_H
#define HEFFTE_BACKEND_CUDA_H

#include "heffte_r2r_executor.h"

#ifdef Heffte_ENABLE_CUDA

#include <cuda_runtime_api.h>
#include <cuda.h>
#include <cufft.h>
#include "heffte_backend_vector.h"

#ifdef Heffte_ENABLE_MAGMA
#include "heffte_magma_helpers.h"
#endif

/*!
 * \ingroup fft3d
 * \addtogroup hefftecuda Backend cufft
 *
 * Wrappers and template specializations related to the cuFFT backend.
 * Requires CMake option:
 * \code
 *  -D Heffte_ENABLE_CUDA=ON
 * \endcode
 *
 * In addition to the cuFFT wrappers, this also includes a series of kernels
 * for packing/unpacking/scaling the data, as well as a simple container
 * that wraps around CUDA arrays for RAII style of resource management.
 */

namespace heffte{

/*!
 * \ingroup hefftecuda
 * \brief Cuda specific methods.
 */
namespace cuda {
    /*!
     * \ingroup hefftecuda
     * \brief Checks the status of a CUDA command and in case of a failure, converts it to a C++ exception.
     */
    inline void check_error(cudaError_t status, const char *function_name){
        if (status != cudaSuccess)
            throw std::runtime_error(std::string(function_name) + " failed with message: " + cudaGetErrorString(status));
    }
    /*!
     * \ingroup hefftecuda
     * \brief Checks the status of a cufft command and in case of a failure, converts it to a C++ exception.
     */
    inline void check_error(cufftResult status, const char *function_name){
        if (status != CUFFT_SUCCESS)
            throw std::runtime_error(std::string(function_name) + " failed with error code: " + std::to_string(status));
    }
    /*!
     * \ingroup hefftecuda
     * \brief Convert real numbers to complex when both are located on the GPU device.
     *
     * Launches a CUDA kernel.
     */
    template<typename precision_type, typename index>
    void convert(cudaStream_t stream, index num_entries, precision_type const source[], std::complex<precision_type> destination[]);
    /*!
     * \ingroup hefftecuda
     * \brief Convert complex numbers to real when both are located on the GPU device.
     *
     * Launches a CUDA kernel.
     */
    template<typename precision_type, typename index>
    void convert(cudaStream_t stream, index num_entries, std::complex<precision_type> const source[], precision_type destination[]);

    /*!
     * \ingroup hefftecuda
     * \brief Scales real data (double or float) by the scaling factor.
     */
    template<typename scalar_type, typename index>
    void scale_data(cudaStream_t stream, index num_entries, scalar_type *data, double scale_factor);

    /*!
     * \ingroup hefftecuda
     * \brief Performs a direct-pack operation for data sitting on the GPU device.
     *
     * Launches a CUDA kernel.
     */
    template<typename scalar_type, typename index>
    void direct_pack(cudaStream_t stream, index nfast, index nmid, index nslow, index line_stride, index plane_stide, scalar_type const source[], scalar_type destination[]);
    /*!
     * \ingroup hefftecuda
     * \brief Performs a direct-unpack operation for data sitting on the GPU device.
     *
     * Launches a CUDA kernel.
     */
    template<typename scalar_type, typename index>
    void direct_unpack(cudaStream_t stream, index nfast, index nmid, index nslow, index line_stride, index plane_stide, scalar_type const source[], scalar_type destination[]);
    /*!
     * \ingroup hefftecuda
     * \brief Performs a transpose-unpack operation for data sitting on the GPU device.
     *
     * Launches a CUDA kernel.
     */
    template<typename scalar_type, typename index>
    void transpose_unpack(cudaStream_t stream, index nfast, index nmid, index nslow, index line_stride, index plane_stide,
                        index buff_line_stride, index buff_plane_stride, int map0, int map1, int map2,
                        scalar_type const source[], scalar_type destination[]);

    /*!
     * \ingroup hefftecuda
     * \brief Implementation of Cosine Transform pre-post processing methods using CUDA.
     */
    struct cos_pre_pos_processor{
        //! \brief Pre-process in the forward transform.
        template<typename precision>
        static void pre_forward(cudaStream_t, int length, precision const input[], precision fft_signal[]);
        //! \brief Post-process in the forward transform.
        template<typename precision>
        static void post_forward(cudaStream_t, int length, std::complex<precision> const fft_result[], precision result[]);
        //! \brief Pre-process in the inverse transform.
        template<typename precision>
        static void pre_backward(cudaStream_t, int length, precision const input[], std::complex<precision> fft_signal[]);
        //! \brief Post-process in the inverse transform.
        template<typename precision>
        static void post_backward(cudaStream_t, int length, precision const fft_result[], precision result[]);
    };
    /*!
     * \ingroup hefftecuda
     * \brief Implementation of Sine Transform pre-post processing methods using CUDA.
     */
    struct sin_pre_pos_processor{
        //! \brief Pre-process in the forward transform.
        template<typename precision>
        static void pre_forward(cudaStream_t, int length, precision const input[], precision fft_signal[]);
        //! \brief Post-process in the forward transform.
        template<typename precision>
        static void post_forward(cudaStream_t, int length, std::complex<precision> const fft_result[], precision result[]);
        //! \brief Pre-process in the inverse transform.
        template<typename precision>
        static void pre_backward(cudaStream_t, int length, precision const input[], std::complex<precision> fft_signal[]);
        //! \brief Post-process in the inverse transform.
        template<typename precision>
        static void post_backward(cudaStream_t, int length, precision const fft_result[], precision result[]);
    };
}

namespace backend{
    /*!
     * \ingroup hefftecuda
     * \brief Indicate that the cuFFT backend has been enabled.
     */
    template<> struct is_enabled<cufft> : std::true_type{};
    /*!
     * \ingroup hefftecuda
     * \brief Indicate that the cuFFT backend has been enabled for Cosine Transform.
     */
    template<> struct is_enabled<cufft_cos> : std::true_type{};
    /*!
     * \ingroup hefftecuda
     * \brief Indicate that the cuFFT backend has been enabled for Sine Transform.
     */
    template<> struct is_enabled<cufft_sin> : std::true_type{};

    /*!
     * \ingroup hefftecuda
     * \brief The CUDA backend uses a CUDA stream.
     */
    template<>
    struct device_instance<tag::gpu>{
        //! \brief Constructor, sets up the stream.
        device_instance(cudaStream_t new_stream = nullptr) : _stream(new_stream){}
        //! \brief Returns the nullptr.
        cudaStream_t stream(){ return _stream; }
        //! \brief Returns the nullptr (const case).
        cudaStream_t stream() const{ return _stream; }
        //! \brief Syncs the execution with the queue, no-op in the CPU case.
        void synchronize_device() const{ cuda::check_error(cudaStreamSynchronize(_stream), "device sync"); }
        //! \brief The CUDA stream to be used in all operations.
        mutable cudaStream_t _stream;
        //! \brief The type for the internal stream.
        using stream_type = cudaStream_t;
    };

    /*!
     * \ingroup hefftecuda
     * \brief In CUDA mode, the default GPU backend is cufft.
     */
    template<> struct default_backend<tag::gpu>{
        //! \brief Set the cufft tag.
        using type = cufft;
    };

    /*!
     * \ingroup hefftecuda
     * \brief Specialization for the data operations in CUDA mode.
     */
    template<> struct data_manipulator<tag::gpu> {
        //! \brief The stream type for the device.
        using stream_type = cudaStream_t;
        //! \brief Defines the backend_device.
        using backend_device = backend::device_instance<tag::gpu>;
        //! \brief Allocate memory.
        template<typename scalar_type>
        static scalar_type* allocate(cudaStream_t stream, size_t num_entries){
            void *new_data;
            if (stream != nullptr) cuda::check_error( cudaStreamSynchronize(stream), "cudaStreamSynchronize()");
            cuda::check_error(cudaMalloc(&new_data, num_entries * sizeof(scalar_type)), "cudaMalloc()");
            return reinterpret_cast<scalar_type*>(new_data);
        }
        //! \brief Free memory.
        template<typename scalar_type>
        static void free(cudaStream_t stream, scalar_type *pntr){
            if (pntr == nullptr) return;
            if (stream != nullptr) cuda::check_error( cudaStreamSynchronize(stream), "cudaStreamSynchronize()");
            cuda::check_error(cudaFree(pntr), "cudaFree()");
        }
        //! \brief Equivalent to std::copy_n() but using CUDA arrays.
        template<typename scalar_type>
        static void copy_n(cudaStream_t stream, scalar_type const source[], size_t num_entries, scalar_type destination[]){
            if (stream == nullptr)
                cuda::check_error(cudaMemcpy(destination, source, num_entries * sizeof(scalar_type), cudaMemcpyDeviceToDevice), "data_manipulator::copy_n()");
            else
                cuda::check_error(cudaMemcpyAsync(destination, source, num_entries * sizeof(scalar_type), cudaMemcpyDeviceToDevice, stream), "data_manipulator::copy_n()");
        }
        //! \brief Copy-convert complex-to-real.
        template<typename scalar_type>
        static void copy_n(cudaStream_t stream, std::complex<scalar_type> const source[], size_t num_entries, scalar_type destination[]){
            cuda::convert(stream, static_cast<long long>(num_entries), source, destination);
        }
        //! \brief Copy-convert real-to-complex.
        template<typename scalar_type>
        static void copy_n(cudaStream_t stream, scalar_type const source[], size_t num_entries, std::complex<scalar_type> destination[]){
            cuda::convert(stream, static_cast<long long>(num_entries), source, destination);
        }
        //! \brief Copy the date from the device to the host.
        template<typename scalar_type>
        static void copy_device_to_host(cudaStream_t stream, scalar_type const source[], size_t num_entries, scalar_type destination[]){
            cuda::check_error(cudaMemcpyAsync(destination, source, num_entries * sizeof(scalar_type), cudaMemcpyDeviceToHost, stream),
                            "device_to_host (cuda)");
        }
        //! \brief Copy the date from the device to the device.
        template<typename scalar_type>
        static void copy_device_to_device(cudaStream_t stream, scalar_type const source[], size_t num_entries, scalar_type destination[]){
                cuda::check_error(cudaMemcpyAsync(destination, source, num_entries * sizeof(scalar_type), cudaMemcpyDeviceToDevice, stream),
                                "device_to_device (cuda)");
        }
        //! \brief Copy the date from the host to the device.
        template<typename scalar_type>
        static void copy_host_to_device(cudaStream_t stream, scalar_type const source[], size_t num_entries, scalar_type destination[]){
            cuda::check_error(cudaMemcpyAsync(destination, source, num_entries * sizeof(scalar_type), cudaMemcpyHostToDevice, stream),
                            "host_to_device (cuda)");
        }
    };

    /*!
     * \ingroup hefftecuda
     * \brief Defines the location type-tag and the cuda container.
     */
    template<>
    struct buffer_traits<cufft>{
        //! \brief The cufft library uses data on the gpu device.
        using location = tag::gpu;
        //! \brief The data is managed by the cuda vector container.
        template<typename T> using container = heffte::gpu::device_vector<T, data_manipulator<tag::gpu>>;
    };
    /*!
     * \ingroup hefftecuda
     * \brief Defines the location type-tag and the cuda container.
     */
    template<>
    struct buffer_traits<cufft_cos>{
        //! \brief The cufft library uses data on the gpu device.
        using location = tag::gpu;
        //! \brief The data is managed by the cuda vector container.
        template<typename T> using container = heffte::gpu::device_vector<T, data_manipulator<tag::gpu>>;
    };
    /*!
     * \ingroup hefftecuda
     * \brief Defines the location type-tag and the cuda container.
     */
    template<>
    struct buffer_traits<cufft_sin>{
        //! \brief The cufft library uses data on the gpu device.
        using location = tag::gpu;
        //! \brief The data is managed by the cuda vector container.
        template<typename T> using container = heffte::gpu::device_vector<T, data_manipulator<tag::gpu>>;
    };
}

/*!
 * \ingroup hefftecuda
 * \brief Recognize the cuFFT single precision complex type.
 */
template<> struct is_ccomplex<cufftComplex> : std::true_type{};
/*!
 * \ingroup hefftecuda
 * \brief Recognize the cuFFT double precision complex type.
 */
template<> struct is_zcomplex<cufftDoubleComplex> : std::true_type{};

/*!
 * \ingroup hefftecuda
 * \brief Wrapper around cufftHandle plans, set for float or double complex.
 *
 *  * \tparam scalar_type must be std::compelx<float> or std::complex<double>
 */
template<typename scalar_type> struct plan_cufft{
    /*!
     * \brief Constructor, takes inputs identical to cufftMakePlanMany().
     *
     * \param stream is the CUDA stream to use for the transform
     * \param size is the number of entries in a 1-D transform
     * \param batch is the number of transforms in the batch
     * \param stride is the distance between entries of the same transform
     * \param dist is the distance between the first entries of consecutive sequences
     */
    plan_cufft(cudaStream_t stream, int size, int batch, int stride, int dist){
        size_t work_size = 0;
        cuda::check_error(cufftCreate(&plan), "plan_cufft::cufftCreate()");
        cuda::check_error(
            cufftMakePlanMany(plan, 1, &size, &size, stride, dist, &size, stride, dist,
                              (std::is_same<scalar_type, std::complex<float>>::value) ? CUFFT_C2C : CUFFT_Z2Z,
                              batch, &work_size),
            "plan_cufft::cufftMakePlanMany() 1D"
        );

        if (stream != nullptr)
            cuda::check_error( cufftSetStream(plan, stream), "cufftSetStream()");
    }
    /*!
     * \brief Constructor, takes inputs identical to cufftMakePlanMany().
     *
     * \param stream is the CUDA stream to use for the transform
     * \param size1 is the number of entries in a 2-D transform (dimension 1)
     * \param size2 is the number of entries in a 2-D transform (dimension 2)
     * \param embed defines the embedding of the 2-D transform into a 3-D box
     * \param batch is the number of transforms in the batch
     * \param stride is the distance between entries of the same transform
     * \param dist is the distance between the first entries of consecutive sequences
     */
    plan_cufft(cudaStream_t stream, int size1, int size2, std::array<int, 2> const &embed, int batch, int stride, int dist){
        size_t work_size = 0;
        int size[2] = {size2, size1};

        cuda::check_error(cufftCreate(&plan), "plan_cufft::cufftCreate()");
        if (embed[0] == 0 and embed[1] == 0){
            cuda::check_error(
                cufftMakePlanMany(plan, 2, size, size, stride, dist, size, stride, dist,
                                (std::is_same<scalar_type, std::complex<float>>::value) ? CUFFT_C2C : CUFFT_Z2Z,
                                batch, &work_size),
                "plan_cufft::cufftMakePlanMany() 2D"
            );
        }else{
            cuda::check_error(
                cufftMakePlanMany(plan, 2, size, const_cast<int*>(embed.data()), stride, dist, const_cast<int*>(embed.data()), stride, dist,
                                (std::is_same<scalar_type, std::complex<float>>::value) ? CUFFT_C2C : CUFFT_Z2Z,
                                batch, &work_size),
                "plan_cufft::cufftMakePlanMany() 2D with embedding"
            );
        }

        if (stream != nullptr)
            cuda::check_error( cufftSetStream(plan, stream), "cufftSetStream()");
    }
    //! \brief Constructor, takes inputs identical to cufftPlan3d()
    plan_cufft(cudaStream_t stream, int size1, int size2, int size3){
        cuda::check_error(
            cufftPlan3d(&plan, size3, size2, size1,
                        (std::is_same<scalar_type, std::complex<float>>::value) ? CUFFT_C2C : CUFFT_Z2Z),
            "plan_cufft::cufftPlan3d()"
        );
        if (stream != nullptr)
            cuda::check_error( cufftSetStream(plan, stream), "cufftSetStream()");
    }
    //! \brief Destructor, deletes the plan.
    ~plan_cufft(){ cufftDestroy(plan); }
    //! \brief Custom conversion to the cufftHandle.
    operator cufftHandle() const{ return plan; }

private:
    //! \brief The cufft opaque structure (pointer to struct).
    cufftHandle plan;
};

/*!
 * \ingroup hefftecuda
 * \brief Wrapper around the cuFFT API.
 *
 * A single class that manages the plans and executions of cuFFT
 * so that a single API is provided for all backends.
 * The executor operates on a box and performs 1-D FFTs
 * for the given dimension.
 * The class silently manages the plans and buffers needed
 * for the different types.
 * All input and output arrays must have size equal to the box.
 */
class cufft_executor : public executor_base{
public:
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::forward;
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::backward;
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::complex_size;
    //! \brief Constructor, specifies the box and dimension.
    template<typename index>
    cufft_executor(cudaStream_t active_stream, box3d<index> const box, int dimension) :
        stream(active_stream),
        size(box.size[dimension]), size2(0),
        howmanyffts(fft1d_get_howmany(box, dimension)),
        stride(fft1d_get_stride(box, dimension)),
        dist((dimension == box.order[0]) ? size : 1),
        blocks((dimension == box.order[1]) ? box.osize(2) : 1),
        block_stride(box.osize(0) * box.osize(1)),
        total_size(box.count()),
        embed({0, 0})
    {}
    //! \brief Merges two FFTs into one.
    template<typename index>
    cufft_executor(cudaStream_t active_stream, box3d<index> const box, int dir1, int dir2) :
        stream(active_stream),
        size(box.size[std::min(dir1, dir2)]), size2(box.size[std::max(dir1, dir2)]),
        blocks(1), block_stride(0), total_size(box.count()), embed({0, 0})
    {
        int odir1 = box.find_order(dir1);
        int odir2 = box.find_order(dir2);

        if (std::min(odir1, odir2) == 0 and std::max(odir1, odir2) == 1){
            stride = 1;
            dist = size * size2;
            howmanyffts = box.size[2];
        }else if (std::min(odir1, odir2) == 1 and std::max(odir1, odir2) == 2){
            stride = box.size[0];
            dist = 1;
            howmanyffts = box.size[0];
        }else{ // case of directions (0, 2)
            stride = 1;
            dist = size;
            embed = {static_cast<int>(box.size[2]), static_cast<int>(box.size[1] * box.size[0])};
            howmanyffts = box.size[1];
        }
    }
    //! \brief Merges three FFTs into one.
    template<typename index>
    cufft_executor(cudaStream_t active_stream, box3d<index> const box) :
        stream(active_stream),
        size(box.size[0]), size2(box.size[1]), howmanyffts(box.size[2]),
        stride(0), dist(0),
        blocks(1), block_stride(0),
        total_size(box.count()),
        embed({0, 0})
    {}

    //! \brief Forward fft, float-complex case.
    void forward(std::complex<float> data[], std::complex<float>*) const override{
        make_plan(ccomplex_plan);
        for(int i=0; i<blocks; i++){
            cufftComplex* block_data = reinterpret_cast<cufftComplex*>(data + i * block_stride);
            cuda::check_error(cufftExecC2C(*ccomplex_plan, block_data, block_data, CUFFT_FORWARD), "cufft_executor::cufftExecC2C() forward");
        }
    }
    //! \brief Backward fft, float-complex case.
    void backward(std::complex<float> data[], std::complex<float>*) const override{
        make_plan(ccomplex_plan);
        for(int i=0; i<blocks; i++){
            cufftComplex* block_data = reinterpret_cast<cufftComplex*>(data + i * block_stride);
            cuda::check_error(cufftExecC2C(*ccomplex_plan, block_data, block_data, CUFFT_INVERSE), "cufft_executor::cufftExecC2C() backward");
        }
    }
    //! \brief Forward fft, double-complex case.
    void forward(std::complex<double> data[], std::complex<double>*) const override{
        make_plan(zcomplex_plan);
        for(int i=0; i<blocks; i++){
            cufftDoubleComplex* block_data = reinterpret_cast<cufftDoubleComplex*>(data + i * block_stride);
            cuda::check_error(cufftExecZ2Z(*zcomplex_plan, block_data, block_data, CUFFT_FORWARD), "cufft_executor::cufftExecZ2Z() forward");
        }
    }
    //! \brief Backward fft, double-complex case.
    void backward(std::complex<double> data[], std::complex<double>*) const override{
        make_plan(zcomplex_plan);
        for(int i=0; i<blocks; i++){
            cufftDoubleComplex* block_data = reinterpret_cast<cufftDoubleComplex*>(data + i * block_stride);
            cuda::check_error(cufftExecZ2Z(*zcomplex_plan, block_data, block_data, CUFFT_INVERSE), "cufft_executor::cufftExecZ2Z() backward");
        }
    }

    //! \brief Converts the deal data to complex and performs float-complex forward transform.
    void forward(float const indata[], std::complex<float> outdata[], std::complex<float> *workspace) const override{
        cuda::convert(stream, total_size, indata, outdata);
        forward(outdata, workspace);
    }
    //! \brief Performs backward float-complex transform and truncates the complex part of the result.
    void backward(std::complex<float> indata[], float outdata[], std::complex<float> *workspace) const{
        backward(indata, workspace);
        cuda::convert(stream, total_size, indata, outdata);
    }
    //! \brief Converts the deal data to complex and performs double-complex forward transform.
    void forward(double const indata[], std::complex<double> outdata[], std::complex<double> *workspace) const override{
        cuda::convert(stream, total_size, indata, outdata);
        forward(outdata, workspace);
    }
    //! \brief Performs backward double-complex transform and truncates the complex part of the result.
    void backward(std::complex<double> indata[], double outdata[], std::complex<double> *workspace) const override{
        backward(indata, workspace);
        cuda::convert(stream, total_size, indata, outdata);
    }

    //! \brief Returns the size of the box.
    int box_size() const override{ return total_size; }
    //! \brief Return the size of the needed workspace.
    size_t workspace_size() const override{ return 0; }

private:
    //! \brief Helper template to create the plan.
    template<typename scalar_type>
    void make_plan(std::unique_ptr<plan_cufft<scalar_type>> &plan) const{
        if (not plan){
            if (dist == 0)
                plan = std::unique_ptr<plan_cufft<scalar_type>>(new plan_cufft<scalar_type>(stream, size, size2, howmanyffts));
            else if (size2 == 0)
                plan = std::unique_ptr<plan_cufft<scalar_type>>(new plan_cufft<scalar_type>(stream, size, howmanyffts, stride, dist));
            else
                plan = std::unique_ptr<plan_cufft<scalar_type>>(new plan_cufft<scalar_type>(stream, size, size2, embed, howmanyffts, stride, dist));
        }
    }

    mutable cudaStream_t stream;

    int size, size2, howmanyffts, stride, dist, blocks, block_stride, total_size;
    std::array<int, 2> embed;
    mutable std::unique_ptr<plan_cufft<std::complex<float>>> ccomplex_plan;
    mutable std::unique_ptr<plan_cufft<std::complex<double>>> zcomplex_plan;
};

/*!
 * \ingroup hefftecuda
 * \brief Plan for the r2c single and double precision transform.
 *
 * \tparam scalar_type must be float or double
 */
template<typename scalar_type> struct plan_cufft_r2c{
    /*!
     * \brief Constructor, takes inputs identical to cufftMakePlanMany().
     *
     * \param stream is the CUDA stream to use for the transform
     * \param dir is the direction (forward or backward) for the plan
     * \param size is the number of entries in a 1-D transform
     * \param batch is the number of transforms in the batch
     * \param stride is the distance between entries of the same transform
     * \param rdist is the distance between the first entries of consecutive real sequences
     * \param cdist is the distance between the first entries of consecutive complex sequences
     */
    plan_cufft_r2c(cudaStream_t stream, direction dir, int size, int batch, int stride, int rdist, int cdist){
        static_assert(std::is_same<scalar_type, float>::value or std::is_same<scalar_type, double>::value,
                      "plan_cufft_r2c can be used only with scalar_type float or double.");

        size_t work_size = 0;
        cuda::check_error(cufftCreate(&plan), "plan_cufft_r2c::cufftCreate()");

        cuda::check_error(
                cufftMakePlanMany(plan, 1, &size, &size, stride,
                                  (dir == direction::forward) ? rdist : cdist,
                                  &size, stride,
                                  (dir == direction::forward) ? cdist : rdist,
                                  (std::is_same<scalar_type, float>::value) ?
                                        (dir == direction::forward) ? CUFFT_R2C : CUFFT_C2R
                                      : (dir == direction::forward) ? CUFFT_D2Z : CUFFT_Z2D,
                                  batch, &work_size),
                "plan_cufft_r2c::cufftMakePlanMany() (forward)"
            );
        if (stream != nullptr)
            cuda::check_error( cufftSetStream(plan, stream), "cufftSetStream()");
    }
    //! \brief Destructor, deletes the plan.
    ~plan_cufft_r2c(){ cufftDestroy(plan); }
    //! \brief Custom conversion to the cufftHandle.
    operator cufftHandle() const{ return plan; }

private:
    //! \brief The cufft opaque structure (pointer to struct).
    cufftHandle plan;
};

/*!
 * \ingroup hefftecuda
 * \brief Wrapper to cuFFT API for real-to-complex transform with shortening of the data.
 *
 * Serves the same purpose of heffte::cufft_executor but only real input is accepted
 * and only the unique (non-conjugate) coefficients are computed.
 * All real arrays must have size of real_size() and all complex arrays must have size complex_size().
 */
class cufft_executor_r2c : public executor_base{
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
    cufft_executor_r2c(cudaStream_t active_stream, box3d<index> const box, int dimension) :
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
        csize(box.r2c(dimension).count())
    {}

    //! \brief Forward transform, single precision.
    void forward(float const indata[], std::complex<float> outdata[], std::complex<float>*) const override{
        make_plan(sforward, direction::forward);
        if (blocks == 1 or rblock_stride % 2 == 0){
            for(int i=0; i<blocks; i++){
                cufftReal *rdata = const_cast<cufftReal*>(indata + i * rblock_stride);
                cufftComplex* cdata = reinterpret_cast<cufftComplex*>(outdata + i * cblock_stride);
                cuda::check_error(cufftExecR2C(*sforward, rdata, cdata), "cufft_executor::cufftExecR2C()");
            }
        }else{
            // need to create a temporary copy of the data since cufftExecR2C() requires aligned input
            backend::buffer_traits<backend::cufft>::container<float> rdata(stream, rblock_stride);
            for(int i=0; i<blocks; i++){
                backend::data_manipulator<tag::gpu>::copy_n(stream, indata + i * rblock_stride, rblock_stride, rdata.data());
                cufftComplex* cdata = reinterpret_cast<cufftComplex*>(outdata + i * cblock_stride);
                cuda::check_error(cufftExecR2C(*sforward, rdata.data(), cdata), "cufft_executor::cufftExecR2C()");
            }
        }
    }
    //! \brief Backward transform, single precision.
    void backward(std::complex<float> indata[], float outdata[], std::complex<float>*) const override{
        make_plan(sbackward, direction::backward);
        if (blocks == 1 or rblock_stride % 2 == 0){
            for(int i=0; i<blocks; i++){
                cufftComplex* cdata = const_cast<cufftComplex*>(reinterpret_cast<cufftComplex const*>(indata + i * cblock_stride));
                cuda::check_error(cufftExecC2R(*sbackward, cdata, outdata + i * rblock_stride), "cufft_executor::cufftExecC2R()");
            }
        }else{
            backend::buffer_traits<backend::cufft>::container<float> odata(stream, rblock_stride);
            for(int i=0; i<blocks; i++){
                cufftComplex* cdata = const_cast<cufftComplex*>(reinterpret_cast<cufftComplex const*>(indata + i * cblock_stride));
                cuda::check_error(cufftExecC2R(*sbackward, cdata, odata.data()), "cufft_executor::cufftExecC2R()");
                backend::data_manipulator<tag::gpu>::copy_n(stream, odata.data(), rblock_stride, outdata + i * rblock_stride);
            }
        }
    }
    //! \brief Forward transform, double precision.
    void forward(double const indata[], std::complex<double> outdata[], std::complex<double>*) const override{
        make_plan(dforward, direction::forward);
        if (blocks == 1 or rblock_stride % 2 == 0){
            for(int i=0; i<blocks; i++){
                cufftDoubleReal *rdata = const_cast<cufftDoubleReal*>(indata + i * rblock_stride);
                cufftDoubleComplex* cdata = reinterpret_cast<cufftDoubleComplex*>(outdata + i * cblock_stride);
                cuda::check_error(cufftExecD2Z(*dforward, rdata, cdata), "cufft_executor::cufftExecD2Z()");
            }
        }else{
            backend::buffer_traits<backend::cufft>::container<double> rdata(stream, rblock_stride);
            for(int i=0; i<blocks; i++){
                backend::data_manipulator<tag::gpu>::copy_n(stream, indata + i * rblock_stride, rblock_stride, rdata.data());
                cufftDoubleComplex* cdata = reinterpret_cast<cufftDoubleComplex*>(outdata + i * cblock_stride);
                cuda::check_error(cufftExecD2Z(*dforward, rdata.data(), cdata), "cufft_executor::cufftExecD2Z()");
            }
        }
    }
    //! \brief Backward transform, double precision.
    void backward(std::complex<double> indata[], double outdata[], std::complex<double>*) const override{
        make_plan(dbackward, direction::backward);
        if (blocks == 1 or rblock_stride % 2 == 0){
            for(int i=0; i<blocks; i++){
                cufftDoubleComplex* cdata = const_cast<cufftDoubleComplex*>(reinterpret_cast<cufftDoubleComplex const*>(indata + i * cblock_stride));
                cuda::check_error(cufftExecZ2D(*dbackward, cdata, outdata + i * rblock_stride), "cufft_executor::cufftExecZ2D()");
            }
        }else{
            backend::buffer_traits<backend::cufft>::container<double> odata(stream, rblock_stride);
            for(int i=0; i<blocks; i++){
                cufftDoubleComplex* cdata = const_cast<cufftDoubleComplex*>(reinterpret_cast<cufftDoubleComplex const*>(indata + i * cblock_stride));
                cuda::check_error(cufftExecZ2D(*dbackward, cdata, odata.data()), "cufft_executor::cufftExecZ2D()");
                backend::data_manipulator<tag::gpu>::copy_n(stream, odata.data(), rblock_stride, outdata + i * rblock_stride);
            }
        }
    }

    //! \brief Returns the size of the box with real data.
    int box_size() const override{ return rsize; }
    //! \brief Returns the size of the box with complex coefficients.
    int complex_size() const override{ return csize; }
    //! \brief Return the size of the needed workspace.
    size_t workspace_size() const override{ return 0; }

private:
    //! \brief Helper template to initialize the plan.
    template<typename scalar_type>
    void make_plan(std::unique_ptr<plan_cufft_r2c<scalar_type>> &plan, direction dir) const{
        if (!plan) plan = std::unique_ptr<plan_cufft_r2c<scalar_type>>(new plan_cufft_r2c<scalar_type>(stream, dir, size, howmanyffts, stride, rdist, cdist));
    }

    cudaStream_t stream;

    int size, howmanyffts, stride, blocks;
    int rdist, cdist, rblock_stride, cblock_stride, rsize, csize;
    mutable std::unique_ptr<plan_cufft_r2c<float>> sforward;
    mutable std::unique_ptr<plan_cufft_r2c<double>> dforward;
    mutable std::unique_ptr<plan_cufft_r2c<float>> sbackward;
    mutable std::unique_ptr<plan_cufft_r2c<double>> dbackward;
};

/*!
 * \ingroup hefftecuda
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 *
 * The struct is specialized for each backend.
 */
template<> struct one_dim_backend<backend::cufft>{
    //! \brief Defines the complex-to-complex executor.
    using executor = cufft_executor;
    //! \brief Defines the real-to-complex executor.
    using executor_r2c = cufft_executor_r2c;
};

/*!
 * \ingroup hefftecuda
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 *
 * The struct is specialized for each backend.
 */
template<> struct one_dim_backend<backend::cufft_cos>{
    //! \brief Defines the complex-to-complex executor.
    using executor = real2real_executor<backend::cufft, cuda::cos_pre_pos_processor>;
    //! \brief Defines the real-to-complex executor.
    using executor_r2c = void;
};
/*!
 * \ingroup hefftecuda
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 *
 * The struct is specialized for each backend.
 */
template<> struct one_dim_backend<backend::cufft_sin>{
    //! \brief Defines the complex-to-complex executor.
    using executor = real2real_executor<backend::cufft, cuda::sin_pre_pos_processor>;
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
    void pack(cudaStream_t stream, pack_plan_3d<index> const &plan, scalar_type const data[], scalar_type buffer[]) const{
        cuda::direct_pack(stream, plan.size[0], plan.size[1], plan.size[2], plan.line_stride, plan.plane_stride, data, buffer);
    }
    //! \brief Execute the planned unpack operation.
    template<typename scalar_type, typename index>
    void unpack(cudaStream_t stream, pack_plan_3d<index> const &plan, scalar_type const buffer[], scalar_type data[]) const{
        cuda::direct_unpack(stream, plan.size[0], plan.size[1], plan.size[2], plan.line_stride, plan.plane_stride, buffer, data);
    }
};

/*!
 * \ingroup hefftepacking
 * \brief GPU version of the transpose packer.
 */
template<> struct transpose_packer<tag::gpu>{
    //! \brief Execute the planned pack operation.
    template<typename scalar_type, typename index>
    void pack(cudaStream_t stream, pack_plan_3d<index> const &plan, scalar_type const data[], scalar_type buffer[]) const{
        direct_packer<tag::gpu>().pack(stream, plan, data, buffer); // packing is done the same way as the direct_packer
    }
    //! \brief Execute the planned transpose-unpack operation.
    template<typename scalar_type, typename index>
    void unpack(cudaStream_t stream, pack_plan_3d<index> const &plan, scalar_type const buffer[], scalar_type data[]) const{
        cuda::transpose_unpack<scalar_type>(stream, plan.size[0], plan.size[1], plan.size[2], plan.line_stride, plan.plane_stride,
                                            plan.buff_line_stride, plan.buff_plane_stride, plan.map[0], plan.map[1], plan.map[2], buffer, data);
    }
};

namespace data_scaling {
    /*!
     * \ingroup hefftecuda
     * \brief Simply multiply the \b num_entries in the \b data by the \b scale_factor.
     */
    template<typename scalar_type, typename index>
    void apply(cudaStream_t stream, index num_entries, scalar_type *data, double scale_factor){
        cuda::scale_data(stream, static_cast<long long>(num_entries), data, scale_factor);
    }
    /*!
     * \ingroup hefftecuda
     * \brief Complex by real scaling.
     */
    template<typename precision_type, typename index>
    void apply(cudaStream_t stream, index num_entries, std::complex<precision_type> *data, double scale_factor){
        apply<precision_type>(stream, 2*num_entries, reinterpret_cast<precision_type*>(data), scale_factor);
    }
}

/*!
 * \ingroup hefftecuda
 * \brief Sets the default options for the cufft backend.
 */
template<> struct default_plan_options<backend::cufft>{
    //! \brief The reshape operations will not transpose the data.
    static const bool use_reorder = false;
};
/*!
 * \ingroup hefftecuda
 * \brief Sets the default options for the cufft backend.
 */
template<> struct default_plan_options<backend::cufft_cos>{
    //! \brief The reshape operations will not transpose the data.
    static const bool use_reorder = true;
};
/*!
 * \ingroup hefftecuda
 * \brief Sets the default options for the cufft backend.
 */
template<> struct default_plan_options<backend::cufft_sin>{
    //! \brief The reshape operations will not transpose the data.
    static const bool use_reorder = true;
};

}

#endif

#endif   /* HEFFTE_BACKEND_FFTW_H */
