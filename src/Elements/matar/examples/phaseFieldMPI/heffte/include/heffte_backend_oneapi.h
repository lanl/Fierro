/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_BACKEND_ONEAPI_H
#define HEFFTE_BACKEND_ONEAPI_H

#include "heffte_r2r_executor.h"

#ifdef Heffte_ENABLE_ONEAPI

#include "heffte_backend_vector.h"

#include <CL/sycl.hpp>
#include "oneapi/mkl.hpp"
#include "oneapi/mkl/dfti.hpp"

#ifdef Heffte_ENABLE_MAGMA
// will enable once MAGMA has a DPC++/SYCL backend
//#include "heffte_magma_helpers.h"
#endif

/*!
 * \ingroup fft3d
 * \addtogroup heffteoneapi Backend oneAPI
 *
 * Wrappers and template specializations related to the oneMKL backend.
 * Requires CMake option:
 * \code
 *  -D Heffte_ENABLE_ONEAPI=ON
 * \endcode
 *
 * In addition to the oneMKL wrappers, this also includes a series of kernels
 * for packing/unpacking/scaling the data, as well as a simple container
 * that wraps around SYCL arrays for RAII style of resource management.
 */

namespace heffte{

/*!
 * \ingroup heffteoneapi
 * \brief SYCL/DPC++ specific methods, vector-like container, error checking, etc.
 *
 * The name is chosen distinct from the oneMKL name that use "oneapi".
 */
namespace oapi {
    //! \brief Creates a new SYCL queue, try to use the GPU but if an issue is encountered then default to the CPU.
    inline sycl::queue make_sycl_queue(){
        try{
            return sycl::queue(sycl::gpu_selector());
        }catch(sycl::exception const&){
            return sycl::queue(sycl::cpu_selector());
        }
    }

    /*!
     * \ingroup heffteoneapi
     * \brief Default queue to use in case the user does not provide one.
     *
     * The SYCL/DPC++ standard does not provide a default queue, unlike CUDA and ROCm.
     * For API consistency, heFFTe will use this queue which will be shard by all objects created
     * by the library. However, it is strongly recommended that the user provide their own
     * queue in the constructor of the heffte::fft3d and heffte::fft3d_r2c objects.
     */
    extern sycl::queue internal_sycl_queue;

    /*!
     * \ingroup heffteoneapi
     * \brief Convert real numbers to complex when both are located on the GPU device.
     *
     * Launches a SYCL/DPC++ kernel.
     */
    template<typename precision_type, typename index>
    void convert(sycl::queue &stream, index num_entries, precision_type const source[], std::complex<precision_type> destination[]);
    /*!
     * \ingroup heffteoneapi
     * \brief Convert complex numbers to real when both are located on the GPU device.
     *
     * Launches a SYCL/DPC++ kernel.
     */
    template<typename precision_type, typename index>
    void convert(sycl::queue &stream, index num_entries, std::complex<precision_type> const source[], precision_type destination[]);

    /*!
     * \ingroup heffteoneapi
     * \brief Scales real data (double or float) by the scaling factor.
     */
    template<typename scalar_type, typename index>
    void scale_data(sycl::queue &stream, index num_entries, scalar_type *data, double scale_factor);

    /*!
     * \ingroup heffteoneapi
     * \brief Performs a direct-pack operation for data sitting on the GPU device.
     *
     * Launches a SYCL/DPC++ kernel.
     */
    template<typename scalar_type, typename index>
    void direct_pack(sycl::queue &stream, index nfast, index nmid, index nslow, index line_stride, index plane_stide, scalar_type const source[], scalar_type destination[]);
    /*!
     * \ingroup heffteoneapi
     * \brief Performs a direct-unpack operation for data sitting on the GPU device.
     *
     * Launches a SYCL/DPC++ kernel.
     */
    template<typename scalar_type, typename index>
    void direct_unpack(sycl::queue &stream, index nfast, index nmid, index nslow, index line_stride, index plane_stide, scalar_type const source[], scalar_type destination[]);
    /*!
     * \ingroup heffteoneapi
     * \brief Performs a transpose-unpack operation for data sitting on the GPU device.
     *
     * Launches a SYCL/DPC++ kernel.
     */
    template<typename scalar_type, typename index>
    void transpose_unpack(sycl::queue &stream, index nfast, index nmid, index nslow, index line_stride, index plane_stide,
                        index buff_line_stride, index buff_plane_stride, int map0, int map1, int map2,
                        scalar_type const source[], scalar_type destination[]);

    /*!
     * \ingroup hefftecuda
     * \brief Implementation of Cosine Transform pre-post processing methods using CUDA.
     */
    struct cos_pre_pos_processor{
        //! \brief Pre-process in the forward transform.
        template<typename precision>
        static void pre_forward(sycl::queue&, int length, precision const input[], precision fft_signal[]);
        //! \brief Post-process in the forward transform.
        template<typename precision>
        static void post_forward(sycl::queue&, int length, std::complex<precision> const fft_result[], precision result[]);
        //! \brief Pre-process in the inverse transform.
        template<typename precision>
        static void pre_backward(sycl::queue&, int length, precision const input[], std::complex<precision> fft_signal[]);
        //! \brief Post-process in the inverse transform.
        template<typename precision>
        static void post_backward(sycl::queue&, int length, precision const fft_result[], precision result[]);
    };
    /*!
     * \ingroup hefftecuda
     * \brief Implementation of Cosine Transform pre-post processing methods using CUDA.
     */
    struct sin_pre_pos_processor{
        //! \brief Pre-process in the forward transform.
        template<typename precision>
        static void pre_forward(sycl::queue&, int length, precision const input[], precision fft_signal[]);
        //! \brief Post-process in the forward transform.
        template<typename precision>
        static void post_forward(sycl::queue&, int length, std::complex<precision> const fft_result[], precision result[]);
        //! \brief Pre-process in the inverse transform.
        template<typename precision>
        static void pre_backward(sycl::queue&, int length, precision const input[], std::complex<precision> fft_signal[]);
        //! \brief Post-process in the inverse transform.
        template<typename precision>
        static void post_backward(sycl::queue&, int length, precision const fft_result[], precision result[]);
    };

}

namespace backend{
    /*!
     * \ingroup heffteoneapi
     * \brief Indicate that the oneMKL backend has been enabled.
     */
    template<> struct is_enabled<onemkl> : std::true_type{};
    /*!
     * \ingroup heffteoneapi
     * \brief Indicate that the oneMKL backend has been enabled.
     */
    template<> struct is_enabled<onemkl_cos> : std::true_type{};
    /*!
     * \ingroup heffteoneapi
     * \brief Indicate that the oneMKL backend has been enabled.
     */
    template<> struct is_enabled<onemkl_sin> : std::true_type{};

    /*!
     * \ingroup heffteoneapi
     * \brief Specialization that contains the sycl::queue needed for the DPC++ backend.
     */
    template<>
    struct device_instance<tag::gpu>{
        //! \brief Empty constructor.
        device_instance() : _stream(heffte::oapi::internal_sycl_queue){}
        //! \brief Constructor assigning the queue.
        device_instance(sycl::queue &new_stream) : _stream(new_stream){}
        //! \brief Constructor assigning from an existing wrapper.
        device_instance(std::reference_wrapper<sycl::queue> &new_stream) : _stream(new_stream){}
        //! \brief Returns the nullptr.
        sycl::queue& stream(){ return _stream; }
        //! \brief Returns the nullptr.
        sycl::queue& stream() const{ return _stream; }
        //! \brief Syncs the execution with the queue.
        void synchronize_device() const{ _stream.get().wait(); }
        //! \brief The sycl::queue, either user provided or created by heFFTe.
        std::reference_wrapper<sycl::queue> _stream;
        //! \brief The type for the internal stream.
        using stream_type = std::reference_wrapper<sycl::queue>;
    };
    /*!
     * \ingroup heffterocm
     * \brief In oneAPI mode, the default GPU backend is onemkl.
     */
    template<> struct default_backend<tag::gpu>{
        //! \brief Set the onemkl tag.
        using type = onemkl;
    };

    /*!
     * \ingroup heffterocm
     * \brief Specialization for the data operations in ROCm mode.
     */
    template<> struct data_manipulator<tag::gpu> {
        //! \brief The stream type for the device.
        using stream_type = sycl::queue&;
        //! \brief Defines the backend_device.
        using backend_device = backend::device_instance<tag::gpu>;
        //! \brief Allocate memory.
        template<typename scalar_type>
        static scalar_type* allocate(sycl::queue &stream, size_t num_entries){
            scalar_type* result = sycl::malloc_device<scalar_type>(num_entries, stream);
            stream.wait();
            return result;
        }
        //! \brief Free memory.
        template<typename scalar_type>
        static void free(sycl::queue &stream, scalar_type *pntr){
            if (pntr == nullptr) return;
            sycl::free(pntr, stream);
        }
        //! \brief Equivalent to std::copy_n() but using CUDA arrays.
        template<typename scalar_type>
        static void copy_n(sycl::queue &stream, scalar_type const source[], size_t num_entries, scalar_type destination[]){
            stream.memcpy(destination, source, num_entries * sizeof(scalar_type)).wait();
        }
        //! \brief Copy-convert complex-to-real.
        template<typename scalar_type>
        static void copy_n(sycl::queue &stream, std::complex<scalar_type> const source[], size_t num_entries, scalar_type destination[]){
            oapi::convert(stream, static_cast<long long>(num_entries), source, destination);
        }
        //! \brief Copy-convert real-to-complex.
        template<typename scalar_type>
        static void copy_n(sycl::queue &stream, scalar_type const source[], size_t num_entries, std::complex<scalar_type> destination[]){
            oapi::convert(stream, static_cast<long long>(num_entries), source, destination);
        }
        //! \brief Copy the date from the device to the host.
        template<typename scalar_type>
        static void copy_device_to_host(sycl::queue &stream, scalar_type const source[], size_t num_entries, scalar_type destination[]){
            stream.memcpy(destination, source, num_entries * sizeof(scalar_type)).wait();
        }
        //! \brief Copy the date from the device to the device.
        template<typename scalar_type>
        static void copy_device_to_device(sycl::queue &stream, scalar_type const source[], size_t num_entries, scalar_type destination[]){
            stream.memcpy(destination, source, num_entries * sizeof(scalar_type)).wait();
        }
        //! \brief Copy the date from the host to the device.
        template<typename scalar_type>
        static void copy_host_to_device(sycl::queue &stream, scalar_type const source[], size_t num_entries, scalar_type destination[]){
            stream.memcpy(destination, source, num_entries * sizeof(scalar_type)).wait();
        }
    };

    /*!
     * \ingroup heffteoneapi
     * \brief Defines the location type-tag and the oneAPI container.
     */
    template<>
    struct buffer_traits<onemkl>{
        //! \brief The oneMKL library uses data on the gpu device.
        using location = tag::gpu;
        //! \brief The data is managed by the oneAPI vector container.
        template<typename T> using container = heffte::gpu::device_vector<T, data_manipulator<tag::gpu>>;
    };
    /*!
     * \ingroup heffteoneapi
     * \brief Defines the location type-tag and the oneAPI container.
     */
    template<>
    struct buffer_traits<onemkl_cos>{
        //! \brief The oneMKL library uses data on the gpu device.
        using location = tag::gpu;
        //! \brief The data is managed by the oneAPI vector container.
        template<typename T> using container = heffte::gpu::device_vector<T, data_manipulator<tag::gpu>>;
    };
    /*!
     * \ingroup heffteoneapi
     * \brief Defines the location type-tag and the oneAPI container.
     */
    template<>
    struct buffer_traits<onemkl_sin>{
        //! \brief The oneMKL library uses data on the gpu device.
        using location = tag::gpu;
        //! \brief The data is managed by the oneAPI vector container.
        template<typename T> using container = heffte::gpu::device_vector<T, data_manipulator<tag::gpu>>;
    };
}

/*!
 * \ingroup heffteoneapi
 * \brief Wrapper around the oneMKL API.
 *
 * A single class that manages the plans and executions of oneMKL FFTs.
 * Handles the complex-to-complex cases.
 */
class onemkl_executor : public executor_base{
public:
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::forward;
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::backward;
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::complex_size;
    //! \brief Constructor, specifies the box and dimension.
    template<typename index>
    onemkl_executor(sycl::queue &inq, box3d<index> const box, int dimension) :
        q(inq),
        size(box.size[dimension]), size2(0),
        howmanyffts(fft1d_get_howmany(box, dimension)),
        stride(fft1d_get_stride(box, dimension)),
        dist((dimension == box.order[0]) ? size : 1),
        blocks((dimension == box.order[1]) ? box.osize(2) : 1),
        block_stride(box.osize(0) * box.osize(1)),
        total_size(box.count()),
        embed({0, static_cast<MKL_LONG>(stride), 0}),
        init_cplan(false), init_zplan(false),
        cplan(size), zplan(size)
    {}
    //! \brief Merges two FFTs into one.
    template<typename index>
    onemkl_executor(sycl::queue &inq, box3d<index> const box, int dir1, int dir2) :
        q(inq),
        size(box.size[std::min(dir1, dir2)]), size2(box.size[std::max(dir1, dir2)]),
        blocks(1), block_stride(0), total_size(box.count()),
        init_cplan(false), init_zplan(false),
        cplan({size, size2}), zplan({size, size2})
    {
        int odir1 = box.find_order(dir1);
        int odir2 = box.find_order(dir2);

        if (std::min(odir1, odir2) == 0 and std::max(odir1, odir2) == 1){
            stride = 1;
            dist = size * size2;
            embed = {0, static_cast<MKL_LONG>(stride), static_cast<MKL_LONG>(size)};
            howmanyffts = box.size[2];
        }else if (std::min(odir1, odir2) == 1 and std::max(odir1, odir2) == 2){
            stride = box.size[0];
            dist = 1;
            embed = {0, static_cast<MKL_LONG>(stride), static_cast<MKL_LONG>(size) * static_cast<MKL_LONG>(stride)};
            howmanyffts = box.size[0];
        }else{ // case of directions (0, 2)
            stride = 1;
            dist = size;
            embed = {0, static_cast<MKL_LONG>(stride), static_cast<MKL_LONG>(box.size[1]) * static_cast<MKL_LONG>(box.size[0])};
            howmanyffts = box.size[1];
        }
    }
    //! \brief Merges two FFTs into one.
    template<typename index>
    onemkl_executor(sycl::queue &inq, box3d<index> const box) :
        q(inq),
        size(box.size[0]), size2(box.size[1]), howmanyffts(box.size[2]),
        stride(0), dist(0),
        blocks(1), block_stride(0), total_size(box.count()),
        init_cplan(false), init_zplan(false),
        cplan({howmanyffts, size2, size}), zplan({howmanyffts, size2, size})
    {}

    //! \brief Forward fft, float-complex case.
    void forward(std::complex<float> data[], std::complex<float>*) const override{
        if (not init_cplan) make_plan(cplan);
        for(int i=0; i<blocks; i++)
            oneapi::mkl::dft::compute_forward(cplan, data + i * block_stride);
        q.wait();
    }
    //! \brief Backward fft, float-complex case.
    void backward(std::complex<float> data[], std::complex<float>*) const override{
        if (not init_cplan) make_plan(cplan);
        for(int i=0; i<blocks; i++)
            oneapi::mkl::dft::compute_backward(cplan, data + i * block_stride);
        q.wait();
    }
    //! \brief Forward fft, double-complex case.
    void forward(std::complex<double> data[], std::complex<double>*) const override{
        if (not init_zplan) make_plan(zplan);
        for(int i=0; i<blocks; i++)
            oneapi::mkl::dft::compute_forward(zplan, data + i * block_stride);
        q.wait();
    }
    //! \brief Backward fft, double-complex case.
    void backward(std::complex<double> data[], std::complex<double>*) const override{
        if (not init_zplan) make_plan(zplan);
        for(int i=0; i<blocks; i++)
            oneapi::mkl::dft::compute_backward(zplan, data + i * block_stride);
        q.wait();
    }

    //! \brief Converts the deal data to complex and performs float-complex forward transform.
    void forward(float const indata[], std::complex<float> outdata[], std::complex<float> *workspace) const override{
        for(int i=0; i<total_size; i++) outdata[i] = std::complex<float>(indata[i]);
        forward(outdata, workspace);
    }
    //! \brief Performs backward float-complex transform and truncates the complex part of the result.
    void backward(std::complex<float> indata[], float outdata[], std::complex<float> *workspace) const override{
        backward(indata, workspace);
        for(int i=0; i<total_size; i++) outdata[i] = std::real(indata[i]);
    }
    //! \brief Converts the deal data to complex and performs double-complex forward transform.
    void forward(double const indata[], std::complex<double> outdata[], std::complex<double> *workspace) const override{
        for(int i=0; i<total_size; i++) outdata[i] = std::complex<double>(indata[i]);
        forward(outdata, workspace);
    }
    //! \brief Performs backward double-complex transform and truncates the complex part of the result.
    void backward(std::complex<double> indata[], double outdata[], std::complex<double> *workspace) const override{
        backward(indata, workspace);
        for(int i=0; i<total_size; i++) outdata[i] = std::real(indata[i]);
    }

    //! \brief Returns the size of the box.
    int box_size() const override{ return total_size; }
    //! \brief Return the size of the needed workspace.
    size_t workspace_size() const override{ return 0; }

private:
    //! \brief Helper template to create the plan.
    template<typename onemkl_plan_type>
    void make_plan(onemkl_plan_type &plan) const{
        if (dist == 0){
            plan.set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, 1);
            plan.set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_INPLACE);
        }else if (size2 == 0){
            plan.set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, (MKL_LONG) howmanyffts);
            plan.set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_INPLACE);
            plan.set_value(oneapi::mkl::dft::config_param::INPUT_STRIDES, embed.data());
            plan.set_value(oneapi::mkl::dft::config_param::OUTPUT_STRIDES, embed.data());
            plan.set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, (MKL_LONG) dist);
            plan.set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, (MKL_LONG) dist);
        }else{
            plan.set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, (MKL_LONG) howmanyffts);
            plan.set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_INPLACE);
            plan.set_value(oneapi::mkl::dft::config_param::INPUT_STRIDES, embed.data());
            plan.set_value(oneapi::mkl::dft::config_param::OUTPUT_STRIDES, embed.data());
            plan.set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, (MKL_LONG) dist);
            plan.set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, (MKL_LONG) dist);
        }

        plan.commit(q);
        q.wait();

        if (std::is_same<oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::SINGLE, oneapi::mkl::dft::domain::COMPLEX>, onemkl_plan_type>::value)
            init_cplan = true;
        else
            init_zplan = true;
    }

    sycl::queue &q;
    int size, size2, howmanyffts, stride, dist, blocks, block_stride, total_size;
    std::array<MKL_LONG, 3> embed;

    mutable bool init_cplan, init_zplan;
    mutable oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::SINGLE, oneapi::mkl::dft::domain::COMPLEX> cplan;
    mutable oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::COMPLEX> zplan;
};

/*!
 * \ingroup heffteoneapi
 * \brief Wrapper to oneMKL API for real-to-complex transform with shortening of the data.
 *
 * Serves the same purpose of heffte::onemkl_executor but only real input is accepted
 * and only the unique (non-conjugate) coefficients are computed.
 * All real arrays must have size of real_size() and all complex arrays must have size complex_size().
 */
class onemkl_executor_r2c : public executor_base{
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
    onemkl_executor_r2c(sycl::queue &inq, box3d<index> const box, int dimension) :
        q(inq),
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
        init_splan(false), init_dplan(false),
        splan(size), dplan(size)
    {}

    //! \brief Forward transform, single precision.
    void forward(float const indata[], std::complex<float> outdata[], std::complex<float>*) const override{
        if (not init_splan) make_plan(splan);
        for(int i=0; i<blocks; i++)
            oneapi::mkl::dft::compute_forward(splan, const_cast<float*>(indata + i * rblock_stride), reinterpret_cast<float*>(outdata + i * cblock_stride));
        q.wait();
    }
    //! \brief Backward transform, single precision.
    void backward(std::complex<float> indata[], float outdata[], std::complex<float>*) const override{
        if (not init_splan) make_plan(splan);
        for(int i=0; i<blocks; i++)
            oneapi::mkl::dft::compute_backward(splan, reinterpret_cast<float*>(const_cast<std::complex<float>*>(indata + i * cblock_stride)), outdata + i * rblock_stride);
        q.wait();
    }
    //! \brief Forward transform, double precision.
    void forward(double const indata[], std::complex<double> outdata[], std::complex<double>*) const override{
        if (not init_dplan) make_plan(dplan);
        for(int i=0; i<blocks; i++)
            oneapi::mkl::dft::compute_forward(dplan, const_cast<double*>(indata + i * rblock_stride), reinterpret_cast<double*>(outdata + i * cblock_stride));
        q.wait();
    }
    //! \brief Backward transform, double precision.
    void backward(std::complex<double> indata[], double outdata[], std::complex<double>*) const override{
        if (not init_dplan) make_plan(dplan);
        for(int i=0; i<blocks; i++)
            oneapi::mkl::dft::compute_backward(dplan, reinterpret_cast<double*>(const_cast<std::complex<double>*>(indata + i * cblock_stride)), outdata + i * rblock_stride);
        q.wait();
    }

    //! \brief Returns the size of the box with real data.
    int box_size() const override{ return rsize; }
    //! \brief Returns the size of the box with complex coefficients.
    int complex_size() const override{ return csize; }
    //! \brief Return the size of the needed workspace.
    size_t workspace_size() const override{ return 0; }

private:
    //! \brief Helper template to initialize the plan.
    template<typename onemkl_plan_type>
    void make_plan(onemkl_plan_type &plan) const{
        plan.set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, (MKL_LONG) howmanyffts);
        plan.set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_NOT_INPLACE);
        plan.set_value(oneapi::mkl::dft::config_param::CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
        MKL_LONG slstride[] = {0, static_cast<MKL_LONG>(stride)};
        plan.set_value(oneapi::mkl::dft::config_param::INPUT_STRIDES, slstride);
        plan.set_value(oneapi::mkl::dft::config_param::OUTPUT_STRIDES, slstride);
        plan.set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, (MKL_LONG) rdist);
        plan.set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, (MKL_LONG) cdist);
        plan.commit(q);

        if (std::is_same<oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::SINGLE, oneapi::mkl::dft::domain::REAL>, onemkl_plan_type>::value)
            init_splan = true;
        else
            init_dplan = true;
        q.wait();
    }

    sycl::queue &q;

    int size, howmanyffts, stride, blocks;
    int rdist, cdist, rblock_stride, cblock_stride, rsize, csize;
    mutable bool init_splan, init_dplan;
    mutable oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::SINGLE, oneapi::mkl::dft::domain::REAL> splan;
    mutable oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL> dplan;
};

/*!
 * \ingroup heffteoneapi
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 *
 * The struct is specialized for each backend.
 * Note that the MKL and oneMKL backends use identical executors where oneMKL decides how to handle
 * the device data in the backend.
 */
template<> struct one_dim_backend<backend::onemkl>{
    //! \brief Defines the complex-to-complex executor.
    using executor = onemkl_executor;
    //! \brief Defines the real-to-complex executor.
    using executor_r2c = onemkl_executor_r2c;
};
/*!
 * \ingroup heffteoneapi
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 */
template<> struct one_dim_backend<backend::onemkl_cos>{
    //! \brief Defines the complex-to-complex executor.
    using executor = real2real_executor<backend::onemkl, oapi::cos_pre_pos_processor>;
    //! \brief Defines the real-to-complex executor.
    using executor_r2c = onemkl_executor_r2c;
};
/*!
 * \ingroup heffteoneapi
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 */
template<> struct one_dim_backend<backend::onemkl_sin>{
    //! \brief Defines the complex-to-complex executor.
    using executor = real2real_executor<backend::onemkl, oapi::sin_pre_pos_processor>;
    //! \brief Defines the real-to-complex executor.
    using executor_r2c = onemkl_executor_r2c;
};

/*!
 * \ingroup hefftepacking
 * \brief Simple packer that copies sub-boxes without transposing the order of the indexes.
 */
template<> struct direct_packer<tag::gpu>{
    //! \brief Execute the planned pack operation.
    template<typename scalar_type, typename index>
    void pack(sycl::queue &stream, pack_plan_3d<index> const &plan, scalar_type const data[], scalar_type buffer[]) const{
        oapi::direct_pack(stream, plan.size[0], plan.size[1], plan.size[2], plan.line_stride, plan.plane_stride, data, buffer);
    }
    //! \brief Execute the planned unpack operation.
    template<typename scalar_type, typename index>
    void unpack(sycl::queue &stream, pack_plan_3d<index> const &plan, scalar_type const buffer[], scalar_type data[]) const{
        oapi::direct_unpack(stream, plan.size[0], plan.size[1], plan.size[2], plan.line_stride, plan.plane_stride, buffer, data);
    }
};

/*!
 * \ingroup hefftepacking
 * \brief GPU version of the transpose packer.
 */
template<> struct transpose_packer<tag::gpu>{
    //! \brief Execute the planned pack operation.
    template<typename scalar_type, typename index>
    void pack(sycl::queue &stream, pack_plan_3d<index> const &plan, scalar_type const data[], scalar_type buffer[]) const{
        direct_packer<tag::gpu>().pack(stream, plan, data, buffer); // packing is done the same way as the direct_packer
    }
    //! \brief Execute the planned transpose-unpack operation.
    template<typename scalar_type, typename index>
    void unpack(sycl::queue &stream, pack_plan_3d<index> const &plan, scalar_type const buffer[], scalar_type data[]) const{
        oapi::transpose_unpack<scalar_type>(stream, plan.size[0], plan.size[1], plan.size[2], plan.line_stride, plan.plane_stride,
                                            plan.buff_line_stride, plan.buff_plane_stride, plan.map[0], plan.map[1], plan.map[2], buffer, data);
    }
};

/*!
 * \ingroup heffteoneapi
 * \brief Specialization for the CPU case.
 */
namespace data_scaling {
    /*!
     * \brief Simply multiply the \b num_entries in the \b data by the \b scale_factor.
     */
    template<typename scalar_type, typename index>
    static void apply(sycl::queue &stream, index num_entries, scalar_type *data, double scale_factor){
        oapi::scale_data(stream, static_cast<long long>(num_entries), data, scale_factor);
    }
    /*!
     * \brief Complex by real scaling.
     */
    template<typename precision_type, typename index>
    static void apply(sycl::queue &stream, index num_entries, std::complex<precision_type> *data, double scale_factor){
        apply<precision_type>(stream, 2*num_entries, reinterpret_cast<precision_type*>(data), scale_factor);
    }
};

/*!
 * \ingroup heffteoneapi
 * \brief Sets the default options for the oneMKL backend.
 */
template<> struct default_plan_options<backend::onemkl>{
    //! \brief The reshape operations will not transpose the data.
    static const bool use_reorder = false;
};
/*!
 * \ingroup heffteoneapi
 * \brief Sets the default options for the oneMKL backend.
 */
template<> struct default_plan_options<backend::onemkl_cos>{
    //! \brief The reshape operations will not transpose the data.
    static const bool use_reorder = true;
};
/*!
 * \ingroup heffteoneapi
 * \brief Sets the default options for the oneMKL backend.
 */
template<> struct default_plan_options<backend::onemkl_sin>{
    //! \brief The reshape operations will not transpose the data.
    static const bool use_reorder = true;
};

}

#endif

#endif
