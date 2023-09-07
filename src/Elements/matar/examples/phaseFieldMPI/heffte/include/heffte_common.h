/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFFTE_COMMON_H
#define HEFFFTE_COMMON_H

#include "heffte_geometry.h"
#include "heffte_trace.h"

namespace heffte {

/*!
 * \ingroup fft3d
 * \addtogroup fft3dbackend Backend common wrappers
 *
 * Sub-module that encompasses all backend wrappers and meta data.
 */

/*!
 * \ingroup fft3dbackend
 * \brief Contains internal type-tags.
 *
 * Empty structs do not generate run-time code,
 * but can be used in type checks and overload resolutions at compile time.
 * Such empty classes are called "type-tags".
 */
namespace tag {

/*!
 * \ingroup fft3dbackend
 * \brief Indicates the use of cpu backend and that all input/output data and arrays will be bound to the cpu.
 *
 * Examples of cpu backends are FFTW and MKL.
 */
struct cpu{};
/*!
 * \ingroup fft3dbackend
 * \brief Indicates the use of gpu backend and that all input/output data and arrays will be bound to the gpu device.
 *
 * Example of gpu backend is cuFFT.
 */
struct gpu{};

}

/*!
 * \ingroup fft3dbackend
 * \brief Contains type tags and templates metadata for the various backends.
 */
namespace backend {

    /*!
    * \ingroup fft3dbackend
    * \brief Common data-transfer operations, must be specializes for each location (cpu/gpu).
    */
    template<typename location_tag> struct data_manipulator{};

    /*!
    * \ingroup fft3dbackend
    * \brief Common data-transfer operations on the cpu.
    */
    template<> struct data_manipulator<tag::cpu> {
        //! \brief The stream type for the device.
        using stream_type = void*;
        //! \brief Wrapper around std::copy_n().
        template<typename source_type, typename destination_type>
        static void copy_n(void*, source_type const source[], size_t num_entries, destination_type destination[]){
            std::copy_n(source, num_entries, destination);
        }
        //! \brief Wrapper around std::copy_n().
        template<typename source_type, typename destination_type>
        static void copy_n(source_type const source[], size_t num_entries, destination_type destination[]){
            std::copy_n(source, num_entries, destination);
        }
        //! \brief Wrapper around std::copy_n().
        template<typename source_type, typename destination_type>
        static void copy_device_to_host(void*, source_type const source[], size_t num_entries, destination_type destination[]){
            std::copy_n(source, num_entries, destination);
        }
        //! \brief Wrapper around std::copy_n().
        template<typename source_type, typename destination_type>
        static void copy_device_to_device(void*, source_type const source[], size_t num_entries, destination_type destination[]){
            std::copy_n(source, num_entries, destination);
        }
        //! \brief Wrapper around std::copy_n().
        template<typename source_type, typename destination_type>
        static void copy_host_to_device(void*, source_type const source[], size_t num_entries, destination_type destination[]){
            std::copy_n(source, num_entries, destination);
        }
    };

    /*!
     * \ingroup hefftefftw
     * \brief Type-tag for the FFTW backend
     */
    struct fftw{};
    /*!
     * \ingroup hefftefftw
     * \brief Type-tag for the Cosine Transform using the FFTW backend
     */
    struct fftw_cos{};
    /*!
     * \ingroup hefftefftw
     * \brief Type-tag for the Sine Transform using the FFTW backend
     */
    struct fftw_sin{};
    /*!
     * \ingroup hefftefftw
     * \brief Type-tag for the Cosine Transform type 1 using the FFTW backend
     */
    struct fftw_cos1{};
    /*!
     * \ingroup hefftefftw
     * \brief Type-tag for the Sine Transform type 1 using the FFTW backend
     */
    struct fftw_sin1{};

    /*!
     * \ingroup hefftestock
     * \brief Type-tag for the stock FFT backend
     */
    struct stock{};
    /*!
     * \ingroup hefftestock
     * \brief Type-tag for the Cosine Transform using the stock FFT backend
     */
    struct stock_cos{};
    /*!
     * \ingroup hefftestock
     * \brief Type-tag for the Sine Transform using the stock FFT backend
     */
    struct stock_sin{};

    /*!
     * \ingroup hefftemkl
     * \brief Type-tag for the MKL backend
     */
    struct mkl{};
    /*!
     * \ingroup hefftemkl
     * \brief Type-tag for the Cosine Transform using the MKL FFT backend
     */
    struct mkl_cos{};
    /*!
     * \ingroup hefftemkl
     * \brief Type-tag for the Sine Transform using the MKL FFT backend
     */
    struct mkl_sin{};

    /*!
     * \ingroup hefftecuda
     * \brief Type-tag for the cuFFT backend
     */
    struct cufft{};

    /*!
     * \ingroup hefftecuda
     * \brief Type-tag for the Cosine Transform using the cuFFT backend
     */
    struct cufft_cos{};
    /*!
     * \ingroup hefftecuda
     * \brief Type-tag for the Sine Transform using the cuFFT backend
     */
    struct cufft_sin{};

    /*!
     * \ingroup heffterocm
     * \brief Type-tag for the rocFFT backend
     */
    struct rocfft{};
    /*!
     * \ingroup heffterocm
     * \brief Type-tag for the Cosine Transform using the rocFFT backend
     */
    struct rocfft_cos{};
    /*!
     * \ingroup heffterocm
     * \brief Type-tag for the Sine Transform using the rocFFT backend
     */
    struct rocfft_sin{};

    /*!
     * \ingroup heffteoneapi
     * \brief Type-tag for the oneMKL backend
     */
    struct onemkl{};
    /*!
     * \ingroup heffteoneapi
     * \brief Type-tag for the Cosine Transform using the oneMKL backend
     */
    struct onemkl_cos{};
    /*!
     * \ingroup heffteoneapi
     * \brief Type-tag for the Sine Transform using the oneMKL backend
     */
    struct onemkl_sin{};

    /*!
     * \ingroup fft3dbackend
     * \brief Allows to define whether a specific backend interface has been enabled.
     *
     * Defaults to std::false_type, but specializations for each enabled backend
     * will overwrite this to the std::true_type, i.e., define const static bool value
     * which is set to true.
     */
    template<typename tag>
    struct is_enabled : std::false_type{};

    /*!
     * \ingroup fft3dbackend
     * \brief Defines the container for the temporary buffers.
     *
     * Specialization for each backend will define whether the raw-arrays are associated
     * with the CPU or GPU devices and the type of the container that will hold temporary
     * buffers.
     */
    template<typename backend_tag, typename std::enable_if<is_enabled<backend_tag>::value, void*>::type = nullptr>
    struct buffer_traits{
        //! \brief Tags the raw-array location tag::cpu or tag::gpu, used by the packers.
        using location = tag::cpu;
        //! \brief Defines the container template to use for the temporary buffers in heffte::fft3d.
        template<typename T> using container = std::vector<T>;
    };

    /*!
     * \ingroup fft3dbackend
     * \brief Struct that specializes to true type if the location of the backend is on the gpu (false type otherwise).
     */
    template<typename backend_tag, typename = void>
    struct uses_gpu : std::false_type{};

    /*!
     * \ingroup fft3dbackend
     * \brief Specialization for the on-gpu case.
     */
    template<typename backend_tag>
    struct uses_gpu<backend_tag,
                    typename std::enable_if<std::is_same<typename buffer_traits<backend_tag>::location, tag::gpu>::value, void>::type>
    : std::true_type{};

    /*!
     * \ingroup fft3dbackend
     * \brief Returns the human readable name of the backend.
     */
    template<typename backend_tag>
    inline std::string name(){ return "unknown"; }

    /*!
     * \ingroup hefftefftw
     * \brief Returns the human readable name of the FFTW backend.
     */
    template<> inline std::string name<fftw>(){ return "fftw"; }
    /*!
     * \ingroup hefftefftw
     * \brief Returns the human readable name of the FFTW backend.
     */
    template<> inline std::string name<fftw_cos>(){ return "fftw-cos-type-II"; }
    /*!
     * \ingroup hefftefftw
     * \brief Returns the human readable name of the FFTW backend.
     */
    template<> inline std::string name<fftw_sin>(){ return "fftw-sin-type-II"; }
    /*!
     * \ingroup hefftefftw
     * \brief Returns the human readable name of the FFTW backend.
     */
    template<> inline std::string name<fftw_cos1>(){ return "fftw-cos-type-I"; }
    /*!
     * \ingroup hefftefftw
     * \brief Returns the human readable name of the FFTW backend.
     */
    template<> inline std::string name<fftw_sin1>(){ return "fftw-sin-type-I"; }

    /*!
     * \ingroup hefftestock
     * \brief Returns the human readable name of the stock backend.
     */
    template<> inline std::string name<stock>(){ return "stock"; }
    /*!
     * \ingroup hefftestock
     * \brief Returns the human readable name of the stock backend.
     */
    template<> inline std::string name<stock_cos>(){ return "stock-cos"; }
    /*!
     * \ingroup hefftestock
     * \brief Returns the human readable name of the stock backend.
     */
    template<> inline std::string name<stock_sin>(){ return "stock-sin"; }

    /*!
     * \ingroup hefftemkl
     * \brief Returns the human readable name of the MKL backend.
     */
    template<> inline std::string name<mkl>(){ return "mkl"; }
    /*!
     * \ingroup hefftemkl
     * \brief Returns the human readable name of the MKL backend.
     */
    template<> inline std::string name<mkl_cos>(){ return "mkl-cos"; }
    /*!
     * \ingroup hefftemkl
     * \brief Returns the human readable name of the MKL backend.
     */
    template<> inline std::string name<mkl_sin>(){ return "mkl-sin"; }
    /*!
     * \ingroup hefftecuda
     * \brief Returns the human readable name of the cuFFT backend.
     */
    template<> inline std::string name<cufft>(){ return "cufft"; }
    /*!
     * \ingroup hefftecuda
     * \brief Returns the human readable name of the cuFFT backend.
     */
    template<> inline std::string name<cufft_cos>(){ return "cufft-cos"; }
    /*!
     * \ingroup hefftecuda
     * \brief Returns the human readable name of the cuFFT backend.
     */
    template<> inline std::string name<cufft_sin>(){ return "cufft-sin"; }

    /*!
     * \ingroup heffterocm
     * \brief Returns the human readable name of the rocFFT backend.
     */
    template<> inline std::string name<rocfft>(){ return "rocfft"; }
    /*!
     * \ingroup heffterocm
     * \brief Returns the human readable name of the rocFFT backend.
     */
    template<> inline std::string name<rocfft_cos>(){ return "rocfft-cos"; }
    /*!
     * \ingroup heffterocm
     * \brief Returns the human readable name of the rocFFT backend.
     */
    template<> inline std::string name<rocfft_sin>(){ return "rocfft-sin"; }

    /*!
     * \ingroup heffteoneapi
     * \brief Returns the human readable name of the oneMKL backend.
     */
    template<> inline std::string name<onemkl>(){ return "onemkl"; }
    /*!
     * \ingroup heffteoneapi
     * \brief Returns the human readable name of the oneMKL backend.
     */
    template<> inline std::string name<onemkl_cos>(){ return "onemkl-cos"; }
    /*!
     * \ingroup heffteoneapi
     * \brief Returns the human readable name of the oneMKL backend.
     */
    template<> inline std::string name<onemkl_sin>(){ return "onemkl-sin"; }

    /*!
     * \ingroup fft3dbackend
     * \brief Indicates the name of the location tag.
     */
    template<> inline std::string name<tag::cpu>(){ return "cpu"; }
    /*!
     * \ingroup fft3dbackend
     * \brief Indicates the name of the location tag.
     */
    template<> inline std::string name<tag::gpu>(){ return "gpu"; }

    /*!
     * \ingroup fft3dbackend
     * \brief Holds the auxiliary variables needed by each backend.
     *
     * The idea is similar to <a href="https://en.cppreference.com/w/cpp/language/crtp">CRTP</a>
     * heffte::fft3d and heffte::fft3d_r2c inherit from this class and specializations based
     * on the backend-tag can define a different set of internal variables.
     * Specifically, this is used to store the sycl::queue used by the DPC++ backend.
     */
    template<typename backend_tag>
    struct device_instance{
        //! \brief Empty constructor.
        device_instance(void* = nullptr){}
        //! \brief Default destructor.
        virtual ~device_instance() = default;
        //! \brief Returns the nullptr.
        void* stream(){ return nullptr; }
        //! \brief Returns the nullptr (const case).
        void* stream() const{ return nullptr; }
        //! \brief Syncs the execution with the queue, no-op in the CPU case.
        void synchronize_device() const{}
        //! \brief The type for the internal stream, the cpu uses just a void pointer.
        using stream_type = void*;
    };

    /*!
     * \ingroup fft3dbackend
     * \brief Defines inverse mapping from the location tag to a default backend tag.
     *
     * Defines a default backend for a given location tag, the backend assumes FFT transform,
     * not a sine or cosine variant.
     */
    template<typename location_tag> struct default_backend{
        //! \brief Defaults to the stock backend.
        using type = stock;
    };

    /*!
     * \ingroup fft3dbackend
     * \brief Defines whether the backend accepts the standard FFT real-complex or complex-complex transform.
     */
    template<typename backend_tag> struct uses_fft_types : std::true_type{};

    /*!
     * \ingroup fft3dbackend
     * \brief Set to true/false type depending whether the types are compatible with the backend transform.
     */
    template<typename backend_tag, typename input, typename output, typename = void> struct check_types : std::false_type{};

    /*!
     * \ingroup fft3dbackend
     * \brief Defines the types compatible for a standard FFT transform.
     */
    template<typename backend_tag, typename input, typename output> struct check_types<backend_tag, input, output,
        typename std::enable_if<uses_fft_types<backend_tag>::value and (
                      (std::is_same<input, float>::value and is_ccomplex<output>::value)
                   or (std::is_same<input, double>::value and is_zcomplex<output>::value)
                   or (is_ccomplex<input>::value and is_ccomplex<output>::value)
                   or (is_zcomplex<input>::value and is_zcomplex<output>::value)
                  )>::type> : std::true_type{};

    /*!
     * \ingroup fft3dbackend
     * \brief Sets the cos() transform C++ types.
     */
    template<> struct uses_fft_types<fftw_cos> : std::false_type{};
    /*!
     * \ingroup fft3dbackend
     * \brief Sets the sin() transform C++ types.
     */
    template<> struct uses_fft_types<fftw_sin> : std::false_type{};
    /*!
     * \ingroup fft3dbackend
     * \brief Sets the cos() transform C++ types.
     */
    template<> struct uses_fft_types<fftw_cos1> : std::false_type{};
    /*!
     * \ingroup fft3dbackend
     * \brief Sets the sin() transform C++ types.
     */
    template<> struct uses_fft_types<fftw_sin1> : std::false_type{};
    /*!
     * \ingroup hefftestock
     * \brief Sets the cos() transform types.
     */
    template<> struct uses_fft_types<stock_cos> : std::false_type{};
    /*!
     * \ingroup hefftestock
     * \brief Sets the sin() transform types.
     */
    template<> struct uses_fft_types<stock_sin> : std::false_type{};
    /*!
     * \ingroup hefftemkl
     * \brief Sets the cos() transform types.
     */
    template<> struct uses_fft_types<mkl_cos> : std::false_type{};
    /*!
     * \ingroup hefftemkl
     * \brief Sets the sin() transform types.
     */
    template<> struct uses_fft_types<mkl_sin> : std::false_type{};
    /*!
     * \ingroup hefftemkl
     * \brief Sets the cos() transform types.
     */
    template<> struct uses_fft_types<cufft_cos> : std::false_type{};
    /*!
     * \ingroup hefftemkl
     * \brief Sets the sin() transform types.
     */
    template<> struct uses_fft_types<cufft_sin> : std::false_type{};
    /*!
     * \ingroup hefftemkl
     * \brief Sets the cos() transform types.
     */
    template<> struct uses_fft_types<rocfft_cos> : std::false_type{};
    /*!
     * \ingroup hefftemkl
     * \brief Sets the sin() transform types.
     */
    template<> struct uses_fft_types<rocfft_sin> : std::false_type{};
    /*!
     * \ingroup hefftemkl
     * \brief Sets the cos() transform types.
     */
    template<> struct uses_fft_types<onemkl_cos> : std::false_type{};
    /*!
     * \ingroup hefftemkl
     * \brief Sets the sin() transform types.
     */
    template<> struct uses_fft_types<onemkl_sin> : std::false_type{};

    /*!
     * \ingroup fft3dbackend
     * \brief Defines the types compatible for a cos() transform.
     */
    template<typename backend_tag, typename input, typename output> struct check_types<backend_tag, input, output,
        typename std::enable_if<not uses_fft_types<backend_tag>::value and (
                      (std::is_same<input, float>::value and std::is_same<output, float>::value)
                   or (std::is_same<input, double>::value and std::is_same<output, double>::value)
                  )>::type> : std::true_type{};

}

/*!
 * \ingroup fft3dbackend
 * \brief Base class for all backend executors.
 */
class executor_base{
public:
    //! \brief Virtual destructor.
    virtual ~executor_base() = default;
    //! \brief Forward r2r, single precision.
    virtual void forward(float[], float*) const{}
    //! \brief Forward r2r, double precision.
    virtual void forward(double[], double*) const{}
    //! \brief Backward r2r, single precision.
    virtual void backward(float[], float*) const{}
    //! \brief Backward r2r, double precision.
    virtual void backward(double[], double*) const{}
    //! \brief Forward FFT, single precision.
    virtual void forward(std::complex<float>[], std::complex<float>*) const{}
    //! \brief Forward FFT, double precision.
    virtual void forward(std::complex<double>[], std::complex<double>*) const{}
    //! \brief Backward FFT, single precision.
    virtual void backward(std::complex<float>[], std::complex<float>*) const{}
    //! \brief Backward FFT, double precision.
    virtual void backward(std::complex<double>[], std::complex<double>*) const{}
    //! \brief Forward FFT real-to-complex, single precision.
    virtual void forward(float const[], std::complex<float>[], std::complex<float>*) const{}
    //! \brief Forward FFT real-to-complex, double precision.
    virtual void forward(double const[], std::complex<double>[], std::complex<double>*) const{}
    //! \brief Backward FFT real-to-complex, single precision.
    virtual void backward(std::complex<float>[], float[], std::complex<float>*) const{}
    //! \brief Backward FFT real-to-complex, double precision.
    virtual void backward(std::complex<double>[], double[], std::complex<double>*) const{}
    //! \brief Return the size of the box.
    virtual int box_size() const{ return 0; }
    //! \brief Return the workspace of the size.
    virtual size_t workspace_size() const{ return 0; }
    //! \brief Return the size of the complex-box (r2c executors).
    virtual int complex_size() const{ return box_size(); }
};

/*!
 * \ingroup fft3dbackend
 * \brief Factory method to create new buffer container for the CPU backends.
 */
template<typename scalar_type>
std::vector<scalar_type> make_buffer_container(void*, size_t size){
    return std::vector<scalar_type>(size);
}

/*!
 * \ingroup fft3dmisc
 * \brief Indicates the direction of the FFT (internal use only).
 */
enum class direction {
    //! \brief Forward DFT transform.
    forward,
    //! \brief Inverse DFT transform.
    backward
};

/*!
 * \ingroup fft3dbackend
 * \brief Indicates the structure that will be used by the fft backend.
 */
template<typename> struct one_dim_backend{};

/*!
 * \ingroup fft3dbackend
 * \brief Factory method to construct an executor for the FFT backend.
 */
template<typename backend_tag, typename index>
static std::unique_ptr<typename one_dim_backend<backend_tag>::executor> make_executor(typename backend::device_instance<typename backend::buffer_traits<backend_tag>::location>::stream_type stream,
                                                                                  box3d<index> const box, int dimension){
    return (box.empty()) ?
        std::unique_ptr<typename one_dim_backend<backend_tag>::executor>() :
        std::unique_ptr<typename one_dim_backend<backend_tag>::executor>(new typename one_dim_backend<backend_tag>::executor(stream, box, dimension));
}
/*!
 * \ingroup fft3dbackend
 * \brief Factory method to construct an executor for the FFT backend, 2D variant.
 */
template<typename backend_tag, typename index>
static std::unique_ptr<typename one_dim_backend<backend_tag>::executor> make_executor(typename backend::device_instance<typename backend::buffer_traits<backend_tag>::location>::stream_type stream,
                                                                                  box3d<index> const box, int dir1, int dir2){
    return (box.empty()) ?
        std::unique_ptr<typename one_dim_backend<backend_tag>::executor>() :
        std::unique_ptr<typename one_dim_backend<backend_tag>::executor>(new typename one_dim_backend<backend_tag>::executor(stream, box, dir1, dir2));
}
/*!
 * \ingroup fft3dbackend
 * \brief Factory method to construct an executor for the FFT backend, 3D variant.
 */
template<typename backend_tag, typename index>
static std::unique_ptr<typename one_dim_backend<backend_tag>::executor> make_executor(typename backend::device_instance<typename backend::buffer_traits<backend_tag>::location>::stream_type stream,
                                                                                  box3d<index> const box){
    return (box.empty()) ?
        std::unique_ptr<typename one_dim_backend<backend_tag>::executor>() :
        std::unique_ptr<typename one_dim_backend<backend_tag>::executor>(new typename one_dim_backend<backend_tag>::executor(stream, box));
}
/*!
 * \ingroup fft3dbackend
 * \brief Factory method to construct an executor for the FFT backend, r2c variant.
 */
template<typename backend_tag, typename index>
static std::unique_ptr<typename one_dim_backend<backend_tag>::executor_r2c> make_executor_r2c(typename backend::device_instance<typename backend::buffer_traits<backend_tag>::location>::stream_type stream,
                                                                                          box3d<index> const box, int dimension){
    return (box.empty()) ?
        std::unique_ptr<typename one_dim_backend<backend_tag>::executor_r2c>() :
        std::unique_ptr<typename one_dim_backend<backend_tag>::executor_r2c>(new typename one_dim_backend<backend_tag>::executor_r2c(stream, box, dimension));
}

/*!
 * \ingroup fft3dbackend
 * \brief Defines whether the executor has a 2D version (slabs).
 */
template<typename backend_tag>
constexpr bool has_executor2d(){
    // cosine transform variants don't have a 2D/3D version yet (due to the missing kernels)
    // most backends are OK with the variants for 2D and 3D (stock isn't)
    return not (std::is_same<backend_tag, backend::stock>::value
            or std::is_same<backend_tag, backend::stock_cos>::value
            or std::is_same<backend_tag, backend::mkl_cos>::value
            or std::is_same<backend_tag, backend::cufft_cos>::value
            or std::is_same<backend_tag, backend::rocfft_cos>::value
            or std::is_same<backend_tag, backend::onemkl_cos>::value
            or std::is_same<backend_tag, backend::stock_sin>::value
            or std::is_same<backend_tag, backend::mkl_sin>::value
            or std::is_same<backend_tag, backend::cufft_sin>::value
            or std::is_same<backend_tag, backend::rocfft_sin>::value
            or std::is_same<backend_tag, backend::onemkl_sin>::value
            );
}
/*!
 * \ingroup fft3dbackend
 * \brief Defines whether the executor has a 3D version (single rank).
 */
template<typename backend_tag>
constexpr bool has_executor3d(){
    return not (std::is_same<backend_tag, backend::stock>::value
            or std::is_same<backend_tag, backend::stock_cos>::value
            or std::is_same<backend_tag, backend::mkl_cos>::value
            or std::is_same<backend_tag, backend::cufft_cos>::value
            or std::is_same<backend_tag, backend::rocfft_cos>::value
            or std::is_same<backend_tag, backend::onemkl_cos>::value
            or std::is_same<backend_tag, backend::stock_sin>::value
            or std::is_same<backend_tag, backend::mkl_sin>::value
            or std::is_same<backend_tag, backend::cufft_sin>::value
            or std::is_same<backend_tag, backend::rocfft_sin>::value
            or std::is_same<backend_tag, backend::onemkl_sin>::value
            );
}

/*!
 * \ingroup fft3dbackend
 * \brief Defines a set of default plan options for a given backend.
 */
template<typename> struct default_plan_options{};

}

#endif   //  #ifndef HEFFTE_COMMON_H
