/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_BACKEND_VECTOR_H
#define HEFFTE_BACKEND_VECTOR_H

#ifdef Heffte_ENABLE_GPU

/*!
 * \ingroup fft3d
 * \addtogroup hefftegpu Common templates for the GPU backends
 *
 * Templates used by the different GPU backends, mostly related
 * to memory management.
 */

namespace heffte{

/*!
 * \ingroup hefftegpu
 * \brief GPU specific methods.
 */
namespace gpu {

    /*!
     * \ingroup hefftegpu
     * \brief Container that wraps around a raw device array.
     *
     * Wrapper around a device array that allows for data to be automatically
     * allocated and freed using RAII style of resource management.
     * The \b scalar_type defines the data-type and the \b memory_manager
     * is a class/struct that implements several static members, e.g.,
     * \code
     * struct cuda_memory_manager{
     *     void* allocate(size_t num_bytes);
     *     void free(void *pnrs);
     *     void host_to_device(void const *source, size_t num_bytes, void *destination);
     *     void device_to_device(void const *source, size_t num_bytes, void *destination);
     *     void device_to_host(void const *source, size_t num_bytes, void *destination);
     * };
     * \endcode
     * All methods accept number of bytes, so sized passed have to be multiplied by
     * sizeof(scalar_type). The free method can accept a nullptr. The device/host
     * methods transfer data between the host, device, or device to device.
     *
     * The device_vector class is movable and copiable (using deep-copy)
     * in both the constructor and assignment operator.
     */
    template<typename scalar_type, typename manipulator>
    class device_vector{
    public:
        //! \brief The value of the array, used for static error checking.
        using value_type = scalar_type;
        //! \brief Define the backend device, same as the manipulator.
        using backend_device = typename manipulator::backend_device;
        //! \brief The type for the stream used by the backend device.
        using stream_type = typename backend_device::stream_type;

        //! \brief Allocate a new vector with the given number of entries.
        device_vector(size_t num_entries = 0) :
            device(),
            num(num_entries),
            device_data(manipulator::template allocate<scalar_type>(device.stream(), num_entries))
        {}
        //! \brief Allocate a new vector with the given number of entries.
        device_vector(backend_device const &new_device, size_t num_entries = 0) :
            device(new_device),
            num(num_entries),
            device_data(manipulator::template allocate<scalar_type>(device.stream(), num_entries))
        {}
        //! \brief Copy a range of entries from the device into the vector.
        device_vector(scalar_type const *begin, scalar_type const *end) :
            device_vector(std::distance(begin, end)){
            manipulator::copy_device_to_device(device.stream(), begin, num, device_data);
        }
        //! \brief Copy a range of entries from the device into the vector.
        device_vector(backend_device const &new_device, scalar_type const *begin, scalar_type const *end) :
            device_vector(new_device, std::distance(begin, end)){
            manipulator::copy_device_to_device(device.stream(), begin, num, device_data);
        }

        //! \brief Copy constructor, copy the data from other to this vector.
        device_vector(const device_vector<scalar_type, manipulator>& other) :
            device_vector(other.device, other.num){
            manipulator::copy_device_to_device(device.stream(), other.device_data, num, device_data);
        }
        //! \brief Move constructor, moves the data from \b other into this vector.
        device_vector(device_vector<scalar_type, manipulator> &&other) :
            device(other.device),
            num(c11_exchange(other.num, 0)),
            device_data(c11_exchange(other.device_data, nullptr))
        {}

        //! \brief Captures ownership of the data in the raw-pointer, resets the pointer to null.
        device_vector(scalar_type* &&raw_pointer, size_t num_entries) :
            device(),
            num(num_entries),
            device_data(c11_exchange(raw_pointer, nullptr))
        {}
        //! \brief Captures ownership of the data in the raw-pointer, resets the pointer to null.
        device_vector(backend_device const &new_device, scalar_type* &&raw_pointer, size_t num_entries) :
            device(new_device),
            num(num_entries),
            device_data(c11_exchange(raw_pointer, nullptr))
        {}

        //! \brief Desructor, deletes all data.
        ~device_vector(){ manipulator::free(device.stream(), device_data); }

        //! \brief Copy assignment, copies the data form \b other to this object.
        void operator =(device_vector<scalar_type, manipulator> const &other){
            device_vector<scalar_type, manipulator> temp(other);
            device = temp.device;
            std::swap(num, temp.num);
            std::swap(device_data, temp.device_data);
        }

        //! \brief Move assignment, moves the data form \b other to this object.
        void operator =(device_vector<scalar_type, manipulator>&& other){
            device_vector<scalar_type, manipulator> temp(std::move(other));
            device = temp.device;
            std::swap(num, temp.num);
            std::swap(device_data, temp.device_data);
        }

        //! \brief Give reference to the array, can be passed directly into cuFFT calls or custom kernels.
        scalar_type* data(){ return device_data; }
        //! \brief Give const reference to the array, can be passed directly into cuFFT calls or custom kernels.
        const scalar_type* data() const{ return device_data; }

        //! \brief Return the current size of the array, i.e., the number of elements.
        size_t size() const{ return num; }
        //! \brief Return \b true if the vector has zero size.
        bool empty() const{ return (num == 0); }

        //! \brief Returns the current array and releases ownership.
        scalar_type* release(){
            num = 0;
            return c11_exchange(device_data, nullptr);
        }

        //! \brief Return a reference to the internal device stream.
        stream_type device_stream(){ return device.stream(); }
        //! \brief Return a const reference to the internal device stream.
        stream_type device_stream() const{ return device.stream(); }

    private:
        //! \brief Save an instance of the device.
        backend_device device;
        //! \brief Stores the number of entries in the vector.
        size_t num;
        //! \brief The array with the GPU data.
        scalar_type *device_data;
    };

    /*!
     * \ingroup hefftegpu
     * \brief Wrapper around cudaGetDeviceCount()
     */
    int device_count();

    /*!
     * \ingroup hefftegpu
     * \brief Wrapper around cudaSetDevice()
     *
     * \param active_device is the new active CUDA device for this thread, see the Nvidia documentation for cudaSetDevice()
     */
    void device_set(int active_device);

    /*!
     * \ingroup hefftegpu
     * \brief Wrapper around cudaStreamSynchronize(nullptr).
     */
    void synchronize_default_stream();

}

}

#endif

#endif   /* HEFFTE_BACKEND_VECTOR_H */
