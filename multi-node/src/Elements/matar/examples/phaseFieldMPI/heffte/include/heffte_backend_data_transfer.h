/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_BACKEND_DATA_TRANSFER_H
#define HEFFTE_BACKEND_DATA_TRANSFER_H

#ifdef Heffte_ENABLE_GPU

namespace heffte{

namespace gpu {

    /*!
     * \ingroup hefftegpu
     * \brief Device vector for the GPU backend.
     */
    template<typename scalar_type>
    using vector = device_vector<scalar_type, heffte::backend::data_manipulator<heffte::tag::gpu>>;

    /*!
     * \ingroup hefftegpu
     * \brief Collection of helper methods that transfer data using a memory manager.
     *
     * Helper templates associated with a specific memory_manager, see device_vector.
     */
    template<typename manipulator>
    struct device_transfer{
        //! \brief Define the backend device, same as the manipulator.
        using backend_device = typename manipulator::backend_device;
        //! \brief Copy data from the vector to the pointer, data size is equal to the vector size.
        template<typename scalar_type>
        static void copy(device_vector<scalar_type, backend_device> const &source, scalar_type destination[]){
            manipulator::copy_device_to_host(source.stream(), source.data(), source.size(), destination);
        }
        //! \brief Copy data from the pointer to the vector, data size is equal to the vector size.
        template<typename scalar_type>
        static void copy(scalar_type const source[], device_vector<scalar_type, backend_device> &destination){
            manipulator::copy_device_to_device(destination.stream(), source, destination.size(), destination.data());
        }

        /*!
         * \brief Copy the data from a buffer on the CPU to a cuda::vector.
         *
         * \tparam scalar_type of the vector entries.
         *
         * \param cpu_source is a buffer with size at least \b num_entries that sits in the CPU
         * \param num_entries is the number of entries to load
         *
         * \returns a device_vector with size equal to \b num_entries and a copy of the CPU data
         */
        template<typename scalar_type>
        static device_vector<scalar_type, manipulator> load(typename backend_device::stream_type stream, scalar_type const *cpu_source, size_t num_entries){
            device_vector<scalar_type, manipulator> result(stream, num_entries);
            manipulator::copy_host_to_device(stream, cpu_source, num_entries, result.data());
            return result;
        }
        //! \brief Using only arrays.
        template<typename scalar_type>
        static void load(typename backend_device::stream_type stream, scalar_type const *cpu_source, size_t num_entries, scalar_type *gpu_destination){
            manipulator::copy_host_to_device(stream, cpu_source, num_entries, gpu_destination);
        }
        //! \brief Using only arrays with both vectors on the CPU.
        template<typename scalar_type>
        static void load(void*, scalar_type const *cpu_source, size_t num_entries, scalar_type *gpu_destination){
            std::copy_n(cpu_source, num_entries, gpu_destination);
        }
        //! \brief Using the default stream.
        template<typename scalar_type>
        static device_vector<scalar_type, manipulator> load(scalar_type const *cpu_source, size_t num_entries){
            return load(backend_device().stream(), cpu_source, num_entries);
        }
        //! \brief Using the default stream.
        template<typename scalar_type>
        static device_vector<scalar_type, manipulator> load(void*, scalar_type const*, size_t){
            return device_vector<scalar_type, manipulator>();
        }

        //! \brief Similar to gpu::load() but loads the data from a std::vector
        template<typename scalar_type>
        static device_vector<scalar_type, manipulator> load(std::vector<scalar_type> const &cpu_source){
            return load(cpu_source.data(), cpu_source.size());
        }
        //! \brief Similar to gpu::load() but loads the data from a std::vector into a pointer.
        template<typename scalar_type>
        static void load(typename backend_device::stream_type stream, std::vector<scalar_type> const &cpu_source, scalar_type gpu_destination[]){
            manipulator::copy_host_to_device(stream, cpu_source.data(), cpu_source.size(), gpu_destination);
        }
        //! \brief Similar to gpu::load() but loads the data from a std::vector into a pointer.
        template<typename scalar_type>
        static void load(std::vector<scalar_type> const &cpu_source, scalar_type gpu_destination[]){
            load(backend_device().stream(), cpu_source, gpu_destination);
        }
        //! \brief Overload to handle the lack of constexpr-if.
        template<typename scalar_type>
        static void load(void*, std::vector<scalar_type> const&, scalar_type[]){}

        /*!
         * \brief Load method that copies two std::vectors, used in template general code.
         *
         * This is never executed.
         * Without if-constexpr (introduced in C++ 2017) generic template code must compile
         * even branches in the if-statements that will never be reached.
         */
        template<typename scalar_type>
        static void load(std::vector<scalar_type> const &a, std::vector<scalar_type> &b){ b = a; }
        /*!
         * \brief Unload method that copies two std::vectors, used in template general code.
         *
         * This is never executed.
         * Without if-constexpr (introduced in C++ 2017) generic template code must compile
         * even branches in the if-statements that will never be reached.
         */
        template<typename scalar_type>
        static std::vector<scalar_type> unload(std::vector<scalar_type> const &a){ return a; }
        //! \brief Dummy method for non-constexpr-if scenarios.
        template<typename scalar_type>
        static std::vector<scalar_type> unload(void*, scalar_type const[], size_t){
            return std::vector<scalar_type>();
        }

        //! \brief Copy number of entries from the GPU pointer into the vector.
        template<typename scalar_type>
        static std::vector<scalar_type> unload(typename backend_device::stream_type stream, scalar_type const gpu_source[], size_t num_entries){
            std::vector<scalar_type> result(num_entries);
            manipulator::copy_device_to_host(stream, gpu_source, num_entries, result.data());
            return result;
        }
        //! \brief Copy number of entries from the GPU pointer into the vector.
        template<typename scalar_type>
        static std::vector<scalar_type> unload(scalar_type const gpu_source[], size_t num_entries){
            return unload(backend_device().stream(), gpu_source, num_entries);
        }

        /*!
         * \brief Copy the data from a cuda::vector to a cpu buffer
         *
         * \tparam scalar_type of the vector entries
         *
         * \param gpu_source is the cuda::vector to holding the data to unload
         * \param cpu_result is a buffer with size at least \b gpu_data.size() that sits in the CPU
         */
        template<typename scalar_type>
        static void unload(device_vector<scalar_type, manipulator> const &gpu_source, scalar_type *cpu_result){
            manipulator::copy_device_to_host(gpu_source.device_stream(), gpu_source.data(), gpu_source.size(), cpu_result);
        }
        //! \brief Unload using a raw-array interface.
        template<typename scalar_type>
        static void unload(typename backend_device::stream_type stream, scalar_type const *gpu_source, size_t num_entries, scalar_type *cpu_result){
            manipulator::copy_device_to_host(stream, gpu_source, num_entries, cpu_result);
        }
        //! \brief Another silly template due to lack of constexpr-if
        template<typename scalar_type>
        static void unload(void*, scalar_type const *gpu_source, size_t num_entries, scalar_type *cpu_result){
            std::copy_n(gpu_source, num_entries, cpu_result); // void* stream indicates CPU backend
        }

        //! \brief Similar to unload() but copies the data into a std::vector.
        template<typename scalar_type>
        static std::vector<scalar_type> unload(device_vector<scalar_type, manipulator> const &gpu_source){
            std::vector<scalar_type> result(gpu_source.size());
            unload(gpu_source, result.data());
            return result;
        }
        /*!
         * \brief Captures ownership of the data in the raw-pointer.
         *
         * The advantage of the factory function over using the constructor is the ability
         * to auto-deduce the scalar type.
         */
        template<typename scalar_type>
        static device_vector<scalar_type, manipulator> capture(scalar_type* &&raw_pointer, size_t num_entries){
            return device_vector<scalar_type, manipulator>(std::forward<scalar_type*>(raw_pointer), num_entries);
        }
    };

    /*!
     * \ingroup hefftegpu
     * \brief Transfer helpers for the GPU backend.
     */
    using transfer = device_transfer<heffte::backend::data_manipulator<heffte::tag::gpu>>;

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

/*!
 * \ingroup hefftegpu
 * \brief Factory method to create new buffer container for the GPU backends.
 */
template<typename scalar_type>
gpu::vector<scalar_type> make_buffer_container(typename gpu::vector<scalar_type>::stream_type stream, size_t size){
    return gpu::vector<scalar_type>(stream, size);
}

}

#endif

#endif
