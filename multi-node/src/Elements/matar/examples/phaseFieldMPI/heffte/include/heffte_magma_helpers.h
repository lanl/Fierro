/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_MAGMA_HELPERS_H
#define HEFFTE_MAGMA_HELPERS_H

#include "heffte_pack3d.h"

/*!
 * \ingroup fft3d
 * \addtogroup fft3magma Helper Methods from MAGMA
 *
 * The UTK MAGMA library provides several methods for GPU accelerated linear algebra
 * which can be incorporated in heFFTe.
 */

namespace heffte{

namespace gpu {
    /*!
     * \ingroup fft3magma
     * \brief Wrapper around a MAGMA handle.
     *
     * The generic template performs no actions and is in fact an empty struct.
     * The MAGMA capabilities are implemented in the tag::gpu specialization.
     */
    template<typename> struct magma_handle{
        //! \brief No-op constructor, see the GPU specialization.
        magma_handle(void*){}
        //! \brief Wrapper around MAGMA sscal()/dscal(), but no-op in a CPU context.
        template<typename scalar_type> void scal(int, double, scalar_type*) const{}
    };
    /*!
     * \ingroup fft3magma
     * \brief Wrapper around a MAGMA handle in a GPU context.
     *
     * Class that wraps around MAGMA capabilities.
     */
    template<>
    struct magma_handle<tag::gpu>{
        //! \brief Constructor, calls magma_init() and creates a new queue on the given stream.
        magma_handle(void* gpu_stream);
        //! \brief Destructor, cleans the queue and calls magma_finalize().
        ~magma_handle();
        //! \brief Wrapper around magma_sscal()
        void scal(int num_entries, double scale_factor, float *data) const;
        //! \brief Wrapper around magma_dscal()
        void scal(int num_entries, double scale_factor, double *data) const;
        //! \brief Template for the complex case, uses reinterpret_cast().
        template<typename precision_type>
        void scal(int num_entries, double scale_factor, std::complex<precision_type> *data) const{
            scal(2*num_entries, scale_factor, reinterpret_cast<precision_type*>(data));
        }
        //! \brief Opaque pointer to a magma_queue.
        mutable void *handle;
    };
}

}

#endif   /* HEFFTE_MAGMA_HELPERS_H */
