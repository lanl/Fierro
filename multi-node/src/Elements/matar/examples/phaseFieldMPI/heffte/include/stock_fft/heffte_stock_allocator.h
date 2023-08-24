#ifndef HEFFTE_STOCK_ALLOCATOR_H
#define HEFFTE_STOCK_ALLOCATOR_H

#include <memory>
#include <vector>

#include "heffte_stock_complex.h"

namespace heffte {
namespace stock{

/*! \ingroup stockbackend
 *  \brief Allocator to use with heffte::stock::Complex types
 *  Class to properly allocate heffte::stock::Complex<F,L> types to ensure
 *  proper alignment of the type when using containers like std::vector
 */
template<typename F>
class complex_allocator_t {
    public:
        //! \brief Mandatory aliases
        typedef F value_type;
        typedef F* pointer;
        typedef const F* const_pointer;
        typedef F& reference;
        typedef const F& const_reference;
        typedef size_t size_type;
        typedef std::ptrdiff_t difference_type;

        //! \brief Defining rebind for the allocator
        template <class U>
        struct rebind {
            typedef complex_allocator_t<U> other;
        };

        //! \brief Get address from a reference
        inline pointer address(reference r) { return &r; }
        //! \brief Get address from a const reference
        inline const_pointer address(const_reference r) const { return &r; }

        //! \brief Define allocation for complex type
        pointer allocate(size_type n, typename std::allocator<void>::const_pointer = nullptr) {
            #ifdef Heffte_ENABLE_AVX
            return reinterpret_cast<pointer>(aligned_alloc(alignof(F), n*sizeof(F)));
            #else
            return reinterpret_cast<pointer>(malloc(n*sizeof(F)));
            #endif
        }
        //! \brief Define deallocation for complex type
        inline void deallocate(pointer p, size_type) {
            free(p);
        };

        //! \brief Copy into pointer
        inline void construct(pointer p, const_reference value) { new (p) value_type(value); }
        //! \brief Destroy a pointer to this type
        inline void destroy(pointer p) { p->~value_type(); }

        //! \brief Define maximum size of an array of this
        inline size_type max_size() const noexcept { return size_type(-1) / sizeof(F); }

        //! \brief Define == operator
        inline bool operator==(const complex_allocator_t&) { return true; }
        //! \brief Define != operator
        inline bool operator!=(const complex_allocator_t& rhs) { return !operator==(rhs); }
};

//! \brief Alias for this complex allocator when used with the Complex type
template<typename F, int L>
using complex_allocator = complex_allocator_t<Complex<F, L>>;

//! \brief Alias for how to use this allocator with a heffte::stock::Complex type
template<typename F, int L>
using complex_vector = std::vector<Complex<F,L>, complex_allocator<F,L>>;

}
}
#endif // End HEFFTE_STOCK_ALLOCATOR.H
