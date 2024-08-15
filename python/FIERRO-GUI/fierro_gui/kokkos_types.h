#ifndef KOKKOS_TYPES_H
#define KOKKOS_TYPES_H
/**********************************************************************************************
 © 2020. Triad National Security, LLC. All rights reserved.
 This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
 National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
 Department of Energy/National Nuclear Security Administration. All rights in the program are
 reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
 Security Administration. The Government is granted for itself and others acting on its behalf a
 nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
 derivative works, distribute copies to the public, perform publicly and display publicly, and
 to permit others to do so.
 This program is open source under the BSD-3 License.
 Redistribution and use in source and binary forms, with or without modification, are permitted
 provided that the following conditions are met:
 
 1.  Redistributions of source code must retain the above copyright notice, this list of
 conditions and the following disclaimer.
 
 2.  Redistributions in binary form must reproduce the above copyright notice, this list of
 conditions and the following disclaimer in the documentation and/or other materials
 provided with the distribution.
 
 3.  Neither the name of the copyright holder nor the names of its contributors may be used
 to endorse or promote products derived from this software without specific prior
 written permission.
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **********************************************************************************************/

#include "host_types.h"


#ifdef HAVE_KOKKOS
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

using HostSpace    = Kokkos::HostSpace;
using MemoryUnmanaged = Kokkos::MemoryUnmanaged;

#ifdef HAVE_CUDA
//using UVMMemSpace     = Kokkos::CudaUVMSpace;
using DefaultMemSpace  = Kokkos::CudaSpace;
using DefaultExecSpace = Kokkos::Cuda;
using DefaultLayout    = Kokkos::LayoutLeft;
#elif HAVE_OPENMP
using DefaultMemSpace  = Kokkos::HostSpace;
using DefaultExecSpace = Kokkos::OpenMP;
using DefaultLayout    = Kokkos::LayoutRight;
#elif HAVE_THREADS
using DefaultMemSpace  = Kokkos::HostSpace;
using DefaultExecSpace = Kokkos::Threads;
using DefaultLayout    = Kokkos::LayoutLeft;
#elif HAVE_HIP
using DefaultMemSpace  = Kokkos::Experimental::HIPSpace;
using DefaultExecSpace = Kokkos::Experimental::HIP;
using DefaultLayout    = Kokkos::LayoutLeft;
#else
using DefaultExecSpace = Kokkos::DefaultExecutionSpace;
using DefaultMemSpace  = Kokkos::DefaultExecutionSpace::memory_space;
using DefaultLayout    = Kokkos::LayoutLeft;
#endif


//MACROS to make the code less scary
#define kmalloc(size) ( Kokkos::kokkos_malloc<DefaultMemSpace>(size) )
#define kfree(pnt)        (  Kokkos::kokkos_free(pnt) )
#define ProfileRegionStart  ( Kokkos::Profiling::pushRegion )
#define ProfileRegionEnd  ( Kokkos::Profiling::popRegion )
#define DEFAULTSTRINGARRAY "array_"
#define DEFAULTSTRINGMATRIX "matrix_"

using policy1D = Kokkos::RangePolicy<DefaultExecSpace>;
using policy2D = Kokkos::MDRangePolicy< Kokkos::Rank<2> >;
using policy3D = Kokkos::MDRangePolicy< Kokkos::Rank<3> >;
using policy4D = Kokkos::MDRangePolicy< Kokkos::Rank<4> >;

using TeamPolicy = Kokkos::TeamPolicy<DefaultExecSpace>;
//using mdrange_policy2 = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
//using mdrange_policy3 = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;

using RMatrix1D    = Kokkos::View<double *,DefaultLayout,DefaultExecSpace>;
using RMatrix2D    = Kokkos::View<double **,DefaultLayout,DefaultExecSpace>;
using RMatrix3D    = Kokkos::View<double ***,DefaultLayout,DefaultExecSpace>;
using RMatrix4D    = Kokkos::View<double ****,DefaultLayout,DefaultExecSpace>;
using RMatrix5D    = Kokkos::View<double *****,DefaultLayout,DefaultExecSpace>;
using IMatrix1D    = Kokkos::View<int *,DefaultLayout,DefaultExecSpace>;
using IMatrix2D    = Kokkos::View<int **,DefaultLayout,DefaultExecSpace>;
using IMatrix3D    = Kokkos::View<int ***,DefaultLayout,DefaultExecSpace>;
using IMatrix4D    = Kokkos::View<int ****,DefaultLayout,DefaultExecSpace>;
using IMatrix5D    = Kokkos::View<int *****,DefaultLayout,DefaultExecSpace>;
using SVar         = Kokkos::View<size_t,DefaultLayout,DefaultExecSpace>;
using SArray2D     = Kokkos::View<size_t **,DefaultLayout,DefaultExecSpace>;
using SArray3D     = Kokkos::View<size_t ***,DefaultLayout,DefaultExecSpace>;
using SArray4D     = Kokkos::View<size_t ****,DefaultLayout,DefaultExecSpace>;
using SArray5D     = Kokkos::View<size_t *****,DefaultLayout,DefaultExecSpace>;

using SHArray1D     = Kokkos::View<size_t *,DefaultLayout,Kokkos::HostSpace>;
#endif

//To disable asserts, uncomment the following line
//#define NDEBUG





#ifdef HAVE_KOKKOS
namespace mtr
{

/*! \brief Kokkos version of the serial FArray class.
 *
 *  This is the Kokkos version of the serial FArray class.
 *  Its usage is analagous to that of the serial FArr5class, and it is to be
 *  used in Kokkos-specific code.
 */


template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace, typename MemoryTraits = void>
class FArrayKokkos {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace, MemoryTraits>;
    
private:
    size_t dims_[7];
    size_t order_;
    size_t length_;
    TArray1D this_array_;

public:

    /*!
     * \brief Default constructor
     */
    FArrayKokkos();

    /*!
     * \brief An overloaded constructor used to construct an 1D FArrayKokkos
              object.

        \param dim0 the length of the first dimension
     */
    FArrayKokkos(size_t dim0, const std::string& tag_string = DEFAULTSTRINGARRAY);

    /*!
     * \brief An overloaded constructor used to construct a 2D FArrayKokkos
              object.

        \param dim0 the length of the first dimension
        \param dim1 the length of the second dimension
     */
    FArrayKokkos(size_t dim0, size_t dim1, const std::string& tag_string = DEFAULTSTRINGARRAY);

    /*!
     * \brief An overloaded constructor used to construct a 3D FArrayKokkos
              object.

        \param dim0 the length of the first dimension
        \param dim1 the length of the second dimension
        \param dim2 the length of the third dimension
     */
    FArrayKokkos(size_t dim0, size_t dim1, size_t dim2, const std::string& tag_string = DEFAULTSTRINGARRAY);

    FArrayKokkos(size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, const std::string& tag_string = DEFAULTSTRINGARRAY);

    FArrayKokkos(size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, const std::string& tag_string = DEFAULTSTRINGARRAY);

    FArrayKokkos(size_t dim0, size_t sone_dim2, size_t dim2,
                 size_t dim3, size_t dim4, size_t dim5, const std::string& tag_string = DEFAULTSTRINGARRAY);

    FArrayKokkos(size_t dim0, size_t sone_dim2, size_t dim2,
                 size_t dim3, size_t dim4, size_t dim5,
                 size_t dim6, const std::string& tag_string = DEFAULTSTRINGARRAY);
    
    // Overload operator() to acces data
    // from 1D to 6D
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i) const;
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_INLINE_FUNCTION
    T& operator() (size_t i, size_t j, size_t k) const;

    KOKKOS_INLINE_FUNCTION
    T& operator() (size_t i, size_t j, size_t k,
                   size_t l) const;

    KOKKOS_INLINE_FUNCTION
    T& operator() (size_t i, size_t j, size_t k,
                   size_t l, size_t m) const;

    KOKKOS_INLINE_FUNCTION
    T& operator() (size_t i, size_t j, size_t k,
                   size_t l, size_t m, size_t n) const;

    KOKKOS_INLINE_FUNCTION
    T& operator() (size_t i, size_t j, size_t k,
                   size_t l, size_t m, size_t n, size_t o) const;

    // Overload = operator
    KOKKOS_INLINE_FUNCTION
    FArrayKokkos& operator= (const FArrayKokkos<T,Layout,ExecSpace,MemoryTraits> &temp);

    KOKKOS_INLINE_FUNCTION
    size_t size() const;
    
    KOKKOS_INLINE_FUNCTION
    size_t extent() const;
    
    KOKKOS_INLINE_FUNCTION
    size_t dims(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    size_t order() const;
   
    KOKKOS_INLINE_FUNCTION
    T* pointer() const;
    
    //return kokkos view
    KOKKOS_INLINE_FUNCTION
    TArray1D get_kokkos_view() const;

    // Destructor
    KOKKOS_INLINE_FUNCTION
    ~FArrayKokkos();

}; //end of FArrayKokkos declarations

// Default constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::FArrayKokkos() {
    length_ = order_ = 0;
    for (int i = 0; i < 7; i++) {
        dims_[i] = 0;
    }
}

// Overloaded 1D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::FArrayKokkos(size_t dim0, const std::string& tag_string){
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    dims_[0] = dim0;
    order_ = 1;
    length_ = dim0;
    this_array_ = TArray1D(tag_string, length_);
}

// Overloaded 2D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::FArrayKokkos(size_t dim0, size_t dim1, const std::string& tag_string) {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    order_ = 2;
    length_ = (dim0 * dim1);
    this_array_ = TArray1D(tag_string, length_);
}

// Overloaded 3D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::FArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, const std::string& tag_string) {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    order_ = 3;
    length_ = (dim0 * dim1 * dim2);
    this_array_ = TArray1D(tag_string, length_);
}

// Overloaded 4D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::FArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3, const std::string& tag_string) {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    order_ = 4;
    length_ = (dim0 * dim1 * dim2 * dim3);
    this_array_ = TArray1D(tag_string, length_);
}

// Overloaded 5D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::FArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3,
                              size_t dim4, const std::string& tag_string) {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    order_ = 5;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4);
    this_array_ = TArray1D(tag_string, length_);
}

// Overloaded 6D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::FArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3,
                              size_t dim4, size_t dim5, const std::string& tag_string) {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    order_ = 6;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4 * dim5);
    this_array_ = TArray1D(tag_string, length_);
}

// Overloaded 7D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::FArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3,
                              size_t dim4, size_t dim5,
                              size_t dim6, const std::string& tag_string) {
    
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    dims_[6] = dim6;
    order_ = 7;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4 * dim5 * dim6);
    this_array_ = TArray1D(tag_string, length_);
}

// Definitions of overload operator()
// for 1D to 7D
// Note: the indices for array all start at 0

// 1D
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()( size_t i) const {
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in FArrayKokkos 1D!");
    assert( i >= 0 && i < dims_[0] && "i is out of bounds in FArrayKokkos 1D!");
    return this_array_(i);
}

// 2D
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j) const {
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in FArrayKokkos 2D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in FArrayKokkos 2D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in FArrayKokkos 2D!");
    return this_array_(i + (j * dims_[0]));
}

// 3D
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k) const {
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in FArrayKokkos 3D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in FArrayKokkos 3D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in FArrayKokkos 3D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in FArrayKokkos 3D!");
    return this_array_(i + (j * dims_[0])
                         + (k * dims_[0] * dims_[1]));
}

// 4D
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in FArrayKokkos 4D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in FArrayKokkos 4D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in FArrayKokkos 4D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in FArrayKokkos 4D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in FArrayKokkos 4D!");
    return this_array_(i + (j * dims_[0])
                         + (k * dims_[0] * dims_[1])
                         + (l * dims_[0] * dims_[1] * dims_[2]));
}

// 5D
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m) const {
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in FArrayKokkos 5D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in FArrayKokkos 5D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in FArrayKokkos 5D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in FArrayKokkos 5D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in FArrayKokkos 5D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in FArrayKokkos 5D!");
    return this_array_(i + (j * dims_[0])
                         + (k * dims_[0] * dims_[1])
                         + (l * dims_[0] * dims_[1] * dims_[2])
                         + (m * dims_[0] * dims_[1] * dims_[2] * dims_[3]));
}

// 6D
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n) const {
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in FArrayKokkos 6D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in FArrayKokkos 6D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in FArrayKokkos 6D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in FArrayKokkos 6D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in FArrayKokkos 6D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in FArrayKokkos 6D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in FArrayKokkos 6D!");
    return this_array_(i + (j * dims_[0])
                         + (k * dims_[0] * dims_[1])
                         + (l * dims_[0] * dims_[1] * dims_[2])
                         + (m * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                         + (n * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4]));
}

// 7D
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n, size_t o) const {
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in FArrayKokkos 7D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in FArrayKokkos 7D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in FArrayKokkos 7D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in FArrayKokkos 7D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in FArrayKokkos 7D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in FArrayKokkos 7D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in FArrayKokkos 7D!");
    assert(o >= 0 && o < dims_[6] && "o is out of bounds in FArrayKokkos 7D!");
    return this_array_(i + (j * dims_[0])
                         + (k * dims_[0] * dims_[1])
                         + (l * dims_[0] * dims_[1] * dims_[2])
                         + (m * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                         + (n * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4])
                         + (o * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4] * dims_[5]));
}

// Overload = operator
// for object assingment THIS = FArrayKokkos<> TEMP(n,m,,,,)
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>& FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator= (const FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>& temp) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    if (this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for

        order_ = temp.order_;
        length_ = temp.length_;
        this_array_ = temp.this_array_;
    }
    return *this;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::size() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::extent() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::dims(size_t i) const {
    assert(i < order_ && "FArrayKokkos order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to FArrayKokkos dims is out of bounds!");
    return dims_[i];
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::order() const {
    return order_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::pointer() const {
    return this_array_.data();
}

//return the stored Kokkos view
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
Kokkos::View<T*, Layout, ExecSpace, MemoryTraits> FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::get_kokkos_view() const {
    return this_array_;
}

// Destructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::~FArrayKokkos() {}

////////////////////////////////////////////////////////////////////////////////
// End of FArrayKokkos
////////////////////////////////////////////////////////////////////////////////

/*! \brief Kokkos version of the serial ViewFArray class.
 *
 */
template <typename T>
class ViewFArrayKokkos {

private:
    size_t dims_[7];
    size_t order_;
    size_t length_;
    T* this_array_;

public:
    KOKKOS_INLINE_FUNCTION
    ViewFArrayKokkos();

    KOKKOS_INLINE_FUNCTION
    ViewFArrayKokkos(T* some_array, size_t dim0);
    
    KOKKOS_INLINE_FUNCTION
    ViewFArrayKokkos(T* some_array, size_t dim0, size_t dim1);

    KOKKOS_INLINE_FUNCTION
    ViewFArrayKokkos(T* some_array, size_t dim0, size_t dim1, size_t dim2);
    
    KOKKOS_INLINE_FUNCTION
    ViewFArrayKokkos(T* some_array, size_t dim0, size_t dim1, size_t dim2,
                     size_t dim3);
    
    KOKKOS_INLINE_FUNCTION
    ViewFArrayKokkos(T* some_array, size_t dim0, size_t dim1, size_t dim2,
                     size_t dim3, size_t dim4);
    
    KOKKOS_INLINE_FUNCTION
    ViewFArrayKokkos(T* some_array, size_t dim0, size_t dim1, size_t dim2,
                     size_t dim3, size_t dim4, size_t dim5);
    
    KOKKOS_INLINE_FUNCTION
    ViewFArrayKokkos(T* some_array, size_t dim0, size_t dim1, size_t dim2,
                     size_t dim3, size_t dim4, size_t dim5, size_t dim6);
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k,
                  size_t l) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k,
                  size_t l, size_t m) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k,
                  size_t l, size_t m, size_t n) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k,
                  size_t l, size_t m, size_t n, size_t o) const;

    
    KOKKOS_INLINE_FUNCTION
    size_t size() const;
    
    KOKKOS_INLINE_FUNCTION
    size_t extent() const;

    KOKKOS_INLINE_FUNCTION
    size_t dims(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    size_t order() const;
 
    KOKKOS_INLINE_FUNCTION
    T* pointer() const;

    KOKKOS_INLINE_FUNCTION
    ~ViewFArrayKokkos();

}; // End of ViewFArrayKokkos declarations

// Default constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewFArrayKokkos<T>::ViewFArrayKokkos() {
    length_ = order_ = 0;
    for (int i = 0; i < 7; i++) {
        dims_[i] = 0;
    }
}

// Overloaded 1D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim0) {
    dims_[0] = dim0;
    order_ = 1;
    length_ = dim0;
    this_array_ = some_array;
}

// Overloaded 2D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim0, size_t dim1) {
    dims_[0] = dim0;
    dims_[1] = dim1;
    order_ = 2;
    length_ = (dim0 * dim1);
    this_array_ = some_array;
}

// Overloaded 3D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim0, size_t dim1,
                                      size_t dim2) {
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    order_ = 3;
    length_ = (dim0 * dim1 * dim2);
    this_array_ = some_array;
}

// Overloaded 4D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim0, size_t dim1,
                                      size_t dim2, size_t dim3) {
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    order_ = 4;
    length_ = (dim0 * dim1 * dim2 * dim3);
    this_array_ = some_array;
}

// Overloaded 5D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim0, size_t dim1,
                                      size_t dim2, size_t dim3, size_t dim4) {
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    order_ = 5;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4);
    this_array_ = some_array;
}

// Overloaded 6D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim0, size_t dim1,
                                      size_t dim2, size_t dim3, size_t dim4,
                                      size_t dim5) {
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    order_ = 6;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4 * dim5);
    this_array_ = some_array;
}

// Overloaded 7D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim0, size_t dim1,
                                      size_t dim2, size_t dim3, size_t dim4,
                                      size_t dim5, size_t dim6) {
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    dims_[6] = dim6;
    order_ = 7;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4 * dim5 * dim6);
    this_array_ = some_array;
}

// Overloaded operator() for 1D array access
template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i) const {
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in ViewFArrayKokkos 1D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewFArrayKokkos 1D!");
    return this_array_[i];
}

//2D
template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i, size_t j) const {
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in ViewFArrayKokkos 2D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewFArrayKokkos 2D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewFArrayKokkos 2D!");
    return this_array_[i + (j * dims_[0])];
}

//3D
template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i, size_t j, size_t k) const {
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in ViewFArrayKokkos 3D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewFArrayKokkos 3D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewFArrayKokkos 3D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewFArrayKokkos 3D!");
    return this_array_[i + (j * dims_[0])
                         + (k * dims_[0] * dims_[1])];
}

//4D
template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i, size_t j, size_t k,
                                   size_t l) const {
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in ViewFArrayKokkos 4D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewFArrayKokkos 4D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewFArrayKokkos 4D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewFArrayKokkos 4D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in ViewFArrayKokkos 4D!");
    return this_array_[i + (j * dims_[0])
                         + (k * dims_[0] * dims_[1])
                         + (l * dims_[0] * dims_[1] *dims_[2])];
}

//5D
template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i, size_t j, size_t k,
                                   size_t l, size_t m) const {
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in ViewFArrayKokkos 5D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewFArrayKokkos 5D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewFArrayKokkos 5D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewFArrayKokkos 5D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in ViewFArrayKokkos 5D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in ViewFArrayKokkos 5D!");
    return this_array_[i + (j * dims_[0])
                         + (k * dims_[0] * dims_[1])
                         + (l * dims_[0] * dims_[1] * dims_[2])
                         + (m * dims_[0] * dims_[1] * dims_[2] * dims_[3])];
}

//6D
template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i, size_t j, size_t k,
                                   size_t l, size_t m, size_t n) const {
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in ViewFArrayKokkos 6D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewFArrayKokkos 6D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewFArrayKokkos 6D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewFArrayKokkos 6D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in ViewFArrayKokkos 6D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in ViewFArrayKokkos 6D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in ViewFArrayKokkos 6D!");
    return this_array_[i + (j * dims_[0])
                         + (k * dims_[0] * dims_[1])
                         + (l * dims_[0] * dims_[1] * dims_[2])
                         + (m * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                         + (n * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4])];
}

//7D
template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i, size_t j, size_t k,
                                   size_t l, size_t m, size_t n,
                                   size_t o) const {
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in ViewFArrayKokkos 7D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewFArrayKokkos 7D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewFArrayKokkos 7D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewFArrayKokkos 7D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in ViewFArrayKokkos 7D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in ViewFArrayKokkos 7D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in ViewFArrayKokkos 7D!");
    assert(o >= 0 && o < dims_[6] && "o is out of bounds in ViewFArrayKokkos 7D!");
    return this_array_[i + (j * dims_[0])
                         + (k * dims_[0] * dims_[1])
                         + (l * dims_[0] * dims_[1] * dims_[2])
                         + (m * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                         + (n * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4])
                         + (o * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4] * dims_[5])];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
size_t ViewFArrayKokkos<T>::size() const {
    return length_;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
size_t ViewFArrayKokkos<T>::extent() const {
    return length_;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
size_t ViewFArrayKokkos<T>::dims(size_t i) const {
    assert(i < order_ && "ViewFArrayKokkos order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to ViewFArrayKokkos dims is out of bounds!");
    return dims_[i];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
size_t ViewFArrayKokkos<T>::order() const {
    return order_;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T* ViewFArrayKokkos<T>::pointer() const {
    return this_array_;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
ViewFArrayKokkos<T>::~ViewFArrayKokkos() {}

////////////////////////////////////////////////////////////////////////////////
// End of ViewFArrayKokkos
////////////////////////////////////////////////////////////////////////////////

/*! \brief Kokkos version of the serial FMatrix class.
 *
 */
template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace, typename MemoryTraits = void>
class FMatrixKokkos {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace, MemoryTraits>;
    
private:

    size_t dims_[7];
    size_t order_;
    size_t length_;
    TArray1D this_matrix_;

public:
    FMatrixKokkos();

    FMatrixKokkos(size_t dim1, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    FMatrixKokkos(size_t dim1, size_t dim2, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    FMatrixKokkos(size_t dim1, size_t dim2, size_t dim3, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    FMatrixKokkos(size_t dim1, size_t dim2, size_t dim3,
                  size_t dim4, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    FMatrixKokkos(size_t dim1, size_t dim2, size_t dim3,
                  size_t dim4, size_t dim5, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    FMatrixKokkos(size_t dim1, size_t dim2, size_t dim3,
                  size_t dim4, size_t dim5, size_t dim6, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    FMatrixKokkos(size_t dim1, size_t dim2, size_t dim3,
                  size_t dim4, size_t dim5, size_t dim6,
                  size_t dim7, const std::string& tag_string = DEFAULTSTRINGMATRIX);
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n, size_t o) const;
    
    KOKKOS_INLINE_FUNCTION
    FMatrixKokkos& operator=(const FMatrixKokkos& temp);

    KOKKOS_INLINE_FUNCTION
    size_t size() const;
    
    KOKKOS_INLINE_FUNCTION
    size_t extent() const;

    KOKKOS_INLINE_FUNCTION
    size_t dims(size_t i) const;
 
    KOKKOS_INLINE_FUNCTION
    size_t order() const;
   
    KOKKOS_INLINE_FUNCTION
    T* pointer() const;
    
    //return kokkos view
    KOKKOS_INLINE_FUNCTION
    TArray1D get_kokkos_view() const;

    KOKKOS_INLINE_FUNCTION
    ~FMatrixKokkos();

}; // End of FMatrixKokkos

// Default constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::FMatrixKokkos() {
    length_ = order_ = 0;
    for (int i = 0; i < 7; i++) {
        dims_[i] = 0;
    }
}

// Overloaded 1D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::FMatrixKokkos(size_t dim1, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim1;
    order_ = 1;
    length_ = dim1;
    this_matrix_ = TArray1D(tag_string, length_);
}

// Overloaded 2D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::FMatrixKokkos(size_t dim1, size_t dim2, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    order_ = 2;
    length_ = (dim1 * dim2);
    this_matrix_ = TArray1D(tag_string, length_);
}

// Overloaded 3D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::FMatrixKokkos(size_t dim1, size_t dim2,
                                size_t dim3, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    order_ = 3;
    length_ = (dim1 * dim2 * dim3);
    this_matrix_ = TArray1D(tag_string, length_);
}

// Overloaded 4D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::FMatrixKokkos(size_t dim1, size_t dim2,
                                size_t dim3, size_t dim4, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    order_ = 4;
    length_ = (dim1 * dim2 * dim3 * dim4);
    this_matrix_ = TArray1D(tag_string, length_);
}

// Overloaded 5D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::FMatrixKokkos(size_t dim1, size_t dim2,
                                size_t dim3, size_t dim4,
                                size_t dim5, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    order_ = 5;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5);
    this_matrix_ = TArray1D(tag_string, length_);
}

// Overloaded 5D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::FMatrixKokkos(size_t dim1, size_t dim2,
                                size_t dim3, size_t dim4,
                                size_t dim5, size_t dim6, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    order_ = 6;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5 * dim6);
    this_matrix_ = TArray1D(tag_string, length_);
}

// Overloaded 5D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::FMatrixKokkos(size_t dim1, size_t dim2,
                                size_t dim3, size_t dim4,
                                size_t dim5, size_t dim6,
                                size_t dim7, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    dims_[6] = dim7;
    order_ = 7;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5 * dim6 * dim7);
    this_matrix_ = TArray1D(tag_string, length_);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i) const {
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in FMatrixKokkos 1D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in FMatrixKokkos in 1D!");
    return this_matrix_((i - 1));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j) const {
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in FMatrixKokkos 2D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in FMatrixKokkos in 2D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in FMatrixKokkos in 2D!");
    return this_matrix_((i - 1) + ((j - 1) * dims_[0]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k) const {
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in FMatrixKokkos 3D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in FMatrixKokkos in 3D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in FMatrixKokkos in 3D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in FMatrixKokkos in 3D!");
    return this_matrix_((i - 1) + ((j - 1) * dims_[0])
                                + ((k - 1) * dims_[0] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in FMatrixKokkos 4D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in FMatrixKokkos in 4D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in FMatrixKokkos in 4D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in FMatrixKokkos in 4D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in FMatrixKokkos in 4D!");
    return this_matrix_((i - 1) + ((j - 1) * dims_[0])
                                + ((k - 1) * dims_[0] * dims_[1])
                                + ((l - 1) * dims_[0] * dims_[1] * dims_[2]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                                size_t m) const {
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in FMatrixKokkos 5D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in FMatrixKokkos in 5D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in FMatrixKokkos in 5D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in FMatrixKokkos in 5D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in FMatrixKokkos in 5D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in FMatrixKokkos in 5D!");
    return this_matrix_((i - 1) + ((j - 1) * dims_[0])
                                + ((k - 1) * dims_[0] * dims_[1])
                                + ((l - 1) * dims_[0] * dims_[1] * dims_[2])
                                + ((m - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                                size_t m, size_t n) const {
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in FMatrixKokkos 6D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in FMatrixKokkos in 6D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in FMatrixKokkos in 6D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in FMatrixKokkos in 6D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in FMatrixKokkos in 6D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in FMatrixKokkos in 6D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in FMatrixKokkos in 6D!");
    return this_matrix_[(i - 1) + ((j - 1) * dims_[0])
                                + ((k - 1) * dims_[0] * dims_[1])
                                + ((l - 1) * dims_[0] * dims_[1] * dims_[2])
                                + ((m - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                                + ((n - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4])];
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                                size_t m, size_t n, size_t o) const {
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in FMatrixKokkos 7D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in FMatrixKokkos in 7D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in FMatrixKokkos in 7D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in FMatrixKokkos in 7D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in FMatrixKokkos in 7D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in FMatrixKokkos in 7D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in FMatrixKokkos in 7D!");
    assert(o >= 1 && o <= dims_[6] && "o is out of bounds in FMatrixKokkos in 7D!");
    return this_matrix_[(i - 1) + ((j - 1) * dims_[0])
                                + ((k - 1) * dims_[0] * dims_[1])
                                + ((l - 1) * dims_[0] * dims_[1] * dims_[2])
                                + ((m - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                                + ((n - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4])
                                + ((o - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4] * dims_[5])];
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>& FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator=(const FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>& temp) {
    // Do nothing if the assignment is of the form x = x
    if (this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for

        order_ = temp.order_;
        length_ = temp.length_;
        this_matrix_ = temp.this_matrix_;
    }
    return *this;
}



template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::size() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::extent() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::dims(size_t i) const {
    i--;
    assert(i < order_ && "FMatrixKokkos order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to FMatrixKokkos dims is out of bounds!");
    return dims_[i];
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::order() const {
    return order_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::pointer() const {
    return this_matrix_.data();
}

//return the stored Kokkos view
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
Kokkos::View<T*, Layout, ExecSpace, MemoryTraits> FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::get_kokkos_view() const {
    return this_matrix_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::~FMatrixKokkos() {}

////////////////////////////////////////////////////////////////////////////////
// End of FMatrixKokkos
////////////////////////////////////////////////////////////////////////////////

/*! \brief Kokkos version of the serial ViewFMatrix class.
 *
 */
template <typename T>
class ViewFMatrixKokkos {

private:

    size_t dims_[7];
    size_t order_;
    size_t length_;
    T* this_matrix_;
    
public:

    KOKKOS_INLINE_FUNCTION
    ViewFMatrixKokkos();
    
    KOKKOS_INLINE_FUNCTION
    ViewFMatrixKokkos(T* some_matrix, size_t dim1);
    
    KOKKOS_INLINE_FUNCTION
    ViewFMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2);
    
    KOKKOS_INLINE_FUNCTION
    ViewFMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2,
                      size_t dim3);
    
    KOKKOS_INLINE_FUNCTION
    ViewFMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2,
                      size_t dim3, size_t dim4);
    
    KOKKOS_INLINE_FUNCTION
    ViewFMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2,
                      size_t dim3, size_t dim4, size_t dim5);
    
    KOKKOS_INLINE_FUNCTION
    ViewFMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2,
                      size_t dim3, size_t dim4, size_t dim5,
                      size_t dim6);
    
    KOKKOS_INLINE_FUNCTION
    ViewFMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2,
                      size_t dim3, size_t dim4, size_t dim5,
                      size_t dim6, size_t dim7);
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i) const;
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;
        
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n) const;
 
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n, size_t o) const;
    
    KOKKOS_INLINE_FUNCTION
    size_t size() const;
    
    KOKKOS_INLINE_FUNCTION
    size_t extent() const;

    KOKKOS_INLINE_FUNCTION
    size_t dims(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    size_t order() const;
 
    KOKKOS_INLINE_FUNCTION
    T* pointer() const;

    KOKKOS_INLINE_FUNCTION
    ~ViewFMatrixKokkos();
    
}; // end of ViewFMatrixKokkos

// Default constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewFMatrixKokkos<T>::ViewFMatrixKokkos() {
    length_ = order_ = 0;
    for (int i = 0; i < 7; i++) {
        dims_[i] = 0;
    }
}

// Overloaded 1D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewFMatrixKokkos<T>::ViewFMatrixKokkos(T* some_matrix, size_t dim1) {
    dims_[0] = dim1;
    order_ = 1;
    length_ = dim1;
    this_matrix_ = some_matrix;
}

// Overloaded 2D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewFMatrixKokkos<T>::ViewFMatrixKokkos(T* some_matrix, size_t dim1,
                                        size_t dim2) {
    dims_[0] = dim1;
    dims_[1] = dim2;
    order_ = 2;
    length_ = (dim1 * dim2);
    this_matrix_ = some_matrix;
}

// Overloaded 3D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewFMatrixKokkos<T>::ViewFMatrixKokkos(T* some_matrix, size_t dim1,
                                        size_t dim2, size_t dim3) {
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    order_ = 3;
    length_ = (dim1 * dim2 * dim3);
    this_matrix_ = some_matrix;
}

// Overloaded 4D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewFMatrixKokkos<T>::ViewFMatrixKokkos(T* some_matrix, size_t dim1,
                                        size_t dim2, size_t dim3,
                                        size_t dim4) {
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    order_ = 4;
    length_ = (dim1 * dim2 * dim3 * dim4);
    this_matrix_ = some_matrix;
}

// Overloaded 5D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewFMatrixKokkos<T>::ViewFMatrixKokkos(T* some_matrix, size_t dim1,
                                        size_t dim2, size_t dim3,
                                        size_t dim4, size_t dim5) {
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    order_ = 5;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5);
    this_matrix_ = some_matrix;
}

// Overloaded 6D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewFMatrixKokkos<T>::ViewFMatrixKokkos(T* some_matrix, size_t dim1,
                                        size_t dim2, size_t dim3,
                                        size_t dim4, size_t dim5,
                                        size_t dim6) {
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    order_ = 6;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5 * dim6);
    this_matrix_ = some_matrix;
}

// Overloaded 6D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewFMatrixKokkos<T>::ViewFMatrixKokkos(T* some_matrix, size_t dim1,
                                        size_t dim2, size_t dim3,
                                        size_t dim4, size_t dim5,
                                        size_t dim6, size_t dim7) {
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    dims_[6] = dim7;
    order_ = 7;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5 * dim6 * dim7);
    this_matrix_ = some_matrix;
}


template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i) const {
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in ViewFMatrixKokkos 1D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewFMatrixKokkos 1D!");
    return this_matrix_[(i - 1)];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i, size_t j) const {
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in ViewFMatrixKokkos 2D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewFMatrixKokkos 2D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewFMatrixKokkos 2D!");
    return this_matrix_[(i - 1) + ((j - 1) * dims_[0])];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k) const
{
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in ViewFMatrixKokkos 3D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewFMatrixKokkos 3D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewFMatrixKokkos 3D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in ViewFMatrixKokkos 3D!");
    
    return this_matrix_[(i - 1) + ((j - 1) * dims_[0])
                                + ((k - 1) * dims_[0] * dims_[1])];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k,
                                    size_t l) const {
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in ViewFMatrixKokkos 4D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewFMatrixKokkos 4D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewFMatrixKokkos 4D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in ViewFMatrixKokkos 4D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in ViewFMatrixKokkos 4D!");
    return this_matrix_[(i - 1) + ((j - 1) * dims_[0])
                                + ((k - 1) * dims_[0] * dims_[1])
                                + ((l - 1) * dims_[0] * dims_[1] * dims_[2])];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l,
                                    size_t m) const {
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in ViewFMatrixKokkos 5D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewFMatrixKokkos 5D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewFMatrixKokkos 5D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in ViewFMatrixKokkos 5D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in ViewFMatrixKokkos 5D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in ViewFMatrixKokkos 5D!");
    return this_matrix_[(i - 1) + ((j - 1) * dims_[0])
                                + ((k - 1) * dims_[0] * dims_[1])
                                + ((l - 1) * dims_[0] * dims_[1] * dims_[2])
                                + ((m - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3])];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l,
                                    size_t m, size_t n) const
{
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in ViewFMatrixKokkos 6D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewFMatrixKokkos 6D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewFMatrixKokkos 6D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in ViewFMatrixKokkos 6D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in ViewFMatrixKokkos 6D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in ViewFMatrixKokkos 6D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in ViewFMatrixKokkos 6D!");
    return this_matrix_[(i - 1) + ((j - 1) * dims_[0])
                                + ((k - 1) * dims_[0] * dims_[1])
                                + ((l - 1) * dims_[0] * dims_[1] * dims_[2])
                                + ((m - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                                + ((n - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4])];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l,
                                    size_t m, size_t n, size_t o) const
{
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in ViewFMatrixKokkos 7D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewFMatrixKokkos 7D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewFMatrixKokkos 7D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in ViewFMatrixKokkos 7D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in ViewFMatrixKokkos 7D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in ViewFMatrixKokkos 7D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in ViewFMatrixKokkos 7D!");
    assert(o >= 1 && o <= dims_[6] && "o is out of bounds in ViewFMatrixKokkos 7D!");
    return this_matrix_[(i - 1) + ((j - 1) * dims_[0])
                                + ((k - 1) * dims_[0] * dims_[1])
                                + ((l - 1) * dims_[0] * dims_[1] * dims_[2])
                                + ((m - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                                + ((n - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4])
                                + ((o - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4] * dims_[5])];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
size_t ViewFMatrixKokkos<T>::size() const {
    return length_;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
size_t ViewFMatrixKokkos<T>::extent() const {
    return length_;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
size_t ViewFMatrixKokkos<T>::dims(size_t i) const {
    i--;
    assert(i < order_ && "ViewFMatrixKokkos order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to ViewFMatrixKokkos dims is out of bounds!");
    return dims_[i];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
size_t ViewFMatrixKokkos<T>::order() const {
    return order_;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T* ViewFMatrixKokkos<T>::pointer() const {
    return this_matrix_;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
ViewFMatrixKokkos<T>::~ViewFMatrixKokkos() {}

////////////////////////////////////////////////////////////////////////////////
// End of ViewFMatrixKokkos
////////////////////////////////////////////////////////////////////////////////


/////////////////////////
// DFArrayKokkos:  Dual type for managing data on both CPU and GPU.
/////////////////////////
template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace, typename MemoryTraits = void>
class DFArrayKokkos {

    // this is manage
    using TArray1D = Kokkos::DualView<T*, Layout, ExecSpace, MemoryTraits>;
    
private:
    size_t dims_[7];
    size_t length_;
    size_t order_;  // tensor order (rank)
    TArray1D this_array_;

public:
    // Data member to access host view
    ViewFArray <T> host;

    DFArrayKokkos();
    
    DFArrayKokkos(size_t dim0, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DFArrayKokkos(size_t dim0, size_t dim1, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DFArrayKokkos (size_t dim0, size_t dim1, size_t dim2, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DFArrayKokkos(size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DFArrayKokkos(size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DFArrayKokkos(size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, size_t dim5, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DFArrayKokkos(size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, size_t dim5,
                 size_t dim6, const std::string& tag_string = DEFAULTSTRINGARRAY);
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n, size_t o) const;
    
    KOKKOS_INLINE_FUNCTION
    DFArrayKokkos& operator=(const DFArrayKokkos& temp);

    // GPU Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t size() const;

    // Host Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t extent() const;

    KOKKOS_INLINE_FUNCTION
    size_t dims(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    size_t order() const;
 
    // Method returns the raw device pointer of the Kokkos DualView
    KOKKOS_INLINE_FUNCTION
    T* device_pointer() const;

    // Method returns the raw host pointer of the Kokkos DualView
    KOKKOS_INLINE_FUNCTION
    T* host_pointer() const;

    // Method that update host view
    void update_host();

    // Method that update device view
    void update_device();

    // Deconstructor
    KOKKOS_INLINE_FUNCTION
    ~DFArrayKokkos ();

}; // End of DFArrayKokkos declarations

// Default constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DFArrayKokkos() {
    length_ = order_ = 0;
    for (int i = 0; i < 7; i++) {
        dims_[i] = 0;
    }
}

// Overloaded 1D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DFArrayKokkos(size_t dim0, const std::string& tag_string) {
    
    dims_[0] = dim0;
    order_ = 1;
    length_ = dim0;
    this_array_ = TArray1D(tag_string, length_);
    // Create host ViewFArray
    host = ViewFArray <T> (this_array_.h_view.data(), dim0);
}

// Overloaded 2D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DFArrayKokkos(size_t dim0, size_t dim1, const std::string& tag_string) {
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    order_ = 2;
    length_ = (dim0 * dim1);
    this_array_ = TArray1D(tag_string, length_);
    // Create host ViewFArray
    host = ViewFArray <T> (this_array_.h_view.data(), dim0, dim1);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DFArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, const std::string& tag_string) {
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    order_ = 3;
    length_ = (dim0 * dim1 * dim2);
    this_array_ = TArray1D(tag_string, length_);
    // Create host ViewFArray
    host = ViewFArray <T> (this_array_.h_view.data(), dim0, dim1, dim2);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DFArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3, const std::string& tag_string) {
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    order_ = 4;
    length_ = (dim0 * dim1 * dim2 * dim3);
    this_array_ = TArray1D(tag_string, length_);
    // Create host ViewFArray
    host = ViewFArray <T> (this_array_.h_view.data(), dim0, dim1, dim2, dim3);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DFArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3,
                              size_t dim4, const std::string& tag_string) {
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    order_ = 5;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4);
    this_array_ = TArray1D(tag_string, length_);
    // Create host ViewFArray
    host = ViewFArray <T> (this_array_.h_view.data(), dim0, dim1, dim2, dim3, dim4);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DFArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3,
                              size_t dim4, size_t dim5, const std::string& tag_string) {
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    order_ = 6;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4 * dim5);
    this_array_ = TArray1D(tag_string, length_);
    // Create host ViewFArray
    host = ViewFArray <T> (this_array_.h_view.data(), dim0, dim1, dim2, dim3, dim4, dim5);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DFArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3,
                              size_t dim4, size_t dim5,
                              size_t dim6, const std::string& tag_string) {
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    dims_[6] = dim6;
    order_ = 7;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4 * dim5 * dim6);
    this_array_ = TArray1D(tag_string, length_);
    // Create host ViewFArray
    host = ViewFArray <T> (this_array_.h_view.data(), dim0, dim1, dim2, dim3, dim4, dim5, dim6);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i) const {
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in DFArrayKokkos 1D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DFArrayKokkos 1D!");
    return this_array_.d_view(i);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j) const {
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in DFArrayKokkos 2D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DFArrayKokkos 2D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DFArrayKokkos 2D!");
    return this_array_.d_view(i + (j * dims_[0]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k) const {
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in DFArrayKokkos 3D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DFArrayKokkos 3D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DFArrayKokkos 3D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DFArrayKokkos 3D!");
    return this_array_.d_view(i + (j * dims_[0])
                                + (k * dims_[0] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in DFArrayKokkos 4D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DFArrayKokkos 4D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DFArrayKokkos 4D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DFArrayKokkos 4D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in DFArrayKokkos 4D!");
    return this_array_.d_view(i + (j * dims_[0])
                                + (k * dims_[0] * dims_[1])
                                + (l * dims_[0] * dims_[1] * dims_[2]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m) const {
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in DFArrayKokkos 5D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DFArrayKokkos 5D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DFArrayKokkos 5D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DFArrayKokkos 5D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in DFArrayKokkos 5D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in DFArrayKokkos 5D!");
    return this_array_.d_view(i + (j * dims_[0])
                                + (k * dims_[0] * dims_[1])
                                + (l * dims_[0] * dims_[1] * dims_[2])
                                + (m * dims_[0] * dims_[1] * dims_[2] * dims_[3]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n) const {
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in DFArrayKokkos 6D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DFArrayKokkos 6D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DFArrayKokkos 6D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DFArrayKokkos 6D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in DFArrayKokkos 6D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in DFArrayKokkos 6D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in DFArrayKokkos 6D!");
    return this_array_.d_view(i + (j * dims_[0])
                                + (k * dims_[0] * dims_[1])
                                + (l * dims_[0] * dims_[1] * dims_[2])
                                + (m * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                                + (n * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n, size_t o) const {
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in DFArrayKokkos 7D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DFArrayKokkos 7D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DFArrayKokkos 7D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DFArrayKokkos 7D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in DFArrayKokkos 7D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in DFArrayKokkos 7D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in DFArrayKokkos 7D!");
    assert(o >= 0 && o < dims_[6] && "o is out of bounds in DFArrayKokkos 7D!");
    return this_array_.d_view(i + (j * dims_[0])
                                + (k * dims_[0] * dims_[1])
                                + (l * dims_[0] * dims_[1] * dims_[2])
                                + (m * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                                + (n * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4])
                                + (o * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4] * dims_[5]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>& DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator= (const DFArrayKokkos& temp) {
    
    // Do nothing if the assignment is of the form x = x
    if (this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for

        order_ = temp.order_;
        length_ = temp.length_;
        this_array_ = temp.this_array_;
    host = temp.host;
    }
    
    return *this;
}

// Return size
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::size() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::extent() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::dims(size_t i) const {
    assert(i < order_ && "DFArrayKokkos order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to DFArrayKokkos dims is out of bounds!");
    return dims_[i];
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::order() const {
    return order_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::device_pointer() const {
    return this_array_.d_view.data();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::host_pointer() const {
    return this_array_.h_view.data();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
void DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::update_host() {

    this_array_.template modify<typename TArray1D::execution_space>();
    this_array_.template sync<typename TArray1D::host_mirror_space>();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
void DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::update_device() {

    this_array_.template modify<typename TArray1D::host_mirror_space>();
    this_array_.template sync<typename TArray1D::execution_space>();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::~DFArrayKokkos() {}
// End DFArrayKokkos


/////////////////////////
// DViewFArrayKokkos:  The DView means dual view of the data, where data is on both CPU and GPU.
//
// This MATAR type is for accepting a pointer to data on the CPU via the constructor and then it copies the data
// data to the GPU where the member functions and overloads access the data on the GPU. The corresponding
// FArrayKokkos type creates memory on the GPU; likewise, the viewFArrayKokkos accesses data already on the GPU.
// To emphasize, the data must be on the CPU prior to calling the constructor for the DView data type.
/////////////////////////
template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace, typename MemoryTraits = void>
class DViewFArrayKokkos {

    // this is always unmanaged
    using TArray1DHost = Kokkos::View<T*, Layout, HostSpace, MemoryUnmanaged>;
    // this is manage
    using TArray1D     = Kokkos::View<T*, Layout, ExecSpace, MemoryTraits>;
    typename ExecSpace::memory_space memspace;
    
private:
    size_t dims_[7];
    size_t length_;
    size_t order_;  // tensor order (rank)
    TArray1D this_array_;
    TArray1DHost this_array_host_;
    T * temp_inp_array_;

public:
    DViewFArrayKokkos();
    
    DViewFArrayKokkos(T * inp_array, size_t dim0, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, size_t dim5, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, size_t dim5,
                 size_t dim6, const std::string& tag_string = DEFAULTSTRINGARRAY);
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n, size_t o) const;
    
    KOKKOS_INLINE_FUNCTION
    DViewFArrayKokkos& operator=(const DViewFArrayKokkos& temp);

    // GPU Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t size() const;

    // Host Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t extent() const;

    KOKKOS_INLINE_FUNCTION
    size_t dims(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    size_t order() const;
 
    // Method returns the raw device pointer of the Kokkos View
    KOKKOS_INLINE_FUNCTION
    T* device_pointer() const;

    // Method returns the raw host pointer of the Kokkos View
    KOKKOS_INLINE_FUNCTION
    T* host_pointer() const;

    // Data member to access host view
    ViewFArray <T> host;

    // Method that update host view
    void update_host();

    // Method that update device view
    void update_device();

    // Deconstructor
    KOKKOS_INLINE_FUNCTION
    ~DViewFArrayKokkos ();
}; // End of DViewFArrayKokkos


// Default constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFArrayKokkos() {
    length_ = order_ = 0;
    temp_inp_array_ = NULL;
    for (int i = 0; i < 7; i++) {
        dims_[i] = 0;
    }
}

// Overloaded 1D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFArrayKokkos(T * inp_array, size_t dim0, const std::string& tag_string) {
    //using TArray1DHost = Kokkos::View<T*, Layout, HostSpace, MemoryUnmanaged>;
    //using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim0;
    order_ = 1;
    length_ = dim0;
    // Create a 1D host view of the external allocation
    this_array_host_ = TArray1DHost(inp_array, length_);
    // Assign temp point to inp_array pointer that is passed in
    temp_inp_array_ = inp_array;
    // Create a device copy of that host view
    this_array_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_array_host_);
    // Create host ViewFArray. Note: inp_array and this_array_host_.data() are the same pointer
    host = ViewFArray <T> (inp_array, dim0);
}

// Overloaded 2D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1, const std::string& tag_string) {
    //using TArray1DHost = Kokkos::View<T*, Layout, HostSpace, MemoryUnmanaged>;
    //using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    //using TArray1Dtemp = TArray1D::HostMirror;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    order_ = 2;
    length_ = (dim0 * dim1);
    // Create a 1D host view of the external allocation
    this_array_host_ = TArray1DHost(inp_array, length_);
    // Assign temp point to inp_array pointer that is passed in
    temp_inp_array_ = inp_array;
    // Create a device copy of that host view
    this_array_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_array_host_);
    // Create host ViewFArray
    host = ViewFArray <T> (inp_array, dim0, dim1);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2, const std::string& tag_string) {
    //using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    order_ = 3;
    length_ = (dim0 * dim1 * dim2);
    // Create a 1D host view of the external allocation
    this_array_host_ = TArray1DHost(inp_array, length_);
    // Assign temp point to inp_array pointer that is passed in
    temp_inp_array_ = inp_array;
    // Create a device copy of that host view
    this_array_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_array_host_);
    // Create host ViewFArray
    host = ViewFArray <T> (inp_array, dim0, dim1, dim2);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3, const std::string& tag_string) {
    //using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    order_ = 4;
    length_ = (dim0 * dim1 * dim2 * dim3);
    // Create a 1D host view of the external allocation
    this_array_host_ = TArray1DHost(inp_array, length_);
    // Assign temp point to inp_array pointer that is passed in
    temp_inp_array_ = inp_array;
    // Create a device copy of that host view
    this_array_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_array_host_);
    // Create host ViewFArray
    host = ViewFArray <T> (inp_array, dim0, dim1, dim2, dim3);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3, 
                              size_t dim4, const std::string& tag_string) {

    //using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    order_ = 5;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4);
    // Create a 1D host view of the external allocation
    this_array_host_ = TArray1DHost(inp_array, length_);
    // Assign temp point to inp_array pointer that is passed in
    temp_inp_array_ = inp_array;
    // Create a device copy of that host view
    this_array_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_array_host_);
    // Create host ViewFArray
    host = ViewFArray <T> (inp_array, dim0, dim1, dim2, dim3, dim4);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3, 
                              size_t dim4, size_t dim5, const std::string& tag_string) {
    //using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    order_ = 6;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4 * dim5);
    // Create a 1D host view of the external allocation
    this_array_host_ = TArray1DHost(inp_array, length_);
    // Assign temp point to inp_array pointer that is passed in
    temp_inp_array_ = inp_array;
    // Create a device copy of that host view
    this_array_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_array_host_);
    // Create host ViewFArray
    host = ViewFArray <T> (inp_array, dim0, dim1, dim2, dim3, dim4, dim5);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3,
                              size_t dim4, size_t dim5,
                              size_t dim6, const std::string& tag_string) {
    //using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    dims_[6] = dim6;
    order_ = 7;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4 * dim5 * dim6);
    // Create a 1D host view of the external allocation
    this_array_host_ = TArray1DHost(inp_array, length_);
    // Assign temp point to inp_array pointer that is passed in
    temp_inp_array_ = inp_array;
    // Create a device copy of that host view
    this_array_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_array_host_);
    // Create host ViewFArray
    host = ViewFArray <T> (inp_array, dim0, dim1, dim2, dim3, dim4, dim5, dim6);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i) const {
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in DViewFArrayKokkos 1D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DViewFArrayKokkos 1D!");
    return this_array_(i);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j) const {
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in DViewFArrayKokkos 2D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DViewFArrayKokkos 2D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DViewFArrayKokkos 2D!");
    return this_array_(i + (j * dims_[0]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k) const {
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in DViewFArrayKokkos 3D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DViewFArrayKokkos 3D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DViewFArrayKokkos 3D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DViewFArrayKokkos 3D!");
    return this_array_(i + (j * dims_[0])
                         + (k * dims_[0] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in DViewFArrayKokkos 4D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DViewFArrayKokkos 4D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DViewFArrayKokkos 4D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DViewFArrayKokkos 4D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in DViewFArrayKokkos 4D!");
    return this_array_(i + (j * dims_[0])
                         + (k * dims_[0] * dims_[1])
                         + (l * dims_[0] * dims_[1] * dims_[2]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m) const {
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in DViewFArrayKokkos 5D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DViewFArrayKokkos 5D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DViewFArrayKokkos 5D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DViewFArrayKokkos 5D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in DViewFArrayKokkos 5D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in DViewFArrayKokkos 5D!");
    return this_array_(i + (j * dims_[0])
                         + (k * dims_[0] * dims_[1])
                         + (l * dims_[0] * dims_[1] * dims_[2])
                         + (m * dims_[0] * dims_[1] * dims_[2] * dims_[3]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n) const {
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in DViewFArrayKokkos 6D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DViewFArrayKokkos 6D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DViewFArrayKokkos 6D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DViewFArrayKokkos 6D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in DViewFArrayKokkos 6D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in DViewFArrayKokkos 6D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in DViewFArrayKokkos 6D!");
    return this_array_(i + (j * dims_[0])
                         + (k * dims_[0] * dims_[1])
                         + (l * dims_[0] * dims_[1] * dims_[2])
                         + (m * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                         + (n * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n, size_t o) const {
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in DViewFArrayKokkos 7D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DViewFArrayKokkos 7D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DViewFArrayKokkos 7D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DViewFArrayKokkos 7D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in DViewFArrayKokkos 7D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in DViewFArrayKokkos 7D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in DViewFArrayKokkos 7D!");
    assert(o >= 0 && o < dims_[6] && "o is out of bounds in DViewFArrayKokkos 7D!");
    return this_array_(i + (j * dims_[0])
                         + (k * dims_[0] * dims_[1])
                         + (l * dims_[0] * dims_[1] * dims_[2])
                         + (m * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                         + (n * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4])
                         + (o * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4] * dims_[5]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>& DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator= (const DViewFArrayKokkos& temp) {
    
    // Do nothing if the assignment is of the form x = x
    if (this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for

        order_ = temp.order_;
        length_ = temp.length_;
        temp_inp_array_ = temp.temp_inp_array_;
        this_array_host_ = temp.this_array_host_;
        this_array_ = temp.this_array_;
        host = temp.host;
    }
    
    return *this;
}

// Return size
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::size() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::extent() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::dims(size_t i) const {
    assert(i < order_ && "DViewFArrayKokkos order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to DViewFArrayKokkos dims is out of bounds!");
    return dims_[i];
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::order() const {
    return order_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::device_pointer() const {
    return this_array_.data();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::host_pointer() const {
    return this_array_host_.data();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
void DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::update_host() {
    // Deep copy of device view to host view
    deep_copy(this_array_host_, this_array_);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
void DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::update_device() {
    // Deep copy of host view to device view
    deep_copy(this_array_, this_array_host_);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::~DViewFArrayKokkos() {}
// End DViewFArrayKokkos


/////////////////////////
// DFMatrixKokkos
/////////////////////////
template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace, typename MemoryTraits = void>
class DFMatrixKokkos {

    // this is manage
    using TArray1D = Kokkos::DualView<T*, Layout, ExecSpace, MemoryTraits>;
    
private:
    size_t dims_[7];
    size_t length_;
    size_t order_;  // tensor order (rank)
    TArray1D this_matrix_;

public:
    DFMatrixKokkos();
    
    DFMatrixKokkos(size_t dim1, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DFMatrixKokkos(size_t dim1, size_t dim2, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DFMatrixKokkos (size_t dim1, size_t dim2, size_t dim3, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DFMatrixKokkos(size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DFMatrixKokkos(size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, size_t dim5, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DFMatrixKokkos(size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, size_t dim5, size_t dim6, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DFMatrixKokkos(size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, size_t dim5, size_t dim6,
                 size_t dim7, const std::string& tag_string = DEFAULTSTRINGMATRIX);
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n, size_t o) const;
    
    KOKKOS_INLINE_FUNCTION
    DFMatrixKokkos& operator=(const DFMatrixKokkos& temp);

    // GPU Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t size() const;

    // Host Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t extent() const;

    KOKKOS_INLINE_FUNCTION
    size_t dims(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    size_t order() const;
 
    // Method returns the raw device pointer of the Kokkos DualView
    KOKKOS_INLINE_FUNCTION
    T* device_pointer() const;

    // Method returns the raw host pointer of the Kokkos DualView
    KOKKOS_INLINE_FUNCTION
    T* host_pointer() const;

    // Data member to access host view
    ViewFMatrix <T> host;

    // Method that update host view
    void update_host();

    // Method that update device view
    void update_device();

    // Deconstructor
    KOKKOS_INLINE_FUNCTION
    ~DFMatrixKokkos ();
}; // End of DFMatrixKokkos declarations

// Default constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DFMatrixKokkos() {
    length_ = order_ = 0;
    for (int i = 0; i < 7; i++) {
        dims_[i] = 0;
    }
}

// Overloaded 1D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DFMatrixKokkos(size_t dim1, const std::string& tag_string) {
    
    dims_[0] = dim1;
    order_ = 1;
    length_ = dim1;
    this_matrix_ = TArray1D(tag_string, length_);
    // Create host ViewFMatrix
    host = ViewFMatrix <T> (this_matrix_.h_view.data(), dim1);
}

// Overloaded 2D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DFMatrixKokkos(size_t dim1, size_t dim2, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    order_ = 2;
    length_ = (dim1 * dim2);
    this_matrix_ = TArray1D(tag_string, length_);
    // Create host ViewFMatrix
    host = ViewFMatrix <T> (this_matrix_.h_view.data(), dim1, dim2);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DFMatrixKokkos(size_t dim1, size_t dim2,
                              size_t dim3, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    order_ = 3;
    length_ = (dim1 * dim2 * dim3);
    this_matrix_ = TArray1D(tag_string, length_);
    // Create host ViewFMatrix
    host = ViewFMatrix <T> (this_matrix_.h_view.data(), dim1, dim2, dim3);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DFMatrixKokkos(size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    order_ = 4;
    length_ = (dim1 * dim2 * dim3 * dim4);
    this_matrix_ = TArray1D(tag_string, length_);
    // Create host ViewFMatrix
    host = ViewFMatrix <T> (this_matrix_.h_view.data(), dim1, dim2, dim3, dim4);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DFMatrixKokkos(size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4,
                              size_t dim5, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    order_ = 5;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5);
    this_matrix_ = TArray1D(tag_string, length_);
    // Create host ViewFMatrix
    host = ViewFMatrix <T> (this_matrix_.h_view.data(), dim1, dim2, dim3, dim4, dim5);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DFMatrixKokkos(size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4,
                              size_t dim5, size_t dim6, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    order_ = 6;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5 * dim6);
    this_matrix_ = TArray1D(tag_string, length_);
    // Create host ViewFMatrix
    host = ViewFMatrix <T> (this_matrix_.h_view.data(), dim1, dim2, dim3, dim4, dim5, dim6);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DFMatrixKokkos(size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4,
                              size_t dim5, size_t dim6,
                              size_t dim7, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    dims_[6] = dim7;
    order_ = 7;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5 * dim6 * dim7);
    this_matrix_ = TArray1D(tag_string, length_);
    // Create host ViewFMatrix
    host = ViewFMatrix <T> (this_matrix_.h_view.data(), dim1, dim2, dim3, dim4, dim5, dim6, dim7);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i) const {
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in DFMatrixKokkos 1D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DFMatrixKokkos 1D!");
    return this_matrix_.d_view((i - 1));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j) const {
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in DFMatrixKokkos 2D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DFMatrixKokkos 2D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DFMatrixKokkos 2D!");
    return this_matrix_.d_view((i - 1) + ((j - 1) * dims_[0]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k) const {
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in DFMatrixKokkos 3D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DFMatrixKokkos 3D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DFMatrixKokkos 3D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DFMatrixKokkos 3D!");
    return this_matrix_.d_view((i - 1) + ((j - 1) * dims_[0])
                                       + ((k - 1) * dims_[0] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in DFMatrixKokkos 4D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DFMatrixKokkos 4D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DFMatrixKokkos 4D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DFMatrixKokkos 4D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in DFMatrixKokkos 4D!");
    return this_matrix_.d_view((i - 1) + ((j - 1) * dims_[0])
                                       + ((k - 1) * dims_[0] * dims_[1])
                                       + ((l - 1) * dims_[0] * dims_[1] * dims_[2]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m) const {
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in DFMatrixKokkos 5D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DFMatrixKokkos 5D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DFMatrixKokkos 5D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DFMatrixKokkos 5D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in DFMatrixKokkos 5D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in DFMatrixKokkos 5D!");
    return this_matrix_.d_view((i - 1) + ((j - 1) * dims_[0])
                                       + ((k - 1) * dims_[0] * dims_[1])
                                       + ((l - 1) * dims_[0] * dims_[1] * dims_[2])
                                       + ((m - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n) const {
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in DFMatrixKokkos 6D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DFMatrixKokkos 6D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DFMatrixKokkos 6D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DFMatrixKokkos 6D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in DFMatrixKokkos 6D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in DFMatrixKokkos 6D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in DFMatrixKokkos 6D!");
    return this_matrix_.d_view((i - 1) + ((j - 1) * dims_[0])
                                       + ((k - 1) * dims_[0] * dims_[1])
                                       + ((l - 1) * dims_[0] * dims_[1] * dims_[2])
                                       + ((m - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                                       + ((n - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n, size_t o) const {
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in DFMatrixKokkos 7D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DFMatrixKokkos 7D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DFMatrixKokkos 7D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DFMatrixKokkos 7D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in DFMatrixKokkos 7D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in DFMatrixKokkos 7D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in DFMatrixKokkos 7D!");
    assert(o >= 1 && o <= dims_[6] && "o is out of bounds in DFMatrixKokkos 7D!");
    return this_matrix_.d_view((i - 1) + ((j - 1) * dims_[0])
                                       + ((k - 1) * dims_[0] * dims_[1])
                                       + ((l - 1) * dims_[0] * dims_[1] * dims_[2])
                                       + ((m - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                                       + ((n - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4])
                                       + ((o - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4] * dims_[5]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>& DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator= (const DFMatrixKokkos& temp) {
    
    // Do nothing if the assignment is of the form x = x
    if (this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for

        order_ = temp.order_;
        length_ = temp.length_;
        this_matrix_ = temp.this_matrix_;
    host = temp.host;
    }
    
    return *this;
}

// Return size
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::size() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::extent() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::dims(size_t i) const {
    i--;
    assert(i < order_ && "DFMatrixKokkos order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to DFMatrixKokkos dims is out of bounds!");
    return dims_[i];
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::order() const {
    return order_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::device_pointer() const {
    return this_matrix_.d_view.data();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::host_pointer() const {
    return this_matrix_.h_view.data();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
void DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::update_host() {

    this_matrix_.template modify<typename TArray1D::execution_space>();
    this_matrix_.template sync<typename TArray1D::host_mirror_space>();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
void DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::update_device() {

    this_matrix_.template modify<typename TArray1D::host_mirror_space>();
    this_matrix_.template sync<typename TArray1D::execution_space>();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::~DFMatrixKokkos() {}
// End DFMatrixKokkos


/////////////////////////
// DViewFMatrixKokkos
/////////////////////////
template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace, typename MemoryTraits = void>
class DViewFMatrixKokkos {

    // this is always unmanaged
    using TArray1DHost = Kokkos::View<T*, Layout, HostSpace, MemoryUnmanaged>;
    // this is manage
    using TArray1D     = Kokkos::View<T*, Layout, ExecSpace, MemoryTraits>;
    typename ExecSpace::memory_space memspace;
    
private:
    size_t dims_[7];
    size_t length_;
    size_t order_;  // tensor order (rank)
    TArray1D this_matrix_;
    TArray1DHost this_matrix_host_;
    T * temp_inp_matrix_;

public:
    DViewFMatrixKokkos();
    
    DViewFMatrixKokkos(T * inp_matrix, size_t dim1, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, size_t dim5, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, size_t dim5, size_t dim6, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, size_t dim5, size_t dim6,
                 size_t dim7, const std::string& tag_string = DEFAULTSTRINGMATRIX);
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n, size_t o) const;
    
    KOKKOS_INLINE_FUNCTION
    DViewFMatrixKokkos& operator=(const DViewFMatrixKokkos& temp);

    // GPU Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t size() const;

    // Host Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t extent() const;

    KOKKOS_INLINE_FUNCTION
    size_t dims(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    size_t order() const;
 
    // Method returns the raw device pointer of the Kokkos View
    KOKKOS_INLINE_FUNCTION
    T* device_pointer() const;

    // Method returns the raw host pointer of the Kokkos View
    KOKKOS_INLINE_FUNCTION
    T* host_pointer() const;

    // Data member to access host view
    ViewFMatrix <T> host;

    // Method that update host view
    void update_host();

    // Method that update device view
    void update_device();

    // Deconstructor
    KOKKOS_INLINE_FUNCTION
    ~DViewFMatrixKokkos ();
}; // End of DViewFMatrixKokkos


// Default constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFMatrixKokkos() {
    length_ = order_ = 0;
    temp_inp_matrix_ = NULL;
    for (int i = 0; i < 7; i++) {
        dims_[i] = 0;
    }
}

// Overloaded 1D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFMatrixKokkos(T * inp_matrix, size_t dim1, const std::string& tag_string) {
    
    dims_[0] = dim1;
    order_ = 1;
    length_ = dim1;
    // Create a 1D host view of the external allocation
    this_matrix_host_ = TArray1DHost(inp_matrix, length_);
    // Assign temp point to inp_matrix pointer that is passed in
    temp_inp_matrix_ = inp_matrix;
    // Create a device copy of that host view
    this_matrix_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_matrix_host_);
    // Create host ViewFMatrix. Note: inp_matrix and this_matrix_host_.data() are the same pointer
    host = ViewFMatrix <T> (inp_matrix, dim1);
}

// Overloaded 2D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    order_ = 2;
    length_ = (dim1 * dim2);
    // Create a 1D host view of the external allocation
    this_matrix_host_ = TArray1DHost(inp_matrix, length_);
    // Assign temp point to inp_matrix pointer that is passed in
    temp_inp_matrix_ = inp_matrix;
    // Create a device copy of that host view
    this_matrix_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_matrix_host_);
    // Create host ViewFMatrix
    host = ViewFMatrix <T> (inp_matrix, dim1, dim2);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    order_ = 3;
    length_ = (dim1 * dim2 * dim3);
    // Create a 1D host view of the external allocation
    this_matrix_host_ = TArray1DHost(inp_matrix, length_);
    // Assign temp point to inp_matrix pointer that is passed in
    temp_inp_matrix_ = inp_matrix;
    // Create a device copy of that host view
    this_matrix_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_matrix_host_);
    // Create host ViewFMatrix
    host = ViewFMatrix <T> (inp_matrix, dim1, dim2, dim3);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    order_ = 4;
    length_ = (dim1 * dim2 * dim3 * dim4);
    // Create a 1D host view of the external allocation
    this_matrix_host_ = TArray1DHost(inp_matrix, length_);
    // Assign temp point to inp_matrix pointer that is passed in
    temp_inp_matrix_ = inp_matrix;
    // Create a device copy of that host view
    this_matrix_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_matrix_host_);
    // Create host ViewFMatrix
    host = ViewFMatrix <T> (inp_matrix, dim1, dim2, dim3, dim4);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4, 
                              size_t dim5, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    order_ = 5;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5);
    // Create a 1D host view of the external allocation
    this_matrix_host_ = TArray1DHost(inp_matrix, length_);
    // Assign temp point to inp_matrix pointer that is passed in
    temp_inp_matrix_ = inp_matrix;
    // Create a device copy of that host view
    this_matrix_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_matrix_host_);
    // Create host ViewFMatrix
    host = ViewFMatrix <T> (inp_matrix, dim1, dim2, dim3, dim4, dim5);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4, 
                              size_t dim5, size_t dim6, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    order_ = 6;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5 * dim6);
    // Create a 1D host view of the external allocation
    this_matrix_host_ = TArray1DHost(inp_matrix, length_);
    // Assign temp point to inp_matrix pointer that is passed in
    temp_inp_matrix_ = inp_matrix;
    // Create a device copy of that host view
    this_matrix_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_matrix_host_);
    // Create host ViewFMatrix
    host = ViewFMatrix <T> (inp_matrix, dim1, dim2, dim3, dim4, dim5, dim6);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4,
                              size_t dim5, size_t dim6,
                              size_t dim7, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    dims_[6] = dim7;
    order_ = 7;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5 * dim6 * dim7);
    // Create a 1D host view of the external allocation
    this_matrix_host_ = TArray1DHost(inp_matrix, length_);
    // Assign temp point to inp_matrix pointer that is passed in
    temp_inp_matrix_ = inp_matrix;
    // Create a device copy of that host view
    this_matrix_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_matrix_host_);
    // Create host ViewFMatrix
    host = ViewFMatrix <T> (inp_matrix, dim1, dim2, dim3, dim4, dim5, dim6, dim7);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i) const {
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in DViewFMatrixKokkos 1D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DViewFMatrixKokkos 1D!");
    return this_matrix_((i - 1));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j) const {
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in DViewFMatrixKokkos 2D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DViewFMatrixKokkos 2D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DViewFMatrixKokkos 2D!");
    return this_matrix_((i - 1) + ((j - 1) * dims_[0]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k) const {
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in DViewFMatrixKokkos 3D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DViewFMatrixKokkos 3D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DViewFMatrixKokkos 3D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DViewFMatrixKokkos 3D!");
    return this_matrix_((i - 1) + ((j - 1) * dims_[0])
                                + ((k - 1) * dims_[0] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in DViewFMatrixKokkos 4D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DViewFMatrixKokkos 4D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DViewFMatrixKokkos 4D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DViewFMatrixKokkos 4D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in DViewFMatrixKokkos 4D!");
    return this_matrix_((i - 1) + ((j - 1) * dims_[0])
                                + ((k - 1) * dims_[0] * dims_[1])
                                + ((l - 1) * dims_[0] * dims_[1] * dims_[2]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m) const {
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in DViewFMatrixKokkos 5D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DViewFMatrixKokkos 5D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DViewFMatrixKokkos 5D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DViewFMatrixKokkos 5D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in DViewFMatrixKokkos 5D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in DViewFMatrixKokkos 5D!");
    return this_matrix_((i - 1) + ((j - 1) * dims_[0])
                                + ((k - 1) * dims_[0] * dims_[1])
                                + ((l - 1) * dims_[0] * dims_[1] * dims_[2])
                                + ((m - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n) const {
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in DViewFMatrixKokkos 6D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DViewFMatrixKokkos 6D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DViewFMatrixKokkos 6D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DViewFMatrixKokkos 6D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in DViewFMatrixKokkos 6D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in DViewFMatrixKokkos 6D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in DViewFMatrixKokkos 6D!");
    return this_matrix_((i - 1) + ((j - 1) * dims_[0])
                                + ((k - 1) * dims_[0] * dims_[1])
                                + ((l - 1) * dims_[0] * dims_[1] * dims_[2])
                                + ((m - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                                + ((n - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n, size_t o) const {
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in DViewFMatrixKokkos 7D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DViewFMatrixKokkos 7D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DViewFMatrixKokkos 7D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DViewFMatrixKokkos 7D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in DViewFMatrixKokkos 7D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in DViewFMatrixKokkos 7D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in DViewFMatrixKokkos 7D!");
    assert(o >= 1 && o <= dims_[6] && "o is out of bounds in DViewFMatrixKokkos 7D!");
    return this_matrix_((i - 1) + ((j - 1) * dims_[0])
                                + ((k - 1) * dims_[0] * dims_[1])
                                + ((l - 1) * dims_[0] * dims_[1] * dims_[2])
                                + ((m - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                                + ((n - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4])
                                + ((o - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4] * dims_[5]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>& DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator= (const DViewFMatrixKokkos& temp) {
    
    // Do nothing if the assignment is of the form x = x
    if (this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for

        order_ = temp.order_;
        length_ = temp.length_;
        temp_inp_matrix_ = temp.temp_inp_matrix_;
        this_matrix_host_ = temp.this_matrix_host_;
        this_matrix_ = temp.this_matrix_;
	host = temp.host;
    }
    
    return *this;
}

// Return size
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::size() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::extent() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::dims(size_t i) const {
    i--;
    assert(i < order_ && "DViewFMatrixKokkos order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to DViewFMatrixKokkos dims is out of bounds!");
    return dims_[i];
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::order() const {
    return order_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::device_pointer() const {
    return this_matrix_.data();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::host_pointer() const {
    return this_matrix_host_.data();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
void DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::update_host() {
    // Deep copy of device view to host view
    deep_copy(this_matrix_host_, this_matrix_);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
void DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::update_device() {
    // Deep copy of host view to device view
    deep_copy(this_matrix_, this_matrix_host_);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::~DViewFMatrixKokkos() {}
// End DViewFMatrixKokkos


/*! \brief Kokkos version of the serial CArray class.
 *
 */
template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace, typename MemoryTraits = void>
class CArrayKokkos {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace, MemoryTraits>;
    
private:
    size_t dims_[7];
    size_t order_;
    size_t length_;
    TArray1D this_array_;

public:
    CArrayKokkos();
    
    CArrayKokkos(size_t dim0, const std::string& tag_string = DEFAULTSTRINGARRAY);

    CArrayKokkos(size_t dim0, size_t dim1, const std::string& tag_string = DEFAULTSTRINGARRAY);

    CArrayKokkos (size_t dim0, size_t dim1, size_t dim2, const std::string& tag_string = DEFAULTSTRINGARRAY);

    CArrayKokkos(size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, const std::string& tag_string = DEFAULTSTRINGARRAY);

    CArrayKokkos(size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, const std::string& tag_string = DEFAULTSTRINGARRAY);

    CArrayKokkos(size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, size_t dim5, const std::string& tag_string = DEFAULTSTRINGARRAY);

    CArrayKokkos(size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, size_t dim5,
                 size_t dim6, const std::string& tag_string = DEFAULTSTRINGARRAY);
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n, size_t o) const;
    
    KOKKOS_INLINE_FUNCTION
    CArrayKokkos& operator=(const CArrayKokkos& temp);

    // GPU Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t size() const;

    // Host Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t extent() const;

    KOKKOS_INLINE_FUNCTION
    size_t dims(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    size_t order() const;
 
    // Methods returns the raw pointer (most likely GPU) of the Kokkos View
    KOKKOS_INLINE_FUNCTION
    T* pointer() const;
    
    //return the view
    KOKKOS_INLINE_FUNCTION
    TArray1D get_kokkos_view() const;

    // Deconstructor
    KOKKOS_INLINE_FUNCTION
    ~CArrayKokkos ();
}; // End of CArrayKokkos

// Default constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::CArrayKokkos() {
    length_ = order_ = 0;
    for (int i = 0; i < 7; i++) {
        dims_[i] = 0;
    }
}

// Overloaded 1D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::CArrayKokkos(size_t dim0, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim0;
    order_ = 1;
    length_ = dim0;
    this_array_ = TArray1D(tag_string, length_);
}

// Overloaded 2D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::CArrayKokkos(size_t dim0, size_t dim1, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    order_ = 2;
    length_ = (dim0 * dim1);
    this_array_ = TArray1D(tag_string, length_);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::CArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    order_ = 3;
    length_ = (dim0 * dim1 * dim2);
    this_array_ = TArray1D(tag_string, length_);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::CArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    order_ = 4;
    length_ = (dim0 * dim1 * dim2 * dim3);
    this_array_ = TArray1D(tag_string, length_);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::CArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3,
                              size_t dim4, const std::string& tag_string) {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    order_ = 5;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4);
    this_array_ = TArray1D(tag_string, length_);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::CArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3,
                              size_t dim4, size_t dim5, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    order_ = 6;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4 * dim5);
    this_array_ = TArray1D(tag_string, length_);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::CArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3,
                              size_t dim4, size_t dim5,
                              size_t dim6, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    dims_[6] = dim6;
    order_ = 7;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4 * dim5 * dim6);
    this_array_ = TArray1D(tag_string, length_);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i) const {
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in CArrayKokkos 1D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in CArrayKokkos 1D!");
    return this_array_(i);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j) const {
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in CArrayKokkos 2D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in CArrayKokkos 2D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in CArrayKokkos 2D!");
    return this_array_(j + (i * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k) const {
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in CArrayKokkos 3D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in CArrayKokkos 3D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in CArrayKokkos 3D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in CArrayKokkos 3D!");
    return this_array_(k + (j * dims_[2])
                         + (i * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in CArrayKokkos 4D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in CArrayKokkos 4D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in CArrayKokkos 4D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in CArrayKokkos 4D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in CArrayKokkos 4D!");
    return this_array_(l + (k * dims_[3])
                         + (j * dims_[3] * dims_[2])
                         + (i * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m) const {
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in CArrayKokkos 5D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in CArrayKokkos 5D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in CArrayKokkos 5D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in CArrayKokkos 5D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in CArrayKokkos 5D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in CArrayKokkos 5D!");
    return this_array_(m + (l * dims_[4])
                         + (k * dims_[4] * dims_[3])
                         + (j * dims_[4] * dims_[3] * dims_[2])
                         + (i * dims_[4] * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n) const {
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in CArrayKokkos 6D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in CArrayKokkos 6D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in CArrayKokkos 6D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in CArrayKokkos 6D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in CArrayKokkos 6D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in CArrayKokkos 6D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in CArrayKokkos 6D!");
    return this_array_(n + (m * dims_[5])
                         + (l * dims_[5] * dims_[4])
                         + (k * dims_[5] * dims_[4] * dims_[3])
                         + (j * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                         + (i * dims_[5] * dims_[4] * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n, size_t o) const {
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in CArrayKokkos 7D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in CArrayKokkos 7D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in CArrayKokkos 7D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in CArrayKokkos 7D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in CArrayKokkos 7D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in CArrayKokkos 7D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in CArrayKokkos 7D!");
    assert(o >= 0 && o < dims_[6] && "o is out of bounds in CArrayKokkos 7D!");
    return this_array_(o + (n * dims_[6])
                         + (m * dims_[6] * dims_[5])
                         + (l * dims_[6] * dims_[5] * dims_[4])
                         + (k * dims_[6] * dims_[5] * dims_[4] * dims_[3])
                         + (j * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                         + (i * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>& CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator= (const CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>& temp) {
    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    // Do nothing if the assignment is of the form x = x
    if (this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for

        order_ = temp.order_;
        length_ = temp.length_;
        this_array_ = temp.this_array_;
    }
    
    return *this;
}

// Return size
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::size() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::extent() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::dims(size_t i) const {
    assert(i < order_ && "CArrayKokkos order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to CArrayKokkos dims is out of bounds!");
    return dims_[i];
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::order() const {
    return order_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::pointer() const {
    return this_array_.data();
}

//return the stored Kokkos view
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
Kokkos::View<T*, Layout, ExecSpace, MemoryTraits> CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::get_kokkos_view() const {
    return this_array_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::~CArrayKokkos() {}

////////////////////////////////////////////////////////////////////////////////
// End of CArrayKokkos
////////////////////////////////////////////////////////////////////////////////

/*! \brief Kokkos version of the serial ViewCArray class.
 *
 */
template <typename T>
class ViewCArrayKokkos {

private:
    size_t dims_[7];
    size_t order_;
    size_t length_;  // Length of 1D array
    T* this_array_;
    
public:
    KOKKOS_INLINE_FUNCTION
    ViewCArrayKokkos();
    
    KOKKOS_INLINE_FUNCTION
    ViewCArrayKokkos(T* some_array, size_t dim0);

    KOKKOS_INLINE_FUNCTION
    ViewCArrayKokkos(T* some_array, size_t dim0, size_t dim1);

    KOKKOS_INLINE_FUNCTION
    ViewCArrayKokkos(T* some_array, size_t dim0, size_t dim1,
                     size_t dim2);

    KOKKOS_INLINE_FUNCTION
    ViewCArrayKokkos(T* some_array, size_t dim0, size_t dim1,
                     size_t dim2, size_t dim3);

    KOKKOS_INLINE_FUNCTION
    ViewCArrayKokkos(T* some_array, size_t dim0, size_t dim1,
                     size_t dim2, size_t dim3, size_t dim4);

    KOKKOS_INLINE_FUNCTION
    ViewCArrayKokkos(T* some_array, size_t dim0, size_t dim1,
                     size_t dim2, size_t dim3, size_t dim4,
                     size_t dim5);
    
    KOKKOS_INLINE_FUNCTION
    ViewCArrayKokkos(T* some_array, size_t dim0, size_t dim1,
                     size_t dim2, size_t dim3, size_t dim4,
                     size_t dim5, size_t dim6);;
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i) const;
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;
        
    KOKKOS_INLINE_FUNCTION
    T& operator() (size_t i, size_t j, size_t k, size_t l, size_t m) const;
        
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n) const;
 
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n, size_t o) const;
    
    KOKKOS_INLINE_FUNCTION
    size_t size() const;
    
    KOKKOS_INLINE_FUNCTION
    size_t extent() const;

    KOKKOS_INLINE_FUNCTION
    size_t dims(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    size_t order() const;
 
    KOKKOS_INLINE_FUNCTION
    T* pointer() const;

    KOKKOS_INLINE_FUNCTION
    ~ViewCArrayKokkos();
    
}; // end of ViewCArrayKokkos

// Default constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewCArrayKokkos<T>::ViewCArrayKokkos() {
    length_ = order_ = 0;
    this_array_ = NULL;
    for (int i = 0; i < 7; i++) {
        dims_[i] = 0;
    }
}

// Overloaded 1D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewCArrayKokkos<T>::ViewCArrayKokkos(T* some_array, size_t dim0) {
    dims_[0] = dim0;
    order_ = 1;
    length_ = dim0;
    this_array_ = some_array;
}

// Overloaded 2D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewCArrayKokkos<T>::ViewCArrayKokkos(T* some_array, size_t dim0,
                                      size_t dim1) {
    dims_[0] = dim0;
    dims_[1] = dim1;
    order_ = 2;
    length_ = (dim0 * dim1);
    this_array_ = some_array;
}

// Overloaded 3D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewCArrayKokkos<T>::ViewCArrayKokkos(T* some_array, size_t dim0,
                                      size_t dim1, size_t dim2) {
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    order_ = 3;
    length_ = (dim0 * dim1 * dim2);
    this_array_ = some_array;
}

// Overloaded 4D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewCArrayKokkos<T>::ViewCArrayKokkos(T* some_array, size_t dim0,
                                      size_t dim1, size_t dim2,
                                      size_t dim3) {
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    order_ = 4;
    length_ = (dim0 * dim1 * dim2 * dim3);
    this_array_ = some_array;
}

// Overloaded 5D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewCArrayKokkos<T>::ViewCArrayKokkos(T* some_array, size_t dim0,
                                      size_t dim1, size_t dim2,
                                      size_t dim3, size_t dim4) {
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    order_ = 5;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4);
    this_array_ = some_array;
}

// Overloaded 6D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewCArrayKokkos<T>::ViewCArrayKokkos(T* some_array, size_t dim0,
                                      size_t dim1, size_t dim2,
                                      size_t dim3, size_t dim4,
                                      size_t dim5) {
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    order_ = 6;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4 * dim5);
    this_array_ = some_array;
}

// Overloaded 7D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewCArrayKokkos<T>::ViewCArrayKokkos(T* some_array, size_t dim0,
                                      size_t dim1, size_t dim2,
                                      size_t dim3, size_t dim4,
                                      size_t dim5, size_t dim6) {
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    dims_[6] = dim6;
    order_ = 7;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4 * dim5 * dim6);
    this_array_ = some_array;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i) const {
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in ViewCArrayKokkos 1D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewCArrayKokkos 1D!");
    return this_array_[i];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i, size_t j) const {
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in ViewCArrayKokkos 2D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewCArrayKokkos 2D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewCArrayKokkos 2D!");
    return this_array_[j + (i * dims_[1])];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i, size_t j, size_t k) const {
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in ViewCArrayKokkos 3D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewCArrayKokkos 3D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewCArrayKokkos 3D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewCArrayKokkos 3D!");
    return this_array_[k + (j * dims_[2])
                         + (i * dims_[2] * dims_[1])];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i, size_t j, size_t k,
                                   size_t l) const {
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in ViewCArrayKokkos 4D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewCArrayKokkos 4D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewCArrayKokkos 4D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewCArrayKokkos 4D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in ViewCArrayKokkos 4D!");
    return this_array_[l + (k * dims_[3])
                         + (j * dims_[3] * dims_[2])
                         + (i * dims_[3] * dims_[2] * dims_[1])];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l,
                                   size_t m) const {
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in ViewCArrayKokkos 5D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewCArrayKokkos 5D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewCArrayKokkos 5D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewCArrayKokkos 5D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in ViewCArrayKokkos 5D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in ViewCArrayKokkos 5D!");
    return this_array_[m + (l * dims_[4])
                         + (k * dims_[4] * dims_[3])
                         + (j * dims_[4] * dims_[3] * dims_[2])
                         + (i * dims_[4] * dims_[3] * dims_[2] * dims_[1])];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l,
                                   size_t m, size_t n) const {
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in ViewCArrayKokkos 6D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewCArrayKokkos 6D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewCArrayKokkos 6D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewCArrayKokkos 6D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in ViewCArrayKokkos 6D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in ViewCArrayKokkos 6D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in ViewCArrayKokkos 6D!");
    return this_array_[n + (m * dims_[5])
                         + (l * dims_[5] * dims_[4])
                         + (k * dims_[5] * dims_[4] * dims_[3])
                         + (j * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                         + (i * dims_[5] * dims_[4] * dims_[3] * dims_[2] * dims_[1])];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l,
                                   size_t m, size_t n, size_t o) const {
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in ViewCArrayKokkos 7D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewCArrayKokkos 7D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewCArrayKokkos 7D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewCArrayKokkos 7D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in ViewCArrayKokkos 7D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in ViewCArrayKokkos 7D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in ViewCArrayKokkos 7D!");
    assert(o >= 0 && o < dims_[6] && "o is out of bounds in ViewCArrayKokkos 7D!");
    return this_array_[o + (n * dims_[6])
                         + (m * dims_[6] * dims_[5])
                         + (l * dims_[6] * dims_[5] * dims_[4])
                         + (k * dims_[6] * dims_[5] * dims_[4] * dims_[3])
                         + (j * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                         + (i * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2] * dims_[1])];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
size_t ViewCArrayKokkos<T>::size() const {
    return length_;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
size_t ViewCArrayKokkos<T>::extent() const {
    return length_;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
size_t ViewCArrayKokkos<T>::dims(size_t i) const {
    assert(i < order_ && "ViewCArrayKokkos order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to ViewCArrayKokkos dims is out of bounds!");
    return dims_[i];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
size_t ViewCArrayKokkos<T>::order() const {
    return order_;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T* ViewCArrayKokkos<T>::pointer() const {
    return this_array_;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
ViewCArrayKokkos<T>::~ViewCArrayKokkos() {}

////////////////////////////////////////////////////////////////////////////////
// End of ViewCArrayKokkos
////////////////////////////////////////////////////////////////////////////////

/*! \brief Kokkos version of the serial CMatrix class.
 *
 */
template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace, typename MemoryTraits = void>
class CMatrixKokkos {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace, MemoryTraits>;
    
private:
    size_t dims_[7];
    size_t order_;
    size_t length_;
    TArray1D this_matrix_;

public:
    CMatrixKokkos();

    CMatrixKokkos(size_t dim1, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    CMatrixKokkos(size_t dim1, size_t dim2, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    CMatrixKokkos(size_t dim1, size_t dim2, size_t dim3, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    CMatrixKokkos(size_t dim1, size_t dim2, size_t dim3,
                  size_t dim4, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    CMatrixKokkos(size_t dim1, size_t dim2, size_t dim3,
                  size_t dim4, size_t dim5, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    CMatrixKokkos(size_t dim1, size_t dim2, size_t dim3,
                  size_t dim4, size_t dim5, size_t dim6, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    CMatrixKokkos(size_t dim1, size_t dim2, size_t dim3,
                  size_t dim4, size_t dim5, size_t dim6,
                  size_t dim7, const std::string& tag_string = DEFAULTSTRINGMATRIX);
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n, size_t o) const;
    
    KOKKOS_INLINE_FUNCTION
    CMatrixKokkos& operator=(const CMatrixKokkos &temp);

    KOKKOS_INLINE_FUNCTION
    size_t size() const;
    
    KOKKOS_INLINE_FUNCTION
    size_t extent() const;

    KOKKOS_INLINE_FUNCTION
    size_t dims(size_t i) const;
 
    KOKKOS_INLINE_FUNCTION
    size_t order() const;
   
    KOKKOS_INLINE_FUNCTION
    T* pointer() const;

    //return the view
    KOKKOS_INLINE_FUNCTION
    TArray1D get_kokkos_view() const;

    KOKKOS_INLINE_FUNCTION
    ~CMatrixKokkos();

}; // End of CMatrixKokkos

// Default constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::CMatrixKokkos() {
    length_ = order_ = 0;
    for (int i = 0; i < 7; i++) {
        dims_[i] = 0;
    }
}

// Overloaded 1D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::CMatrixKokkos(size_t dim1, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim1;
    order_ = 1;
    length_ = dim1;
    this_matrix_ = TArray1D(tag_string, length_);
}

// Overloaded 2D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::CMatrixKokkos(size_t dim1, size_t dim2, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    order_ = 2;
    length_ = (dim1 * dim2);
    this_matrix_ = TArray1D(tag_string, length_);
}

// Overloaded 3D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::CMatrixKokkos(size_t dim1, size_t dim2,
                                size_t dim3, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    order_ = 3;
    length_ = (dim1 * dim2 * dim3);
    this_matrix_ = TArray1D(tag_string, length_);
}

// Overloaded 4D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::CMatrixKokkos(size_t dim1, size_t dim2,
                                size_t dim3, size_t dim4, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    order_ = 4;
    length_ = (dim1 * dim2 * dim3 * dim4);
    this_matrix_ = TArray1D(tag_string, length_);
}

// Overloaded 5D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::CMatrixKokkos(size_t dim1, size_t dim2,
                                size_t dim3, size_t dim4,
                                size_t dim5, const std::string& tag_string) {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    order_ = 5;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5);
    this_matrix_ = TArray1D(tag_string, length_);
}

// Overloaded 6D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::CMatrixKokkos(size_t dim1, size_t dim2,
                                size_t dim3, size_t dim4,
                                size_t dim5, size_t dim6, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    order_ = 6;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5 * dim6);
    this_matrix_ = TArray1D(tag_string, length_);
}

// Overloaded 7D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::CMatrixKokkos(size_t dim1, size_t dim2,
                                size_t dim3, size_t dim4,
                                size_t dim5, size_t dim6,
                                size_t dim7, const std::string& tag_string) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    dims_[6] = dim7;
    order_ = 7;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5 * dim6 * dim7);
    this_matrix_ = TArray1D(tag_string, length_);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i) const {
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in CMatrixKokkos 1D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in CMatrixKokkos 1D!");
    return this_matrix_((i - 1));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j) const {
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in CMatrixKokkos 2D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in CMatrixKokkos 2D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in CMatrixKokkos 2D!");
    return this_matrix_((j - 1) + ((i - 1) * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k) const {
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in CMatrixKokkos 3D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in CMatrixKokkos 3D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in CMatrixKokkos 3D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in CMatrixKokkos 3D!");
    return this_matrix_((k - 1) + ((j - 1) * dims_[2])
                                + ((i - 1) * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in CMatrixKokkos 4D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in CMatrixKokkos 4D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in CMatrixKokkos 4D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in CMatrixKokkos 4D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in CMatrixKokkos 4D!");
    return this_matrix_((l - 1) + ((k - 1) * dims_[3])
                                + ((j - 1) * dims_[3] * dims_[2])
                                + ((i - 1) * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                                size_t m) const {
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in CMatrixKokkos 5D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in CMatrixKokkos 5D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in CMatrixKokkos 5D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in CMatrixKokkos 5D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in CMatrixKokkos 5D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in CMatrixKokkos 5D!");
    return this_matrix_((m - 1) + ((l - 1) * dims_[4])
                                + ((k - 1) * dims_[4] * dims_[3])
                                + ((j - 1) * dims_[4] * dims_[3] * dims_[2])
                                + ((i - 1) * dims_[4] * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                                size_t m, size_t n) const {
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in CMatrixKokkos 6D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in CMatrixKokkos 6D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in CMatrixKokkos 6D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in CMatrixKokkos 6D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in CMatrixKokkos 6D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in CMatrixKokkos 6D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in CMatrixKokkos 6D!");
    return this_matrix_((n - 1) + ((m - 1) * dims_[5])
                                + ((l - 1) * dims_[5] * dims_[4])
                                + ((k - 1) * dims_[5] * dims_[4] * dims_[3])
                                + ((j - 1) * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                                + ((i - 1) * dims_[5] * dims_[4] * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                                size_t m, size_t n, size_t o) const {
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in CMatrixKokkos 7D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in CMatrixKokkos 7D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in CMatrixKokkos 7D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in CMatrixKokkos 7D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in CMatrixKokkos 7D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in CMatrixKokkos 7D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in CMatrixKokkos 7D!");
    assert(o >= 1 && o <= dims_[6] && "o is out of bounds in CMatrixKokkos 7D!");
    return this_matrix_((o-1) + ((n - 1) * dims_[6])
                              + ((m - 1) * dims_[6] * dims_[5])
                              + ((l - 1) * dims_[6] * dims_[5] * dims_[4])
                              + ((k - 1) * dims_[6] * dims_[5] * dims_[4] * dims_[3])
                              + ((j - 1) * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                              + ((i - 1) * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2] * dims_[1]));
}

// Overload = operator
// for object assignment THIS = CMatrixKokkos <> temp
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits> & CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator=(const CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits> &temp) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;

    if( this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for

        order_ = temp.order_;
        length_ = temp.length_;
        this_matrix_ = temp.this_matrix_;
    }
    
    return *this;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::size() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::extent() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::dims(size_t i) const {
    i--;
    assert(i < order_ && "CMatrixKokkos order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to CMatrixKokkos dims is out of bounds!");
    return dims_[i];
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::order() const {
    return order_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::pointer() const {
    return this_matrix_.data();
}

//return the stored Kokkos view
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
Kokkos::View<T*, Layout, ExecSpace, MemoryTraits> CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::get_kokkos_view() const {
    return this_matrix_;
}

// Deconstructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::~CMatrixKokkos() {}

////////////////////////////////////////////////////////////////////////////////
// End of CMatrixKokkos
////////////////////////////////////////////////////////////////////////////////

/*! \brief Kokkos version of the serial ViewCMatrix class.
 *
 */
template <typename T>
class ViewCMatrixKokkos {

private:
    size_t dims_[7];
    size_t order_;
    size_t length_;
    T* this_matrix_;

public:
    KOKKOS_INLINE_FUNCTION
    ViewCMatrixKokkos();

    KOKKOS_INLINE_FUNCTION
    ViewCMatrixKokkos(T* some_matrix, size_t dim1);

    KOKKOS_INLINE_FUNCTION
    ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2);

    KOKKOS_INLINE_FUNCTION
    ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2, size_t dim3);

    KOKKOS_INLINE_FUNCTION
    ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2, size_t dim3,
                      size_t dim4);

    KOKKOS_INLINE_FUNCTION
    ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2, size_t dim3,
                      size_t dim4, size_t dim5);

    KOKKOS_INLINE_FUNCTION
    ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2, size_t dim3,
                      size_t dim4, size_t dim5, size_t dim6);

    KOKKOS_INLINE_FUNCTION
    ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2, size_t dim3,
                      size_t dim4, size_t dim5, size_t dim6, size_t dim7);
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j , size_t k) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k , size_t l) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t o) const;

    KOKKOS_INLINE_FUNCTION
    size_t size() const;

    KOKKOS_INLINE_FUNCTION
    size_t extent() const;

    KOKKOS_INLINE_FUNCTION
    size_t dims(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    size_t order() const;
 
    KOKKOS_INLINE_FUNCTION
    T* pointer() const;

    KOKKOS_INLINE_FUNCTION
    ~ViewCMatrixKokkos();

}; // End of ViewCMatrixKokkos

// Default constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(){
    length_ = order_ = 0;
    this_matrix_ = NULL;
    for (int i = 0; i < 7; i++) {
        dims_[i] = 0;
    }
}

// Overloaded 1D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T* some_matrix, size_t dim1) {
    dims_[0] = dim1;
    order_ = 1;
    length_ = dim1;
    this_matrix_ = some_matrix;
}

// Overloaded 2D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T* some_matrix, size_t dim1,
                                        size_t dim2) {
    dims_[0] = dim1;
    dims_[1] = dim2;
    order_ = 2;
    length_ = (dim1 * dim2);
    this_matrix_ = some_matrix;
}

// Overloaded 3D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2,
                                        size_t dim3) {
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    order_ = 3;
    length_ = (dim1 * dim2 * dim3);
    this_matrix_ = some_matrix;
}

// Overloaded 4D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2,
                                        size_t dim3, size_t dim4) {
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    order_ = 4;
    length_ = (dim1 * dim2 * dim3 * dim4);
    this_matrix_ = some_matrix;
}

// Overloaded 5D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2,
                                        size_t dim3, size_t dim4, size_t dim5) {
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    order_ = 5;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5);
    this_matrix_ = some_matrix;
}

// Overloaded 6D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2,
                                        size_t dim3, size_t dim4, size_t dim5,
                                        size_t dim6) {
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    order_ = 6;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5 * dim6);
    this_matrix_ = some_matrix;
}

// Overloaded 7D constructor
template <typename T>
KOKKOS_INLINE_FUNCTION
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2,
                                        size_t dim3, size_t dim4, size_t dim5,
                                        size_t dim6, size_t dim7) {
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    dims_[6] = dim7;
    order_ = 7;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5 * dim6 * dim7);
    this_matrix_ = some_matrix;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i) const {
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in ViewCMatrixKokkos 1D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewCMatrixKokkos 1D!");
    return this_matrix_[(i - 1)];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i, size_t j) const {
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in ViewCMatrixKokkos 2D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewCMatrixKokkos 2D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewCMatrixKokkos 2D!");
    return this_matrix_[(j - 1) + ((i - 1) * dims_[1])];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k) const {
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in ViewCMatrixKokkos 3D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewCMatrixKokkos 3D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewCMatrixKokkos 3D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in ViewCMatrixKokkos 3D!");
    return this_matrix_[(k - 1) + ((j - 1) * dims_[2])
                                + ((i - 1) * dims_[2] * dims_[1])];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i, size_t j , size_t k, size_t l) const {
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in ViewCMatrixKokkos 4D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewCMatrixKokkos 4D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewCMatrixKokkos 4D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in ViewCMatrixKokkos 4D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in ViewCMatrixKokkos 4D!");
    return this_matrix_[(l - 1) + ((k - 1) * dims_[3])
                                + ((j - 1) * dims_[3] * dims_[2])
                                + ((i - 1) * dims_[3] * dims_[2] * dims_[1])];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l,
                                    size_t m) const {
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in ViewCMatrixKokkos 5D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds for ViewCMatrixKokkos 5D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds for ViewCMatrixKokkos 5D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds for ViewCMatrixKokkos 5D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds for ViewCMatrixKokkos 5D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds for ViewCMatrixKokkos 5D!");
    return this_matrix_[(m - 1) + ((l - 1) * dims_[4])
                                + ((k - 1) * dims_[4] * dims_[3])
                                + ((j - 1) * dims_[4] * dims_[3] * dims_[2])
                                + ((i - 1) * dims_[4] * dims_[3] * dims_[2] * dims_[1])];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l,
                                    size_t m, size_t n) const {
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in ViewCMatrixKokkos 6D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds for ViewCMatrixKokkos 6D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds for ViewCMatrixKokkos 6D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds for ViewCMatrixKokkos 6D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds for ViewCMatrixKokkos 6D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds for ViewCMatrixKokkos 6D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds for ViewCMatrixKokkos 6D!");
    return this_matrix_[(n - 1) + ((m - 1) * dims_[5])
                                + ((l - 1) * dims_[5] * dims_[4])
                                + ((k - 1) * dims_[5] * dims_[4] * dims_[3])
                                + ((j - 1) * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                                + ((i - 1) * dims_[5] * dims_[4] * dims_[3] * dims_[2] * dims_[1])];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l,
                                    size_t m, size_t n, size_t o) const {
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in ViewCMatrixKokkos 7D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds for ViewCMatrixKokkos 7D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds for ViewCMatrixKokkos 7D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds for ViewCMatrixKokkos 7D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds for ViewCMatrixKokkos 7D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds for ViewCMatrixKokkos 7D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds for ViewCMatrixKokkos 7D!");
    assert(o >= 1 && o <= dims_[6] && "o is out of bounds for ViewCMatrixKokkos 7D!");
    return this_matrix_[o + ((n - 1) * dims_[6])
                          + ((m - 1) * dims_[6] * dims_[5])
                          + ((l - 1) * dims_[6] * dims_[5] * dims_[4])
                          + ((k - 1) * dims_[6] * dims_[5] * dims_[4] * dims_[3])
                          + ((j - 1) * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                          + ((i - 1) * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2] * dims_[1])];
}


template <typename T>
KOKKOS_INLINE_FUNCTION
size_t ViewCMatrixKokkos<T>::size() const {
    return length_;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
size_t ViewCMatrixKokkos<T>::extent() const {
    return length_;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
size_t ViewCMatrixKokkos<T>::dims(size_t i) const {
    i--;
    assert(i < order_ && "ViewCMatrixKokkos order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to ViewCMatrixKokkos dims is out of bounds!");
    return dims_[i];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
size_t ViewCMatrixKokkos<T>::order() const {
    return order_;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T* ViewCMatrixKokkos<T>::pointer() const {
    return this_matrix_;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
ViewCMatrixKokkos<T>::~ViewCMatrixKokkos() {}

////////////////////////////////////////////////////////////////////////////////
// End of ViewCMatrixKokkos
////////////////////////////////////////////////////////////////////////////////

/////////////////////////
// DCArrayKokkos:  Dual type for managing data on both CPU and GPU.
/////////////////////////
template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace, typename MemoryTraits = void>
class DCArrayKokkos {

    // this is manage
    using TArray1D = Kokkos::DualView <T*, Layout, ExecSpace, MemoryTraits>;
    
private:
    size_t dims_[7];
    size_t length_;
    size_t order_;  // tensor order (rank)
    TArray1D this_array_;

public:
    // Data member to access host view
    ViewCArray <T> host;

    DCArrayKokkos();
    
    DCArrayKokkos(size_t dim0, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DCArrayKokkos(size_t dim0, size_t dim1, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DCArrayKokkos (size_t dim0, size_t dim1, size_t dim2, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DCArrayKokkos(size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DCArrayKokkos(size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DCArrayKokkos(size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, size_t dim5, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DCArrayKokkos(size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, size_t dim5,
                 size_t dim6, const std::string& tag_string = DEFAULTSTRINGARRAY);
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n, size_t o) const;
    
    KOKKOS_INLINE_FUNCTION
    DCArrayKokkos& operator=(const DCArrayKokkos& temp);

    // GPU Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t size() const;

    // Host Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t extent() const;

    KOKKOS_INLINE_FUNCTION
    size_t dims(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    size_t order() const;
 
    // Method returns the raw device pointer of the Kokkos DualView
    KOKKOS_INLINE_FUNCTION
    T* device_pointer() const;

    // Method returns the raw host pointer of the Kokkos DualView
    KOKKOS_INLINE_FUNCTION
    T* host_pointer() const;

    // Method returns kokkos dual view
    KOKKOS_INLINE_FUNCTION
    TArray1D get_kokkos_dual_view() const;

    // Method that update host view
    void update_host();

    // Method that update device view
    void update_device();

    // Deconstructor
    KOKKOS_INLINE_FUNCTION
    ~DCArrayKokkos ();
}; // End of DCArrayKokkos


// Default constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DCArrayKokkos() {
    length_ = order_ = 0;
    for (int i = 0; i < 7; i++) {
        dims_[i] = 0;
    }
}

// Overloaded 1D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DCArrayKokkos(size_t dim0, const std::string& tag_string) {
    
    dims_[0] = dim0;
    order_ = 1;
    length_ = dim0;
    this_array_ = TArray1D(tag_string, length_);
    // Create host ViewCArray
    host = ViewCArray <T> (this_array_.h_view.data(), dim0);
}

// Overloaded 2D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DCArrayKokkos(size_t dim0, size_t dim1, const std::string& tag_string) {
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    order_ = 2;
    length_ = (dim0 * dim1);
    this_array_ = TArray1D(tag_string, length_);
    // Create host ViewCArray
    host = ViewCArray <T> (this_array_.h_view.data(), dim0, dim1);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DCArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, const std::string& tag_string) {
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    order_ = 3;
    length_ = (dim0 * dim1 * dim2);
    this_array_ = TArray1D(tag_string, length_);
    // Create host ViewCArray
    host = ViewCArray <T> (this_array_.h_view.data(), dim0, dim1, dim2);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DCArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3, const std::string& tag_string) {
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    order_ = 4;
    length_ = (dim0 * dim1 * dim2 * dim3);
    this_array_ = TArray1D(tag_string, length_);
    // Create host ViewCArray
    host = ViewCArray <T> (this_array_.h_view.data(), dim0, dim1, dim2, dim3);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DCArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3,
                              size_t dim4, const std::string& tag_string) {
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    order_ = 5;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4);
    this_array_ = TArray1D(tag_string, length_);
    // Create host ViewCArray
    host = ViewCArray <T> (this_array_.h_view.data(), dim0, dim1, dim2, dim3, dim4);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DCArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3,
                              size_t dim4, size_t dim5, const std::string& tag_string) {
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    order_ = 6;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4 * dim5);
    this_array_ = TArray1D(tag_string, length_);
    // Create host ViewCArray
    host = ViewCArray <T> (this_array_.h_view.data(), dim0, dim1, dim2, dim3, dim4, dim5);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DCArrayKokkos(size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3,
                              size_t dim4, size_t dim5,
                              size_t dim6, const std::string& tag_string) {
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    dims_[6] = dim6;
    order_ = 7;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4 * dim5 * dim6);
    this_array_ = TArray1D(tag_string, length_);
    // Create host ViewCArray
    host = ViewCArray <T> (this_array_.h_view.data(), dim0, dim1, dim2, dim3, dim4, dim5, dim6);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i) const {
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in DCArrayKokkos 1D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DCArrayKokkos 1D!");
    return this_array_.d_view(i);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j) const {
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in DCArrayKokkos 2D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DCArrayKokkos 2D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DCArrayKokkos 2D!");
    return this_array_.d_view(j + (i * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k) const {
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in DCArrayKokkos 3D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DCArrayKokkos 3D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DCArrayKokkos 3D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DCArrayKokkos 3D!");
    return this_array_.d_view(k + (j * dims_[2])
                                + (i * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in DCArrayKokkos 4D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DCArrayKokkos 4D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DCArrayKokkos 4D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DCArrayKokkos 4D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in DCArrayKokkos 4D!");
    return this_array_.d_view(l + (k * dims_[3])
                                + (j * dims_[3] * dims_[2])
                                + (i * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m) const {
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in DCArrayKokkos 5D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DCArrayKokkos 5D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DCArrayKokkos 5D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DCArrayKokkos 5D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in DCArrayKokkos 5D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in DCArrayKokkos 5D!");
    return this_array_.d_view(m + (l * dims_[4])
                                + (k * dims_[4] * dims_[3])
                                + (j * dims_[4] * dims_[3] * dims_[2])
                                + (i * dims_[4] * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n) const {
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in DCArrayKokkos 6D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DCArrayKokkos 6D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DCArrayKokkos 6D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DCArrayKokkos 6D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in DCArrayKokkos 6D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in DCArrayKokkos 6D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in DCArrayKokkos 6D!");
    return this_array_.d_view(n + (m * dims_[5])
                                + (l * dims_[5] * dims_[4])
                                + (k * dims_[5] * dims_[4] * dims_[3])
                                + (j * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                                + (i * dims_[5] * dims_[4] * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n, size_t o) const {
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in DCArrayKokkos 7D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DCArrayKokkos 7D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DCArrayKokkos 7D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DCArrayKokkos 7D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in DCArrayKokkos 7D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in DCArrayKokkos 7D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in DCArrayKokkos 7D!");
    assert(o >= 0 && o < dims_[6] && "o is out of bounds in DCArrayKokkos 7D!");
    return this_array_.d_view(o + (n * dims_[6])
                                + (m * dims_[6] * dims_[5])
                                + (l * dims_[6] * dims_[5] * dims_[4])
                                + (k * dims_[6] * dims_[5] * dims_[4] * dims_[3])
                                + (j * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                                + (i * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>& DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator= (const DCArrayKokkos& temp) {
    
    // Do nothing if the assignment is of the form x = x
    if (this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for

        order_ = temp.order_;
        length_ = temp.length_;
        this_array_ = temp.this_array_;
        host = temp.host;
    }
    
    return *this;
}

// Return size
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::size() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::extent() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::dims(size_t i) const {
    assert(i < order_ && "DCArrayKokkos order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to DCArrayKokkos dims is out of bounds!");
    return dims_[i];
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::order() const {
    return order_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::device_pointer() const {
    return this_array_.d_view.data();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::host_pointer() const {
    return this_array_.h_view.data();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
Kokkos::DualView <T*, Layout, ExecSpace, MemoryTraits> DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::get_kokkos_dual_view() const {
  return this_array_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
void DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::update_host() {

    this_array_.template modify<typename TArray1D::execution_space>();
    this_array_.template sync<typename TArray1D::host_mirror_space>();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
void DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::update_device() {

    this_array_.template modify<typename TArray1D::host_mirror_space>();
    this_array_.template sync<typename TArray1D::execution_space>();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::~DCArrayKokkos() {}
// End DCArrayKokkos


/////////////////////////
// DViewCArrayKokkos
/////////////////////////
template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace, typename MemoryTraits = void>
class DViewCArrayKokkos {

    // this is always unmanaged
    using TArray1DHost = Kokkos::View<T*, Layout, HostSpace, MemoryUnmanaged>;
    // this is manage
    using TArray1D     = Kokkos::View<T*, Layout, ExecSpace, MemoryTraits>;
    typename ExecSpace::memory_space memspace;
    
private:
    size_t dims_[7];
    size_t length_;
    size_t order_;  // tensor order (rank)
    TArray1D this_array_;
    TArray1DHost this_array_host_;
    T * temp_inp_array_;
    //typename Kokkos::View<T*, Layout, ExecSpace>::HostMirror  h_this_array_;

public:
    DViewCArrayKokkos();
    
    DViewCArrayKokkos(T * inp_array, size_t dim0, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, size_t dim5, const std::string& tag_string = DEFAULTSTRINGARRAY);

    DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, size_t dim5,
                 size_t dim6, const std::string& tag_string = DEFAULTSTRINGARRAY);
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, 
                  size_t n) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n, size_t o) const;
    
    KOKKOS_INLINE_FUNCTION
    DViewCArrayKokkos& operator=(const DViewCArrayKokkos& temp);

    // GPU Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t size() const;

    // Host Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t extent() const;

    KOKKOS_INLINE_FUNCTION
    size_t dims(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    size_t order() const;
 
    // Method returns the raw device pointer of the Kokkos View
    KOKKOS_INLINE_FUNCTION
    T* device_pointer() const;

    // Method returns the raw host pointer of the Kokkos View
    KOKKOS_INLINE_FUNCTION
    T* host_pointer() const;

    // Data member to access host view
    ViewCArray <T> host;

    // Method that update host view
    void update_host();

    // Method that update device view
    void update_device();

    // Deconstructor
    KOKKOS_INLINE_FUNCTION
    ~DViewCArrayKokkos ();
}; // End of DViewCArrayKokkos


// Default constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCArrayKokkos() {
    length_ = order_ = 0;
    temp_inp_array_ = NULL;
    for (int i = 0; i < 7; i++) {
        dims_[i] = 0;
    }
}

// Overloaded 1D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCArrayKokkos(T * inp_array, size_t dim0, const std::string& tag_string) {
    //using TArray1DHost = Kokkos::View<T*, Layout, HostSpace, MemoryUnmanaged>;
    //using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim0;
    order_ = 1;
    length_ = dim0;
    // Create a 1D host view of the external allocation
    this_array_host_ = TArray1DHost(inp_array, length_);
    // Assign temp point to inp_array pointer that is passed in
    temp_inp_array_ = inp_array;
    // Create a device copy of that host view
    //this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    this_array_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_array_host_);

    // Create host ViewCArray. Note: inp_array and this_array_host_.data() are the same pointer 
    host = ViewCArray <T> (inp_array, dim0);
}

// Overloaded 2D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1, const std::string& tag_string) {
    //using TArray1DHost = Kokkos::View<T*, Layout, HostSpace, MemoryUnmanaged>;
    //using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    //using TArray1Dtemp = TArray1D::HostMirror;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    order_ = 2;
    length_ = (dim0 * dim1);
    // Create a 1D host view of the external allocation
    this_array_host_ = TArray1DHost(inp_array, length_);
    // Assign temp point to inp_array pointer that is passed in
    temp_inp_array_ = inp_array;
    //this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    this_array_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_array_host_);
    // Create host ViewCArray
    host = ViewCArray <T> (inp_array, dim0, dim1);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2, const std::string& tag_string) {
    //using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    order_ = 3;
    length_ = (dim0 * dim1 * dim2);
    // Create a 1D host view of the external allocation
    this_array_host_ = TArray1DHost(inp_array, length_);
    // Assign temp point to inp_array pointer that is passed in
    temp_inp_array_ = inp_array;
    // Create a device copy of that host view
    //this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    this_array_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_array_host_);
    // Create host ViewCArray
    host = ViewCArray <T> (inp_array, dim0, dim1, dim2);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3, const std::string& tag_string) {
    //using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    order_ = 4;
    length_ = (dim0 * dim1 * dim2 * dim3);
    // Create a 1D host view of the external allocation
    this_array_host_ = TArray1DHost(inp_array, length_);
    // Assign temp point to inp_array pointer that is passed in
    temp_inp_array_ = inp_array;
    // Create a device copy of that host view
    //this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    this_array_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_array_host_);
    // Create host ViewCArray
    host = ViewCArray <T> (inp_array, dim0, dim1, dim2, dim3);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3, 
                              size_t dim4, const std::string& tag_string) {

    //using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    order_ = 5;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4);
    // Create a 1D host view of the external allocation
    this_array_host_ = TArray1DHost(inp_array, length_);
    // Assign temp point to inp_array pointer that is passed in
    temp_inp_array_ = inp_array;
    // Create a device copy of that host view
    //this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    this_array_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_array_host_);
    // Create host ViewCArray
    host = ViewCArray <T> (inp_array, dim0, dim1, dim2, dim3, dim4);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3, 
                              size_t dim4, size_t dim5, const std::string& tag_string) {
    //using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    order_ = 6;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4 * dim5);
    // Create a 1D host view of the external allocation
    this_array_host_ = TArray1DHost(inp_array, length_);
    // Assign temp point to inp_array pointer that is passed in
    temp_inp_array_ = inp_array;
    // Create a device copy of that host view
    //this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    this_array_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_array_host_);
    // Create host ViewCArray
    host = ViewCArray <T> (inp_array, dim0, dim1, dim2, dim3, dim4, dim5);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3,
                              size_t dim4, size_t dim5,
                              size_t dim6, const std::string& tag_string) {
    //using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    dims_[6] = dim6;
    order_ = 7;
    length_ = (dim0 * dim1 * dim2 * dim3 * dim4 * dim5 * dim6);
    // Create a 1D host view of the external allocation
    this_array_host_ = TArray1DHost(inp_array, length_);
    // Assign temp point to inp_array pointer that is passed in
    temp_inp_array_ = inp_array;
    //this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    this_array_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_array_host_);
    // Create host ViewCArray
    host = ViewCArray <T> (inp_array, dim0, dim1, dim2, dim3, dim4, dim5, dim6);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i) const {
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in DViewCArrayKokkos 1D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DViewCArrayKokkos 1D!");
    return this_array_(i);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j) const {
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in DViewCArrayKokkos 2D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DViewCArrayKokkos 2D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DViewCArrayKokkos 2D!");
    return this_array_(j + (i * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k) const {
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in DViewCArrayKokkos 3D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DViewCArrayKokkos 3D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DViewCArrayKokkos 3D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DViewCArrayKokkos 3D!");
    return this_array_(k + (j * dims_[2])
                         + (i * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in DViewCArrayKokkos 4D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DViewCArrayKokkos 4D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DViewCArrayKokkos 4D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DViewCArrayKokkos 4D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in DViewCArrayKokkos 4D!");
    return this_array_(l + (k * dims_[3])
                         + (j * dims_[3] * dims_[2])
                         + (i * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m) const {
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in DViewCArrayKokkos 5D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DViewCArrayKokkos 5D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DViewCArrayKokkos 5D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DViewCArrayKokkos 5D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in DViewCArrayKokkos 5D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in DViewCArrayKokkos 5D!");
    return this_array_(m + (l * dims_[4])
                         + (k * dims_[4] * dims_[3])
                         + (j * dims_[4] * dims_[3] * dims_[2])
                         + (i * dims_[4] * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n) const {
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in DViewCArrayKokkos 6D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DViewCArrayKokkos 6D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DViewCArrayKokkos 6D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DViewCArrayKokkos 6D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in DViewCArrayKokkos 6D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in DViewCArrayKokkos 6D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in DViewCArrayKokkos 6D!");
    return this_array_(n + (m * dims_[5])
                         + (l * dims_[5] * dims_[4])
                         + (k * dims_[5] * dims_[4] * dims_[3])
                         + (j * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                         + (i * dims_[5] * dims_[4] * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n, size_t o) const {
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in DViewCArrayKokkos 7D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in DViewCArrayKokkos 7D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in DViewCArrayKokkos 7D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in DViewCArrayKokkos 7D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in DViewCArrayKokkos 7D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in DViewCArrayKokkos 7D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in DViewCArrayKokkos 7D!");
    assert(o >= 0 && o < dims_[6] && "o is out of bounds in DViewCArrayKokkos 7D!");
    return this_array_(o + (n * dims_[6])
                         + (m * dims_[6] * dims_[5])
                         + (l * dims_[6] * dims_[5] * dims_[4])
                         + (k * dims_[6] * dims_[5] * dims_[4] * dims_[3])
                         + (j * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                         + (i * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>& DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator= (const DViewCArrayKokkos& temp) {
    //using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    // Do nothing if the assignment is of the form x = x
    if (this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for

        order_ = temp.order_;
        length_ = temp.length_;
        temp_inp_array_ = temp.temp_inp_array_;
        this_array_host_ = temp.this_array_host_;
        this_array_ = temp.this_array_;
        host = temp.host;
    }
    
    return *this;
}

// Return size
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::size() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::extent() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::dims(size_t i) const {
    assert(i < order_ && "DViewCArrayKokkos order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to DViewCArrayKokkos dims is out of bounds!");
    return dims_[i];
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::order() const {
    return order_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::device_pointer() const {
    return this_array_.data();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::host_pointer() const {
    return this_array_host_.data();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
void DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::update_host() {
    // Deep copy of device view to host view
    deep_copy(this_array_host_, this_array_);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
void DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::update_device() {
    // Deep copy of host view to device view
    deep_copy(this_array_, this_array_host_);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::~DViewCArrayKokkos() {}
// End DViewCArrayKokkos


/////////////////////////
// DCMatrixKokkos
/////////////////////////
template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace, typename MemoryTraits = void>
class DCMatrixKokkos {

    // this is manage
    using TArray1D = Kokkos::DualView<T*, Layout, ExecSpace, MemoryTraits>;
    
private:
    size_t dims_[7];
    size_t length_;
    size_t order_;  // tensor order (rank)
    TArray1D this_matrix_;

public:
    // Data member to access host view
    ViewCMatrix <T> host;

    DCMatrixKokkos();
    
    DCMatrixKokkos(size_t dim1, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DCMatrixKokkos(size_t dim1, size_t dim2, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DCMatrixKokkos (size_t dim1, size_t dim2, size_t dim3, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DCMatrixKokkos(size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DCMatrixKokkos(size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, size_t dim5, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DCMatrixKokkos(size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, size_t dim5, size_t dim6, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DCMatrixKokkos(size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, size_t dim5, size_t dim6,
                 size_t dim7, const std::string& tag_string = DEFAULTSTRINGMATRIX);
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n, size_t o) const;
    
    KOKKOS_INLINE_FUNCTION
    DCMatrixKokkos& operator=(const DCMatrixKokkos& temp);

    // GPU Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t size() const;

    // Host Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t extent() const;

    KOKKOS_INLINE_FUNCTION
    size_t dims(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    size_t order() const;
 
    // Method returns the raw device pointer of the Kokkos DualView
    KOKKOS_INLINE_FUNCTION
    T* device_pointer() const;

    // Method returns the raw host pointer of the Kokkos DualView
    KOKKOS_INLINE_FUNCTION
    T* host_pointer() const;

    // Method that update host view
    void update_host();

    // Method that update device view
    void update_device();

    // Deconstructor
    KOKKOS_INLINE_FUNCTION
    ~DCMatrixKokkos ();

}; // End of DCMatrixKokkos declarations

// Default constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DCMatrixKokkos() {
    length_ = order_ = 0;
    for (int i = 0; i < 7; i++) {
        dims_[i] = 0;
    }
}

// Overloaded 1D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DCMatrixKokkos(size_t dim1, const std::string& tag_string) {
    
    dims_[0] = dim1;
    order_ = 1;
    length_ = dim1;
    this_matrix_ = TArray1D(tag_string, length_);
    // Create host ViewCMatrix
    host = ViewCMatrix <T> (this_matrix_.h_view.data(), dim1);
}

// Overloaded 2D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DCMatrixKokkos(size_t dim1, size_t dim2, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    order_ = 2;
    length_ = (dim1 * dim2);
    this_matrix_ = TArray1D(tag_string, length_);
    // Create host ViewCMatrix
    host = ViewCMatrix <T> (this_matrix_.h_view.data(), dim1, dim2);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DCMatrixKokkos(size_t dim1, size_t dim2,
                              size_t dim3, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    order_ = 3;
    length_ = (dim1 * dim2 * dim3);
    this_matrix_ = TArray1D(tag_string, length_);
    // Create host ViewCMatrix
    host = ViewCMatrix <T> (this_matrix_.h_view.data(), dim1, dim2, dim3);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DCMatrixKokkos(size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    order_ = 4;
    length_ = (dim1 * dim2 * dim3 * dim4);
    this_matrix_ = TArray1D(tag_string, length_);
    // Create host ViewCMatrix
    host = ViewCMatrix <T> (this_matrix_.h_view.data(), dim1, dim2, dim3, dim4);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DCMatrixKokkos(size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4,
                              size_t dim5, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    order_ = 5;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5);
    this_matrix_ = TArray1D(tag_string, length_);
    // Create host ViewCMatrix
    host = ViewCMatrix <T> (this_matrix_.h_view.data(), dim1, dim2, dim3, dim4, dim5);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DCMatrixKokkos(size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4,
                              size_t dim5, size_t dim6, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    order_ = 6;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5 * dim6);
    this_matrix_ = TArray1D(tag_string, length_);
    // Create host ViewCMatrix
    host = ViewCMatrix <T> (this_matrix_.h_view.data(), dim1, dim2, dim3, dim4, dim5, dim6);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DCMatrixKokkos(size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4,
                              size_t dim5, size_t dim6,
                              size_t dim7, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    dims_[6] = dim7;
    order_ = 7;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5 * dim6 * dim7);
    this_matrix_ = TArray1D(tag_string, length_);
    // Create host ViewCMatrix
    host = ViewCMatrix <T> (this_matrix_.h_view.data(), dim1, dim2, dim3, dim4, dim5, dim6, dim7);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i) const {
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in DCMatrixKokkos 1D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DCMatrixKokkos 1D!");
    return this_matrix_.d_view((i - 1));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j) const {
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in DCMatrixKokkos 2D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DCMatrixKokkos 2D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DCMatrixKokkos 2D!");
    return this_matrix_.d_view((j - 1) + ((i - 1) * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k) const {
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in DCMatrixKokkos 3D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DCMatrixKokkos 3D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DCMatrixKokkos 3D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DCMatrixKokkos 3D!");
    return this_matrix_.d_view((k - 1) + ((j - 1) * dims_[2])
                                       + ((i - 1) * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in DCMatrixKokkos 4D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DCMatrixKokkos 4D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DCMatrixKokkos 4D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DCMatrixKokkos 4D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in DCMatrixKokkos 4D!");
    return this_matrix_.d_view((l - 1) + ((k - 1) * dims_[3])
                                       + ((j - 1) * dims_[3] * dims_[2])
                                       + ((i - 1) * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m) const {
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in DCMatrixKokkos 5D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DCMatrixKokkos 5D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DCMatrixKokkos 5D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DCMatrixKokkos 5D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in DCMatrixKokkos 5D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in DCMatrixKokkos 5D!");
    return this_matrix_.d_view((m - 1) + ((l - 1) * dims_[4])
                                       + ((k - 1) * dims_[4] * dims_[3])
                                       + ((j - 1) * dims_[4] * dims_[3] * dims_[2])
                                       + ((i - 1) * dims_[4] * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n) const {
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in DCMatrixKokkos 6D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DCMatrixKokkos 6D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DCMatrixKokkos 6D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DCMatrixKokkos 6D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in DCMatrixKokkos 6D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in DCMatrixKokkos 6D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in DCMatrixKokkos 6D!");
    return this_matrix_.d_view((n - 1) + ((m - 1) * dims_[5])
                                       + ((l - 1) * dims_[5] * dims_[4])
                                       + ((k - 1) * dims_[5] * dims_[4] * dims_[3])
                                       + ((j - 1) * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                                       + ((i - 1) * dims_[5] * dims_[4] * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n, size_t o) const {
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in DCMatrixKokkos 7D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DCMatrixKokkos 7D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DCMatrixKokkos 7D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DCMatrixKokkos 7D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in DCMatrixKokkos 7D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in DCMatrixKokkos 7D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in DCMatrixKokkos 7D!");
    assert(o >= 1 && o <= dims_[6] && "o is out of bounds in DCMatrixKokkos 7D!");
    return this_matrix_.d_view((o-1) + ((n - 1) * dims_[6])
                                     + ((m - 1) * dims_[6] * dims_[5])
                                     + ((l - 1) * dims_[6] * dims_[5] * dims_[4])
                                     + ((k - 1) * dims_[6] * dims_[5] * dims_[4] * dims_[3])
                                     + ((j - 1) * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                                     + ((i - 1) * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>& DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator= (const DCMatrixKokkos& temp) {
    
    // Do nothing if the assignment is of the form x = x
    if (this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for

        order_ = temp.order_;
        length_ = temp.length_;
        this_matrix_ = temp.this_matrix_;
        host = temp.host;
    }
    
    return *this;
}

// Return size
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::size() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::extent() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::dims(size_t i) const {
    i--;
    assert(i < order_ && "DCMatrixKokkos order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to DCMatrixKokkos dims is out of bounds!");
    return dims_[i];
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::order() const {
    return order_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::device_pointer() const {
    return this_matrix_.d_view.data();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::host_pointer() const {
    return this_matrix_.h_view.data();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
void DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::update_host() {

    this_matrix_.template modify<typename TArray1D::execution_space>();
    this_matrix_.template sync<typename TArray1D::host_mirror_space>();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
void DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::update_device() {

    this_matrix_.template modify<typename TArray1D::host_mirror_space>();
    this_matrix_.template sync<typename TArray1D::execution_space>();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::~DCMatrixKokkos() {}
// End DCMatrixKokkos


/////////////////////////
// DViewCMatrixKokkos
/////////////////////////
template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace, typename MemoryTraits = void>
class DViewCMatrixKokkos {

    // this is always unmanaged
    using TArray1DHost = Kokkos::View<T*, Layout, HostSpace, MemoryUnmanaged>;
    // this is manage
    using TArray1D     = Kokkos::View<T*, Layout, ExecSpace, MemoryTraits>;
    typename ExecSpace::memory_space memspace;
    
private:
    size_t dims_[7];
    size_t length_;
    size_t order_;  // tensor order (rank)
    TArray1D this_matrix_;
    TArray1DHost this_matrix_host_;
    T * temp_inp_matrix_;

public:
    DViewCMatrixKokkos();
    
    DViewCMatrixKokkos(T * inp_matrix, size_t dim1, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, size_t dim5, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, size_t dim5, size_t dim6, const std::string& tag_string = DEFAULTSTRINGMATRIX);

    DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, size_t dim5, size_t dim6,
                 size_t dim7, const std::string& tag_string = DEFAULTSTRINGMATRIX);
    
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n) const;

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m,
                  size_t n, size_t o) const;
    
    KOKKOS_INLINE_FUNCTION
    DViewCMatrixKokkos& operator=(const DViewCMatrixKokkos& temp);

    // GPU Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t size() const;

    // Host Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t extent() const;

    KOKKOS_INLINE_FUNCTION
    size_t dims(size_t i) const;

    KOKKOS_INLINE_FUNCTION
    size_t order() const;
 
    // Method returns the raw device pointer of the Kokkos View
    KOKKOS_INLINE_FUNCTION
    T* device_pointer() const;

    // Method returns the raw host pointer of the Kokkos View
    KOKKOS_INLINE_FUNCTION
    T* host_pointer() const;

    // Data member to access host view
    ViewCMatrix <T> host;

    // Method that update host view
    void update_host();

    // Method that update device view
    void update_device();

    // Deconstructor
    KOKKOS_INLINE_FUNCTION
    ~DViewCMatrixKokkos ();
}; // End of DViewCMatrixKokkos


// Default constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCMatrixKokkos() {
    length_ = order_ = 0;
    temp_inp_matrix_ = NULL;
    for (int i = 0; i < 7; i++) {
        dims_[i] = 0;
    }
}

// Overloaded 1D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCMatrixKokkos(T * inp_matrix, size_t dim1, const std::string& tag_string) {
    
    dims_[0] = dim1;
    order_ = 1;
    length_ = dim1;
    // Create a 1D host view of the external allocation
    this_matrix_host_ = TArray1DHost(inp_matrix, length_);
    // Assign temp point to inp_matrix pointer that is passed in
    temp_inp_matrix_ = inp_matrix;
    // Create a device copy of that host view
    this_matrix_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_matrix_host_);
    // Create host ViewCMatrix. Note: inp_matrix and this_matrix_host_.data() are the same pointer
    host = ViewCMatrix <T> (inp_matrix, dim1);
}

// Overloaded 2D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    order_ = 2;
    length_ = (dim1 * dim2);
    // Create a 1D host view of the external allocation
    this_matrix_host_ = TArray1DHost(inp_matrix, length_);
    // Assign temp point to inp_matrix pointer that is passed in
    temp_inp_matrix_ = inp_matrix;
    // Create a device copy of that host view
    this_matrix_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_matrix_host_);
    // Create host ViewCMatrix
    host = ViewCMatrix <T> (inp_matrix, dim1, dim2);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    order_ = 3;
    length_ = (dim1 * dim2 * dim3);
    // Create a 1D host view of the external allocation
    this_matrix_host_ = TArray1DHost(inp_matrix, length_);
    // Assign temp point to inp_matrix pointer that is passed in
    temp_inp_matrix_ = inp_matrix;
    // Create a device copy of that host view
    this_matrix_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_matrix_host_);
    // Create host ViewCMatrix
    host = ViewCMatrix <T> (inp_matrix, dim1, dim2, dim3);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    order_ = 4;
    length_ = (dim1 * dim2 * dim3 * dim4);
    // Create a 1D host view of the external allocation
    this_matrix_host_ = TArray1DHost(inp_matrix, length_);
    // Assign temp point to inp_matrix pointer that is passed in
    temp_inp_matrix_ = inp_matrix;
    // Create a device copy of that host view
    this_matrix_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_matrix_host_);
    // Create host ViewCMatrix
    host = ViewCMatrix <T> (inp_matrix, dim1, dim2, dim3, dim4);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4, 
                              size_t dim5, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    order_ = 5;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5);
    // Create a 1D host view of the external allocation
    this_matrix_host_ = TArray1DHost(inp_matrix, length_);
    // Assign temp point to inp_matrix pointer that is passed in
    temp_inp_matrix_ = inp_matrix;
    // Create a device copy of that host view
    this_matrix_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_matrix_host_);
    // Create host ViewCMatrix
    host = ViewCMatrix <T> (inp_matrix, dim1, dim2, dim3, dim4, dim5);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4, 
                              size_t dim5, size_t dim6, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    order_ = 6;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5 * dim6);
    // Create a 1D host view of the external allocation
    this_matrix_host_ = TArray1DHost(inp_matrix, length_);
    // Assign temp point to inp_matrix pointer that is passed in
    temp_inp_matrix_ = inp_matrix;
    // Create a device copy of that host view
    this_matrix_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_matrix_host_);
    // Create host ViewCMatrix
    host = ViewCMatrix <T> (inp_matrix, dim1, dim2, dim3, dim4, dim5, dim6);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4,
                              size_t dim5, size_t dim6,
                              size_t dim7, const std::string& tag_string) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    dims_[6] = dim7;
    order_ = 7;
    length_ = (dim1 * dim2 * dim3 * dim4 * dim5 * dim6 * dim7);
    // Create a 1D host view of the external allocation
    this_matrix_host_ = TArray1DHost(inp_matrix, length_);
    // Assign temp point to inp_matrix pointer that is passed in
    temp_inp_matrix_ = inp_matrix;
    // Create a device copy of that host view
    this_matrix_ = create_mirror_view_and_copy(Kokkos::view_alloc(memspace, tag_string), this_matrix_host_);
    // Create host ViewCMatrix
    host = ViewCMatrix <T> (inp_matrix, dim1, dim2, dim3, dim4, dim5, dim6, dim7);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i) const {
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in DViewCMatrixKokkos 1D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DViewCMatrixKokkos 1D!");
    return this_matrix_((i - 1));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j) const {
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in DViewCMatrixKokkos 2D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DViewCMatrixKokkos 2D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DViewCMatrixKokkos 2D!");
    return this_matrix_((j - 1) + ((i - 1) * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k) const {
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in DViewCMatrixKokkos 3D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DViewCMatrixKokkos 3D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DViewCMatrixKokkos 3D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DViewCMatrixKokkos 3D!");
    return this_matrix_((k - 1) + ((j - 1) * dims_[2])
                                + ((i - 1) * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in DViewCMatrixKokkos 4D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DViewCMatrixKokkos 4D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DViewCMatrixKokkos 4D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DViewCMatrixKokkos 4D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in DViewCMatrixKokkos 4D!");
    return this_matrix_((l - 1) + ((k - 1) * dims_[3])
                                + ((j - 1) * dims_[3] * dims_[2])
                                + ((i - 1) * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m) const {
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in DViewCMatrixKokkos 5D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DViewCMatrixKokkos 5D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DViewCMatrixKokkos 5D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DViewCMatrixKokkos 5D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in DViewCMatrixKokkos 5D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in DViewCMatrixKokkos 5D!");
    return this_matrix_((m - 1) + ((l - 1) * dims_[4])
                                + ((k - 1) * dims_[4] * dims_[3])
                                + ((j - 1) * dims_[4] * dims_[3] * dims_[2])
                                + ((i - 1) * dims_[4] * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n) const {
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in DViewCMatrixKokkos 6D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DViewCMatrixKokkos 6D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DViewCMatrixKokkos 6D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DViewCMatrixKokkos 6D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in DViewCMatrixKokkos 6D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in DViewCMatrixKokkos 6D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in DViewCMatrixKokkos 6D!");
    return this_matrix_((n - 1) + ((m - 1) * dims_[5])
                                + ((l - 1) * dims_[5] * dims_[4])
                                + ((k - 1) * dims_[5] * dims_[4] * dims_[3])
                                + ((j - 1) * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                                + ((i - 1) * dims_[5] * dims_[4] * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n, size_t o) const {
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in DViewCMatrixKokkos 7D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in DViewCMatrixKokkos 7D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in DViewCMatrixKokkos 7D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in DViewCMatrixKokkos 7D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in DViewCMatrixKokkos 7D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in DViewCMatrixKokkos 7D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in DViewCMatrixKokkos 7D!");
    assert(o >= 1 && o <= dims_[6] && "o is out of bounds in DViewCMatrixKokkos 7D!");
    return this_matrix_(o + ((n - 1) * dims_[6])
                          + ((m - 1) * dims_[6] * dims_[5])
                          + ((l - 1) * dims_[6] * dims_[5] * dims_[4])
                          + ((k - 1) * dims_[6] * dims_[5] * dims_[4] * dims_[3])
                          + ((j - 1) * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                          + ((i - 1) * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2] * dims_[1]));
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>& DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::operator= (const DViewCMatrixKokkos& temp) {
    
    // Do nothing if the assignment is of the form x = x
    if (this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for

        order_ = temp.order_;
        length_ = temp.length_;
        temp_inp_matrix_ = temp.temp_inp_matrix_;
        this_matrix_host_ = temp.this_matrix_host_;
        this_matrix_ = temp.this_matrix_;
        host = temp.host;
    }
    
    return *this;
}

// Return size
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::size() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::extent() const {
    return length_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::dims(size_t i) const {
    i--;
    assert(i < order_ && "DViewCMatrixKokkos order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to DViewCMatrixKokkos dims is out of bounds!");
    return dims_[i];
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::order() const {
    return order_;
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::device_pointer() const {
    return this_matrix_.data();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::host_pointer() const {
    return this_matrix_host_.data();
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
void DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::update_host() {
    // Deep copy of device view to host view
    deep_copy(this_matrix_host_, this_matrix_);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
void DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::update_device() {
    // Deep copy of host view to device view
    deep_copy(this_matrix_, this_matrix_host_);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::~DViewCMatrixKokkos() {}
// End DViewCMatrixKokkos


/*! \brief Kokkos version of the serial RaggedRightArray class.
 *
 */
template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace,
 typename MemoryTraits = void, typename ILayout = Layout>
class RaggedRightArrayKokkos {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace, MemoryTraits>;
    using SArray1D = Kokkos::View<size_t *,Layout, ExecSpace, MemoryTraits>;
    using Strides1D = Kokkos::View<size_t *,ILayout, ExecSpace, MemoryTraits>;
    
private:
    TArray1D array_;
    
    size_t dim1_;
    size_t length_;
    
public:
    // Default constructor
    RaggedRightArrayKokkos();
    
    //--- 2D array access of a ragged right array ---
    
    // Overload constructor for a CArrayKokkos
    RaggedRightArrayKokkos(CArrayKokkos<size_t,ILayout,ExecSpace,MemoryTraits> &strides_array, const std::string& tag_string = DEFAULTSTRINGARRAY);

    // Overload constructor for a DCArrayKokkos
    RaggedRightArrayKokkos(DCArrayKokkos<size_t,ILayout,ExecSpace,MemoryTraits> &strides_array, const std::string& tag_string = DEFAULTSTRINGARRAY);
    
    // Overload constructor for a ViewCArray
    RaggedRightArrayKokkos(ViewCArray<size_t> &strides_array, const std::string& tag_string = DEFAULTSTRINGARRAY);
    
    // Overloaded constructor for a traditional array
    RaggedRightArrayKokkos(size_t* strides_array, size_t some_dim1, const std::string& tag_string = DEFAULTSTRINGARRAY);
    
    // A method to return the stride size
    KOKKOS_INLINE_FUNCTION
    size_t stride(size_t i) const;

    // Host method to return the stride size
    size_t stride_host(size_t i) const;
    
    // A method to increase the number of column entries, i.e.,
    // the stride size. Used with the constructor for building
    // the stride_array dynamically.
    // DO NOT USE with the constructures with a strides_array
    KOKKOS_INLINE_FUNCTION
    size_t& build_stride(const size_t i) const;
    
    KOKKOS_INLINE_FUNCTION
    void stride_finalize() const;
    
    // Overload operator() to access data as array(i,j)
    // where i=[0:N-1], j=[stride(i)]
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;

    // method to return total size
    KOKKOS_INLINE_FUNCTION
    size_t size(){
      return length_;
    }
    
    //setup start indices
    void data_setup(const std::string& tag_string);
    
    KOKKOS_INLINE_FUNCTION
    T* pointer();

    //return the view
    KOKKOS_INLINE_FUNCTION
    TArray1D get_kokkos_view();

    // Kokkos views of strides and start indices
    Strides1D mystrides_;
    SArray1D start_index_;
    
    KOKKOS_INLINE_FUNCTION
    RaggedRightArrayKokkos& operator= (const RaggedRightArrayKokkos &temp);

    //initialize start indices view
    class init_start_indices_functor{
      public:
      SArray1D mystart_index_;
      init_start_indices_functor(SArray1D tempstart_index_){
        mystart_index_ = tempstart_index_;
      }
      KOKKOS_INLINE_FUNCTION void operator()(const int index) const {
        mystart_index_(index) = 0;
      }
    };

    //setup start indices view
    class setup_start_indices_functor{
        public:
        SArray1D mystart_index_;
        Strides1D mytemp_strides_;
        setup_start_indices_functor(SArray1D tempstart_index_, Strides1D temp_strides_){
          mystart_index_ = tempstart_index_;
          mytemp_strides_ = temp_strides_;
        }
        KOKKOS_INLINE_FUNCTION void operator()(const int index, int& update, bool final) const {
          // Load old value in case we update it before accumulating
            const size_t count = mytemp_strides_(index);
            update += count;
            if (final) {
                mystart_index_((index+1)) = update;
            }
        }
    };

    //setup length of view
    class setup_length_functor{
        public:
        //kokkos needs this typedef named
        typedef size_t value_type;
        // This is helpful for determining the right index type,
        // especially if you expect to need a 64-bit index.
        //typedef Kokkos::View<size_t*>::size_type size_type;
        Strides1D mytemp_strides_;
        setup_length_functor(Strides1D temp_strides_){
          mytemp_strides_ = temp_strides_;
        }
        KOKKOS_INLINE_FUNCTION void operator()(const int index, size_t& update) const {
            //const size_t count = mytemp_strides_(index);
            update += mytemp_strides_(index);
        }
    };

    //sets final 1D array size
    class finalize_stride_functor{
        public:
        SArray1D mystart_index_;
        finalize_stride_functor(SArray1D tempstart_index_){
          mystart_index_ = tempstart_index_;
        }
        KOKKOS_INLINE_FUNCTION void operator()(const int index, int& update, bool final) const {
          // Load old value in case we update it before accumulating
            const size_t count = mystart_index_(index+1);
            update += count;
            if (final) {
                mystart_index_((index+1)) = update;
            }
        }
    };

    // Destructor
    KOKKOS_INLINE_FUNCTION
    ~RaggedRightArrayKokkos ( );
}; // End of RaggedRightArray

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
RaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::RaggedRightArrayKokkos() {
    dim1_ = length_ = 0;
}

// Overloaded constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
RaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::RaggedRightArrayKokkos(CArrayKokkos<size_t,ILayout,ExecSpace,MemoryTraits> &strides_array,
                                                                                        const std::string& tag_string) {
    mystrides_ = strides_array.get_kokkos_view();
    dim1_ = strides_array.extent();
    data_setup(tag_string);
} // End constructor

// Overloaded constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
RaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::RaggedRightArrayKokkos(DCArrayKokkos<size_t,ILayout,ExecSpace,MemoryTraits> &strides_array,
                                                                                        const std::string& tag_string) {
    mystrides_ = strides_array.get_kokkos_dual_view().d_view;
    dim1_ = strides_array.extent();
    data_setup(tag_string);
} // End constructor

// Overloaded constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
RaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::RaggedRightArrayKokkos(ViewCArray<size_t> &strides_array,
                                                                                         const std::string& tag_string) {
} // End constructor

// Overloaded constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
RaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::RaggedRightArrayKokkos(size_t* strides_array,  size_t some_dim1,
                                                                                        const std::string& tag_string) {
    mystrides_.assign_data(strides_array);
    dim1_ = some_dim1;
    data_setup(tag_string);
} // End constructor

//setup start indices
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
void RaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::data_setup(const std::string& tag_string) {
    //allocate start indices
    std::string append_indices_string("start_indices");
    std::string append_array_string("array");
    std::string temp_copy_string = tag_string;
    std::string start_index_tag_string = temp_copy_string.append(append_indices_string);
    temp_copy_string = tag_string;
    std::string array_tag_string = temp_copy_string.append(append_array_string);

    start_index_ = SArray1D(start_index_tag_string,dim1_ + 1);
    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_for("StartValuesInit", dim1_+1, KOKKOS_CLASS_LAMBDA(const int i) {
      start_index_(i) = 0;
    });
    #else
    init_start_indices_functor execution_functor(start_index_);
    Kokkos::parallel_for("StartValuesInit", dim1_+1,execution_functor);
    #endif

    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_scan("StartValuesSetup", dim1_, KOKKOS_CLASS_LAMBDA(const int i, int& update, const bool final) {
            // Load old value in case we update it before accumulating
            const size_t count = mystrides_(i);
            update += count;
            if (final) {
                start_index_((i+1)) = update;
            }

        });
    #else
    setup_start_indices_functor setup_execution_functor(start_index_, mystrides_);
    Kokkos::parallel_scan("StartValuesSetup", dim1_,setup_execution_functor);
    #endif

    //compute length of the storage
    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_reduce("LengthSetup", dim1_, KOKKOS_CLASS_LAMBDA(const int i, int& update) {
            // Load old value in case we update it before accumulating
            update += mystrides_(i);
        }, length_);
    #else
    setup_length_functor length_functor(mystrides_);
    Kokkos::parallel_reduce("LengthSetup", dim1_, length_functor, length_);
    #endif

    //allocate view
    array_ = TArray1D(array_tag_string, length_);
}

// A method to return the stride size
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
size_t RaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::stride(size_t i) const {
    // Ensure that i is within bounds
    assert(i < (dim1_) && "i is greater than dim1_ in RaggedRightArray");
    return mystrides_(i);
}

// Method to build the stride (non-Kokkos push back)
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
size_t& RaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::build_stride(const size_t i) const {
    return start_index_(i+1);
}

// Method to finalize stride
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
void RaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::stride_finalize() const {
    
    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_scan("StartValues", dim1_, KOKKOS_CLASS_LAMBDA(const int i, int& update, const bool final) {
            // Load old value in case we update it before accumulating
            const size_t count = start_index_(i+1);
            update += count;
            if (final) {
                start_index_((i+1)) = update;
            }

        });
    #else
    finalize_stride_functor execution_functor(start_index_);
    Kokkos::parallel_scan("StartValues", dim1_,execution_functor);
    #endif
    Kokkos::fence();
}


// Overload operator() to access data as array(i,j)
// where i=[0:N-1], j=[0:stride(i)]
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
T& RaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::operator()(size_t i, size_t j) const {
    // Get the 1D array index
    size_t start = start_index_(i);
    
    // asserts
    assert(i < dim1_ && "i is out of dim1 bounds in RaggedRightArrayKokkos");  // die if >= dim1
    assert(j < stride(i) && "j is out of stride bounds in RaggedRightArrayKokkos");  // die if >= stride
    
    return array_(j + start);
} // End operator()

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
T* RaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::pointer() {
    return array_.data();
}


template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
RaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout> & RaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::
  operator= (const RaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout> &temp) {

  if (this != &temp) {
      /*
    SArray1D tempdim = SArray1D("tempdim", 1);
    auto h_tempdim = HostMirror(tempdim);
    Kokkos::parallel_for("StrideDim", 1, KOKKOS_CLASS_LAMBDA(const int&) {
            tempdim(0)  = strides_array.size();
            //dim1_  = strides_array.size();
        });
    Kokkos::fence();
    deep_copy(h_tempdim, tempdim);
    dim1_ = h_tempdim(0);
    */
    dim1_ = temp.dim1_;
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = temp.start_index_;
    //start_index_(0) = 0; // the 1D array starts at 0

    /*
    size_t * h_start_index = new size_t [dim1_+1];
    h_start_index[0] = 0;
    size_t * herenow = new size_t [2];
    herenow[0] = 1;
    herenow[1] = 2;
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += herenow[i];
        h_start_index[(i + 1)] = count;
        printf("%d) Start check %ld\n", i, h_start_index[i]);
    } // end for i
    */
    /*
    SArray1D templen = SArray1D("templen", 1);
    auto h_templen = Kokkos::create_mirror_view(templen);
    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_for("ArrayLength", 1, KOKKOS_CLASS_LAMBDA(const int&) {
            templen(0) = start_index_(dim1_);
            //length_ = start_index_(dim1_);
        });
    #else
    templen_functor templen_execution_functor(templen);
    Kokkos::parallel_for("ArrayLength", 1, templen_execution_functor);
    #endif
    Kokkos::fence();
    Kokkos::deep_copy(h_templen, templen);
    if (h_templen(0) != 0)
        length_ = h_templen(0);
    else
    */
    length_ = temp.length_;


    //printf("Length %ld\n", length_);

    //Kokkos::parallel_for("StartCheck", dim1_+1, KOKKOS_CLASS_LAMBDA(const int i) {
    //        printf("%d) Start %ld\n", i, start_index_(i));
    //    });
    //Kokkos::fence();
    
    array_ = temp.array_;
    mystrides_ = temp.mystrides_;

    /*
        dim1_ = temp.dim1_;
        length_ = temp.length_;
        start_index_ = SArray1D("start_index_", dim1_ + 1);
        Kokkos::parallel_for("EqualOperator", dim1_+1, KOKKOS_CLASS_LAMBDA(const int j) {
                start_index_(j) = temp.start_index_(j);
            });
        //for (int j = 0; j < dim1_; j++) {
        //    start_index_(j) = temp.start_index_(j);
        //}
        array_ = TArray1D("array_", length_);
    */
  }
    
    return *this;
}

//return the stored Kokkos view
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
Kokkos::View<T*, Layout, ExecSpace, MemoryTraits> RaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::get_kokkos_view() {
    return array_;
}

// Destructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
RaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::~RaggedRightArrayKokkos() { }

////////////////////////////////////////////////////////////////////////////////
// End of RaggedRightArrayKokkos
////////////////////////////////////////////////////////////////////////////////

/*! \brief Kokkos version of the RaggedRightArray class.
 *
 */
template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace, typename MemoryTraits = void, typename ILayout = Layout>
class RaggedRightArrayofVectorsKokkos {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace, MemoryTraits>;
    using SArray1D = Kokkos::View<size_t *,Layout, ExecSpace, MemoryTraits>;
    using Strides1D = Kokkos::View<size_t *,ILayout, ExecSpace, MemoryTraits>;
    
private:
    TArray1D array_;
    
    size_t dim1_, vector_dim_;
    size_t length_;
    
public:
    // Default constructor
    RaggedRightArrayofVectorsKokkos();
    
    //--- 2D array access of a ragged right array ---
    
    // Overload constructor for a CArrayKokkos
    RaggedRightArrayofVectorsKokkos(CArrayKokkos<size_t,ILayout,ExecSpace,MemoryTraits> &strides_array, size_t vector_dim,
                                    const std::string& tag_string = DEFAULTSTRINGARRAY );
    
    // Overload constructor for a ViewCArray
    RaggedRightArrayofVectorsKokkos(ViewCArray<size_t> &strides_array, size_t vector_dim, const std::string& tag_string = DEFAULTSTRINGARRAY);
    
    // Overloaded constructor for a traditional array
    RaggedRightArrayofVectorsKokkos(size_t* strides_array, size_t some_dim1, size_t vector_dim, const std::string& tag_string = DEFAULTSTRINGARRAY);
    
    // A method to return the stride size
    KOKKOS_INLINE_FUNCTION
    size_t stride(size_t i) const;
    
    // A method to increase the number of column entries, i.e.,
    // the stride size. Used with the constructor for building
    // the stride_array dynamically.
    // DO NOT USE with the constructures with a strides_array
    KOKKOS_INLINE_FUNCTION
    size_t& build_stride(const size_t i) const;
    
    KOKKOS_INLINE_FUNCTION
    void stride_finalize() const;
    
    // Overload operator() to access data as array(i,j)
    // where i=[0:N-1], j=[stride(i)]
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    // method to return total size
    KOKKOS_INLINE_FUNCTION
    size_t size(){
      return length_;
    }
    
    //setup start indices
    void data_setup(const std::string& tag_string);
    
    KOKKOS_INLINE_FUNCTION
    T* pointer();

    //return the view
    KOKKOS_INLINE_FUNCTION
    TArray1D get_kokkos_view();

    // Kokkos views of strides and start indices
    Strides1D mystrides_;
    SArray1D start_index_;
    
    KOKKOS_INLINE_FUNCTION
    RaggedRightArrayofVectorsKokkos& operator= (const RaggedRightArrayofVectorsKokkos &temp);

    //functors for kokkos execution policies
    // initialize start indices view
    class init_start_indices_functor{
      public:
      SArray1D mystart_index_;
      init_start_indices_functor(SArray1D tempstart_index_){
        mystart_index_ = tempstart_index_;
      }
      KOKKOS_INLINE_FUNCTION void operator()(const int index) const {
        mystart_index_(index) = 0;
      }
    };

    //setup start indices view
    class setup_start_indices_functor{
        public:
        SArray1D mystart_index_;
        Strides1D mytemp_strides_;
        size_t myvector_dim_;
        setup_start_indices_functor(SArray1D tempstart_index_, Strides1D temp_strides_, size_t myvector_dim){
          mystart_index_ = tempstart_index_;
          mytemp_strides_ = temp_strides_;
          myvector_dim_ = myvector_dim;
        }
        KOKKOS_INLINE_FUNCTION void operator()(const int index, int& update, bool final) const {
          // Load old value in case we update it before accumulating
            const size_t count = mytemp_strides_(index)*myvector_dim_;
            update += count;
            if (final) {
                mystart_index_((index+1)) = update;
            }
        }
    };

    //setup length of view
    class setup_length_functor{
        public:
        //kokkos needs this typedef named
        typedef size_t value_type;
        // This is helpful for determining the right index type,
        // especially if you expect to need a 64-bit index.
        //typedef Kokkos::View<size_t*>::size_type size_type;

        Strides1D mytemp_strides_;
        size_t myvector_dim_;

        setup_length_functor(Strides1D temp_strides_, size_t myvector_dim){
          mytemp_strides_ = temp_strides_;
          myvector_dim_ = myvector_dim;
        }
        KOKKOS_INLINE_FUNCTION void operator()(const int index, size_t& update) const {
            //const size_t count = mytemp_strides_(index)*myvector_dim_;
            update += mytemp_strides_(index)*myvector_dim_;;
        }
    };

    //sets final 1D array size
    class finalize_stride_functor{
        public:
        SArray1D mystart_index_;
        finalize_stride_functor(SArray1D tempstart_index_){
          mystart_index_ = tempstart_index_;
        }
        KOKKOS_INLINE_FUNCTION void operator()(const int index, int& update, bool final) const {
          // Load old value in case we update it before accumulating
            const size_t count = mystart_index_(index+1);
            update += count;
            if (final) {
                mystart_index_((index+1)) = update;
            }
        }
    };

    // Destructor
    KOKKOS_INLINE_FUNCTION
    ~RaggedRightArrayofVectorsKokkos ( );
}; // End of RaggedRightArrayofVectorsKokkos

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
RaggedRightArrayofVectorsKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::RaggedRightArrayofVectorsKokkos() {
    dim1_ = length_ = vector_dim_ = 0;
}

// Overloaded constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
RaggedRightArrayofVectorsKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::RaggedRightArrayofVectorsKokkos(CArrayKokkos<size_t,ILayout,ExecSpace,MemoryTraits>
                                                                                                          &strides_array, size_t vector_dim,
                                                                                                          const std::string& tag_string) {
    //mystrides_.assign_data(strides_array.pointer());
    vector_dim_ = vector_dim;
    mystrides_ = strides_array.get_kokkos_view();
    dim1_ = strides_array.extent();
    data_setup(tag_string);
} // End constructor

/*
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
RaggedRightArrayofVectorsKokkos<T,Layout,ExecSpace,MemoryTraits>::RaggedRightArrayofVectorsKokkos(CArrayKokkos<size_t,Kokkos::LayoutLeft,ExecSpace,MemoryTraits>
                                                                                                  &strides_array, size_t vector_dim) {
    //mystrides_.assign_data(strides_array.pointer());
    vector_dim_ = vector_dim;
    mystrides_ = strides_array;
    dim1_ = strides_array.extent();
} // End constructor
*/

// Overloaded constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
RaggedRightArrayofVectorsKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::RaggedRightArrayofVectorsKokkos(ViewCArray<size_t> &strides_array, size_t vector_dim,
                                                                                                          const std::string& tag_string) {
} // End constructor

// Overloaded constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
RaggedRightArrayofVectorsKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::RaggedRightArrayofVectorsKokkos(size_t* strides_array, size_t some_dim1, size_t vector_dim,
                                                                                                          const std::string& tag_string) {
    vector_dim_ = vector_dim;
    mystrides_.assign_data(strides_array);
    dim1_ = some_dim1;
    data_setup(tag_string);
} // End constructor

//setup start indices
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
void RaggedRightArrayofVectorsKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::data_setup(const std::string& tag_string) {

    //allocate start indices
    std::string append_indices_string("start_indices");
    std::string append_array_string("array");
    std::string temp_copy_string = tag_string;
    std::string start_index_tag_string = temp_copy_string.append(append_indices_string);
    temp_copy_string = tag_string;
    std::string array_tag_string = temp_copy_string.append(append_array_string);

    start_index_ = SArray1D(start_index_tag_string,dim1_ + 1);
    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_for("StartValuesInit", dim1_+1, KOKKOS_CLASS_LAMBDA(const int i) {
      start_index_(i) = 0;
    });
    #else
    init_start_indices_functor execution_functor(start_index_);
    Kokkos::parallel_for("StartValuesInit", dim1_+1,execution_functor);
    #endif

    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_scan("StartValuesSetup", dim1_, KOKKOS_CLASS_LAMBDA(const int i, int& update, const bool final) {
            // Load old value in case we update it before accumulating
            const size_t count = mystrides_(i)*vector_dim_;
            update += count;
            if (final) {
                start_index_((i+1)) = update;
            }

        });
    #else
    setup_start_indices_functor setup_execution_functor(start_index_, mystrides_, vector_dim_);
    Kokkos::parallel_scan("StartValuesSetup", dim1_,setup_execution_functor);
    #endif

    //compute length of the storage
    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_reduce("LengthSetup", dim1_, KOKKOS_CLASS_LAMBDA(const int i, int& update) {
            // Load old value in case we update it before accumulating
            update += mystrides_(i)*vector_dim_;
        }, length_);
    #else
    setup_length_functor length_functor(mystrides_, vector_dim_);
    Kokkos::parallel_reduce("LengthSetup", dim1_, length_functor,length_);
    #endif

    //allocate view
    array_ = TArray1D(array_tag_string, length_);
}

// A method to return the stride size
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
size_t RaggedRightArrayofVectorsKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::stride(size_t i) const {
    // Ensure that i is within bounds
    assert(i < (dim1_) && "i is greater than dim1_ in RaggedRightArray");
    return mystrides_(i);
}

// Method to build the stride (non-Kokkos push back)
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
size_t& RaggedRightArrayofVectorsKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::build_stride(const size_t i) const {
    return start_index_(i+1);
}

// Method to finalize stride
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
void RaggedRightArrayofVectorsKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::stride_finalize() const {
    
    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_scan("StartValues", dim1_, KOKKOS_CLASS_LAMBDA(const int i, int& update, const bool final) {
            // Load old value in case we update it before accumulating
            const size_t count = start_index_(i+1);
            update += count;
            if (final) {
                start_index_((i+1)) = update;
            }

        });
    #else
    finalize_stride_functor execution_functor(start_index_);
    Kokkos::parallel_scan("StartValues", dim1_,execution_functor);
    #endif
    Kokkos::fence();
}


// Overload operator() to access data as array(i,j)
// where i=[0:N-1], j=[0:stride(i)]
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
T& RaggedRightArrayofVectorsKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::operator()(size_t i, size_t j, size_t k) const {
    // Get the 1D array index
    size_t start = start_index_(i);
    
    // asserts
    assert(i < dim1_ && "i is out of dim1 bounds in RaggedRightArrayKokkos");  // die if >= dim1
    assert(j < stride(i) && "j is out of stride bounds in RaggedRightArrayKokkos");  // die if >= stride
    assert(k < vector_dim_ && "k is out of vector_dim bounds in RaggedRightArrayKokkos");  // die if >= vector_dim
    
    return array_(j*vector_dim_ + start + k);
} // End operator()

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
T* RaggedRightArrayofVectorsKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::pointer() {
    return array_.data();
}


template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
RaggedRightArrayofVectorsKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout> & RaggedRightArrayofVectorsKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::
  operator= (const RaggedRightArrayofVectorsKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout> &temp) {

  if (this != &temp) {
    dim1_ = temp.dim1_;
    vector_dim_ = temp.vector_dim_;
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = temp.start_index_;
    length_ = temp.length_;
    
    array_ = temp.array_;
    mystrides_ = temp.mystrides_;
  }
    
    return *this;
}

//return the stored Kokkos view
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
Kokkos::View<T*, Layout, ExecSpace, MemoryTraits> RaggedRightArrayofVectorsKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::get_kokkos_view() {
    return array_;
}

// Destructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
RaggedRightArrayofVectorsKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::~RaggedRightArrayofVectorsKokkos() { }

////////////////////////////////////////////////////////////////////////////////
// End of RaggedRightArrayofVectorsKokkos
////////////////////////////////////////////////////////////////////////////////

/*! \brief Kokkos version of the serial RaggedDownArray class.
 *
 */
template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace,
          typename MemoryTraits = void, typename ILayout = Layout>
class RaggedDownArrayKokkos {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace, MemoryTraits>;
    using SArray1D = Kokkos::View<size_t *, Layout, ExecSpace, MemoryTraits>;
    using Strides1D = Kokkos::View<size_t *, ILayout, ExecSpace, MemoryTraits>;
    
private:
    TArray1D array_;
    
    size_t dim2_;
    size_t length_;
    
public:
    // Default constructor
    RaggedDownArrayKokkos();
    
    //--- 2D array access of a ragged right array ---
    
    // Overload constructor for a CArray
    RaggedDownArrayKokkos(CArrayKokkos<size_t, Layout, ExecSpace, MemoryTraits> &strides_array, const std::string& tag_string = DEFAULTSTRINGARRAY);
    
    // Overload constructor for a ViewCArray
    RaggedDownArrayKokkos(ViewCArray<size_t> &strides_array, const std::string& tag_string = DEFAULTSTRINGARRAY);
    
    // Overloaded constructor for a traditional array
    RaggedDownArrayKokkos(size_t* strides_array, size_t some_dim2, const std::string& tag_string = DEFAULTSTRINGARRAY);

    // A method to return the stride size
    KOKKOS_INLINE_FUNCTION
    size_t stride(size_t j) const;

    //setup start indices
    void data_setup(const std::string& tag_string);
    
    // Overload operator() to access data as array(i,j)
    // where i=[0:N-1], j=[stride(i)]
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;
    
    KOKKOS_INLINE_FUNCTION
    T* pointer();

    //return the view
    KOKKOS_INLINE_FUNCTION
    TArray1D get_kokkos_view();
    
    KOKKOS_INLINE_FUNCTION
    RaggedDownArrayKokkos& operator= (const RaggedDownArrayKokkos &temp);

    // Kokkos views of strides and start indices
    Strides1D mystrides_;
    SArray1D start_index_;
    
    //functors for kokkos execution policies
    //initialize start indices view
    class init_start_indices_functor{
      public:
      SArray1D mystart_index_;
      init_start_indices_functor(SArray1D tempstart_index_){
        mystart_index_ = tempstart_index_;
      }
      KOKKOS_INLINE_FUNCTION void operator()(const int index) const {
        mystart_index_(index) = 0;
      }
    };

    //setup start indices view
    class setup_start_indices_functor{
        public:
        SArray1D mystart_index_;
        Strides1D mytemp_strides_;
        setup_start_indices_functor(SArray1D tempstart_index_, Strides1D temp_strides_){
          mystart_index_ = tempstart_index_;
          mytemp_strides_ = temp_strides_;
        }
        KOKKOS_INLINE_FUNCTION void operator()(const int index, int& update, bool final) const {
          // Load old value in case we update it before accumulating
            const size_t count = mytemp_strides_(index);
            update += count;
            if (final) {
                mystart_index_((index+1)) = update;
            }
        }
    };

    //setup length of view
    class setup_length_functor{
        public:
        //kokkos needs this typedef named
        typedef size_t value_type;
        // This is helpful for determining the right index type,
        // especially if you expect to need a 64-bit index.
        //typedef Kokkos::View<size_t*>::size_type size_type;
        Strides1D mytemp_strides_;
        setup_length_functor(Strides1D temp_strides_){
          mytemp_strides_ = temp_strides_;
        }
        KOKKOS_INLINE_FUNCTION void operator()(const int index, size_t& update) const {
            //const size_t count = mytemp_strides_(index);
            update += mytemp_strides_(index);
        }
    };

    // Destructor
    KOKKOS_INLINE_FUNCTION
    ~RaggedDownArrayKokkos ( );
}; // End of RaggedDownArray

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
RaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::RaggedDownArrayKokkos() {
    dim2_ = length_ = 0;
}

// Overloaded constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
RaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::RaggedDownArrayKokkos(CArrayKokkos<size_t, Layout, ExecSpace, MemoryTraits> &strides_array,
                                                                              const std::string& tag_string) {
    mystrides_ = strides_array.get_kokkos_view();
    dim2_ = strides_array.extent();
    data_setup(tag_string);
} // End constructor

// Overloaded constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
RaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::RaggedDownArrayKokkos(ViewCArray<size_t> &strides_array, const std::string& tag_string) {
} // End constructor

// Overloaded constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
RaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::RaggedDownArrayKokkos(size_t* strides_array, size_t some_dim2,
                                                                              const std::string& tag_string) {
    mystrides_.assign_data(strides_array);
    dim2_ = some_dim2;
    data_setup(tag_string);
} // End constructor


//setup start indices
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
void RaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::data_setup(const std::string& tag_string) {
    //allocate start indices
    std::string append_indices_string("start_indices");
    std::string append_array_string("array");
    std::string temp_copy_string = tag_string;
    std::string start_index_tag_string = temp_copy_string.append(append_indices_string);
    temp_copy_string = tag_string;
    std::string array_tag_string = temp_copy_string.append(append_array_string);

    start_index_ = SArray1D(start_index_tag_string,dim2_ + 1);
    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_for("StartValuesInit", dim2_+1, KOKKOS_CLASS_LAMBDA(const int i) {
      start_index_(i) = 0;
    });
    #else
    init_start_indices_functor execution_functor(start_index_);
    Kokkos::parallel_for("StartValuesInit", dim2_+1,execution_functor);
    #endif

    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_scan("StartValuesSetup", dim2_, KOKKOS_CLASS_LAMBDA(const int i, int& update, const bool final) {
            // Load old value in case we update it before accumulating
            const size_t count = mystrides_(i);
            update += count;
            if (final) {
                start_index_((i+1)) = update;
            }

        });
    #else
    setup_start_indices_functor setup_execution_functor(start_index_, mystrides_);
    Kokkos::parallel_scan("StartValuesSetup", dim2_,setup_execution_functor);
    #endif

    //compute length of the storage
    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_reduce("LengthSetup", dim2_, KOKKOS_CLASS_LAMBDA(const int i, int& update) {
            // Load old value in case we update it before accumulating
            update += mystrides_(i);
        }, length_);
    #else
    setup_length_functor length_functor(mystrides_);
    Kokkos::parallel_reduce("LengthSetup", dim2_, length_functor, length_);
    #endif

    //allocate view
    array_ = TArray1D(array_tag_string, length_);
}

// A method to return the stride size
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
size_t RaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::stride(size_t j) const {
    // Ensure that j is within bounds
    assert(j < (dim2_) && "j is greater than dim1_ in RaggedDownArray");

    return mystrides_(j);
}

// Overload operator() to access data as array(i,j)
// where i=[0:N-1], j=[0:stride(i)]
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
T& RaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::operator()(size_t i, size_t j) const {
    // Get the 1D array index
    size_t start = start_index_(j);
    
    // asserts
    assert(i < stride(j) && "i is out of stride bounds in RaggedDownArrayKokkos");  // die if >= stride
    assert(j < dim2_ && "j is out of dim1 bounds in RaggedDownArrayKokkos");  // die if >= dim1
    
    return array_(i + start);
} // End operator()

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
RaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>& RaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::
operator= (const RaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout> &temp) {

  if (this != &temp) {
      /*
    SArray1D tempdim = SArray1D("tempdim", 1);
    auto h_tempdim = HostMirror(tempdim);
    Kokkos::parallel_for("StrideDim", 1, KOKKOS_CLASS_LAMBDA(const int&) {
            tempdim(0)  = strides_array.size();
            //dim1_  = strides_array.size();
        });
    Kokkos::fence();
    deep_copy(h_tempdim, tempdim);
    dim1_ = h_tempdim(0);
    */
    dim2_ = temp.dim2_;
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = temp.start_index_;
    /*
    //start_index_(0) = 0; // the 1D array starts at 0
    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_for("StartFirst", 1, KOKKOS_CLASS_LAMBDA(const int&) {
            start_index_(0) = 0;
        });
    #else
    assignment_init_functor init_execution_functor;
    Kokkos::parallel_for("StartFirst", 1, init_execution_functor);
    #endif
    Kokkos::fence();
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_scan("StartValues", dim2_, KOKKOS_CLASS_LAMBDA(const int j, double& update, const bool final) {
            // Load old value in case we update it before accumulating
            const size_t count = temp.mystrides_[j];
            update += count;
            if (final) {
                start_index_((j+1)) = update;
            }

        });
    #else
    assignment_scan_functor scan_execution_functor(temp);
    Kokkos::parallel_scan("StartValues", dim2_, scan_execution_functor);
    #endif
    Kokkos::fence();
    */
    /*
    size_t * h_start_index = new size_t [dim1_+1];
    h_start_index[0] = 0;
    size_t * herenow = new size_t [2];
    herenow[0] = 1;
    herenow[1] = 2;
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += herenow[i];
        h_start_index[(i + 1)] = count;
        printf("%d) Start check %ld\n", i, h_start_index[i]);
    } // end for i
    */
    /*
    SArray1D templen = SArray1D("templen", 1);
    auto h_templen = Kokkos::create_mirror_view(templen);
    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_for("ArrayLength", 1, KOKKOS_CLASS_LAMBDA(const int&) {
            templen(0) = start_index_(dim2_);
            //length_ = start_index_(dim2_);
        });
    #else
    templen_functor templen_execution_functor(templen);
    Kokkos::parallel_for("ArrayLength", 1, templen_execution_functor);
    #endif
    Kokkos::fence();
    deep_copy(h_templen, templen);
    length_ = h_templen(0);

    printf("Length %ld\n", length_);
    
    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_for("StartCheck", dim2_+1, KOKKOS_CLASS_LAMBDA(const int j) {
            printf("%d) Start %ld\n", j, start_index_(j));
        });
    #else
    stride_check_functor check_execution_functor;
    Kokkos::parallel_for("StartCheck", dim2_+1, check_execution_functor);
    #endif
    Kokkos::fence();
    */
    length_ = temp.length_;
    array_ = temp.length_;
    mystrides_ = temp.mystrides_;

    /*
        dim1_ = temp.dim1_;
        length_ = temp.length_;
        start_index_ = SArray1D("start_index_", dim1_ + 1);
        Kokkos::parallel_for("EqualOperator", dim1_+1, KOKKOS_CLASS_LAMBDA(const int j) {
                start_index_(j) = temp.start_index_(j);
            });
        //for (int j = 0; j < dim1_; j++) {
        //    start_index_(j) = temp.start_index_(j);
        //}
        array_ = TArray1D("array_", length_);
    */
  }
    
    return *this;
}

//return the stored Kokkos view
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
Kokkos::View<T*, Layout, ExecSpace, MemoryTraits> RaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::get_kokkos_view() {
    return array_;
}

// Destructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits, typename ILayout>
KOKKOS_INLINE_FUNCTION
RaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::~RaggedDownArrayKokkos() { }

////////////////////////////////////////////////////////////////////////////////
// End of RaggedDownArrayKokkos
////////////////////////////////////////////////////////////////////////////////

//11. DynamicRaggedRightArray
template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace, typename MemoryTraits = void>
class DynamicRaggedRightArrayKokkos {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace, MemoryTraits>;
    using SArray1D = Kokkos::View<size_t *,Layout, ExecSpace, MemoryTraits>;
    
private:
    // THIS WILL BE A GPU POINTER!
    SArray1D stride_;
    TArray1D array_;
    
    size_t dim1_;
    size_t dim2_;
    size_t length_;
    
public:
    // Default constructor
    DynamicRaggedRightArrayKokkos ();
    
    //--- 2D array access of a ragged right array ---
    
    // overload constructor
    DynamicRaggedRightArrayKokkos (size_t dim1, size_t dim2, const std::string& tag_string = DEFAULTSTRINGARRAY);
    
    // A method to return or set the stride size
    KOKKOS_INLINE_FUNCTION
    size_t& stride(size_t i) const;
    
    // A method to return the size
    KOKKOS_INLINE_FUNCTION
    size_t size() const;

    //return the view
    KOKKOS_INLINE_FUNCTION
    TArray1D get_kokkos_view();
    
    // Overload operator() to access data as array(i,j),
    // where i=[0:N-1], j=[stride(i)]
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;
    
    // Overload copy assignment operator
    KOKKOS_INLINE_FUNCTION
    DynamicRaggedRightArrayKokkos& operator= (const DynamicRaggedRightArrayKokkos &temp);
    
    //kokkos policy functors

    //functors for kokkos execution policies
    //set strides to a constant value
    class set_strides_functor{
      public:
      SArray1D functor_strides_;
      size_t init_stride_;
      set_strides_functor(size_t init_stride, SArray1D temp_strides_){
        init_stride_ = init_stride;
        functor_strides_ = temp_strides_;
      }
      KOKKOS_INLINE_FUNCTION void operator()(const int index) const {
        functor_strides_(index) = init_stride_;
      }
    };
    
    // Destructor
    KOKKOS_INLINE_FUNCTION
    ~DynamicRaggedRightArrayKokkos ();
};

//nothing
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DynamicRaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DynamicRaggedRightArrayKokkos () {
    dim1_ = dim2_ = length_ = 0;
}

// Overloaded constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DynamicRaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DynamicRaggedRightArrayKokkos (size_t dim1, size_t dim2, const std::string& tag_string) {
    // The dimensions of the array;
    dim1_  = dim1;
    dim2_  = dim2;
    length_ = dim1*dim2;
    
    std::string append_stride_string("strides");
    std::string append_array_string("array");
    std::string temp_copy_string = tag_string;
    std::string strides_tag_string = temp_copy_string.append(append_stride_string);
    temp_copy_string = tag_string;
    std::string array_tag_string = temp_copy_string.append(append_array_string);

    stride_ = SArray1D(strides_tag_string, dim1_);
    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_for("StridesInit", dim1_, KOKKOS_CLASS_LAMBDA(const int i) {
      strides_(i) = 0;
    });
    #else
    set_strides_functor execution_functor(0, stride_);
    Kokkos::parallel_for("StridesInit", dim1_,execution_functor);
    #endif

    //allocate view
    array_ = TArray1D(array_tag_string, length_);
}

// A method to set the stride size for row i
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t& DynamicRaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::stride(size_t i) const {
    return stride_(i);
}

//return size
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DynamicRaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::size() const{
    return length_;
}

// Overload operator() to access data as array(i,j),
// where i=[0:N-1], j=[0:stride(i)]
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DynamicRaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j) const {
    // Asserts
    assert(i < dim1_ && "i is out of dim1 bounds in DynamicRaggedRight");  // die if >= dim1
    assert(j < stride_(i) && "j is out of stride bounds in DynamicRaggedRight");  // die if >= dim2
    // Cannot assert on Kokkos View
    //assert(j < stride_[i] && "j is out of stride bounds in DynamicRaggedRight");  // die if >= stride
    
    return array_(j + i*dim2_);
}

//overload = operator
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DynamicRaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits>&
       DynamicRaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator= (const DynamicRaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits> &temp)
{
    
    if( this != &temp) {
        dim1_ = temp.dim1_;
        dim2_ = temp.dim2_;
        length_ = temp.length_;
        stride_ = temp.stride_;
        array_ = temp.array_;
        /*
        #ifdef HAVE_CLASS_LAMBDA
        Kokkos::parallel_for("StrideZeroOut", dim1_, KOKKOS_CLASS_LAMBDA(const int i) {
            stride_(i) = 0;
        });
        #else
        stride_zero_functor execution_functor;
        Kokkos::parallel_for("StrideZeroOut", dim1_, execution_functor);
        #endif
        */
    }
    
    return *this;
}

//return the stored Kokkos view
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
Kokkos::View<T*, Layout, ExecSpace, MemoryTraits> DynamicRaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::get_kokkos_view() {
    return array_;
}

// Destructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DynamicRaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::~DynamicRaggedRightArrayKokkos() {
}




//----end DynamicRaggedRightArray class definitions----


//12. DynamicRaggedDownArray

template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace, typename MemoryTraits = void>
class DynamicRaggedDownArrayKokkos {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace, MemoryTraits>;
    using SArray1D = Kokkos::View<size_t *,Layout, ExecSpace, MemoryTraits>;

private:
    SArray1D stride_;
    TArray1D array_;
    
    size_t dim1_;
    size_t dim2_;
    size_t length_;
    
public:
    // Default constructor
    DynamicRaggedDownArrayKokkos ();
    
    //--- 2D array access of a ragged right array ---
    
    // overload constructor
    DynamicRaggedDownArrayKokkos (size_t dim1, size_t dim2, const std::string& tag_string = DEFAULTSTRINGARRAY);
    
    // A method to return or set the stride size
    KOKKOS_INLINE_FUNCTION
    size_t& stride(size_t j) const;
    
    // A method to return the size
    KOKKOS_INLINE_FUNCTION
    size_t size() const;

    //return the view
    KOKKOS_INLINE_FUNCTION
    TArray1D get_kokkos_view();
    
    // Overload operator() to access data as array(i,j),
    // where i=[stride(j)], j=[0:N-1]
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;
    
    // Overload copy assignment operator
    KOKKOS_INLINE_FUNCTION
    DynamicRaggedDownArrayKokkos& operator= (const DynamicRaggedDownArrayKokkos &temp);

    //kokkos policy functors
    //set strides to 0 functor
    //set strides to a constant value
    class set_strides_functor{
      public:
      SArray1D functor_strides_;
      size_t init_stride_;
      set_strides_functor(size_t init_stride, SArray1D temp_strides_){
        init_stride_ = init_stride;
        functor_strides_ = temp_strides_;
      }
      KOKKOS_INLINE_FUNCTION void operator()(const int index) const {
        functor_strides_(index) = init_stride_;
      }
    };
    
    // Destructor
    KOKKOS_INLINE_FUNCTION
    ~DynamicRaggedDownArrayKokkos ();
};

//nothing
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DynamicRaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DynamicRaggedDownArrayKokkos () {
    dim1_ = dim2_ = length_ = 0;
}

// Overloaded constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DynamicRaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DynamicRaggedDownArrayKokkos (size_t dim1, size_t dim2, const std::string& tag_string) {
    // The dimensions of the array;
    dim1_  = dim1;
    dim2_  = dim2;
    length_ = dim1*dim2;

    std::string append_stride_string("strides");
    std::string append_array_string("array");
    std::string temp_copy_string = tag_string;
    std::string strides_tag_string = temp_copy_string.append(append_stride_string);
    temp_copy_string = tag_string;
    std::string array_tag_string = temp_copy_string.append(append_array_string);

    stride_ = SArray1D(strides_tag_string, dim2_);
    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_for("StridesInit", dim2_, KOKKOS_CLASS_LAMBDA(const int i) {
      strides_(i) = 0;
    });
    #else
    set_strides_functor execution_functor(0, stride_);
    Kokkos::parallel_for("StridesInit", dim2_,execution_functor);
    #endif

    //allocate view
    array_ = TArray1D(array_tag_string, length_);
}

// A method to set the stride size for column j
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t& DynamicRaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::stride(size_t j) const {
    return stride_(j);
}

//return size
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t DynamicRaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::size() const{
    return length_;
}

// overload operator () to access data as an array(i,j)
// Note: i = 0:stride(j), j = 0:N-1

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& DynamicRaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator()(size_t i, size_t j) const {
    // Asserts
    assert(j < dim2_ && "j is out of dim2 bounds in DynamicRaggedDownArrayKokkos");  // die if >= dim2
    assert(i < stride(j) && "i is out of stride bounds in DynamicRaggedDownArrayKokkos");  // die if >= stride(j)
    // Can't do this assert with a Kokkos View
    //assert(i < stride_[j] && "i is out of stride bounds in DynamicRaggedDownArrayKokkos");  // die if >= stride
    
    return array_(i + j*dim1_);
}

//overload = operator
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DynamicRaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits>&
  DynamicRaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator= (const DynamicRaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits> &temp)
{
    
    if( this != &temp) {
        dim1_ = temp.dim1_;
        dim2_ = temp.dim2_;
        length_ = temp.length_;
        stride_ = temp.stride_;
        array_ = temp.array_;
        /*
        #ifdef HAVE_CLASS_LAMBDA
        Kokkos::parallel_for("StrideZeroOut", dim2_, KOKKOS_CLASS_LAMBDA(const int j) {
            stride_(j) = 0;
        });
        #else
        stride_zero_functor execution_functor;
        Kokkos::parallel_for("StrideZeroOut", dim2_, execution_functor);
        #endif
        */
    }
    
    return *this;
}

//return the stored Kokkos view
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
Kokkos::View<T*, Layout, ExecSpace, MemoryTraits> DynamicRaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::get_kokkos_view() {
    return array_;
}

// Destructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
DynamicRaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::~DynamicRaggedDownArrayKokkos() {
}

////// END DynamicRaggedDownArrayKokkos

// KokkosCSRArray
template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace, typename MemoryTraits = void>
class CSRArrayKokkos {
   
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace, MemoryTraits>;
    using SArray1D = Kokkos::View<size_t*, Layout, ExecSpace, MemoryTraits>;

  private: // What ought to be private ?
    size_t dim1_, dim2_;
    size_t nnz_;
    TArray1D array_;
    SArray1D column_index_;
    SArray1D start_index_;
    TArray1D miss_;
  public:

    /**
     * @brief Construct a new Sparse Row Array Kokkos object
     *
     */
    CSRArrayKokkos();
    //CSRArray(CArray<T> data, CArray<T> col_ptrs, CArray<T> row_ptrs, size_t rows, size_t cols);

   CSRArrayKokkos(
               CArrayKokkos<T, Layout, ExecSpace, MemoryTraits> &array,
               CArrayKokkos<size_t, Layout, ExecSpace, MemoryTraits> &start_index,
               CArrayKokkos<size_t, Layout, ExecSpace, MemoryTraits> &colum_index,
               size_t dim1, size_t dim2, const std::string & tag_string = DEFAULTSTRINGARRAY);


    /**
     * @brief Constructor takes in dense matrix
     */
    //KOKKOS_INLINE_FUNCTION
    CSRArrayKokkos(const CArrayKokkos<T, Layout, ExecSpace, MemoryTraits> &dense, const size_t dim1, const  size_t dim2);
    
    void data_setup(const std::string& tag_string);
    /**
     * @brief Access method to A(i,j) returns a dummy value of 0 if value is not allocated
     *
     * @param i row
     * @param j column
     * @return KOKKOS_INLINE_FUNCTION&
     */
    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t j) const;
    
    /**
     * @brief Same functionality as operator. Included for compatibility with other matar types
     *
     * @param i
     * @param j
     * @return KOKKOS_INLINE_FUNCTION&
     */
    KOKKOS_INLINE_FUNCTION
    T& value(size_t i, size_t j) const;

    /**
     * @brief Copy operator
     *
     * @param temp
     */
    KOKKOS_INLINE_FUNCTION
    CSRArrayKokkos& operator=(const CSRArrayKokkos &temp);
    
  
    /**
     * @brief Pointer to start of array_ data
     *
     */
    KOKKOS_INLINE_FUNCTION
    T* pointer() const;
    
    /**
     * @brief Get the beginning of the start_index_ array
     *
     */
    KOKKOS_INLINE_FUNCTION
    size_t* get_starts() const;

     
    /**
     * @brief Number of columns
     *
     * @return KOKKOS_INLINE_FUNCTION
     */
    KOKKOS_INLINE_FUNCTION
    size_t dim2() const ;
    
    /**
     * @brief Number of rows
     *
     * @return KOKKOS_INLINE_FUNCTION
     */
    KOKKOS_INLINE_FUNCTION
    size_t dim1() const;
    
    /**
     * @brief iterator notation to access the non zero elements of row i. Returns pointer to first element in row i
     *
     * @param i row
     * @return KOKKOS_INLINE_FUNCTION*
     */
    KOKKOS_INLINE_FUNCTION
    T* begin(size_t i);
    
    /**
     * @brief Iteator notation to access the non zero elements of row i. Returns pointer first element of the next row
     *
     * @param i
     * @return KOKKOS_INLINE_FUNCTION*
     */
    KOKKOS_INLINE_FUNCTION
    T* end(size_t i);
    
    /**
     * @brief Get the size of row i. Same functionality as nnz(i) but included for compatiblity.
     *
     * @param i
     * @return KOKKOS_INLINE_FUNCTION
     */
    KOKKOS_INLINE_FUNCTION
    size_t stride(size_t i) const;


    /*
     * @brief get values from dense array
     */
    void from_dense(CArrayKokkos<T, Layout, ExecSpace, MemoryTraits> &starts,
                    CArrayKokkos<T, Layout, ExecSpace, MemoryTraits> &columns,
                    CArrayKokkos<T, Layout, ExecSpace, MemoryTraits> &array);

    /*
     * iterator for the raw data at row i
     * i.e. return the index each element is the index in the 1 array
     *  This as the use of providing a reasonable way to get the column
     * index and data value in the case you need both
    */
    KOKKOS_INLINE_FUNCTION
    size_t begin_index(size_t i) const;
    KOKKOS_INLINE_FUNCTION
    size_t end_index(size_t i) const;

    /**
     * @brief get the number of non zero elements in row i
     */
    KOKKOS_INLINE_FUNCTION
    size_t nnz(size_t i);
    
    /**
     * @brief get the total number of non zero elements
     */
    KOKKOS_INLINE_FUNCTION
    size_t nnz() const ;
   
    // Use the index into the 1d array to get what value is stored there and what is the corresponding row
    KOKKOS_INLINE_FUNCTION
    T& get_val_flat(size_t k) const;
    KOKKOS_INLINE_FUNCTION
    size_t get_col_flat(size_t k) const;
    // reverse map function from A(i,j) to what element of data/col_pt_ it corersponds to
    int flat_index(size_t i, size_t j);
    // Convertor
    
    // int toCSC(CArray<T> &data, CArray<size_t> &col_ptrs, CArray<size_t> &row_ptrs);

    void to_dense(CArrayKokkos<T,Layout, ExecSpace, MemoryTraits>& A);
    
    class init_start_indices_functor{
       public:
              SArray1D mystart_index_;
              init_start_indices_functor(SArray1D tempstart_index_){
                     mystart_index_ = tempstart_index_;
              }
              KOKKOS_INLINE_FUNCTION void operator()(const int index) const {
                      mystart_index_(index) = 0;
             }
    };
    //destructor
    KOKKOS_INLINE_FUNCTION
    ~CSRArrayKokkos();

   
};

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CSRArrayKokkos<T, Layout,ExecSpace, MemoryTraits>::CSRArrayKokkos() {
    dim1_ = dim2_ = nnz_ = 0;
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CSRArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::CSRArrayKokkos(
               CArrayKokkos<T, Layout, ExecSpace, MemoryTraits> &array,
               CArrayKokkos<size_t, Layout, ExecSpace, MemoryTraits> &start_index,
               CArrayKokkos<size_t, Layout, ExecSpace, MemoryTraits> &colum_index,
               size_t dim1, size_t dim2, const std::string & tag_string){
    dim1_ = dim1;
    dim2_ = dim2;
    start_index_ = start_index.get_kokkos_view();
    array_ = array.get_kokkos_view();
    column_index_ = colum_index.get_kokkos_view();
    nnz_ = colum_index.extent();
    miss_ = TArray1D("miss", 1);
}

/*
template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
CSRArrayKokkos<T,Layout, ExecSpace,MemoryTraits>::CSRArrayKokkos(const CArrayKokkos<T, Layout, ExecSpace, MemoryTraits> &dense, const size_t dim1, const size_t dim2){
    dim1_ = dim1;
    dim2_ = dim2;
    miss_ = TArray1D("miss",1);
    start_index_ = Kokkos::View<size_t*>("start indices", dim1 + 1);
    nnz_ = 0;
    

    start_index_(0) = 0;
    // TODO MAKE parallel
    for(size_t i = 0; i < dim1_; i++){
        start_index_(i+1) = 0;
        for(size_t j =0; j < dim2_; j++){
                if(dense(i,j) != 0){
                        start_index_(i+1) ++;
                        nnz_++;
                }
        }
   }

    
    for(size_t i = 1; i < dim1_ + 1; i++){
            start_index_(i) = start_index_[i] + start_index_[i-1];
    }
    
    column_index_ = Kokkos::View<size_t*>("column Indices", nnz_);
    array_ = Kokkos::View<T*>("array elements", nnz_);
    size_t next = 0 ;
    for(size_t i = 0; i < dim1_; i++){
            for(size_t j =0 ; j < dim2_; j++){
                    if(dense(i,j) != 0){
         //               column_index_(next) = j;
         //               array_(next) = dense(i,j);
                        next++;
                    }
            }
    }
}
*/

//setup start indices
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
void CSRArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::data_setup(const std::string& tag_string) {
    //allocate start indices
    std::string append_indices_string("start_indices");
    std::string temp_copy_string = tag_string;
    std::string start_index_tag_string = temp_copy_string.append(append_indices_string);
    temp_copy_string = tag_string;

    start_index_ = SArray1D(start_index_tag_string, dim1_ + 1);
    #ifdef HAVE_CLASS_LAMBDA
    Kokkos::parallel_for("StartValuesInit", dim1_+1, KOKKOS_CLASS_LAMBDA(const int i) {
      start_index_(i) = 0;
    });
    #else
    init_start_indices_functor execution_functor(start_index_);
    Kokkos::parallel_for("StartValuesInit", dim1_+1,execution_functor);
    #endif

}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CSRArrayKokkos<T, Layout, ExecSpace, MemoryTraits>::operator()(size_t i, size_t j) const {
    size_t row_start = start_index_[i];
    size_t row_end = start_index_[i+1];
    size_t k;
    for(k = 0; k < row_end - row_start; k++){
        if(column_index_[row_start + k] == j){
            return array_.data()[row_start + k];
        }
    }
    miss_[0] = (T) NULL;
    return miss_[0];
}


template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CSRArrayKokkos<T, Layout, ExecSpace, MemoryTraits>::value(size_t i, size_t j) const {
    size_t row_start = start_index_[i];
    size_t row_end = start_index_[i+1];
    size_t k;
    for(k = 0; k < row_end - row_start; k++){
        if(column_index_[row_start + k] == j){
            return array_.data()[row_start + k];
        }
    }
    miss_[0] = (T) NULL;
    return miss_[0];
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* CSRArrayKokkos<T, Layout, ExecSpace, MemoryTraits>::pointer() const{
    return array_.data();
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t* CSRArrayKokkos<T, Layout, ExecSpace, MemoryTraits>::get_starts() const {
    return start_index_.data();
}

template<typename T,typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
CSRArrayKokkos<T,Layout, ExecSpace, MemoryTraits>& CSRArrayKokkos<T, Layout, ExecSpace, MemoryTraits>::operator=(const CSRArrayKokkos<T, Layout,ExecSpace,MemoryTraits> &temp){
    if(this != temp) {
        nnz_ = temp.nnz_;
        dim1_ = temp.dim1_;
        dim2_ = temp.dim2_;
        
        start_index_ = temp.start_index_;
        column_index_ = temp.column_index_;
        array_ = temp.array_;
    }
    return *this;
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
void CSRArrayKokkos<T,Layout, ExecSpace, MemoryTraits>::to_dense(CArrayKokkos<T,Layout, ExecSpace, MemoryTraits>& A){
    size_t i,j;
    for(i = 0; i < dim1_; i++){
        for(j = 0; j < dim2_; j++){
            A(i,j) = (*this)(i,j);
        }
    }
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CSRArrayKokkos<T, Layout, ExecSpace, MemoryTraits>::stride(size_t i) const {
   assert(i <= dim1_ && "Index i out of bounds in CSRArray.stride()");
   return start_index_.data()[i+i] - start_index_.data()[i];
}


template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CSRArrayKokkos<T, Layout, ExecSpace, MemoryTraits>::dim2() const {
    return dim2_;
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CSRArrayKokkos<T, Layout, ExecSpace, MemoryTraits>::dim1() const{
    return dim1_;
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* CSRArrayKokkos<T, Layout, ExecSpace, MemoryTraits>::begin(size_t i){
    assert(i <= dim1_ && "i is out of bounds in CSRArray.begin()");
    size_t row_start = start_index_.data()[i];
    return &array_.data()[row_start];
}
template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* CSRArrayKokkos<T, Layout, ExecSpace, MemoryTraits>::end(size_t i){
    assert(i <= dim1_ && "i is out of bounds in CSRArray.end()");
    size_t row_start = start_index_.data()[i+1];
    return &array_.data()[row_start];
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CSRArrayKokkos<T, Layout, ExecSpace,MemoryTraits>::begin_index(size_t i) const{
    assert(i <= dim1_ && "i is out of bounds in CSRArray.begin_index()");
    return start_index_.data()[i];
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CSRArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::end_index(size_t i) const{
    assert(i <= dim1_ && "i is out of bounds in CSRArray.begin_index()");
    return start_index_.data()[i+1];
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CSRArrayKokkos<T, Layout, ExecSpace, MemoryTraits>::nnz() const{
    return nnz_;
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CSRArrayKokkos<T,Layout, ExecSpace,MemoryTraits>::nnz(size_t i){
    assert(i <= dim1_ && "Index i out of bounds in CSRArray.stride()");
    return start_index_.data()[i+1] - start_index_.data()[i];
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CSRArrayKokkos<T,Layout,ExecSpace, MemoryTraits>::get_val_flat(size_t k) const{
   assert(k < nnz_ && "Index k is out of bounds in CSRArray.get_val_flat()");
   return array_.data()[k];
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CSRArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::get_col_flat(size_t k) const{
    assert(k < nnz_ && "Index k is out of bounds in CSRArray.get_col_lat()");
    return column_index_.data()[k];
}


template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
int CSRArrayKokkos<T,Layout, ExecSpace, MemoryTraits>::flat_index(size_t i, size_t j){
    size_t k;
    size_t row_start = start_index_.data()[i];
    size_t row_end = start_index_.data()[i+1];
    for(k = 0; k < row_end - row_start; k++){
        if(column_index_.data()[row_start+k] == j){
            return row_start+k;
        }
    }
    return  -1;
}

//template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
//void CSRArrrayKokkos<T,Layout,ExecSpace, MemoryTraits>::from_dense(CArrayKokkos<T, Layout, ExecSpace, MemoryTraits> &starts,
//                    CArrayKokkos<T, Layout, ExecSpace, MemoryTraits> &columns,
//                    CArrayKokkos<T, Layout, ExecSpace, MemoryTraits> &array);

                    /*
int CSRArray<T>::toCSC(CArray<T> &data, CArray<size_t> &col_ptrs, CArray<size_t> &row_ptrs ){
    int nnz_cols[ncols_ + 1];
    int col_counts[ncols_];
    int i = 0;
    // How many elements are each column
    for(i =0 ; i < ncols_; i++){
        nnz_cols[i] = 0;
        col_counts[i] = 0;
    }
    nnz_cols[ncols_] = 0;
    col_ptrs(0) = 0;
    for(i =0; i < nnz_; i++){
        nnz_cols[column_index_[i] + 1] += 1;
    }
    // What we actually care about is how many elements are
    // in all the  columns preceeding this column.
    for(i = 1; i <= ncols_; i++){
        nnz_cols[i] = nnz_cols[i-1] + nnz_cols[i];
        col_ptrs(i) = nnz_cols[i];
    }
    size_t row = 1;
    // if b is at A(i,j)  stored in csr format
    // it needs to go where the where the ith column starts
    // + how many things we have put in the "window"
    // we allocated for this column already
    // For row we simply keep track of what row we are currently in
    // as we scan through the 1d array of data.
    for(i = 0; i < nnz_; i++){
        if(i >= start_index_[row]){
            row++;
        }
        int idx = nnz_cols[column_index_[i]] + col_counts[column_index_[i]];
        col_counts[column_index_[i]] += 1;
        data(idx) = array_[i];
        row_ptrs(idx) = row - 1;
    }
    // I return an int because I thought I might need to return an error code
    // Not sure that is true
    return 0;
}
*/

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CSRArrayKokkos<T, Layout, ExecSpace, MemoryTraits>::~CSRArrayKokkos() {}

// 16 CSCArrayKokkos
template <typename T, typename Layout = DefaultLayout, typename ExecSpace = DefaultExecSpace, typename MemoryTraits = void>
class CSCArrayKokkos
{

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace, MemoryTraits>;
    using SArray1D = Kokkos::View<size_t*, Layout, ExecSpace, MemoryTraits>;
private: // What ought to be private ?
    size_t dim1_, dim2_;
    size_t nnz_;
    TArray1D array_;
    TArray1D miss_;
    SArray1D start_index_;
    SArray1D row_index_;
    
  public:

      /**
       * @brief Construct a new empty Sparse Col Array object
       *
       */
      CSCArrayKokkos();

      /**
      * @brief Construct a new Sparse Col Array object
      *
      * @param array: 1d array of data values in order as read top to bottom, left to right
      * @param row_index: 1d array that marks what row each element is in
      * @param start_index: 1d array that marks where the first element of each column starts
      * @param dim1: number of rows the matrix should have
      * @param dim2: number of columns the matrix should have
      */
     CSCArrayKokkos(
               CArrayKokkos<T, Layout, ExecSpace, MemoryTraits> &array,
               CArrayKokkos<size_t, Layout, ExecSpace, MemoryTraits> &start_index,
               CArrayKokkos<size_t, Layout, ExecSpace, MemoryTraits> &row_index,
               size_t dim1, size_t dim2, const std::string & tag_string = DEFAULTSTRINGARRAY);


      /**
       * @brief Access A(i,j). Returns a dummy address with value 0 if A(i,j) is not allocated
       *
       * @param i : row
       * @param j : column
       * @return T& : address of array_ that corresponds to A(i,j)
       */
      KOKKOS_INLINE_FUNCTION
      T &operator()(size_t i, size_t j) const;

      /**
       * @brief Overloaded copy operator
       *
       * @param temp : Array to copy
       * @return CSCArray&
       */
      KOKKOS_INLINE_FUNCTION
      CSCArrayKokkos &operator=(const CSCArrayKokkos &temp);

      /**
       * @brief  returns pointer to array_
       *
       */
      KOKKOS_INLINE_FUNCTION
      T *pointer() const;

      /**
       * @brief Same functionality as nnz(i) included for compatibility with the rest of matar
       *
       * @param i : row
       * @return size_t
       */
      KOKKOS_INLINE_FUNCTION
      size_t stride(size_t i) const;
      
      /**
       * @brief Same functionality as operator() included for compatibility with the rest of matar
       *
       * @param i: row
       * @param j: column
       * @return T&
       */
      KOKKOS_INLINE_FUNCTION
      T &value(size_t i, size_t j) const;

      /**
       * @brief Get the start_index array
       *
       * @return size_t* : returns start_index_
       */
      KOKKOS_INLINE_FUNCTION
      size_t *get_starts() const;

      /**
       * @brief Get number of rows
       *
       * @return size_t  number of rows
       */
      KOKKOS_INLINE_FUNCTION
      size_t dim1() const;
      
      /**
       * @brief Get number of columns
       *
       * @return size_t number of columns
       */
      KOKKOS_INLINE_FUNCTION
      size_t dim2() const;

      /**
       * @brief iterator notation for iterating through the non zeros values of row i.
       *
       * @param i : row
       * @return T*
       */
      KOKKOS_INLINE_FUNCTION
      T *begin(size_t i);
 
       /**
       * @brief iterator notation for iterating through the non zeros values of row i.
       *
       * @param i : row
       * @return T*
       */
      KOKKOS_INLINE_FUNCTION
      T *end(size_t i);

      // iterator for the raw data at row i
      // i.e. return the index each element is the index in the 1 array
      // This as the use of providing a reasonable way to get the column
      // index and data value in the case you need both
      KOKKOS_INLINE_FUNCTION
      size_t begin_index(size_t i);

      KOKKOS_INLINE_FUNCTION
      size_t end_index(size_t i);

      /**
       * @brief Get the number of non zero elements in row i
       *
       * @param i : row to get
       * @return size_t  : size of row
       */
      KOKKOS_INLINE_FUNCTION
      size_t nnz(size_t i);
      
      /**
       * @brief Get number of non zero elements total in array
       *
       * @return size_t
       */
      KOKKOS_INLINE_FUNCTION
      size_t nnz() const;

      // Use the index into the 1d array to get what value is stored there and what is the corresponding row
      KOKKOS_INLINE_FUNCTION
      T &get_val_flat(size_t k);

      KOKKOS_INLINE_FUNCTION
      size_t get_row_flat(size_t k);
      
      // reverse map function from A(i,j) to what element of data/col_pt_ it corersponds to
      KOKKOS_INLINE_FUNCTION
      int flat_index(size_t i, size_t j);
      
      // Convertor
      //int toCSR(CArray<T> &data, CArray<size_t> &row_ptrs, CArray<size_t> &col_ptrs);
      //void to_dense(FArray<T> &A);
      // destructor
      KOKKOS_INLINE_FUNCTION
      ~CSCArrayKokkos();
};

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CSCArrayKokkos<T, Layout, ExecSpace, MemoryTraits>::CSCArrayKokkos() {
    dim1_ = dim2_ = nnz_ = 0;
}

 
template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CSCArrayKokkos<T, Layout, ExecSpace, MemoryTraits>::CSCArrayKokkos(
               CArrayKokkos<T, Layout, ExecSpace, MemoryTraits> &array,
               CArrayKokkos<size_t, Layout, ExecSpace, MemoryTraits> &start_index,
               CArrayKokkos<size_t, Layout, ExecSpace, MemoryTraits> &row_index,
               size_t dim1, size_t dim2, const std::string & tag_string){

    dim1_ = dim1;
    dim2_ = dim2;
    start_index_ = start_index.get_kokkos_view();
    array_ = array.get_kokkos_view();
    row_index_ = row_index.get_kokkos_view();
    nnz_ = row_index.extent();
    miss_ = TArray1D("miss", 1);
}


template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CSCArrayKokkos<T, Layout, ExecSpace, MemoryTraits>::operator()(size_t i, size_t j) const {
    size_t col_start = start_index_[j];
    size_t col_end = start_index_[j + 1];
    size_t k;
    for(k =0; k < col_end - col_start;k++){
        if(row_index_[col_start + k] == i){
                return array_.data()[col_start + k];
        }
    }
    //return array_.data()[nnz_];
    miss_[0] = (T) NULL;
    return miss_[0];
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* CSCArrayKokkos<T,Layout, ExecSpace, MemoryTraits>::pointer() const {
    return array_.data();
}


template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CSCArrayKokkos<T,Layout, ExecSpace, MemoryTraits>::value(size_t i, size_t j) const {
    size_t col_start = start_index_.data()[j];
    size_t col_end = start_index_.data()[j + 1];
    size_t k;
    for(k =0; k < col_end - col_start;k++){
        if(row_index_.data()[col_start + k] == i){
                return array_.data()[col_start + k];
        }
    }
    miss_[0] = (T) NULL;
    return miss_[0];
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t* CSCArrayKokkos<T,Layout, ExecSpace, MemoryTraits>::get_starts() const{
    return &start_index_.data()[0];
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
CSCArrayKokkos<T,Layout, ExecSpace, MemoryTraits>& CSCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::operator=(const CSCArrayKokkos<T,Layout,ExecSpace,MemoryTraits> &temp){
    if(this != temp) {
        nnz_ = temp.nnz_;
        dim2_ = temp.dim2_;
        dim1_ = temp.dim1_;;
        
        start_index_ = temp.start_index_;
        row_index_ = temp.row_row_index_;
        array_ = temp.array_;
    }
    return *this;
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CSCArrayKokkos<T,Layout, ExecSpace, MemoryTraits>::stride(size_t i) const{
    assert(i < dim2_ && "i is out of bounds in CSCArray.stride()");
    return start_index_.data()[i+1] - start_index_.data()[i];
}


template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CSCArrayKokkos<T,Layout, ExecSpace, MemoryTraits>::dim1() const {
    return dim1_;
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CSCArrayKokkos<T,Layout, ExecSpace, MemoryTraits>::dim2() const{
    return dim2_;
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* CSCArrayKokkos<T,Layout, ExecSpace, MemoryTraits>::begin(size_t i){
    assert(i <= dim2_ && "index i out of bounds at CSCArray.begin()");
    size_t col_start = start_index_.data()[i];
    return &array_.data()[col_start];
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T* CSCArrayKokkos<T,Layout, ExecSpace, MemoryTraits>::end(size_t i){
    assert(i <= dim2_ && "index i out of bounds at CSCArray.endt()");
    size_t col_start = start_index_.data()[i+1];
    return &array_.data()[col_start];
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CSCArrayKokkos<T,Layout, ExecSpace, MemoryTraits>::begin_index(size_t i){
    assert(i <= dim2_ && "index i out of bounds at CSCArray.begin_index()");
    return start_index_.data()[i];
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CSCArrayKokkos<T,Layout, ExecSpace, MemoryTraits>::end_index(size_t i){
    assert(i <= dim2_ && "index i out of bounds at CSCArray.end_index()");
    return start_index_.data()[i + 1];
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CSCArrayKokkos<T,Layout, ExecSpace, MemoryTraits>::nnz() const{
    return nnz_;
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CSCArrayKokkos<T,Layout, ExecSpace, MemoryTraits>::nnz(size_t i){
    return start_index_.data()[i+1] - start_index_.data()[i];
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
T& CSCArrayKokkos<T,Layout, ExecSpace, MemoryTraits>::get_val_flat(size_t k){
    return array_.data()[k];
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
size_t CSCArrayKokkos<T,Layout, ExecSpace, MemoryTraits>::get_row_flat(size_t k){
    return row_index_.data()[k];
}

template<typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
KOKKOS_INLINE_FUNCTION
int CSCArrayKokkos<T,Layout, ExecSpace, MemoryTraits>::flat_index(size_t i, size_t j){
    size_t col_start = start_index_.data()[j];
    size_t col_end = start_index_.data()[j+1];
    size_t k;
    for (k = 0; k < col_end - col_start; k++)
    {
        if(row_index_.data()[col_start + k] == i){
                return col_start + k;
        }
    }
    return  -1;
}

// Assumes that data, col_ptrs, and row_ptrs
// have been allocated size already before this call
// Returns the data in this csr format but as represented as the appropriatte vectors
// for a csc format
/*template<typename T>
int CSCArray<T>::toCSR(CArray<T> &data, CArray<size_t> &col_ptrs, CArray<size_t> &row_ptrs ){
    int nnz_rows[dim1_ + 1];
    int row_counts[dim1_];
    int i = 0;
    // How many elements are each column
    for(i =0 ; i < dim1_; i++){
        nnz_rows[i] = 0;
        row_counts[i] = 0;
    }
    nnz_rows[dim1_] = 0;
    row_ptrs(i) = 0 ;
    for(i =0; i < nnz_; i++){
        nnz_rows[row_index_[i] + 1] += 1;
    }
    // What we actually care about is how many elements are
    // in all the columns preceeding this column.
    for(i = 1; i < dim1_; i++){
        nnz_rows[i] = nnz_rows[i-1] + nnz_rows[i];
        row_ptrs(i) = nnz_rows[i];
    }
    size_t col = 1;
    // if b is at A(i,j)  stored in csr format
    // it needs to go where the where the ith column starts
    // + how many things we have put in the "window"
    // we allocated for this column already
    // For row we simply keep track of what row we are currently in
    // as we scan through the 1d array of data.
    for(i = 0; i < nnz_; i++){
        if(i >= start_index_[col]){
            col++;
        }
        int idx = nnz_rows[row_index_[i]] + row_counts[row_index_[i]];
        row_counts[row_index_[i]] += 1;
        data(idx) = array_[i];
        col_ptrs(idx) = col - 1;
    }
    // I return an int because I thought I might need to return an error code
    // Not sure that is true
    return 0;
}
*/

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
CSCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::~CSCArrayKokkos() {}

//////////////////////////
// Inherited Class Array
//////////////////////////

/*
//template<class T, class Layout, class ExecSpace>
template<typename T>
class InheritedArray2L {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;

private:
    size_t dim1_, length_;

public:
    TArray1D this_array_;
    typename Kokkos::View<T*, Layout, ExecSpace>::HostMirror  h_this_array_;

    InheritedArray2L();
    
    InheritedArray2L(size_t some_dim1);

    KOKKOS_INLINE_FUNCTION
    T& operator()(size_t i, size_t dest) const;

    template <typename U>
    void AllocateHost(size_t size, U *obj);

    void AllocateGPU();

    template <typename U, typename V>
    void InitModels(U *obj, V input);

    template <typename U>
    void ClearModels(U obj);

    InheritedArray2L& operator=(const InheritedArray2L& temp);

    // GPU Method
    // Method that returns size
    KOKKOS_INLINE_FUNCTION
    size_t size();

    // Host Method
    // Method that returns size
    size_t extent();

    // Methods returns the raw pointer (most likely GPU) of the Kokkos View
    T* pointer();

    // Deconstructor
    KOKKOS_INLINE_FUNCTION
    ~InheritedArray2L ();
}; // End of InheritedArray2L

// Default constructor
template <typename T>
InheritedArray2L<T>::InheritedArray2L() {}

// Overloaded 1D constructor
template <typename T>
InheritedArray2L<T>::InheritedArray2L(size_t some_dim1) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    length_ = dim1_;
    this_array_ = TArray1D("this_array_", length_);
    h_this_array_ = Kokkos::create_mirror_view(this_array_);
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& InheritedArray2L<T>::operator()(size_t i, size_t dest) const {
    assert(i < dim1_ && "i is out of bounds in InheritedArray2L 1D!");
    assert(dest < 2 && "dest is out of bounds in InheritedArray2L 1D!");
    if (dest == 0)
        return h_this_array_(i);
    else
        return this_array_(i);
}

template <typename T>
template <typename U>
void InheritedArray2L<T>::AllocateHost(size_t size, U *obj) {
    obj = (U *) kmalloc(size);
}

template <typename T>
void InheritedArray2L<T>::AllocateGPU() {
    Kokkos::deep_copy(this_array_, h_this_array_);
}

template <typename T>
template <typename U, typename V>
void InheritedArray2L<T>::InitModels(U *obj, V input) {
    Kokkos::parallel_for(
            "CreateObjects", 1, KOKKOS_CLASS_LAMBDA(const int&) {
                new ((V *)obj) V{input};
            });
}

template <typename T>
template <typename U>
void InheritedArray2L<T>::ClearModels(U obj) {
    Kokkos::parallel_for(
            "DestroyObjects", 1, KOKKOS_LAMBDA(const int&) {
              this_array_(0).obj->~U();
              this_array_(1).obj->~U();
            });
}

template <typename T>
InheritedArray2L<T>& InheritedArray2L<T>::operator= (const InheritedArray2L& temp) {
    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    // Do nothing if the assignment is of the form x = x
    if (this != &temp) {
        dim1_ = temp.dim1_;
        length_ = temp.length_;
        this_array_ = TArray1D("this_array_", length_);
    }
    
    return *this;
}

// Return size
template <typename T>
KOKKOS_INLINE_FUNCTION
size_t InheritedArray2L<T>::size() {
    return length_;
}

template <typename T>
size_t InheritedArray2L<T>::extent() {
    return length_;
}

template <typename T>
T* InheritedArray2L<T>::pointer() {
    return this_array_.data();
}

template <typename T>
KOKKOS_INLINE_FUNCTION
InheritedArray2L<T>::~InheritedArray2L() {}
*/

////////////////////////////////////////////////////////////////////////////////
// End of InheritedArray2L
////////////////////////////////////////////////////////////////////////////////


} // end namespace

#endif // end if have Kokkos


#endif // KOKKOS_TYPES_H
