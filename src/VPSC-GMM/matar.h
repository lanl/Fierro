#ifndef MATAR_H
#define MATAR_H
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

// Order
//
//  Standard (non-Kokkos data structures)
//   1. FArray
//   2. ViewFArray
//   3. FMatrix
//   4. ViewFMatrix
//   5. CArray
//   6. ViewCArray
//   7. CMatrix
//   8. ViewCMatrix
//   9. RaggedRightArray
//   10. RaggedDownArray
//   11. DynamicRaggedRightArray
//   12. DynamicRaggedDownArray
//   13. SparseRowArray
//   14. SparseColArray
//
//   Kokkos Data structures
//   15. FArrayKokkos
//   16. ViewFArrayKokkos
//   17. FMatrixKokkos
//   18. ViewFMatrixKokkos
//   19. CArrayKokkos
//   20. ViewCArrayKokkos
//   21. CMatrixKokkos
//   22. ViewCMatrixKokkos
//   23. RaggedRightArrayKokkos
//   24. RaggedDownArrayKokkos
//   25. DynamicRaggedRightArrayKokkos
//   26. DynamicRaggedDownArrayKokkos
//   27. SparseRowArrayKokkos
//   28. SparseColArrayKokkos


#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <assert.h>
#include <memory> // for shared_ptr
#include "macros.h"

using real_t = double;
using u_int  = unsigned int;


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
using DefaultMemSpace  = Kokkos::Serial;
using DefaultExecSpace = Kokkos::Serial;
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

using RMatrix1D    = Kokkos::View<real_t *,DefaultLayout,DefaultExecSpace>;
using RMatrix2D    = Kokkos::View<real_t **,DefaultLayout,DefaultExecSpace>;
using RMatrix3D    = Kokkos::View<real_t ***,DefaultLayout,DefaultExecSpace>;
using RMatrix4D    = Kokkos::View<real_t ****,DefaultLayout,DefaultExecSpace>;
using RMatrix5D    = Kokkos::View<real_t *****,DefaultLayout,DefaultExecSpace>;
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


//---Begin Standard Data Structures---

//1. FArray
// indicies are [0:N-1]
template <typename T>
class FArray {
    
private:
    size_t dims_[7];
    size_t length_;
    size_t order_;  // tensor order (rank)
    std::shared_ptr <T []> array_;
    
public:
    
    // default constructor
   FArray ();
   
    //overload constructors from 1D to 7D
     
   FArray(size_t dim0);
    
   FArray(size_t dim0,
          size_t dim1);
    
   FArray(size_t dim0,
          size_t dim1,
          size_t dim2);
    
   FArray(size_t dim0,
          size_t dim1,
          size_t dim2,
          size_t dim3);
    
   FArray(size_t dim0,
          size_t dim1,
          size_t dim2,
          size_t dim3,
          size_t dim4);

   FArray(size_t dim0,
          size_t dim1,
          size_t dim2,
          size_t dim3,
          size_t dim4,
          size_t dim5);

   FArray(size_t dim0,
          size_t dim1,
          size_t dim2,
          size_t dim3,
          size_t dim4,
          size_t dim5,
          size_t dim6);

    FArray (const FArray& temp);
    
    // overload operator() to access data as array(i,....,n);
    T& operator()(size_t i) const;
    
    T& operator()(size_t i,
                  size_t j) const;
    
    T& operator()(size_t i,
                  size_t j,
                  size_t k) const;
    
    T& operator()(size_t i,
                  size_t j,
                  size_t k,
                  size_t l) const;
    
    T& operator()(size_t i,
                  size_t j,
                  size_t k,
                  size_t l,
                  size_t m) const;
    
    T& operator()(size_t i,
                  size_t j,
                  size_t k,
                  size_t l,
                  size_t m,
                  size_t n) const;
    
    T& operator()(size_t i,
                  size_t j,
                  size_t k,
                  size_t l,
                  size_t m,
                  size_t n,
                  size_t o) const;
    
    //overload = operator
    FArray& operator=(const FArray& temp);
    
    //return array size
    size_t size() const;

    // return array dims
    size_t dims(size_t i) const;
    
    // return array order (rank)
    size_t order() const;
    
    //return pointer
    T* pointer() const;
    
    // deconstructor
    ~FArray ();
    
}; // end of f_array_t

//---FArray class definnitions----

//constructors
template <typename T>
FArray<T>::FArray(){
    array_ = NULL;
    length_ = 0;
}

//1D
template <typename T>
FArray<T>::FArray(size_t dim0)
{
    dims_[0] = dim0;
    length_ = dim0;
    order_ = 1;
    array_ = std::shared_ptr <T []> (new T[length_]);
}

template <typename T>
FArray<T>::FArray(size_t dim0,
                  size_t dim1)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    order_ = 2;
    length_ = dim0*dim1;
    array_ = std::shared_ptr <T []> (new T[length_]);
}

//3D
template <typename T>
FArray<T>::FArray(size_t dim0,
                  size_t dim1,
                  size_t dim2)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    order_ = 3;
    length_ = dim0*dim1*dim2;
    array_ = std::shared_ptr <T []> (new T[length_]);
}

//4D
template <typename T>
FArray<T>::FArray(size_t dim0,
                  size_t dim1,
                  size_t dim2,
                  size_t dim3)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    order_ = 4;
    length_ = dim0*dim1*dim2*dim3;
    array_ = std::shared_ptr <T []> (new T[length_]);
}

//5D
template <typename T>
FArray<T>::FArray(size_t dim0,
                  size_t dim1,
                  size_t dim2,
                  size_t dim3,
                  size_t dim4)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    order_ = 5;
    length_ = dim0*dim1*dim2*dim3*dim4;
    array_ = std::shared_ptr <T []> (new T[length_]);
}

//6D
template <typename T>
FArray<T>::FArray(size_t dim0,
                  size_t dim1,
                  size_t dim2,
                  size_t dim3,
                  size_t dim4,
                  size_t dim5)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    order_ = 6;
    length_ = dim0*dim1*dim2*dim3*dim4*dim5;
    array_ = std::shared_ptr <T []> (new T[length_]);
}


//7D
template <typename T>
FArray<T>::FArray(size_t dim0,
                  size_t dim1,
                  size_t dim2,
                  size_t dim3,
                  size_t dim4,
                  size_t dim5,
                  size_t dim6)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    dims_[6] = dim6;
    order_ = 7;
    length_ = dim0*dim1*dim2*dim3*dim4*dim5*dim6;
    array_ = std::shared_ptr <T []> (new T[length_]);
        
}

//Copy constructor

template <typename T>
FArray<T>::FArray(const FArray& temp) {
    
    // Do nothing if the assignment is of the form x = x
    
    if (this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for
        
        order_  = temp.order_;
        length_ = temp.length_;       
        array_ = temp.array_;
    } // end if
    
} // end constructor

//overload operator () for 1D to 7D
//indices are from [0:N-1]

//1D
template <typename T>
T& FArray<T>::operator()(size_t i) const
{
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in FArray 1D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in FArray 1D!");
    return array_[i];
}

//2D
template <typename T>
T& FArray<T>::operator()(size_t i,
                         size_t j) const
{
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in FArray 2D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in FArray 2D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in FArray 2D!");
    return array_[i + j*dims_[0]];
}

//3D
template <typename T>
T& FArray<T>::operator()(size_t i,
                         size_t j,
                         size_t k) const
{
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in FArray 3D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in FArray 3D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in Farray 3D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in FArray 3D!");
    return array_[i + j*dims_[0]
                    + k*dims_[0]*dims_[1]];
}

//4D
template <typename T>
T& FArray<T>::operator()(size_t i,
                         size_t j,
                         size_t k,
                         size_t l) const
{
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in FArray 4D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in FArray 4D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in FArray 4D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in FArray 4D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in FArray 4D!");
    return array_[i + j*dims_[0]
                    + k*dims_[0]*dims_[1]
                    + l*dims_[0]*dims_[1]*dims_[2]];
}

//5D
template <typename T>
T& FArray<T>::operator()(size_t i,
                         size_t j,
                         size_t k,
                         size_t l,
                         size_t m) const
{
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in FArray 5D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in FArray 5D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in FArray 5D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in FArray 5D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in FArray 5D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in FArray 5D!");
    return array_[i + j*dims_[0]
                    + k*dims_[0]*dims_[1]
                    + l*dims_[0]*dims_[1]*dims_[2]
                    + m*dims_[0]*dims_[1]*dims_[2]*dims_[3]];
}

//6D
template <typename T>
T& FArray<T>::operator()(size_t i,
                         size_t j,
                         size_t k,
                         size_t l,
                         size_t m,
                         size_t n) const
{
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in FArray 6D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in FArray 6D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in FArray 6D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in FArray 6D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in FArray 6D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in FArray 6D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in FArray 6D!");
    return array_[i + j*dims_[0]
                    + k*dims_[0]*dims_[1]
                    + l*dims_[0]*dims_[1]*dims_[2]
                    + m*dims_[0]*dims_[1]*dims_[2]*dims_[3]
                    + n*dims_[0]*dims_[1]*dims_[2]*dims_[3]*dims_[4]];
}

//7D
template <typename T>
T& FArray<T>::operator()(size_t i,
                         size_t j,
                         size_t k,
                         size_t l,
                         size_t m,
                         size_t n,
                         size_t o) const
{
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in FArray 7D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in FArray 7D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in FArray 7D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in FArray 7D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in FArray 7D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in FArray 7D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in FArray 7D!");
    assert(o >= 0 && o < dims_[6] && "o is out of bounds in FArray 7D!");
    return array_[i + j*dims_[0]
                    + k*dims_[0]*dims_[1]
                    + l*dims_[0]*dims_[1]*dims_[2]
                    + m*dims_[0]*dims_[1]*dims_[2]*dims_[3]
                    + n*dims_[0]*dims_[1]*dims_[2]*dims_[3]*dims_[4]
                    + o*dims_[0]*dims_[1]*dims_[2]*dims_[3]*dims_[4]*dims_[5]];
}
    
// = operator
//THIS = FArray <> TEMP(n,m,...)
template <typename T>
FArray<T>& FArray<T>::operator= (const FArray& temp)
{
    if(this != & temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for

        order_  = temp.order_;
        length_ = temp.length_;
        array_  = temp.array_;
    }
    return *this;
}

template <typename T>
inline size_t FArray<T>::size() const {
    return length_;
}

template <typename T>
inline size_t FArray<T>::dims(size_t i) const {
    assert(i < order_ && "FArray order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to FArray dims is out of bounds!");
    return dims_[i];
}

template <typename T>
inline size_t FArray<T>::order() const {
    return order_;
}


template <typename T>
inline T* FArray<T>::pointer() const {
    return array_.get();
}

//delete FArray
template <typename T>
FArray<T>::~FArray(){}

//---end of FArray class definitions----


//2. ViewFArray
// indicies are [0:N-1]
template <typename T>
class ViewFArray {

private:
    size_t dims_[7];
    size_t length_; // Length of 1D array
    size_t order_;  // tensor order (rank)
    T * array_;
    
public:
    
    // default constructor
    ViewFArray ();

    //---1D to 7D array ---
    ViewFArray(T *array,
               size_t dim0);
    
    ViewFArray (T *array,
                size_t dim0,
                size_t dim1);

    ViewFArray (T *array,
                size_t dim0,
                size_t dim1,
                size_t dim2);

    ViewFArray (T *array,
                size_t dim0,
                size_t dim1,
                size_t dim2,
                size_t dim3);
    
    ViewFArray (T *array,
                size_t dim0,
                size_t dim1,
                size_t dim2,
                size_t dim3,
                size_t dim4);

    ViewFArray (T *array,
                size_t dim0,
                size_t dim1,
                size_t dim2,
                size_t dim3,
                size_t dim4,
                size_t dim5);
    
    ViewFArray (T *array,
                size_t dim0,
                size_t dim1,
                size_t dim2,
                size_t dim3,
                size_t dim4,
                size_t dim5,
                size_t dim6);
    
    T& operator()(size_t i) const;
    
    T& operator()(size_t i,
                  size_t j) const;
    
    T& operator()(size_t i,
                  size_t j,
                  size_t k) const;
    
    T& operator()(size_t i,
                  size_t j,
                  size_t k,
                  size_t l) const;
    
    T& operator()(size_t i,
                  size_t j,
                  size_t k,
                  size_t l,
                  size_t m) const;
    
    T& operator()(size_t i,
                  size_t j,
                  size_t k,
                  size_t l,
                  size_t m,
                  size_t n) const;
    
    T& operator()(size_t i,
                  size_t j,
                  size_t k,
                  size_t l,
                  size_t m,
                  size_t n,
                  size_t o) const;
    
    // calculate C = math(A,B)
    template <typename M>
    void operator=(M do_this_math);
    
    //return array size
    size_t size() const;
    
    //return array dims
    size_t dims(size_t i) const;
    
    // return array order (rank)
    size_t order() const;

    // return pointer
    T* pointer() const;
    
}; // end of viewFArray

//class definitions for viewFArray

//~~~~constructors for viewFArray for 1D to 7D~~~~~~~

//no dimension
template <typename T>
ViewFArray<T>::ViewFArray(){
  array_ = NULL;
  length_ = 0;
}

//1D
template <typename T>
ViewFArray<T>::ViewFArray(T *array,
                          size_t dim0)
{
    dims_[0] = dim0;
    order_ = 1;
    length_ = dim0;
    array_  = array;
}

//2D
template <typename T>
ViewFArray<T>::ViewFArray(T *array,
                          size_t dim0,
                          size_t dim1)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    order_ = 2;
    length_ = dim0*dim1;
    array_  = array;
}

//3D
template <typename T>
ViewFArray<T>::ViewFArray(T *array,
                          size_t dim0,
                          size_t dim1,
                          size_t dim2)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    order_ = 3;
    length_ = dim0*dim1*dim2;
    array_  = array;
}

//4D
template <typename T>
ViewFArray<T>::ViewFArray(T *array,
                          size_t dim0,
                          size_t dim1,
                          size_t dim2,
                          size_t dim3)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    order_ = 4;
    length_ = dim0*dim1*dim2*dim3;
    array_  = array;
}

//5D
template <typename T>
ViewFArray<T>::ViewFArray(T *array,
                          size_t dim0,
                          size_t dim1,
                          size_t dim2,
                          size_t dim3,
                          size_t dim4)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    order_ = 5;
    length_ = dim0*dim1*dim2*dim3*dim4;
    array_  = array;
}

//6D
template <typename T>
ViewFArray<T>::ViewFArray(T *array,
                          size_t dim0,
                          size_t dim1,
                          size_t dim2,
                          size_t dim3,
                          size_t dim4,
                          size_t dim5)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    order_ = 6;
    length_ = dim0*dim1*dim2*dim3*dim4*dim5;
    array_  = array;
}

//7D
template <typename T>
ViewFArray<T>::ViewFArray(T *array,
                          size_t dim0,
                          size_t dim1,
                          size_t dim2,
                          size_t dim3,
                          size_t dim4,
                          size_t dim5,
                          size_t dim6)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    dims_[6] = dim6;
    order_ = 7;
    length_ = dim0*dim1*dim2*dim3*dim4*dim5*dim6;
    array_  = array;
}

//~~~~~~operator () overload 
//for dimensions 1D to 7D
//indices for array are from 0...N-1

//1D
template <typename T>
T& ViewFArray<T>::operator()(size_t i) const
{
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in ViewFArray 1D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewFArray 1D!");
    return array_[i];
}

//2D
template <typename T>
T& ViewFArray<T>::operator()(size_t i,
                             size_t j) const
{
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in ViewFArray 2D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewFArray 2D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewFArray 2D!");
    return array_[i + j*dims_[0]];
}

//3D
template <typename T>
T& ViewFArray<T>::operator()(size_t i,
                             size_t j,
                             size_t k) const
{
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in ViewFArray 3D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewFArray 3D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewFArray 3D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewFArray 3D!");
    return array_[i + j*dims_[0]
                    + k*dims_[0]*dims_[1]];
}

//4D
template <typename T>
T& ViewFArray<T>::operator()(size_t i,
                             size_t j,
                             size_t k,
                             size_t l) const
{
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in ViewFArray 4D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewFArray 4D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewFArray 4D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewFArray 4D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in ViewFArray 4D!");
    return array_[i + j*dims_[0]
                    + k*dims_[0]*dims_[1]
                    + l*dims_[0]*dims_[1]*dims_[2]];
}

//5D
template <typename T>
T& ViewFArray<T>::operator()(size_t i,
                             size_t j,
                             size_t k,
                             size_t l,
                             size_t m) const
{
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in ViewFArray 5D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewFArray 5D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewFArray 5D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewFArray 5D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in ViewFArray 5D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in ViewFArray 5D!");
    return array_[i + j*dims_[0]
                    + k*dims_[0]*dims_[1]
                    + l*dims_[0]*dims_[1]*dims_[2]
                    + m*dims_[0]*dims_[1]*dims_[2]*dims_[3]];
}

//6D
template <typename T>
T& ViewFArray<T>:: operator()(size_t i,
                              size_t j,
                              size_t k,
                              size_t l,
                              size_t m,
                              size_t n) const
{
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in ViewFArray 6D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewFArray 6D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewFArray 6D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewFArray 6D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in ViewFArray 6D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in ViewFArray 6D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in ViewFArray 6D!");
    return array_[i + j*dims_[0]
                    + k*dims_[0]*dims_[1]
                    + l*dims_[0]*dims_[1]*dims_[2]
                    + m*dims_[0]*dims_[1]*dims_[2]*dims_[3]
                    + n*dims_[0]*dims_[1]*dims_[2]*dims_[3]*dims_[4]];
}

//7D
template <typename T>
T& ViewFArray<T>:: operator()(size_t i,
                              size_t j,
                              size_t k,
                              size_t l,
                              size_t m,
                              size_t n,
                              size_t o) const
{
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in ViewFArray 7D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewFArray 7D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewFArray 7D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewFArray 7D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in ViewFArray 7D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in ViewFArray 7D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in ViewFArray 7D!");
    assert(o >= 0 && o < dims_[6] && "n is out of bounds in ViewFArray 7D!");
    return array_[i + j*dims_[0]
                    + k*dims_[0]*dims_[1]
                    + l*dims_[0]*dims_[1]*dims_[2]
                    + m*dims_[0]*dims_[1]*dims_[2]*dims_[3]
                    + n*dims_[0]*dims_[1]*dims_[2]*dims_[3]*dims_[4]
                    + o*dims_[0]*dims_[1]*dims_[2]*dims_[3]*dims_[4]*dims_[5]];
}

// calculate this ViewFArray object = math(A,B)
template <typename T>
template <typename M>
void ViewFArray<T>::operator=(M do_this_math){
    do_this_math(*this); // pass in this ViewFArray object
}// end of math opperation

template <typename T>
inline size_t ViewFArray<T>::dims(size_t i) const {
    assert(i < order_ && "ViewFArray order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to ViewFArray dims is out of bounds!");
    return dims_[i];
}

template <typename T>
inline size_t ViewFArray<T>::order() const {
    return order_;
}

template <typename T>
inline size_t ViewFArray<T>::size() const {
    return length_;
}

template <typename T>
inline T* ViewFArray<T>::pointer() const {
    return array_;
}

//---end of ViewFArray class definitions---


//3. FMatrix
// indicies are [1:N]
template <typename T>
class FMatrix {
private:
    size_t dims_[7];
    size_t length_; // Length of 1D array
    size_t order_;  // tensor order (rank)
    std::shared_ptr <T []> matrix_;

public:
    // Default constructor
    FMatrix ();

    //---1D to 7D matrix ---
    FMatrix (size_t dim1);

    FMatrix (size_t dim1,
             size_t dim2);

    FMatrix (size_t dim1,
             size_t dim2,
             size_t dim3);

    FMatrix (size_t dim1,
             size_t dim2,
             size_t dim3,
             size_t dim4);

    FMatrix (size_t dim1,
             size_t dim2,
             size_t dim3,
             size_t dim4,
             size_t dim5);

    FMatrix (size_t dim1,
             size_t dim2,
             size_t dim3,
             size_t dim4,
             size_t dim5,
             size_t dim6);

    FMatrix (size_t dim1,
             size_t dim2,
             size_t dim3,
             size_t dim4,
             size_t dim5,
             size_t dim6,
             size_t dim7);
    
    FMatrix (const FMatrix& temp);
    
    T& operator() (size_t i) const;
    
    T& operator() (size_t i,
                   size_t j) const;
    
    T& operator() (size_t i,
                   size_t j,
                   size_t k) const;
    
    T& operator() (size_t i,
                   size_t j,
                   size_t k,
                   size_t l) const;
    
    T& operator() (size_t i,
                   size_t j,
                   size_t k,
                   size_t l,
                   size_t m) const;
    
    T& operator() (size_t i,
                   size_t j,
                   size_t k,
                   size_t l,
                   size_t m,
                   size_t n) const;

    T& operator() (size_t i,
                   size_t j,
                   size_t k,
                   size_t l,
                   size_t m,
                   size_t n,
                   size_t o) const;
    
    
    // Overload copy assignment operator
    FMatrix& operator=(const FMatrix& temp);

    // the length of the 1D storage array
    size_t size() const;

    // matrix dims
    size_t dims(size_t i) const;
    
    // return matrix order (rank)
    size_t order() const;
    
    //return pointer
    T* pointer() const;

    // Deconstructor
    ~FMatrix ();

}; // End of FMatrix

//---FMatrix class definitions---

//constructors
template <typename T>
FMatrix<T>::FMatrix(){
    matrix_ = NULL;
    length_ = 0;
}

//1D
template <typename T>
FMatrix<T>::FMatrix(size_t dim1)
{
    dims_[0] = dim1;
    order_ = 1;
    length_ = dim1;
    matrix_ = std::shared_ptr <T []> (new T[length_]);
}

//2D
template <typename T>
FMatrix<T>::FMatrix(size_t dim1,
                    size_t dim2)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    order_ = 2;
    length_ = dim1 * dim2;
    matrix_ = std::shared_ptr <T []> (new T[length_]);
}

//3D
template <typename T>
FMatrix<T>::FMatrix(size_t dim1,
                    size_t dim2,
                    size_t dim3)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    order_ = 3;
    length_ = dim1 * dim2 * dim3;
    matrix_ = std::shared_ptr <T []> (new T[length_]);
}

//4D
template <typename T>
FMatrix<T>::FMatrix(size_t dim1,
                    size_t dim2,
                    size_t dim3,
                    size_t dim4)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    order_ = 4;
    length_ = dim1 * dim2 * dim3 * dim4;
    matrix_ = std::shared_ptr <T []> (new T[length_]);
}

//5D
template <typename T>
FMatrix<T>::FMatrix(size_t dim1,
                    size_t dim2,
                    size_t dim3,
                    size_t dim4,
                    size_t dim5)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    order_ = 5;
    length_ = dim1 * dim2 * dim3 * dim4 * dim5;
    matrix_ = std::shared_ptr <T []> (new T[length_]);
}

//6D
template <typename T>
FMatrix<T>::FMatrix(size_t dim1,
                    size_t dim2,
                    size_t dim3,
                    size_t dim4,
                    size_t dim5,
                    size_t dim6)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    order_ = 6;
    length_ = dim1 * dim2 * dim3 * dim4 * dim5 * dim6;
    matrix_ = std::shared_ptr <T []> (new T[length_]);

}

template <typename T>
FMatrix<T>::FMatrix(size_t dim1,
                    size_t dim2,
                    size_t dim3,
                    size_t dim4,
                    size_t dim5,
                    size_t dim6,
                    size_t dim7)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    dims_[6] = dim7;
    order_ = 7;
    length_ = dim1 * dim2 * dim3 * dim4 * dim5 * dim6 * dim7;
    matrix_ = std::shared_ptr <T []> (new T[length_]);
    
}

template <typename T>
FMatrix<T>::FMatrix(const FMatrix& temp) {
    
    // Do nothing if the assignment is of the form x = x
    
    if (this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for
        
        order_  = temp.order_;
        length_ = temp.length_;
        matrix_ = temp.matrix_;
    } // end if
    
} // end constructor


//overload operators

//1D
template <typename T>
inline T& FMatrix<T>::operator() (size_t i) const
{
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in FMatrix 1D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in FMatrix 1D!");
    return matrix_[i - 1];
}

//2D
template <typename T>
inline T& FMatrix<T>::operator() (size_t i,
                                  size_t j) const
{
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in FMatrix 2D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in FMatrix 2D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in FMatrix 2D!");
    return matrix_[(i - 1) + ((j - 1) * dims_[0])];
}

//3D
template <typename T>
inline T& FMatrix<T>::operator() (size_t i,
                                  size_t j,
                                  size_t k) const
{
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in FMatrix 3D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in FMatrix 3D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in FMatrix 3D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in FMatrix 3D!");
    return matrix_[(i - 1) + ((j - 1) * dims_[0])
                           + ((k - 1) * dims_[0] * dims_[1])];
}

//4D
template <typename T>
inline T& FMatrix<T>::operator() (size_t i,
                                  size_t j,
                                  size_t k,
                                  size_t l) const
{
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in FMatrix 4D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in FMatrix 4D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in FMatrix 4D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in FMatrix 4D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in FMatrix 4D!");
    return matrix_[(i - 1) + ((j - 1) * dims_[0])
                           + ((k - 1) * dims_[0] * dims_[1])
                           + ((l - 1) * dims_[0] * dims_[1] * dims_[2])];
}

//5D
template <typename T>
inline T& FMatrix<T>::operator() (size_t i,
                                  size_t j,
                                  size_t k,
                                  size_t l,
                                  size_t m) const
{
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in FMatrix 5D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in FMatrix 5D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in FMatrix 5D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in FMatrix 5D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in FMatrix 5D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in FMatrix 5D!");
    return matrix_[(i - 1) + ((j - 1) * dims_[0])
                           + ((k - 1) * dims_[0] * dims_[1])
                           + ((l - 1) * dims_[0] * dims_[1] * dims_[2])
                           + ((m - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3])];
}

//6D
template <typename T>
inline T& FMatrix<T>::operator() (size_t i,
                                  size_t j,
                                  size_t k,
                                  size_t l,
                                  size_t m,
                                  size_t n) const
{
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in FMatrix 6D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in FMatrix 6D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in FMatrix 6D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in FMatrix 6D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in FMatrix 6D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in FMatrix 6D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in FMatrix 6D!");
    return matrix_[(i - 1) + ((j - 1) * dims_[0])
                           + ((k - 1) * dims_[0] * dims_[1])
                           + ((l - 1) * dims_[0] * dims_[1] * dims_[2])
                           + ((m - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                           + ((n - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4])];
}

//7D
template <typename T>
inline T& FMatrix<T>::operator() (size_t i,
                                  size_t j,
                                  size_t k,
                                  size_t l,
                                  size_t m,
                                  size_t n,
                                  size_t o) const
{
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in FMatrix 7D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in FMatrix 7D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in FMatrix 7D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in FMatrix 7D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in FMatrix 7D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in FMatrix 7D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in FMatrix 7D!");
    assert(o >= 1 && o <= dims_[6] && "o is out of bounds in FMatrix 7D!");
    return matrix_[(i - 1) + ((j - 1) * dims_[0])
                           + ((k - 1) * dims_[0] * dims_[1])
                           + ((l - 1) * dims_[0] * dims_[1] * dims_[2])
                           + ((m - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                           + ((n - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4])
                           + ((o - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4] * dims_[5])];
}


template <typename T>
inline FMatrix<T>& FMatrix<T>::operator= (const FMatrix& temp)
{
    // Do nothing if assignment is of the form x = x
    if (this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for

        order_  = temp.order_;
        length_ = temp.length_;
	matrix_ = temp.matrix_;
    }
    
    return *this;
}

template <typename T>
inline size_t FMatrix<T>::size() const {
    return length_;
}

template <typename T>
inline size_t FMatrix<T>::dims(size_t i) const {
    i--; // i starts at 1
    assert(i < order_ && "FMatrix order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to FMatrix dims is out of bounds!");
    return dims_[i];
}

template <typename T>
inline size_t FMatrix<T>::order() const {
    return order_;
}

template <typename T>
inline T* FMatrix<T>::pointer() const{
    return matrix_.get();
}

template <typename T>
FMatrix<T>::~FMatrix() {}

//----end of FMatrix class definitions----


//4. ViewFMatrix
//  indices are [1:N]
template <typename T>
class ViewFMatrix {

private:
    size_t dims_[7];
    size_t length_; // Length of 1D array
    size_t order_;  // tensor order (rank)
    T * matrix_;
    
public:
    
    // Default constructor
    ViewFMatrix ();
    
    //--- 1D to 7D matrix ---

    ViewFMatrix(T *matrix,
                size_t dim1);
    
    ViewFMatrix(T *some_matrix,
                size_t dim1,
                size_t dim2);
    
    ViewFMatrix(T *matrix,
                size_t dim1,
                size_t dim2,
                size_t dim3);
    
    ViewFMatrix(T *matrix,
                size_t dim1,
                size_t dim2,
                size_t dim3,
                size_t dim4);
    
    ViewFMatrix (T *matrix,
                 size_t dim1,
                 size_t dim2,
                 size_t dim3,
                 size_t dim4,
                 size_t dim5);
    
    ViewFMatrix (T *matrix,
                 size_t dim1,
                 size_t dim2,
                 size_t dim3,
                 size_t dim4,
                 size_t dim5,
                 size_t dim6);
    
    ViewFMatrix (T *matrix,
                 size_t dim1,
                 size_t dim2,
                 size_t dim3,
                 size_t dim4,
                 size_t dim5,
                 size_t dim6,
                 size_t dim7);
    
    T& operator()(size_t i) const;
    
    T& operator()(size_t i,
                  size_t j) const;
    
    T& operator()(size_t i,
                  size_t j,
                  size_t k) const;
    
    T& operator()(size_t i,
                  size_t j,
                  size_t k,
                  size_t l) const;
    
    T& operator() (size_t i,
                   size_t j,
                   size_t k,
                   size_t l,
                   size_t m) const;
    
    T& operator()(size_t i, 
                  size_t j, 
                  size_t k, 
                  size_t l, 
                  size_t m, 
                  size_t n) const;

    T& operator()(size_t i,
                  size_t j,
                  size_t k,
                  size_t l,
                  size_t m,
                  size_t n,
                  size_t o) const;
    
    // calculate C = math(A,B)
    template <typename M>
    void operator=(M do_this_math);
    
    // length of 1D array
    size_t size() const;
    
    // matrix dims
    size_t dims(size_t i) const;
    
    // return matrix order (rank)
    size_t order() const;

    // return pointer
    T* pointer() const;
    
}; // end of ViewFMatrix

//constructors

//no dimension
template <typename T>
ViewFMatrix<T>::ViewFMatrix() {
  matrix_ = NULL;
  length_ = 0;
}

//1D
template <typename T>
ViewFMatrix<T>::ViewFMatrix(T *matrix,
                            size_t dim1)
{
    dims_[0] = dim1;
    order_ = 1;
    length_ = dim1;
    matrix_ = matrix;
}

//2D
template <typename T>
ViewFMatrix<T>::ViewFMatrix(T *matrix,
                            size_t dim1,
                            size_t dim2)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    order_ = 2;
    length_ = dim1 * dim2;
    matrix_ = matrix;
}

//3D
template <typename T>
ViewFMatrix<T>::ViewFMatrix (T *matrix,
                             size_t dim1,
                             size_t dim2,
                             size_t dim3)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    order_ = 3;
    length_ = dim1 * dim2 * dim3;
    matrix_ = matrix;
}

//4D
template <typename T>
ViewFMatrix<T>::ViewFMatrix(T *matrix,
                            size_t dim1,
                            size_t dim2,
                            size_t dim3,
                            size_t dim4)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    order_ = 4;
    length_ = dim1 * dim2 * dim3 * dim4;
    matrix_ = matrix;
}

//5D
template <typename T>
ViewFMatrix<T>::ViewFMatrix(T *matrix,
                            size_t dim1,
                            size_t dim2,
                            size_t dim3,
                            size_t dim4,
                            size_t dim5)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    order_ = 5;
    length_ = dim1 * dim2 * dim3 * dim4 * dim5;
    matrix_ = matrix;
}

//6D
template <typename T>
ViewFMatrix<T>::ViewFMatrix(T *matrix,
                            size_t dim1,
                            size_t dim2,
                            size_t dim3,
                            size_t dim4,
                            size_t dim5,
                            size_t dim6)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    order_ = 6;
    length_ = dim1 * dim2 * dim3 * dim4 * dim5 * dim6;
    matrix_ = matrix;
}

//6D
template <typename T>
ViewFMatrix<T>::ViewFMatrix(T *matrix,
                            size_t dim1,
                            size_t dim2,
                            size_t dim3,
                            size_t dim4,
                            size_t dim5,
                            size_t dim6,
                            size_t dim7)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    dims_[6] = dim7;
    order_ = 7;
    length_ = dim1 * dim2 * dim3 * dim4 * dim5 * dim6 * dim7;
    matrix_ = matrix;
}


//overload operator ()

//1D
template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i) const
{
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in ViewFMatrix 1D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewFMatrix 1D");  // die if >= dim1
        
    return matrix_[(i - 1)];
}

//2D
template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i, 
                                     size_t j) const
{
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in ViewFMatrix 2D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewFMatrix 2D");  // die if >= dim1
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewFMatrix 2D");  // die if >= dim2
        
    return matrix_[(i - 1) + ((j - 1) * dims_[0])];
}

//3D
template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i, 
                                     size_t j, 
                                     size_t k) const
{
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in ViewFMatrix 3D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewFMatrix 3D");  // die if >= dim1
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewFMatrix 3D");  // die if >= dim2
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in ViewFMatrix 3D");  // die if >= dim3
        
    return matrix_[(i - 1) + ((j - 1) * dims_[0])
                           + ((k - 1) * dims_[0] * dims_[1])];
}

//4D
template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i, 
                                     size_t j, 
                                     size_t k, 
                                     size_t l) const
{
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in ViewFMatrix 4D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewFMatrix 4D");  // die if >= dim1
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewFMatrix 4D");  // die if >= dim2
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in ViewFMatrix 4D");  // die if >= dim3
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in ViewFMatrix 4D");  // die if >= dim4
        
    return matrix_[(i - 1) + ((j - 1) * dims_[0])
                           + ((k - 1) * dims_[0] * dims_[1])
                           + ((l - 1) * dims_[0] * dims_[1] * dims_[2])];
}

//5D
template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i, 
                                     size_t j, 
                                     size_t k, 
                                     size_t l, 
                                     size_t m) const
{
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in ViewFMatrix 5D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewFMatrix 5D");  // die if >= dim1
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewFMatrix 5D");  // die if >= dim2
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in ViewFMatrix 5D");  // die if >= dim3
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in ViewFMatrix 5D");  // die if >= dim4
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in ViewFMatrix 5D");  // die if >= dim5
       
    return matrix_[(i - 1) + ((j - 1) * dims_[0])
                           + ((k - 1) * dims_[0] * dims_[1])
                           + ((l - 1) * dims_[0] * dims_[1] * dims_[2])
                           + ((m - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3])];
}

//6D
template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i,
                                     size_t j,
                                     size_t k,
                                     size_t l,
                                     size_t m,
                                     size_t n) const
{
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in ViewFMatrix 6D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewFMatrix 6D");  // die if >= dim1
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewFMatrix 6D");  // die if >= dim2
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in ViewFMatrix 6D");  // die if >= dim3
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in ViewFMatrix 6D");  // die if >= dim4
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in ViewFMatrix 6D");  // die if >= dim5
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in ViewFMatrix 6D");  // die if >= dim6
    return matrix_[(i - 1) + ((j - 1) * dims_[0])
                           + ((k - 1) * dims_[0] * dims_[1])
                           + ((l - 1) * dims_[0] * dims_[1] * dims_[2])
                           + ((m - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                           + ((n - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4])];
}

//6D
template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i,
                                     size_t j,
                                     size_t k,
                                     size_t l,
                                     size_t m,
                                     size_t n,
                                     size_t o) const
{
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in ViewFMatrix 7D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewFMatrix 7D");  // die if >= dim1
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewFMatrix 7D");  // die if >= dim2
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in ViewFMatrix 7D");  // die if >= dim3
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in ViewFMatrix 7D");  // die if >= dim4
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in ViewFMatrix 7D");  // die if >= dim5
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in ViewFMatrix 7D");  // die if >= dim6
    assert(o >= 1 && o <= dims_[6] && "o is out of bounds in ViewFMatrix 7D");  // die if >= dim7
    
    return matrix_[(i - 1) + ((j - 1) * dims_[0])
                           + ((k - 1) * dims_[0] * dims_[1])
                           + ((l - 1) * dims_[0] * dims_[1] * dims_[2])
                           + ((m - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3])
                           + ((n - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4])
                           + ((o - 1) * dims_[0] * dims_[1] * dims_[2] * dims_[3] * dims_[4] * dims_[5])];
}

// calculate this ViewFMatrix object = math(A,B)
template <typename T>
template <typename M>
void ViewFMatrix<T>::operator=(M do_this_math){
    do_this_math(*this); // pass in this ViewFArray object
}// end of math opperation

template <typename T>
inline size_t ViewFMatrix<T>::dims(size_t i) const {
    i--; // i starts at 1
    assert(i < order_ && "ViewFMatrix order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to ViewFMatrix dims is out of bounds!");
    return dims_[i];
}

template <typename T>
inline size_t ViewFMatrix<T>::order() const {
    return order_;
}

template <typename T>
inline T* ViewFMatrix<T>::pointer() const {
    return matrix_;
}
//-----end ViewFMatrix-----


//5. CArray
// indicies are [0:N-1]
template <typename T>
class CArray {
    
private:
    size_t dims_[7];
    size_t length_; // Length of 1D array
    size_t order_;  // tensor order (rank)
    std::shared_ptr <T []> array_;

public:
    // Default constructor
    CArray ();

    // --- 1D to 7D array ---
    
    CArray (size_t dim0);

    CArray (size_t dim0,
            size_t dim1);

    CArray (size_t dim0,
            size_t dim1,
            size_t dim2);

    CArray (size_t dim0,
            size_t dim1,
            size_t dim2,
            size_t dim3);

    CArray (size_t dim0,
            size_t dim1,
            size_t dim2,
            size_t dim3,
            size_t dim4);

    CArray (size_t dim0,
            size_t dim1,
            size_t dim2,
            size_t dim3,
            size_t dim4,
            size_t dim5);

    CArray (size_t dim0,
            size_t dim1,
            size_t dim2,
            size_t dim3,
            size_t dim4,
            size_t dim5,
            size_t dim6);
    
    CArray (const CArray& temp);
    
    // Overload operator()
    T& operator() (size_t i) const;
    
    T& operator() (size_t i,
                   size_t j) const;
    
    T& operator() (size_t i,
                   size_t j,
                   size_t k) const;
    
    T& operator() (size_t i,
                   size_t j,
                   size_t k,
                   size_t l) const;
    
    T& operator() (size_t i,
                   size_t j,
                   size_t k,
                   size_t l,
                   size_t m) const;
    
    T& operator() (size_t i,
                   size_t j,
                   size_t k,
                   size_t l,
                   size_t m,
                   size_t n) const;
    
    T& operator() (size_t i,
                   size_t j,
                   size_t k,
                   size_t l,
                   size_t m,
                   size_t n,
                   size_t o) const;
    
    // Overload copy assignment operator
    CArray& operator= (const CArray& temp); 

     //return array size
    size_t size() const;

    // return array dims
    size_t dims(size_t i) const;
    
    // return array order (rank)
    size_t order() const;
    
    //return pointer
    T* pointer() const;

    // Deconstructor
    ~CArray ();

}; // End of CArray

//---carray class declarations---

//constructors

//no dim
template <typename T>
CArray<T>::CArray() {
    array_ = NULL;
    length_ = order_ = 0;
}

//1D
template <typename T>
CArray<T>::CArray(size_t dim0)
{
    dims_[0] = dim0;
    order_ = 1;
    length_ = dim0;
    array_ = std::shared_ptr <T[]> (new T[length_]);
}

//2D
template <typename T>
CArray<T>::CArray(size_t dim0,
                  size_t dim1)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    order_ = 2;
    length_ = dim0 * dim1;
    array_ = std::shared_ptr <T[]> (new T[length_]);
}

//3D
template <typename T>
CArray<T>::CArray(size_t dim0,
                  size_t dim1,
                  size_t dim2)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    order_ = 3;
    length_ = dim0 * dim1 * dim2;
    array_ = std::shared_ptr <T[]> (new T[length_]);
}

//4D
template <typename T>
CArray<T>::CArray(size_t dim0,
                  size_t dim1,
                  size_t dim2,
                  size_t dim3)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    order_ = 4;
    length_ = dim0 * dim1 * dim2 * dim3;
    array_ = std::shared_ptr <T[]> (new T[length_]);
}

//5D
template <typename T>
CArray<T>::CArray(size_t dim0,
                  size_t dim1,
                  size_t dim2,
                  size_t dim3,
                  size_t dim4) {
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    order_ = 5;
    length_ = dim0 * dim1 * dim2 * dim3 * dim4;
    array_ = std::shared_ptr <T[]> (new T[length_]);
}

//6D
template <typename T>
CArray<T>::CArray(size_t dim0,
                  size_t dim1,
                  size_t dim2,
                  size_t dim3,
                  size_t dim4,
                  size_t dim5) {
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    order_ = 6;
    length_ = dim0 * dim1 * dim2 * dim3 * dim4 * dim5;
    array_ = std::shared_ptr <T[]> (new T[length_]);
}

//7D
template <typename T>
CArray<T>::CArray(size_t dim0,
                  size_t dim1,
                  size_t dim2,
                  size_t dim3,
                  size_t dim4,
                  size_t dim5,
                  size_t dim6) {
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    dims_[6] = dim6;
    order_ = 7;
    length_ = dim0 * dim1 * dim2 * dim3 * dim4 * dim5 * dim6;
    array_ = std::shared_ptr <T[]> (new T[length_]);
}

//Copy constructor

template <typename T>
CArray<T>::CArray(const CArray& temp) {
    
    // Do nothing if the assignment is of the form x = x
    
    if (this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for
        
        order_  = temp.order_;
        length_ = temp.length_;
        array_ = temp.array_;
    } // end if
    
} // end constructor


//overload () operator

//1D
template <typename T>
inline T& CArray<T>::operator() (size_t i) const
{
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in CArray 1D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in CArray 1D!");

    return array_[i];
}

//2D
template <typename T>
inline T& CArray<T>::operator() (size_t i,
                                 size_t j) const
{
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in CArray 2D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in CArray 2D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in CArray 2D!");
    
    return array_[j + (i *  dims_[1])];
}

//3D
template <typename T>
inline T& CArray<T>::operator() (size_t i,
                                 size_t j,
                                 size_t k) const
{
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in CArray 3D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in CArray 3D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in Carray 3D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in CArray 3D!");
    
    return array_[k + (j * dims_[2])
                    + (i * dims_[2] *  dims_[1])];
}

//4D
template <typename T>
inline T& CArray<T>::operator() (size_t i,
                                 size_t j,
                                 size_t k,
                                 size_t l) const
{
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in CArray 4D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in CArray 4D");  // die if >= dim0
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in CArray 4D");  // die if >= dim1
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in CArray 4D");  // die if >= dim2
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in CArray 4D");  // die if >= dim3

    return array_[l + (k * dims_[3])
                    + (j * dims_[3] * dims_[2])
                    + (i * dims_[3] * dims_[2] *  dims_[1])];
}

//5D
template <typename T>
inline T& CArray<T>::operator() (size_t i,
                                 size_t j,
                                 size_t k,
                                 size_t l,
                                 size_t m) const
{
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in CArray 5D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in CArray 5D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in CArray 5D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in CArray 5D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in CArray 5D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in CArray 5D!");
    
    return array_[m + (l * dims_[4])
                    + (k * dims_[4] * dims_[3])
                    + (j * dims_[4] * dims_[3] * dims_[2])
                    + (i * dims_[4] * dims_[3] * dims_[2] *  dims_[1])];
}

//6D
template <typename T>
inline T& CArray<T>::operator() (size_t i,
                                 size_t j,
                                 size_t k,
                                 size_t l,
                                 size_t m,
                                 size_t n) const
{
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in CArray 6D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in CArray 6D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in CArray 6D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in CArray 6D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in CArray 6D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in CArray 6D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in CArray 6D!");
    
    return array_[n + (m * dims_[5])
                    + (l * dims_[5] * dims_[4])
                    + (k * dims_[5] * dims_[4] * dims_[3])
                    + (j * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                    + (i * dims_[5] * dims_[4] * dims_[3] * dims_[2] *  dims_[1])];
}

//7D
template <typename T>
inline T& CArray<T>::operator() (size_t i,
                                 size_t j,
                                 size_t k,
                                 size_t l,
                                 size_t m,
                                 size_t n,
                                 size_t o) const
{
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in CArray 7D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in CArray 7D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in CArray 7D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in CArray 7D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in CArray 7D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in CArray 7D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in CArray 7D!");
    assert(o >= 0 && o < dims_[6] && "o is out of bounds in CArray 7D!");
    
    return array_[o + (n * dims_[6])
                    + (m * dims_[6] * dims_[5])
                    + (l * dims_[6] * dims_[5] * dims_[4])
                    + (k * dims_[6] * dims_[5] * dims_[4] * dims_[3])
                    + (j * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                    + (i * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2] *  dims_[1])];
    
}


//overload = operator
template <typename T>
inline CArray<T>& CArray<T>::operator= (const CArray& temp)
{
    
    // Do nothing if the assignment is of the form x = x
    if (this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for

        order_  = temp.order_;
        length_ = temp.length_;
        array_  = temp.array_;
    }
    return *this;
}



//return size
template <typename T>
inline size_t CArray<T>::size() const {
    return length_;
}

template <typename T>
inline size_t CArray<T>::dims(size_t i) const {
    assert(i < order_ && "CArray order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to CArray dims is out of bounds!");
    return dims_[i];
}

template <typename T>
inline size_t CArray<T>::order() const {
    return order_;
}


template <typename T>
inline T* CArray<T>::pointer() const{
    return array_.get();
}

//destructor
template <typename T>
CArray<T>::~CArray() {}

//----endof carray class definitions----


//6. ViewCArray
// indicies are [0:N-1]
template <typename T>
class ViewCArray {

private:
    size_t dims_[7];
    size_t length_; // Length of 1D array
    size_t order_;  // tensor order (rank)
    T * array_;
    
public:
    
    // Default constructor
    ViewCArray ();
    
    //--- 1D to 7D array ---
    ViewCArray(T *array,
               size_t dim0);

    ViewCArray(T *array,
               size_t dim0,
               size_t dim1);
    
    ViewCArray(T *some_array,
               size_t dim0,
               size_t dim1,
               size_t dim2);
    
    ViewCArray(T *some_array,
               size_t dim0,
               size_t dim1,
               size_t dim2,
               size_t dim3);
    
    ViewCArray (T *some_array,
                size_t dim0,
                size_t dim1,
                size_t dim2,
                size_t dim3,
                size_t dim4);

    ViewCArray (T *some_array,
                size_t dim0,
                size_t dim1,
                size_t dim2,
                size_t dim3,
                size_t dim4,
                size_t dim5);
 
    ViewCArray (T *some_array,
                size_t dim0,
                size_t dim1,
                size_t dim2,
                size_t dim3,
                size_t dim4,
                size_t dim5,
                size_t dim6);
    
    T& operator()(size_t i) const;
    
    T& operator()(size_t i,
                  size_t j) const;
    
    T& operator()(size_t i,
                  size_t j,
                  size_t k) const;
    
    T& operator()(size_t i,
                  size_t j,
                  size_t k,
                  size_t l) const;
    T& operator()(size_t i,
                  size_t j,
                  size_t k,
                  size_t l,
                  size_t m) const;
    
    T& operator()(size_t i,
                  size_t j,
                  size_t k,
                  size_t l,
                  size_t m,
                  size_t n) const;
    
    T& operator()(size_t i,
                  size_t j,
                  size_t k,
                  size_t l,
                  size_t m,
                  size_t n,
                  size_t o) const;

    // calculate C = math(A,B)
    template <typename M>
    void operator=(M do_this_math);
    
    //return array size
    size_t size() const;
    
    // return array dims
    size_t dims(size_t i) const;
    
    // return array order (rank)
    size_t order() const;

    // return pointer
    T* pointer() const;
    
}; // end of ViewCArray

//class definitions

//constructors

//no dim
template <typename T>
ViewCArray<T>::ViewCArray() {
  array_ = NULL;
  length_ = order_ = 0;
}

//1D
template <typename T>
ViewCArray<T>::ViewCArray(T *array,
                          size_t dim0)
{
    dims_[0] = dim0;
    order_ = 1;
    length_ = dim0;
    array_ = array;
}

//2D
template <typename T>
ViewCArray<T>::ViewCArray(T *array,
                          size_t dim0,
                          size_t dim1)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    order_ = 2;
    length_ = dim0 * dim1;
    array_ = array;
}

//3D
template <typename T>
ViewCArray<T>::ViewCArray (T *array,
                           size_t dim0,
                           size_t dim1,
                           size_t dim2)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    order_ = 3;
    length_ = dim0 * dim1 * dim2;
    array_ = array;
}

//4D
template <typename T>
ViewCArray<T>::ViewCArray(T *array,
                          size_t dim0,
                          size_t dim1,
                          size_t dim2,
                          size_t dim3)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    order_ = 4;
    length_ = dim0 * dim1 * dim2 * dim3;
    array_ = array;
}

//5D
template <typename T>
ViewCArray<T>::ViewCArray(T *array,
                          size_t dim0,
                          size_t dim1,
                          size_t dim2,
                          size_t dim3,
                          size_t dim4)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    order_ = 5;
    length_ = dim0 * dim1 * dim2 * dim3 * dim4;
    array_ = array;
}

//6D
template <typename T>
ViewCArray<T>::ViewCArray(T *array,
                          size_t dim0,
                          size_t dim1,
                          size_t dim2,
                          size_t dim3,
                          size_t dim4,
                          size_t dim5)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    order_ = 6;
    length_ = dim0 * dim1 * dim2 * dim3 * dim4 * dim5;
    array_ = array;
}

//7D
template <typename T>
ViewCArray<T>::ViewCArray(T *array,
                          size_t dim0,
                          size_t dim1,
                          size_t dim2,
                          size_t dim3,
                          size_t dim4,
                          size_t dim5,
                          size_t dim6)
{
    dims_[0] = dim0;
    dims_[1] = dim1;
    dims_[2] = dim2;
    dims_[3] = dim3;
    dims_[4] = dim4;
    dims_[5] = dim5;
    dims_[6] = dim6;
    order_ = 7;
    length_ = dim0 * dim1 * dim2 * dim3 * dim4 * dim5 * dim6;
    array_ = array;
}

//overload () operator

//1D
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i) const
{
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in ViewCArray 1D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewCArray 1D!");
    
    return array_[i];
}

/*
//specification for CArray type
//1D
template <typename T>
inline T& ViewCArray<CArray<T>>::operator()(size_t i) const
{
    assert(i < dim1_ && "i is out of bounds in c_array 1D");  // die if >= dim1
    
    return (*this_array_)(i);
}
*/

//2D
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i, 
                                    size_t j) const
{
   
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in ViewCArray 2D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewCArray 2D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewCArray 2D!");
    
    return array_[j + (i *  dims_[1])];
}

//3D
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i, 
                                    size_t j, 
                                    size_t k) const
{
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in ViewCArray 3D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewCArray 3D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewCarray 3D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewCArray 3D!");
    
    return array_[k + (j * dims_[2])
                    + (i * dims_[2] *  dims_[1])];
}

//4D
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i, 
                                    size_t j, 
                                    size_t k, 
                                    size_t l) const
{
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in ViewCArray 4D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewCArray 4D");  // die if >= dim0
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewCArray 4D");  // die if >= dim1
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewCArray 4D");  // die if >= dim2
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in ViewCArray 4D");  // die if >= dim3
    
    return array_[l + (k * dims_[3])
                    + (j * dims_[3] * dims_[2])
                    + (i * dims_[3] * dims_[2] *  dims_[1])];
}

//5D
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i, 
                                    size_t j, 
                                    size_t k, 
                                    size_t l, 
                                    size_t m) const
{
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in ViewCArray 5D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewCArray 5D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewCArray 5D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewCArray 5D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in ViewCArray 5D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in ViewCArray 5D!");
    
    return array_[m + (l * dims_[4])
                    + (k * dims_[4] * dims_[3])
                    + (j * dims_[4] * dims_[3] * dims_[2])
                    + (i * dims_[4] * dims_[3] * dims_[2] *  dims_[1])];
}

//6D
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i,
                                    size_t j,
                                    size_t k,
                                    size_t l,
                                    size_t m,
                                    size_t n) const
{
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in ViewCArray 6D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewCArray 6D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewCArray 6D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewCArray 6D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in ViewCArray 6D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in ViewCArray 6D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in ViewCArray 6D!");
    
    return array_[n + (m * dims_[5])
                    + (l * dims_[5] * dims_[4])
                    + (k * dims_[5] * dims_[4] * dims_[3])
                    + (j * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                    + (i * dims_[5] * dims_[4] * dims_[3] * dims_[2] *  dims_[1])];
}

//7D
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i,
                                    size_t j,
                                    size_t k,
                                    size_t l,
                                    size_t m,
                                    size_t n,
                                    size_t o) const
{
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in ViewCArray 7D!");
    assert(i >= 0 && i < dims_[0] && "i is out of bounds in ViewCArray 7D!");
    assert(j >= 0 && j < dims_[1] && "j is out of bounds in ViewCArray 7D!");
    assert(k >= 0 && k < dims_[2] && "k is out of bounds in ViewCArray 7D!");
    assert(l >= 0 && l < dims_[3] && "l is out of bounds in ViewCArray 7D!");
    assert(m >= 0 && m < dims_[4] && "m is out of bounds in ViewCArray 7D!");
    assert(n >= 0 && n < dims_[5] && "n is out of bounds in ViewCArray 7D!");
    assert(o >= 0 && o < dims_[6] && "o is out of bounds in ViewCArray 7D!");
    
    return array_[o + (n * dims_[6])
                    + (m * dims_[6] * dims_[5])
                    + (l * dims_[6] * dims_[5] * dims_[4])
                    + (k * dims_[6] * dims_[5] * dims_[4] * dims_[3])
                    + (j * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2])
                    + (i * dims_[6] * dims_[5] * dims_[4] * dims_[3] * dims_[2] *  dims_[1])];
}


// calculate this ViewFArray object = math(A,B)
template <typename T>
template <typename M>
void ViewCArray<T>::operator=(M do_this_math){
    do_this_math(*this); // pass in this ViewFArray object
}// end of math opperation

//return size    
template <typename T>
inline size_t ViewCArray<T>::size() const {
    return length_;
}

template <typename T>
inline size_t ViewCArray<T>::dims(size_t i) const {
    assert(i < order_ && "ViewCArray order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to ViewCArray dims is out of bounds!");
    return dims_[i];
}

template <typename T>
inline size_t ViewCArray<T>::order() const {
    return order_;
}

template <typename T>
inline T* ViewCArray<T>::pointer() const {
    return array_;
}

//---end of ViewCArray class definitions----


//7. CMatrix
template <typename T>
class CMatrix {
        
private:
    size_t dims_[7];
    size_t length_; // Length of 1D array
    size_t order_;  // tensor order (rank)
    std::shared_ptr <T []> matrix_;
            
public:
        
    // default constructor
    CMatrix();

    CMatrix(size_t dim1);

    CMatrix(size_t dim1,
            size_t dim2);

    CMatrix(size_t dim1,
            size_t dim2,
            size_t dim3);

    CMatrix(size_t dim1,
            size_t dim2,
            size_t dim3,
            size_t dim4);

    CMatrix(size_t dim1,
            size_t dim2,
            size_t dim3,
            size_t dim4,
            size_t dim5);

    CMatrix (size_t dim1,
            size_t dim2,
            size_t dim3,
            size_t dim4,
            size_t dim5,
            size_t dim6);

    CMatrix (size_t dim1,
            size_t dim2,
            size_t dim3,
            size_t dim4,
            size_t dim5,
            size_t dim6,
            size_t dim7);

    CMatrix(const CMatrix& temp);
    
    //overload operators to access data
    T& operator()(size_t i) const;

    T& operator()(size_t i,
                    size_t j) const;

    T& operator()(size_t i,
                    size_t j,
                    size_t k) const;

    T& operator()(size_t i,
                    size_t j,
                    size_t k,
                    size_t l) const;

    T& operator()(size_t i,
                    size_t j,
                    size_t k,
                    size_t l,
                    size_t m) const;

    T& operator()(size_t i,
                    size_t j,
                    size_t k,
                    size_t l,
                    size_t m,
                    size_t n) const;

    T& operator()(size_t i,
                    size_t j,
                    size_t k,
                    size_t l,
                    size_t m,
                    size_t n,
                    size_t o) const;

    //overload = operator
    CMatrix& operator= (const CMatrix &temp);

    //return array size
    size_t size() const;
    
    // return array dims
    size_t dims(size_t i) const;
    
    // return array order (rank)
    size_t order() const;

    //return pointer
    T* pointer() const;
    
    // deconstructor
    ~CMatrix( );
        
}; // end of CMatrix

// CMatrix class definitions

//constructors

//no dim

//1D
template <typename T>
CMatrix<T>::CMatrix() {
    matrix_ = NULL;
    length_ = 0;
}

//1D
template <typename T>
CMatrix<T>::CMatrix(size_t dim1)
{
    dims_[0] = dim1;
    order_ = 1;
    length_ = dim1;
    matrix_ = std::shared_ptr <T[]> (new T[length_]);
}

//2D
template <typename T>
CMatrix<T>::CMatrix(size_t dim1,
                    size_t dim2)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    order_ = 2;
    length_ = dim1 * dim2;
    matrix_ = std::shared_ptr <T[]> (new T[length_]);
}

//3D
template <typename T>
CMatrix<T>::CMatrix(size_t dim1,
                    size_t dim2,
                    size_t dim3)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    order_ = 3;
    length_ = dim1 * dim2 * dim3;
    matrix_ = std::shared_ptr <T[]> (new T[length_]);
}

//4D
template <typename T>
CMatrix<T>::CMatrix(size_t dim1,
                    size_t dim2,
                    size_t dim3,
                    size_t dim4)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    order_ = 4;
    length_ = dim1 * dim2 * dim3 * dim4;
    matrix_ = std::shared_ptr <T[]> (new T[length_]);
}   

//5D
template <typename T>
CMatrix<T>::CMatrix(size_t dim1,
                    size_t dim2,
                    size_t dim3,
                    size_t dim4,
                    size_t dim5)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    order_ = 5;
    length_ = dim1 * dim2 * dim3 * dim4 * dim5;
    matrix_ = std::shared_ptr <T[]> (new T[length_]);
}

//6D
template <typename T>
CMatrix<T>::CMatrix(size_t dim1,
                    size_t dim2,
                    size_t dim3,
                    size_t dim4,
                    size_t dim5,
                    size_t dim6)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    order_ = 6;
    length_ = dim1 * dim2 * dim3 * dim4 * dim5 * dim6;
    matrix_ = std::shared_ptr <T[]> (new T[length_]);
}

//7D
template <typename T>
CMatrix<T>::CMatrix(size_t dim1,
                    size_t dim2,
                    size_t dim3,
                    size_t dim4,
                    size_t dim5,
                    size_t dim6,
                    size_t dim7)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    dims_[6] = dim7;
    order_ = 7;
    length_ = dim1 * dim2 * dim3 * dim4 * dim5 * dim6 * dim7;
    matrix_ = std::shared_ptr <T[]> (new T[length_]);
}

template <typename T>
CMatrix<T>::CMatrix(const CMatrix& temp) {
    
    // Do nothing if the assignment is of the form x = x
    
    if (this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for
        
        order_  = temp.order_;
        length_ = temp.length_;
        matrix_ = temp.matrix_;
    } // end if
    
} // end constructor

//overload () operator

//1D
template <typename T>
T& CMatrix<T>::operator()(size_t i) const
{
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in CMatrix 1D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in CMatrix 1D!");
    
    return matrix_[i-1];
}

//2D
template <typename T>
T& CMatrix<T>::operator()(size_t i,
                          size_t j) const
{
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in CMatrix 2D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in CMatrix 2D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in CMatrix 2D!");
    
    return matrix_[(j-1) + (i-1)*dims_[1]];
}

//3D
template <typename T>
T& CMatrix<T>::operator()(size_t i,
                          size_t j,
                          size_t k) const
{
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in CMatrix 3D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in CMatrix 3D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in CMatrix 3D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in CMatrix 3D!");
    
    return matrix_[(k-1) + (j-1)*dims_[2]
                         + (i-1)*dims_[2]*dims_[1]];
}

//4D
template <typename T>
T& CMatrix<T>::operator()(size_t i,
                          size_t j,
                          size_t k,
                          size_t l) const
{
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in CMatrix 4D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in CMatrix 4D");  // die if >= dim0
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in CMatrix 4D");  // die if >= dim1
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in CMatrix 4D");  // die if >= dim2
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in CMatrix 4D");  // die if >= dim3
    
    return matrix_[(l-1) + (k-1)*dims_[3]
                         + (j-1)*dims_[3]*dims_[2]
                         + (i-1)*dims_[3]*dims_[2]*dims_[1]];
}

//5D
template <typename T>
T& CMatrix<T>::operator()(size_t i,
                          size_t j,
                          size_t k,
                          size_t l,
                          size_t m) const
{
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in CMatrix 5D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in CMatrix 5D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in CMatrix 5D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in CMatrix 5D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in CMatrix 5D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in CMatrix 5D!");
    
    return matrix_[(m-1) + (l-1)*dims_[4]
                         + (k-1)*dims_[4]*dims_[3]
                         + (j-1)*dims_[4]*dims_[3]*dims_[2]
                         + (i-1)*dims_[4]*dims_[3]*dims_[2]*dims_[1]];
}

//6D
template <typename T>
T& CMatrix<T>::operator()(size_t i,
                          size_t j,
                          size_t k,
                          size_t l,
                          size_t m,
                          size_t n) const
{
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in CMatrix 6D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in CMatrix 6D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in CMatrix 6D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in CMatrix 6D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in CMatrix 6D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in CMatrix 6D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in CMatrix 6D!");
    
    return matrix_[ (n-1) + (m-1)*dims_[5]
                          + (l-1)*dims_[5]*dims_[4]
                          + (k-1)*dims_[5]*dims_[4]*dims_[3]
                          + (j-1)*dims_[5]*dims_[4]*dims_[3]*dims_[2]
                          + (i-1)*dims_[5]*dims_[4]*dims_[3]*dims_[2]*dims_[1]];
}

//7D
template <typename T>
T& CMatrix<T>::operator()(size_t i,
                          size_t j,
                          size_t k,
                          size_t l,
                          size_t m,
                          size_t n,
                          size_t o) const
{
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in CMatrix 7D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in CMatrix 7D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in CMatrix 7D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in CMatrix 7D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in CMatrix 7D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in CMatrix 7D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in CMatrix 7D!");
    assert(o >= 1 && o <= dims_[6] && "o is out of bounds in CMatrix 7D!");
    
    return matrix_[(o-1) + (n-1)*dims_[6]
                         + (m-1)*dims_[6]*dims_[5]
                         + (l-1)*dims_[6]*dims_[5]*dims_[4]
                         + (k-1)*dims_[6]*dims_[5]*dims_[4]*dims_[3]
                         + (j-1)*dims_[6]*dims_[5]*dims_[4]*dims_[3]*dims_[2]
                         + (i-1)*dims_[6]*dims_[5]*dims_[4]*dims_[3]*dims_[2]*dims_[1]];
}

//overload = operator
//THIS = CMatrix<> temp
template <typename T>
CMatrix<T> &CMatrix<T>::operator= (const CMatrix &temp) {
    if(this != &temp) {
        for (int iter = 0; iter < temp.order_; iter++){
            dims_[iter] = temp.dims_[iter];
        } // end for

        order_  = temp.order_;
        length_ = temp.length_;
        matrix_ = temp.matrix_;
    }
  return *this;
}

template <typename T>
inline size_t CMatrix<T>::size() const {
    return length_;
}

template <typename T>
inline size_t CMatrix<T>::dims(size_t i) const {
    i--; // i starts at 1
    assert(i < order_ && "CMatrix order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to CMatrix dims is out of bounds!");
    return dims_[i];
}

template <typename T>
inline size_t CMatrix<T>::order() const {
    return order_;
}

template <typename T>
inline T* CMatrix<T>::pointer() const{
    return matrix_.get();
}

// Destructor
template <typename T>
CMatrix<T>::~CMatrix(){}

//----end of CMatrix class definitions----


//8. ViewCMatrix
//  indices [1:N]
template <typename T>
class ViewCMatrix {

private:
    size_t dims_[7];
    size_t length_; // Length of 1D array
    size_t order_;  // tensor order (rank)
     T * matrix_;
		    
public:
		    
    // default constructor
    ViewCMatrix();
		    
		    
    //--- 1D array ---	   	    
    // overloaded constructor
    ViewCMatrix (T *matrix,
                 size_t dim1);
    
    ViewCMatrix (T *matrix,
                 size_t dim1,
                 size_t dim2);

    ViewCMatrix (T *matrix,
		size_t dim1,
		size_t dim2,
		size_t dim3);

    ViewCMatrix (T *matrix,
		size_t dim1,
		size_t dim2,
		size_t dim3,
		size_t dim4);

    ViewCMatrix (T *matrix,
		size_t dim1,
		size_t dim2,
		size_t dim3,
		size_t dim4,
		size_t dim5);

    ViewCMatrix (T *matrix,
		   size_t dim1,
		   size_t dim2,
		   size_t dim3,
		   size_t dim4,
		   size_t dim5,
		   size_t dim6);

    ViewCMatrix (T *matrix,
                 size_t dim1,
                 size_t dim2,
                 size_t dim3,
                 size_t dim4,
                 size_t dim5,
                 size_t dim6,
                 size_t dim7);
    
    T& operator() (size_t i) const;
    
    T& operator() (size_t i,
                   size_t j) const;
    
    T& operator() (size_t i,
                   size_t j,
                   size_t k) const;
    
    T& operator() (size_t i,
                   size_t j,
                   size_t k,
                   size_t l) const;
    
    T& operator() (size_t i,
                   size_t j,
                   size_t k,
                   size_t l,
                   size_t m) const;
    
    T& operator() (size_t i,
                   size_t j,
                   size_t k,
                   size_t l,
                   size_t m,
                   size_t n) const;
    T& operator() (size_t i,
                   size_t j,
                   size_t k,
                   size_t l,
                   size_t m,
                   size_t n,
                   size_t o) const;

    // calculate C = math(A,B)
    template <typename M>
    void operator=(M do_this_math);
    
    //return array size
    size_t size() const;
    
    // return array dims
    size_t dims(size_t i) const;
    
    // return array order (rank)
    size_t order() const;

    // return pointer
    T* pointer() const;
    
}; // end of ViewCMatrix

//class definitions

//constructors

//no dim
template <typename T>
ViewCMatrix<T>::ViewCMatrix(){
  matrix_ = NULL;
  length_ = 0;
}

//1D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *matrix,
                            size_t dim1)
{
    dims_[0] = dim1;
    order_ = 1;
    length_ = dim1;
	matrix_ = matrix;
}

//2D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *matrix,
                            size_t dim1,
                            size_t dim2)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    order_ = 2;
    length_ = dim1 * dim2;
	matrix_ = matrix;
}

//3D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *matrix,
                            size_t dim1,
                            size_t dim2,
                            size_t dim3)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    order_ = 3;
    length_ = dim1 * dim2 * dim3;
	matrix_ = matrix;
}

//4D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *matrix,
                            size_t dim1,
                            size_t dim2,
                            size_t dim3,
                            size_t dim4)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    order_ = 4;
    length_ = dim1 * dim2 * dim3 * dim4;
	matrix_ = matrix;
}

//5D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *matrix,
                            size_t dim1,
                            size_t dim2,
                            size_t dim3,
                            size_t dim4,
                            size_t dim5)
{
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    order_ = 5;
    length_ = dim1 * dim2 * dim3 * dim4 * dim5;
	matrix_ = matrix;
}

//6D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *matrix,
                            size_t dim1,
                            size_t dim2,
                            size_t dim3,
                            size_t dim4,
                            size_t dim5,
                            size_t dim6) {
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    order_ = 6;
    length_ = dim1 * dim2 * dim3 * dim4 * dim5 * dim6;
	matrix_ = matrix;
}

//7D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *matrix,
                            size_t dim1,
                            size_t dim2,
                            size_t dim3,
                            size_t dim4,
                            size_t dim5,
                            size_t dim6,
                            size_t dim7) {
    dims_[0] = dim1;
    dims_[1] = dim2;
    dims_[2] = dim3;
    dims_[3] = dim4;
    dims_[4] = dim5;
    dims_[5] = dim6;
    dims_[6] = dim7;
    order_ = 7;
    length_ = dim1 * dim2 * dim3 * dim4 * dim5 * dim6 * dim7;
    matrix = matrix_;
}

//overload () operator

//1D
template <typename T>
T& ViewCMatrix<T>:: operator() (size_t i) const
{
    assert(order_ == 1 && "Tensor order (rank) does not match constructor in ViewCMatrix 1D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewCMatrix 1D!");
    
	return matrix_[i-1];
}

//2D
template <typename T>
T& ViewCMatrix<T>::operator() (size_t i,
                               size_t j) const
{
    assert(order_ == 2 && "Tensor order (rank) does not match constructor in ViewCMatrix 2D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewCMatrix 2D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewCMatrix 2D!");
    
    return matrix_[(j-1) + (i-1)*dims_[1]];
}

//3D
template <typename T>
T& ViewCMatrix<T>::operator () (size_t i,
                                size_t j,
                                size_t k) const
{
    assert(order_ == 3 && "Tensor order (rank) does not match constructor in ViewCMatrix 3D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewCMatrix 3D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewCMatrix 3D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in ViewCMatrix 3D!");
    
    return matrix_[(k-1) + (j-1)*dims_[2]
                         + (i-1)*dims_[2]*dims_[1]];
}

//4D
template <typename T>
T& ViewCMatrix<T>::operator()(size_t i,
                              size_t j,
                              size_t k,
                              size_t l) const
{
    assert(order_ == 4 && "Tensor order (rank) does not match constructor in ViewCMatrix 4D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewCMatrix 4D");  // die if >= dim0
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewCMatrix 4D");  // die if >= dim1
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in ViewCMatrix 4D");  // die if >= dim2
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in ViewCMatrix 4D");  // die if >= dim3
    
    return matrix_[(l-1) + (k-1)*dims_[3]
                         + (j-1)*dims_[3]*dims_[2]
                         + (i-1)*dims_[3]*dims_[2]*dims_[1]];
}

//5D
template <typename T>
T& ViewCMatrix<T>::operator()(size_t i,
                              size_t j,
                              size_t k,
                              size_t l,
                              size_t m) const
{
    assert(order_ == 5 && "Tensor order (rank) does not match constructor in ViewCMatrix 5D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewCMatrix 5D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewCMatrix 5D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in ViewCMatrix 5D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in ViewCMatrix 5D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in ViewCMatrix 5D!");
    
    return matrix_[(m-1) + (l-1)*dims_[4]
                         + (k-1)*dims_[4]*dims_[3]
                         + (j-1)*dims_[4]*dims_[3]*dims_[2]
                         + (i-1)*dims_[4]*dims_[3]*dims_[2]*dims_[1]];
}

//6D
template <typename T>
T& ViewCMatrix<T>::operator()(size_t i,
                              size_t j,
                              size_t k,
                              size_t l,
                              size_t m,
                              size_t n) const
{
    assert(order_ == 6 && "Tensor order (rank) does not match constructor in ViewCMatrix 6D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewCMatrix 6D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewCMatrix 6D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in ViewCMatrix 6D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in ViewCMatrix 6D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in ViewCMatrix 6D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in ViewCMatrix 6D!");
    
    return matrix_[(n-1) + (m-1)*dims_[5]
                         + (l-1)*dims_[5]*dims_[4]
                         + (k-1)*dims_[5]*dims_[4]*dims_[3]
                         + (j-1)*dims_[5]*dims_[4]*dims_[3]*dims_[2]
                         + (i-1)*dims_[5]*dims_[4]*dims_[3]*dims_[2]*dims_[1]];
}

//7D
template <typename T>
T& ViewCMatrix<T>::operator()(size_t i,
                              size_t j,
                              size_t k,
                              size_t l,
                              size_t m,
                              size_t n,
                              size_t o) const
{
    assert(order_ == 7 && "Tensor order (rank) does not match constructor in ViewCMatrix 7D!");
    assert(i >= 1 && i <= dims_[0] && "i is out of bounds in ViewCMatrix 7D!");
    assert(j >= 1 && j <= dims_[1] && "j is out of bounds in ViewCMatrix 7D!");
    assert(k >= 1 && k <= dims_[2] && "k is out of bounds in ViewCMatrix 7D!");
    assert(l >= 1 && l <= dims_[3] && "l is out of bounds in ViewCMatrix 7D!");
    assert(m >= 1 && m <= dims_[4] && "m is out of bounds in ViewCMatrix 7D!");
    assert(n >= 1 && n <= dims_[5] && "n is out of bounds in ViewCMatrix 7D!");
    assert(o >= 1 && o <= dims_[6] && "o is out of bounds in ViewCMatrix 7D!");
    
    return matrix_[(o-1) + (n-1)*dims_[6]
                         + (m-1)*dims_[6]*dims_[5]
                         + (l-1)*dims_[6]*dims_[5]*dims_[4]
                         + (k-1)*dims_[6]*dims_[5]*dims_[4]*dims_[3]
                         + (j-1)*dims_[6]*dims_[5]*dims_[4]*dims_[3]*dims_[2]
                         + (i-1)*dims_[6]*dims_[5]*dims_[4]*dims_[3]*dims_[2]*dims_[1]];
}

// calculate this ViewFArray object = math(A,B)
template <typename T>
template <typename M>
void ViewCMatrix<T>::operator=(M do_this_math){
    do_this_math(*this); // pass in this ViewFArray object
}// end of math opperation

template <typename T>
inline size_t ViewCMatrix<T>::size() const {
    return length_;
}

template <typename T>
inline size_t ViewCMatrix<T>::dims(size_t i) const {
    i--; // i starts at 1
    assert(i < order_ && "ViewCMatrix order (rank) does not match constructor, dim[i] does not exist!");
    assert(i >= 0 && dims_[i]>0 && "Access to ViewCMatrix dims is out of bounds!");
    return dims_[i];
}

template <typename T>
inline size_t ViewCMatrix<T>::order() const {
    return order_;
}

template <typename T>
inline T* ViewCMatrix<T>::pointer() const {
    return matrix_;
}


//----end of ViewCMatrix class definitions----

//9. RaggedRightArray
template <typename T>
class RaggedRightArray {
private:
    size_t *start_index_;
    T * array_;
    
    size_t dim1_, length_;
    size_t num_saved_; // the number saved in the 1D array
    
public:
    // Default constructor
    RaggedRightArray ();
    
    //--- 2D array access of a ragged right array ---
    
    // Overload constructor for a CArray
    RaggedRightArray (CArray<size_t> &strides_array);
    
    // Overload constructor for a ViewCArray
    RaggedRightArray (ViewCArray<size_t> &strides_array);
    
    // Overloaded constructor for a traditional array
    RaggedRightArray (size_t *strides_array, size_t some_dim1);
    
    // Overload constructor for a RaggedRightArray to
    // support a dynamically built stride_array
    RaggedRightArray (size_t some_dim1, size_t buffer);
    
    // A method to return the stride size
    size_t stride(size_t i) const;
    
    // A method to increase the number of column entries, i.e.,
    // the stride size. Used with the constructor for building
    // the stride_array dynamically.
    // DO NOT USE with the constructures with a strides_array
    void push_back(size_t i);
    
    // Overload operator() to access data as array(i,j)
    // where i=[0:N-1], j=[stride(i)]
    T& operator()(size_t i, size_t j) const;

    // method to return total size
    size_t size() const;

    //return pointer
    T* pointer() const;
    
    //get row starts array
    size_t* get_starts() const;

    RaggedRightArray& operator+= (const size_t i);

    RaggedRightArray& operator= (const RaggedRightArray &temp);

    // Destructor
    ~RaggedRightArray ( );
}; // End of RaggedRightArray

// Default constructor
template <typename T>
RaggedRightArray<T>::RaggedRightArray () {
    array_ = NULL;
    start_index_ = NULL;
    length_ = 0;
}


// Overloaded constructor with CArray
template <typename T>
RaggedRightArray<T>::RaggedRightArray (CArray<size_t> &strides_array){
    // The length of the stride array is some_dim1;
    dim1_  = strides_array.size();
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[(dim1_ + 1)];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array(i);
        start_index_[(i + 1)] = count;
    } // end for i
    length_ = count;
    
    array_ = new T[length_];
} // End constructor

// Overloaded constructor with a view c array
template <typename T>
RaggedRightArray<T>::RaggedRightArray (ViewCArray<size_t> &strides_array) {
    // The length of the stride array is some_dim1;
    dim1_  = strides_array.size();
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[(dim1_ + 1)];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array(i);
        start_index_[(i + 1)] = count;
    } // end for i
    length_ = count;
    
    array_ = new T[length_];
} // End constructor

// Overloaded constructor with a regular cpp array
template <typename T>
RaggedRightArray<T>::RaggedRightArray (size_t *strides_array, size_t dim1){
    // The length of the stride array is some_dim1;
    dim1_ = dim1;
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[(dim1_ + 1)];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array[i];
        start_index_[(i + 1)] = count;
    } // end for i
    length_ = count;
    
    array_ = new T[length_];
} // End constructor

// overloaded constructor for a dynamically built strides_array.
// buffer is the max number of columns needed
template <typename T>
RaggedRightArray<T>::RaggedRightArray (size_t some_dim1, size_t buffer){
    
    dim1_ = some_dim1;
    
    // create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[dim1_+1]();  // note the dim1+1
    //start_index_[0] = 0; // the 1D array starts at 0

    num_saved_ = 0;
    
    length_ = some_dim1*buffer;
    array_ = new T[some_dim1*buffer];
    
} // end constructor

// A method to return the stride size
template <typename T>
inline size_t RaggedRightArray<T>::stride(size_t i) const {
    // Ensure that i is within bounds
    assert(i < dim1_ && "i is greater than dim1_ in RaggedRightArray");

    return start_index_[(i + 1)] - start_index_[i];
}

// A method to increase the stride size, in other words,
// this is used to build the stride array dynamically
// DO NOT USE with constructors that are given a stride array
template <typename T>
void RaggedRightArray<T>::push_back(size_t i){
    num_saved_ ++;
    start_index_[i+1] = num_saved_;
}

// Overload operator() to access data as array(i,j)
// where i=[0:N-1], j=[0:stride(i)]
template <typename T>
inline T& RaggedRightArray<T>::operator()(size_t i, size_t j) const {
    // get the 1D array index
    size_t start = start_index_[i];
    
    // asserts
    assert(i < dim1_ && "i is out of dim1 bounds in RaggedRightArray");  // die if >= dim1
    //assert(j < stride(i) && "j is out of stride bounds in RaggedRightArray");  // die if >= stride
    assert(j+start < length_ && "j+start is out of bounds in RaggedRightArray");  // die if >= 1D array length)
    
    return array_[j + start];
} // End operator()

//return size
template <typename T>
size_t RaggedRightArray<T>::size() const {
    return length_;
}

template <typename T>
RaggedRightArray<T> & RaggedRightArray<T>::operator+= (const size_t i) {
    this->num_saved_ ++;
    this->start_index_[i+1] = num_saved_;
    return *this;
}

//overload = operator
template <typename T>
RaggedRightArray<T> & RaggedRightArray<T>::operator= (const RaggedRightArray &temp) {

    if( this != &temp) {
        dim1_ = temp.dim1_;
        length_ = temp.length_;
        num_saved_ = temp.num_saved_;
        if(start_index_!=NULL)
          delete[] start_index_;
        start_index_ = new size_t[dim1_ + 1];
        for (int j = 0; j < dim1_ + 1; j++) {
            start_index_[j] = temp.start_index_[j];  
        }

        if(array_!=NULL)
          delete[] array_;
        array_ = new T[length_];
        //copy contents
        for(int iter = 0; iter < length_; iter++)
          array_[iter] = temp.array_[iter];
    }
	
    return *this;
}

template <typename T>
inline T* RaggedRightArray<T>::pointer() const{
    return array_;
}

template <typename T>
inline size_t* RaggedRightArray<T>::get_starts() const{
    return start_index_;
}

// Destructor
template <typename T>
RaggedRightArray<T>::~RaggedRightArray () {
    if(array_!=NULL)
      delete[] array_;
    if(start_index_!=NULL)
      delete[] start_index_;
}

//----end of RaggedRightArray class definitions----

//9. RaggedRightArrayofVectors
template <typename T>
class RaggedRightArrayofVectors {
private:
    size_t *start_index_;
    T * array_;
    
    size_t dim1_, length_, vector_dim_;
    size_t num_saved_; // the number saved in the 1D array
    
public:
    // Default constructor
    RaggedRightArrayofVectors ();
    
    //--- 3D array access of a ragged right array storing a vector of size vector_dim_ at each (i,j)---
    
    // Overload constructor for a CArray
    RaggedRightArrayofVectors (CArray<size_t> &strides_array, size_t vector_dim);
    
    // Overload constructor for a ViewCArray
    RaggedRightArrayofVectors (ViewCArray<size_t> &strides_array, size_t vector_dim);
    
    // Overloaded constructor for a traditional array
    RaggedRightArrayofVectors (size_t *strides_array, size_t some_dim1, size_t vector_dim);
    
    // Overload constructor for a RaggedRightArray to
    // support a dynamically built stride_array
    RaggedRightArrayofVectors (size_t some_dim1, size_t buffer, size_t vector_dim);
    
    // A method to return the stride size
    size_t stride(size_t i) const;

    // A method to return the vector dim
    size_t vector_dim() const;
    
    // A method to increase the number of column entries, i.e.,
    // the stride size. Used with the constructor for building
    // the stride_array dynamically.
    // DO NOT USE with the constructures with a strides_array
    void push_back(size_t i);
    
    // Overload operator() to access data as array(i,j)
    // where i=[0:N-1], j=[stride(i)], k=[0,vector_dim_]
    T& operator()(size_t i, size_t j, size_t k) const;

    // method to return total size
    size_t size() const;

    //return pointer
    T* pointer() const;
    
    //get row starts array
    size_t* get_starts() const;

    RaggedRightArrayofVectors& operator+= (const size_t i);

    RaggedRightArrayofVectors& operator= (const RaggedRightArrayofVectors &temp);

    // Destructor
    ~RaggedRightArrayofVectors ( );
}; // End of RaggedRightArray

// Default constructor
template <typename T>
RaggedRightArrayofVectors<T>::RaggedRightArrayofVectors () {
    array_ = NULL;
    start_index_ = NULL;
    length_ = 0;
}


// Overloaded constructor with CArray
template <typename T>
RaggedRightArrayofVectors<T>::RaggedRightArrayofVectors (CArray<size_t> &strides_array, size_t vector_dim){
    // The length of the stride array is some_dim1;
    dim1_  = strides_array.size();
    vector_dim_ = vector_dim;
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[(dim1_ + 1)];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array(i)*vector_dim_;
        start_index_[(i + 1)] = count;
    } // end for i
    length_ = count;
    
    array_ = new T[length_];
} // End constructor

// Overloaded constructor with a view c array
template <typename T>
RaggedRightArrayofVectors<T>::RaggedRightArrayofVectors (ViewCArray<size_t> &strides_array, size_t vector_dim) {
    // The length of the stride array is some_dim1;
    dim1_  = strides_array.size();
    vector_dim_ = vector_dim;
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[(dim1_ + 1)];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array(i)*vector_dim_;
        start_index_[(i + 1)] = count;
    } // end for i
    length_ = count;
    
    array_ = new T[length_];
} // End constructor

// Overloaded constructor with a regular cpp array
template <typename T>
RaggedRightArrayofVectors<T>::RaggedRightArrayofVectors (size_t *strides_array, size_t dim1, size_t vector_dim){
    // The length of the stride array is some_dim1;
    dim1_ = dim1;
    vector_dim_ = vector_dim;

    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[(dim1_ + 1)];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array of vectors and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array[i]*vector_dim_;
        start_index_[(i + 1)] = count;
    } // end for i
    length_ = count;
    
    array_ = new T[length_];
} // End constructor

// overloaded constructor for a dynamically built strides_array.
// buffer is the max number of columns needed
template <typename T>
RaggedRightArrayofVectors<T>::RaggedRightArrayofVectors (size_t some_dim1, size_t buffer, size_t vector_dim){
    
    dim1_ = some_dim1;
    vector_dim_ = vector_dim;

    // create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[dim1_+1]();  // note the dim1+1
    //start_index_[0] = 0; // the 1D array starts at 0

    num_saved_ = 0;
    
    length_ = some_dim1*buffer*vector_dim;
    array_ = new T[some_dim1*buffer];
    
} // end constructor

// A method to return the stride size
template <typename T>
inline size_t RaggedRightArrayofVectors<T>::stride(size_t i) const {
    // Ensure that i is within bounds
    assert(i < dim1_ && "i is greater than dim1_ in RaggedRightArray");

    return (start_index_[(i + 1)] - start_index_[i])/vector_dim_;
}

// A method to increase the stride size, in other words,
// this is used to build the stride array dynamically
// DO NOT USE with constructors that are given a stride array
template <typename T>
void RaggedRightArrayofVectors<T>::push_back(size_t i){
    num_saved_ += vector_dim_;
    start_index_[i+1] = num_saved_;
}

// Overload operator() to access data as array(i,j,k)
// where i=[0:N-1], j=[0:stride(i)], k=[0:vector_dim_]
template <typename T>
inline T& RaggedRightArrayofVectors<T>::operator()(size_t i, size_t j, size_t k) const {
    // get the 1D array index
    size_t start = start_index_[i];
    
    // asserts
    assert(i < dim1_ && "i is out of dim1 bounds in RaggedRightArray");  // die if >= dim1
    //assert(j < stride(i) && "j is out of stride bounds in RaggedRightArray");  // die if >= stride
    assert(j*vector_dim_+start + k < length_ && "j+start is out of bounds in RaggedRightArray");  // die if >= 1D array length)
    
    return array_[j*vector_dim_ + start + k];
} // End operator()

//return size
template <typename T>
size_t RaggedRightArrayofVectors<T>::size() const {
    return length_;
}

template <typename T>
RaggedRightArrayofVectors<T> & RaggedRightArrayofVectors<T>::operator+= (const size_t i) {
    this->num_saved_ += vector_dim_;
    this->start_index_[i+1] = num_saved_;
    return *this;
}

//overload = operator
template <typename T>
RaggedRightArrayofVectors<T> & RaggedRightArrayofVectors<T>::operator= (const RaggedRightArrayofVectors &temp) {

    if( this != &temp) {
        dim1_ = temp.dim1_;
        vector_dim_ = temp.vector_dim_;
        length_ = temp.length_;
        num_saved_ = temp.num_saved_;
        if(start_index_!=NULL)
          delete[] start_index_;
        start_index_ = new size_t[dim1_ + 1];
        for (int j = 0; j < dim1_ + 1; j++) {
            start_index_[j] = temp.start_index_[j];  
        }

        if(array_!=NULL)
          delete[] array_;
        array_ = new T[length_];
        //copy contents
        for(int iter = 0; iter < length_; iter++)
          array_[iter] = temp.array_[iter];
    }
	
    return *this;
}

template <typename T>
inline T* RaggedRightArrayofVectors<T>::pointer() const{
    return array_;
}

template <typename T>
inline size_t* RaggedRightArrayofVectors<T>::get_starts() const{
    return start_index_;
}

// Destructor
template <typename T>
RaggedRightArrayofVectors<T>::~RaggedRightArrayofVectors () {
    if(array_!=NULL)
      delete[] array_;
    if(start_index_!=NULL)
      delete[] start_index_;
}

//----end of RaggedRightArrayofVectors class definitions----

//10. RaggedDownArray
template <typename T>
class RaggedDownArray { 
private:
    size_t *start_index_;
	T * array_;

	size_t dim2_;
    size_t length_;
    size_t num_saved_; // the number saved in the 1D array

public:
    //default constructor
    RaggedDownArray() ;

    //~~~~2D`~~~~
	//overload constructor with CArray
	RaggedDownArray(CArray<size_t> &strides_array);

	//overload with ViewCArray
	RaggedDownArray(ViewCArray <size_t> &strides_array);

	//overload with traditional array
	RaggedDownArray(size_t *strides_array, size_t dome_dim1);

    // Overload constructor for a RaggedDownArray to
    // support a dynamically built stride_array
    RaggedDownArray (size_t some_dim2, size_t buffer);
    
	//method to return stride size
	size_t stride(size_t j);

    // A method to increase the number of column entries, i.e.,
    // the stride size. Used with the constructor for building
    // the stride_array dynamically.
    // DO NOT USE with the constructures with a strides_array
    void push_back(size_t j);
    
	//overload () operator to access data as array (i,j)
	T& operator()(size_t i, size_t j);

    // method to return total size
    size_t size();

    //return pointer
    T* pointer() const;
    
    //get row starts array
    size_t* get_starts() const;

    //overload = operator
    RaggedDownArray& operator= (const RaggedDownArray &temp);

    //destructor
    ~RaggedDownArray();

}; //~~~~~end of RaggedDownArray class declarations~~~~~~~~	

//no dims
template <typename T>
RaggedDownArray<T>::RaggedDownArray() {
    array_ = NULL;
    start_index_ = NULL;
    length_ = 0;
}

//overload constructor with CArray 
template <typename T>
RaggedDownArray<T>::RaggedDownArray( CArray <size_t> &strides_array) {
    // Length of stride array
    //dim2_ = strides_array.size();

    // Create and initialize startding indices
    start_index_ = new size_t[dim2_+1]; //theres a plus 1, because 
    start_index_[0] = 0; //1D array starts at 0

		
	//length of strides
	dim2_ = strides_array.size();

    // Loop to find total length of 1D array
    size_t count = 0;
    for(size_t j = 0; j < dim2_ ; j++) { 
        count += strides_array(j);
        start_index_[j+1] = count;
    } 
    length_ = count;

    array_ = new T[length_];

} // End constructor 

// Overload constructor with ViewCArray
template <typename T>
RaggedDownArray<T>::RaggedDownArray( ViewCArray <size_t> &strides_array) {
    // Length of strides
    //dim2_ = strides_array.size();

    //create array for holding start indices
    start_index_ = new size_t[dim2_+1];
    start_index_[0] = 0;

    size_t count = 0;
    // Loop over to get total length of 1D array
    for(size_t j = 0; j < dim2_ ;j++ ) {
        count += strides_array(j);
        start_index_[j+1] = count;
    }
    length_ = count;	
    array_ = new T[length_];

} // End constructor 

// Overload constructor with regualar array
template <typename T>
RaggedDownArray<T>::RaggedDownArray( size_t *strides_array, size_t dim2){
    // Length of stride array
    dim2_ = dim2;

    // Create and initialize starting index of entries
    start_index_ = new size_t[dim2_+1];
    start_index_[0] = 0;

    // Loop over to find length of 1D array
    // Represent ragged down array and set 1D index
    size_t count = 0;
    for(size_t j = 0; j < dim2_; j++) {
        count += strides_array[j];
        start_index_[j+1] = count;
	}

    length_ = count;	
    array_ = new T[length_];

} //end construnctor

// overloaded constructor for a dynamically built strides_array.
// buffer is the max number of columns needed
template <typename T>
RaggedDownArray<T>::RaggedDownArray (size_t some_dim2, size_t buffer){
    
    dim2_ = some_dim2;
    
    // create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[dim2_+1]();  // note the dim2+1
    //start_index_[0] = 0; // the 1D array starts at 0
    
    num_saved_ = 0;
    
    length_ = some_dim2*buffer;
    array_ = new T[some_dim2*buffer];
    
} // end constructor

// Check the stride size
template <typename T>
size_t RaggedDownArray<T>::stride(size_t j) {
    assert(j < dim2_ && "j is greater than dim2_ in RaggedDownArray");

    return start_index_[j+1] - start_index_[j];
}

// A method to increase the stride size, in other words,
// this is used to build the stride array dynamically
// DO NOT USE with constructors that are given a stride array
template <typename T>
void RaggedDownArray<T>::push_back(size_t j){
    num_saved_ ++;
    start_index_[j+1] = num_saved_;
}

//return size
template <typename T>
size_t RaggedDownArray<T>::size() {
    return length_;
}

// overload operator () to access data as an array(i,j)
// Note: i = 0:stride(j), j = 0:N-1
template <typename T>
T& RaggedDownArray<T>::operator()(size_t i, size_t j) {
    // Where is the array starting?
    // look at start index
    size_t start = start_index_[j]; 

    // Make sure we are within array bounds
    assert(i < stride(j) && "i is out of bounds in RaggedDownArray");
    assert(j < dim2_ && "j is out of dim2_ bounds in RaggedDownArray");
    assert(i+start < length_ && "i+start is out of bounds in RaggedDownArray");  // die if >= 1D array length)
    
    return array_[i + start];

} // End () operator

//overload = operator
template <typename T>
RaggedDownArray<T> & RaggedDownArray<T>::operator= (const RaggedDownArray &temp) {

    if( this != &temp) {
        dim2_ = temp.dim2_;
        length_ = temp.length_;
        num_saved_ = temp.num_saved_;
        if(start_index_!=NULL)
          delete[] start_index_;

        start_index_ = new size_t[dim2_ + 1];
        for (int j = 0; j < dim2_ + 1; j++) {
            start_index_[j] = temp.start_index_[j];  
        }

        if(array_!=NULL)
          delete[] array_;
        array_ = new T[length_];
        //copy contents
        for(int iter = 0; iter < length_; iter++)
          array_[iter] = temp.array_[iter];
    }
	
    return *this;
}

template <typename T>
inline T* RaggedDownArray<T>::pointer() const{
    return array_;
}


template <typename T>
inline size_t* RaggedDownArray<T>::get_starts() const{
    return start_index_;
}

// Destructor
template <typename T>
RaggedDownArray<T>::~RaggedDownArray() {
  if(array_!=NULL)
    delete[] array_;
  if(start_index_!=NULL)
    delete[] start_index_;

} // End destructor


//----end of RaggedDownArray----


//11. DynamicRaggedRightArray

template <typename T>
class DynamicRaggedRightArray {
private:
    size_t *stride_;
    T * array_;
    
    size_t dim1_;
    size_t dim2_;
    size_t length_;
    
public:
    // Default constructor
    DynamicRaggedRightArray ();
    
    //--- 2D array access of a ragged right array ---
    
    // overload constructor
    DynamicRaggedRightArray (size_t dim1, size_t dim2);
    
    // A method to return or set the stride size
    size_t& stride(size_t i) const;
    
    // A method to return the size
    size_t size() const;

    //return pointer
    T* pointer() const;
    
    // Overload operator() to access data as array(i,j),
    // where i=[0:N-1], j=[stride(i)]
    T& operator()(size_t i, size_t j) const;
    
    // Overload copy assignment operator
    DynamicRaggedRightArray& operator= (const DynamicRaggedRightArray &temp);
    
    // Destructor
    ~DynamicRaggedRightArray ();
};

//nothing
template <typename T>
DynamicRaggedRightArray<T>::DynamicRaggedRightArray () {
    array_ = NULL;
    stride_ = NULL;
    length_ = 0;
}

// Overloaded constructor
template <typename T>
DynamicRaggedRightArray<T>::DynamicRaggedRightArray (size_t dim1, size_t dim2) {
    // The dimensions of the array;
    dim1_  = dim1;
    dim2_  = dim2;
    length_ = dim1*dim2;
    
    // Create memory on the heap for the values
    array_ = new T[dim1*dim2];
    
    // Create memory for the stride size in each row
    stride_ = new size_t[dim1];
    
    // Initialize the stride
    for (int i=0; i<dim1_; i++){
        stride_[i] = 0;
    }
    
    // Start index is always = j + i*dim2
}

// A method to set the stride size for row i
template <typename T>
size_t& DynamicRaggedRightArray<T>::stride(size_t i) const {
    return stride_[i];
}

//return size
template <typename T>
size_t DynamicRaggedRightArray<T>::size() const{
    return length_;
}

// Overload operator() to access data as array(i,j),
// where i=[0:N-1], j=[0:stride(i)]
template <typename T>
inline T& DynamicRaggedRightArray<T>::operator()(size_t i, size_t j) const {
    // Asserts
    assert(i < dim1_ && "i is out of dim1 bounds in DynamicRaggedRight");  // die if >= dim1
    assert(j < dim2_ && "j is out of dim2 bounds in DynamicRaggedRight");  // die if >= dim2
    assert(j < stride_[i] && "j is out of stride bounds in DynamicRaggedRight");  // die if >= stride
    
    return array_[j + i*dim2_];
}

//overload = operator
template <typename T>
inline DynamicRaggedRightArray<T>& DynamicRaggedRightArray<T>::operator= (const DynamicRaggedRightArray &temp)
{
    
    if( this != &temp) {
        dim1_ = temp.dim1_;
        dim2_ = temp.dim2_;
        length_ = temp.length_;
        if(stride_!=NULL)
          delete[] stride_;
        stride_ = new size_t[dim1_];
        for (int i = 0; i < dim1_; i++) {
            stride_[i] = temp.stride_[i];
        }

        if(array_!=NULL)
          delete[] array_;
        array_ = new T[length_];
        //copy contents
        for(int iter = 0; iter < length_; iter++)
          array_[iter] = temp.array_[iter];
    }
    
    return *this;
}

template <typename T>
inline T* DynamicRaggedRightArray<T>::pointer() const{
    return array_;
}

// Destructor
template <typename T>
DynamicRaggedRightArray<T>::~DynamicRaggedRightArray() {
    if(array_!=NULL)
      delete[] array_;
    if(stride_!=NULL)
      delete[] stride_;
}




//----end DynamicRaggedRightArray class definitions----


//12. DynamicRaggedDownArray

template <typename T>
class DynamicRaggedDownArray {
private:
    size_t *stride_;
    T * array_;
    
    size_t dim1_;
    size_t dim2_;
    size_t length_;
    
public:
    // Default constructor
    DynamicRaggedDownArray ();
    
    //--- 2D array access of a ragged right array ---
    
    // overload constructor
    DynamicRaggedDownArray (size_t dim1, size_t dim2);
    
    // A method to return or set the stride size
    size_t& stride(size_t j) const;
    
    // A method to return the size
    size_t size() const;
    
    // Overload operator() to access data as array(i,j),
    // where i=[stride(j)], j=[0:N-1]
    T& operator()(size_t i, size_t j) const;
    
    // Overload copy assignment operator
    DynamicRaggedDownArray& operator= (const DynamicRaggedDownArray &temp);

    //return pointer
    T* pointer() const;
    
    // Destructor
    ~DynamicRaggedDownArray ();
};

//nothing
template <typename T>
DynamicRaggedDownArray<T>::DynamicRaggedDownArray () {
    array_ = NULL;
    stride_ = NULL;
    length_ = 0;
}

// Overloaded constructor
template <typename T>
DynamicRaggedDownArray<T>::DynamicRaggedDownArray (size_t dim1, size_t dim2) {
    // The dimensions of the array;
    dim1_  = dim1;
    dim2_  = dim2;
    length_ = dim1*dim2;
    
    // Create memory on the heap for the values
    array_ = new T[dim1*dim2];
    
    // Create memory for the stride size in each row
    stride_ = new size_t[dim2];
    
    // Initialize the stride
    for (int j=0; j<dim2_; j++){
        stride_[j] = 0;
    }
    
    // Start index is always = i + j*dim1
}

// A method to set the stride size for column j
template <typename T>
size_t& DynamicRaggedDownArray<T>::stride(size_t j) const {
    return stride_[j];
}

//return size
template <typename T>
size_t DynamicRaggedDownArray<T>::size() const{
    return length_;
}

// overload operator () to access data as an array(i,j)
// Note: i = 0:stride(j), j = 0:N-1

template <typename T>
inline T& DynamicRaggedDownArray<T>::operator()(size_t i, size_t j) const {
    // Asserts
    assert(i < dim1_ && "i is out of dim1 bounds in DynamicRaggedDownArray");  // die if >= dim1
    assert(j < dim2_ && "j is out of dim2 bounds in DynamicRaggedDownArray");  // die if >= dim2
    assert(i < stride_[j] && "i is out of stride bounds in DynamicRaggedDownArray");  // die if >= stride
    
    return array_[i + j*dim1_];
}

//overload = operator
template <typename T>
inline DynamicRaggedDownArray<T>& DynamicRaggedDownArray<T>::operator= (const DynamicRaggedDownArray &temp)
{
    
    if( this != &temp) {
        dim1_ = temp.dim1_;
        dim2_ = temp.dim2_;
        length_ = temp.length_;
        if(stride_!=NULL)
          delete[] stride_;
        stride_ = new size_t[dim1_];
        for (int j = 0; j < dim2_; j++) {
            stride_[j] = temp.stride_[j];
        }
        if(array_!=NULL)
          delete[] array_;
        array_ = new T[length_];
        //copy contents
        for(int iter = 0; iter < length_; iter++)
          array_[iter] = temp.array_[iter];
    }
    
    return *this;
}

template <typename T>
inline T* DynamicRaggedDownArray<T>::pointer() const{
    return array_;
}

// Destructor
template <typename T>
DynamicRaggedDownArray<T>::~DynamicRaggedDownArray() {
    if(array_!=NULL)
      delete[] array_;
    if(stride_!=NULL)
      delete[] stride_;
}

//----end of DynamicRaggedDownArray class definitions-----




//13. SparseRowArray
template <typename T>
class SparseRowArray {
private:
    size_t *start_index_;
    size_t *column_index_;
    
    T * array_;
    
    size_t dim1_, length_;
    
public:
    // Default constructor
    SparseRowArray ();
    
    //--- 2D array access of a ragged right array ---
    
    // Overload constructor for a CArray
    SparseRowArray (CArray<size_t> &strides_array);
    
    // Overload constructor for a ViewCArray
    SparseRowArray (ViewCArray<size_t> &strides_array);
    
    // Overloaded constructor for a traditional array
    SparseRowArray (size_t *strides_array, size_t some_dim1);
    
    // A method to return the stride size
    size_t stride(size_t i) const;
    
    // A method to return the column index as array.column_index(i,j)
    size_t& column_index(size_t i, size_t j) const;
    
    // A method to access data as array.value(i,j),
    // where i=[0:N-1], j=[stride(i)]
    T& value(size_t i, size_t j) const;

    // A method to return the total size of the array
    size_t size() const;

    //return pointer
    T* pointer() const;

    //get row starts array
    size_t* get_starts() const;
    
    // Destructor
    ~SparseRowArray ();
}; 

//Default Constructor
template <typename T>
SparseRowArray<T>::SparseRowArray (){
    array_ = NULL;
    start_index_ = NULL;
    column_index_ = NULL;
    length_ = 0;
}
// Overloaded constructor
template <typename T>
SparseRowArray<T>::SparseRowArray (CArray<size_t> &strides_array) {
    // The length of the stride array is some_dim1;
    dim1_  = strides_array.size();
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[dim1_+1];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array(i);
        start_index_[i+1] = count;
    } // end for i
    
    length_ = count;
    array_ = new T[count];
    column_index_ = new size_t[count];
} 


// Overloaded constructor
template <typename T>
SparseRowArray<T>::SparseRowArray (ViewCArray<size_t> &strides_array) {
    // The length of the stride array is some_dim1;
    dim1_  = strides_array.size();
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[dim1_+1];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array(i);
        start_index_[i+1] = count;
    } // end for i
    
    length_ = count;
    array_ = new T[count];
    column_index_ = new size_t[count];
} 

// Overloaded constructor
template <typename T>
SparseRowArray<T>::SparseRowArray (size_t *strides_array, size_t dim1) {
    // The length of the stride array is some_dim1;
    dim1_ = dim1;
    
    // Create and initialize the starting index of the entries in the 1D array
    start_index_ = new size_t[dim1_+1];  // note the dim1+1
    start_index_[0] = 0; // the 1D array starts at 0
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    size_t count = 0;
    for (size_t i = 0; i < dim1_; i++){
        count += strides_array[i];
        start_index_[i+1] = count;
    } // end for i
    
    length_ = count;
    array_ = new T[count];
    column_index_ = new size_t[count];
} 


// A method to return the stride size
template <typename T>
size_t SparseRowArray<T>::stride(size_t i) const {
    return start_index_[i+1] - start_index_[i];
}

// A method to return the column index
template <typename T>
size_t& SparseRowArray<T>::column_index(size_t i, size_t j) const {
    // Get the 1D array index
    size_t start = start_index_[i];
    
    // Asserts
    assert(i < dim1_ && "i is out of dim1 bounds in SparseRowArray");  // die if >= dim1
    assert(j < stride(i) && "j is out of stride bounds in SparseRowArray");  // die if >= stride
    
    return column_index_[j + start];
}

// Access data as array.value(i,j), 
// where i=[0:N-1], j=[0:stride(i)]
template <typename T>
inline T& SparseRowArray<T>::value(size_t i, size_t j) const {
    // Get the 1D array index
    size_t start = start_index_[i];
    
    // Asserts
    assert(i < dim1_ && "i is out of dim1 bounds in sparseRowArray");  // die if >= dim1
    assert(j < stride(i) && "j is out of stride bounds in sparseRowArray");  // die if >= stride
    
    return array_[j + start];
} 

//return size
template <typename T>
size_t SparseRowArray<T>::size() const{
    return length_;
}

template <typename T>
inline T* SparseRowArray<T>::pointer() const{
    return array_;
}

template <typename T>
inline size_t* SparseRowArray<T>::get_starts() const{
    return start_index_;
}

// Destructor
template <typename T>
SparseRowArray<T>::~SparseRowArray() {
    if(array_!=NULL)
      delete[] array_;
    if(start_index_!=NULL)
      delete[] start_index_;
    if(column_index_!=NULL)
      delete[] column_index_;
}

//---- end of SparseRowArray class definitions-----



//14. SparseColArray
template <typename T>
class SparseColArray {

private:
	size_t *start_index_;
	size_t *row_index_;
	T * array_;

	size_t dim2_, length_;

public:

	//default constructor
	SparseColArray ();

	//constructor with CArray
	SparseColArray(CArray<size_t> &strides_array);

	//constructor with ViewCArray
	SparseColArray(ViewCArray<size_t> &strides_array);

	//constructor with regular array
	SparseColArray(size_t *strides_array, size_t some_dim1);

	//method return stride size
	size_t stride(size_t j) const;

	//method return row index ass array.row_index(i,j)
	size_t& row_index(size_t i, size_t j) const;

	//method access data as an array
	T& value(size_t i, size_t j) const;

    // A method to return the total size of the array
    size_t size() const;

    //return pointer
    T* pointer() const;

    //get row starts array
    size_t* get_starts() const;

	//destructor
	~SparseColArray();
};

//Default Constructor
template <typename T>
SparseColArray<T>::SparseColArray (){
    array_ = NULL;
    start_index_ = NULL;
    row_index_ = NULL;
    length_ = 0;
}
//overload constructor with CArray
template <typename T>
SparseColArray<T>::SparseColArray(CArray<size_t> &strides_array) {

	dim2_ = strides_array.size();

	start_index_ = new size_t[dim2_+1];
	start_index_[0] = 0;

	//loop over to find total length of the 1D array
	size_t count = 0;
	for(size_t j = 0; j < dim2_; j++) {
	  count+= strides_array(j);
	  start_index_[j+1] = count;
	}
    
    length_ = count;
	array_ = new T[count];
	row_index_ = new T[count];

} //end constructor with CArray


//overload constructor with ViewCArray
template <typename T>
SparseColArray<T>::SparseColArray(ViewCArray<size_t> &strides_array) {

	dim2_ = strides_array.size();

	//create and initialize starting index of 1D array
	start_index_ = new size_t[dim2_+1];
	start_index_[0] = 0;

	//loop over to find total length of 1D array
	size_t count = 0;
	for(size_t j = 0; j < dim2_ ; j++) {
	  count += strides_array(j);
	  start_index_[j+1] = count;
	}
    
    length_ = count;
	array_ = new T[count];
	row_index_ = new T[count];

} //end constructor

//overload constructor with traditional array
template <typename T>
SparseColArray<T>::SparseColArray(size_t *strides_array, size_t dim2) {

	dim2_ = dim2;

	//create and initialize the starting index 
	start_index_ = new size_t[dim2_ +1];
	start_index_[0] = 0;

	//loop over to find the total length of the 1D array
	size_t count = 0;
	for(size_t j = 0; j < dim2_; j++) {
	  count += strides_array[j];
	  start_index_[j+1] = count;
	}
    
    length_ = count;
	array_ = new T[count];
	row_index_ = new T[count];

} //end constructor

//method to return stride size
template <typename T>
size_t SparseColArray<T>::stride(size_t j) const{
	return start_index_[j+1] - start_index_[j];
}

//acces data ass arrow.row_index(i,j)
// where i = 0:stride(j), j = 0:N-1
template <typename T>
size_t& SparseColArray<T>::row_index(size_t i, size_t j) const {

	//get 1D array index
	size_t start = start_index_[j];

	//asserts to make sure we are in bounds
	assert(i < stride(j) && "i is out of stride bounnds in SparseColArray!");
	assert(j < dim2_ && "j is out of dim1 bounds in SparseColArray");

	return row_index_[i + start];

} //end row index method	


//access values as array.value(i,j)
// where i = 0:stride(j), j = 0:N-1
template <typename T>
T& SparseColArray<T>::value(size_t i, size_t j) const {

	size_t start = start_index_[j];

	//asserts
	assert(i < stride(j) && "i is out of stride boundns in SparseColArray");
	assert(j < dim2_ && "j is out of dim1 bounds in SparseColArray");

	return array_[i + start];

}

//return size
template <typename T>
size_t SparseColArray<T>::size() const{
    return length_;
}

template <typename T>
inline T* SparseColArray<T>::pointer() const{
    return array_;
}

template <typename T>
inline size_t* SparseColArray<T>::get_starts() const{
    return start_index_;
}

//destructor
template <typename T>
SparseColArray<T>::~SparseColArray() {
    if(array_!=NULL)
	  delete [] array_;
    if(start_index_!=NULL)
	  delete [] start_index_;
    if(row_index_!=NULL)
	  delete [] row_index_;
}

//----end SparseColArray----

//=======================================================================
//	end of standard MATAR data-types
//========================================================================

/*! \brief Kokkos version of the serial FArray class.
 *
 *  This is the Kokkos version of the serial FArray class. 
 *  Its usage is analagous to that of the serial FArray class, and it is to be
 *  used in Kokkos-specific code.
 */
#ifdef HAVE_KOKKOS
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
FArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::FArrayKokkos() {}

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
ViewFArrayKokkos<T>::ViewFArrayKokkos() {}

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
FMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::FMatrixKokkos() {}

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
ViewFMatrixKokkos<T>::ViewFMatrixKokkos() {}

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

    // Data member to access host view
    ViewFArray <T> host;

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
DFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DFArrayKokkos() {}

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
    
private:
    size_t dims_[7];
    size_t length_;
    size_t order_;  // tensor order (rank)
    TArray1D this_array_; 
    TArray1DHost this_array_host_; 
    T * temp_inp_array_;

public:
    DViewFArrayKokkos();
    
    DViewFArrayKokkos(T * inp_array, size_t dim0);

    DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1);

    DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2);

    DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3);

    DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4);

    DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, size_t dim5);

    DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, size_t dim5,
                 size_t dim6);
    
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
DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFArrayKokkos() {}

// Overloaded 1D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFArrayKokkos(T * inp_array, size_t dim0) {
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
    this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    // Create host ViewFArray. Note: inp_array and this_array_host_.data() are the same pointer 
    host = ViewFArray <T> (inp_array, dim0);
}

// Overloaded 2D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1) {
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
    this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    // Create host ViewFArray
    host = ViewFArray <T> (inp_array, dim0, dim1);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2) {
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
    this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    // Create host ViewFArray
    host = ViewFArray <T> (inp_array, dim0, dim1, dim2);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3) {
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
    this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    // Create host ViewFArray
    host = ViewFArray <T> (inp_array, dim0, dim1, dim2, dim3);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3, 
                              size_t dim4) {

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
    this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    // Create host ViewFArray
    host = ViewFArray <T> (inp_array, dim0, dim1, dim2, dim3, dim4);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3, 
                              size_t dim4, size_t dim5) {
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
    this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    // Create host ViewFArray
    host = ViewFArray <T> (inp_array, dim0, dim1, dim2, dim3, dim4, dim5);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3,
                              size_t dim4, size_t dim5,
                              size_t dim6) {
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
    this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
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
DFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DFMatrixKokkos() {}

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
    
private:
    size_t dims_[7];
    size_t length_;
    size_t order_;  // tensor order (rank)
    TArray1D this_matrix_; 
    TArray1DHost this_matrix_host_; 
    T * temp_inp_matrix_;

public:
    DViewFMatrixKokkos();
    
    DViewFMatrixKokkos(T * inp_matrix, size_t dim1);

    DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2);

    DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3);

    DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4);

    DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, size_t dim5);

    DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, size_t dim5, size_t dim6);

    DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, size_t dim5, size_t dim6,
                 size_t dim7);
    
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
DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFMatrixKokkos() {}

// Overloaded 1D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFMatrixKokkos(T * inp_matrix, size_t dim1) {
    
    dims_[0] = dim1;
    order_ = 1;
    length_ = dim1;
    // Create a 1D host view of the external allocation
    this_matrix_host_ = TArray1DHost(inp_matrix, length_);
    // Assign temp point to inp_matrix pointer that is passed in
    temp_inp_matrix_ = inp_matrix;
    // Create a device copy of that host view
    this_matrix_ = create_mirror_view_and_copy(ExecSpace(), this_matrix_host_);
    // Create host ViewFMatrix. Note: inp_matrix and this_matrix_host_.data() are the same pointer 
    host = ViewFMatrix <T> (inp_matrix, dim1);
}

// Overloaded 2D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    order_ = 2;
    length_ = (dim1 * dim2);
    // Create a 1D host view of the external allocation
    this_matrix_host_ = TArray1DHost(inp_matrix, length_);
    // Assign temp point to inp_matrix pointer that is passed in
    temp_inp_matrix_ = inp_matrix;
    // Create a device copy of that host view
    this_matrix_ = create_mirror_view_and_copy(ExecSpace(), this_matrix_host_);
    // Create host ViewFMatrix
    host = ViewFMatrix <T> (inp_matrix, dim1, dim2);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3) {
    
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
    this_matrix_ = create_mirror_view_and_copy(ExecSpace(), this_matrix_host_);
    // Create host ViewFMatrix
    host = ViewFMatrix <T> (inp_matrix, dim1, dim2, dim3);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4) {
    
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
    this_matrix_ = create_mirror_view_and_copy(ExecSpace(), this_matrix_host_);
    // Create host ViewFMatrix
    host = ViewFMatrix <T> (inp_matrix, dim1, dim2, dim3, dim4);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4, 
                              size_t dim5) {
    
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
    this_matrix_ = create_mirror_view_and_copy(ExecSpace(), this_matrix_host_);
    // Create host ViewFMatrix
    host = ViewFMatrix <T> (inp_matrix, dim1, dim2, dim3, dim4, dim5);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4, 
                              size_t dim5, size_t dim6) {
    
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
    this_matrix_ = create_mirror_view_and_copy(ExecSpace(), this_matrix_host_);
    // Create host ViewFMatrix
    host = ViewFMatrix <T> (inp_matrix, dim1, dim2, dim3, dim4, dim5, dim6);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewFMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewFMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4,
                              size_t dim5, size_t dim6,
                              size_t dim7) {
    
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
    this_matrix_ = create_mirror_view_and_copy(ExecSpace(), this_matrix_host_);
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
CArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::CArrayKokkos() {}

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
ViewCArrayKokkos<T>::ViewCArrayKokkos() {}

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
CMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::CMatrixKokkos() {}

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
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(){ }

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

template <typename T>
class ViewCArrayMat {

private:
	size_t dim1_ = 1;
	size_t dim2_ = 1;
	size_t dim3_ = 1;
	size_t dim4_ = 1;
	size_t dim5_ = 1;
	size_t dim6_ = 1;
	size_t length_ = 1;
	T* this_matrix_;

public:

	KOKKOS_FUNCTION
		ViewCArrayMat();

	//copy constructor
	KOKKOS_FUNCTION
		ViewCArrayMat(const ViewCArrayMat& temp);

	//View Constructors
	KOKKOS_FUNCTION
		ViewCArrayMat(T* some_matrix, size_t dim1);

	KOKKOS_FUNCTION
		ViewCArrayMat(T* some_matrix, size_t dim1, size_t dim2);

	KOKKOS_FUNCTION
		ViewCArrayMat(T* some_matrix, size_t dim1, size_t dim2, size_t dim3);

	KOKKOS_FUNCTION
		ViewCArrayMat(T* some_matrix, size_t dim1, size_t dim2, size_t dim3,
			size_t dim4);

	KOKKOS_FUNCTION
		ViewCArrayMat(T* some_matrix, size_t dim1, size_t dim2, size_t dim3,
			size_t dim4, size_t dim5);

	KOKKOS_FUNCTION
		ViewCArrayMat(T* some_matrix, size_t dim1, size_t dim2, size_t dim3,
			size_t dim4, size_t dim5, size_t dim6);

	KOKKOS_FUNCTION
		ViewCArrayMat(T* some_matrix, size_t dim1, size_t dim2, size_t dim3,
			size_t dim4, size_t dim5, size_t dim6, size_t dim7);

	//Reverse Constructors to quickly convert between Array and Matrix
	KOKKOS_FUNCTION
		ViewCArrayMat(size_t length, size_t dim1, size_t dim2, size_t dim3,
			size_t dim4, size_t dim5, size_t dim6, T* some_matrix);


	KOKKOS_FORCEINLINE_FUNCTION
		T& operator()(size_t i) const;

	KOKKOS_FORCEINLINE_FUNCTION
		T& operator()(size_t i, size_t j) const;

	KOKKOS_FORCEINLINE_FUNCTION
		T& operator()(size_t i, size_t j, size_t k) const;

	KOKKOS_FORCEINLINE_FUNCTION
		T& operator()(size_t i, size_t j, size_t k, size_t l) const;

	KOKKOS_FORCEINLINE_FUNCTION
		T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

	KOKKOS_FORCEINLINE_FUNCTION
		T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const;

	KOKKOS_FORCEINLINE_FUNCTION
		T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t o) const;



	//overload = operator
	KOKKOS_FUNCTION
		ViewCArrayMat& operator=    (const ViewCArrayMat& temp);
	KOKKOS_FUNCTION
		ViewCArrayMat& set(const ViewCArrayMat& temp);
	KOKKOS_FUNCTION
		ViewCArrayMat& operator=    (const       T& temp);
	KOKKOS_FUNCTION
		ViewCArrayMat& operator=    (T* temp);

	KOKKOS_FUNCTION
		operator T* ();

	KOKKOS_FUNCTION
		size_t size();

	size_t extent();

	KOKKOS_FUNCTION
		~ViewCArrayMat();

}; // End of ViewCArrayMat

// Default constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayMat<T>::ViewCArrayMat() { }


//Copy Constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayMat<T>::ViewCArrayMat(const ViewCArrayMat& temp)
{
	dim1_ = temp.dim1_;
	dim2_ = temp.dim2_;
	dim3_ = temp.dim3_;
	dim4_ = temp.dim4_;
	dim5_ = temp.dim5_;
	dim6_ = temp.dim6_;
	length_ = temp.length_;

	this_matrix_ = temp.this_matrix_;
	printf("copy constructor has been called but its not really a copy constructor\n");																			
}

// Overloaded 1D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayMat<T>::ViewCArrayMat(T* some_matrix, size_t dim1) {
	length_ = dim1;
	this_matrix_ = some_matrix;
}

// Overloaded 2D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayMat<T>::ViewCArrayMat(T* some_matrix, size_t dim1,
	size_t dim2) {
	dim1_ = dim2;
	length_ = (dim1_ * dim1);
	this_matrix_ = some_matrix;
}

// Overloaded 3D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayMat<T>::ViewCArrayMat(T* some_matrix, size_t dim1, size_t dim2,
	size_t dim3) {
	dim2_ = dim3;
	dim1_ = dim2_ * dim2;
	length_ = (dim1_ * dim1);
	this_matrix_ = some_matrix;
}

// Overloaded 4D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayMat<T>::ViewCArrayMat(T* some_matrix, size_t dim1, size_t dim2,
	size_t dim3, size_t dim4) {
	dim3_ = dim4;
	dim2_ = dim3_ * dim3;
	dim1_ = dim2_ * dim2;
	length_ = (dim1_ * dim1);
	this_matrix_ = some_matrix;
}

// Overloaded 5D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayMat<T>::ViewCArrayMat(T* some_matrix, size_t dim1, size_t dim2,
	size_t dim3, size_t dim4, size_t dim5) {
	dim4_ = dim5;
	dim3_ = dim4_ * dim4;
	dim2_ = dim3_ * dim3;
	dim1_ = dim2_ * dim2;

	length_ = (dim1_ * dim1);
	this_matrix_ = some_matrix;
}

// Overloaded 6D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayMat<T>::ViewCArrayMat(T* some_matrix, size_t dim1, size_t dim2,
	size_t dim3, size_t dim4, size_t dim5,
	size_t dim6) {
	dim5_ = dim6;
	dim4_ = dim5_ * dim5;
	dim3_ = dim4_ * dim4;
	dim2_ = dim3_ * dim3;
	dim1_ = dim2_ * dim2;

	length_ = (dim1_ * dim1);
	this_matrix_ = some_matrix;
}

// Overloaded 7D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayMat<T>::ViewCArrayMat(T* some_matrix, size_t dim1, size_t dim2,
	size_t dim3, size_t dim4, size_t dim5,
	size_t dim6, size_t dim7) {

	dim6_ = dim7;
	dim5_ = dim6_ * dim6;
	dim4_ = dim5_ * dim5;
	dim3_ = dim4_ * dim4;
	dim2_ = dim3_ * dim3;
	dim1_ = dim2_ * dim2;
	length_ = (dim1_ * dim1);
	this_matrix_ = some_matrix;
}

// Reverse Constructor to convert from Array to Matrix quickly
template <typename T>
KOKKOS_FUNCTION
ViewCArrayMat<T>::ViewCArrayMat(size_t length, size_t dim1, size_t dim2, size_t dim3, size_t dim4, size_t dim5, size_t dim6, T* some_matrix) {
	length_ = length;
	dim1_ = dim1;
	dim2_ = dim2;
	dim3_ = dim3;
	dim4_ = dim4;
	dim5_ = dim5;
	dim6_ = dim6;
	this_matrix_ = some_matrix;
}



template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
T& ViewCArrayMat<T>::operator()(size_t i) const {
	//assert(i < length_ && "i is out of bounds in ViewCArrayMat 1D!");
	return this_matrix_[(i)];
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
T& ViewCArrayMat<T>::operator()(size_t i, size_t j) const {
	//assert(i < length_/dim1_ && "i is out of bounds in ViewCArrayMat 2D!");
	//assert(j < dim1_/dim2_ && "j is out of bounds in ViewCArrayMat 2D!");
	return this_matrix_[(i) * dim1_ + (j)];
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
T& ViewCArrayMat<T>::operator()(size_t i, size_t j, size_t k) const {
	//assert(i < length_ / dim1_ && "i is out of bounds in ViewCArrayMat 3D!");
	//assert(j < dim1_ / dim2_ && "j is out of bounds in ViewCArrayMat 3D!");
	//assert(k < dim2_/dim3_ && "k is out of bounds in ViewCArrayMat 3D!");
	return this_matrix_[(i) * dim1_ + (j) * dim2_ + k];
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
T& ViewCArrayMat<T>::operator()(size_t i, size_t j, size_t k, size_t l) const {
	//assert(i < length_/dim1_ && "i is out of bounds in ViewCArrayMat 4D!");
	//assert(j < dim1_ / dim2_ && "j is out of bounds in ViewCArrayMat 4D!");
	//assert(k < dim2_/dim3_ && "k is out of bounds in ViewCArrayMat 4D!");
	//assert(l < dim3_/dim4_ && "l is out of bounds in ViewCArrayMat 4D!");
	return this_matrix_[(i) * dim1_ + (j) * dim2_ + (k) * dim3_ + (l)];
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
T& ViewCArrayMat<T>::operator()(size_t i, size_t j, size_t k, size_t l,
	size_t m) const {
	//assert(i >= 1 && i <= dim1_ && "i is out of bounds for ViewCArrayMat 5D!");
	//assert(j >= 1 && j <= dim2_ && "j is out of bounds for ViewCArrayMat 5D!");
	//assert(k >= 1 && k <= dim3_ && "k is out of bounds for ViewCArrayMat 5D!");
	//assert(l >= 1 && l <= dim4_ && "l is out of bounds for ViewCArrayMat 5D!");
	//assert(m >= 1 && m <= dim5_ && "m is out of bounds for ViewCArrayMat 5D!");
	return this_matrix_[(i) * dim1_ + (j) * dim2_ + (k) * dim3_ + (l) * dim4_ + (m)];
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
T& ViewCArrayMat<T>::operator()(size_t i, size_t j, size_t k, size_t l,
	size_t m, size_t n) const {
	//assert(i >= 1 && i <= dim1_ && "i is out of bounds for ViewCArrayMat 6D!");
	//assert(j >= 1 && j <= dim2_ && "j is out of bounds for ViewCArrayMat 6D!");
	//assert(k >= 1 && k <= dim3_ && "k is out of bounds for ViewCArrayMat 6D!");
	//assert(l >= 1 && l <= dim4_ && "l is out of bounds for ViewCArrayMat 6D!");
	//assert(m >= 1 && m <= dim5_ && "m is out of bounds for ViewCArrayMat 6D!");
	//assert(n >= 1 && n <= dim6_ && "n is out of bounds for ViewCArrayMat 6D!");
	return this_matrix_[(i) * dim1_ + (j) * dim2_ + (k) * dim3_ + (l) * dim4_ + (m) * dim5_ + (n)];
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
T& ViewCArrayMat<T>::operator()(size_t i, size_t j, size_t k, size_t l,
	size_t m, size_t n, size_t o) const {
	//assert(i >= 1 && i <= dim1_ && "i is out of bounds for ViewCArrayMat 7D!");
	//assert(j >= 1 && j <= dim2_ && "j is out of bounds for ViewCArrayMat 7D!");
	//assert(k >= 1 && k <= dim3_ && "k is out of bounds for ViewCArrayMat 7D!");
	//assert(l >= 1 && l <= dim4_ && "l is out of bounds for ViewCArrayMat 7D!");
	//assert(m >= 1 && m <= dim5_ && "m is out of bounds for ViewCArrayMat 7D!");
	//assert(n >= 1 && n <= dim6_ && "n is out of bounds for ViewCArrayMat 7D!");
	//assert(o >= 1 && o <= dim7_ && "o is out of bounds for ViewCArrayMat 7D!");
	return this_matrix_[(i) * dim1_ + (j) * dim2_ + (k) * dim3_ + (l) * dim4_ + (m) * dim5_ + (n) * dim6_ + (o)];
}


template <typename T>
KOKKOS_FUNCTION
size_t ViewCArrayMat<T>::size() {
	return length_;
}

template <typename T>
size_t ViewCArrayMat<T>::extent() {
	return length_;
}


template <typename T>
KOKKOS_FUNCTION ViewCArrayMat<T>& ViewCArrayMat<T>::operator=  (const ViewCArrayMat <T>& temp)
{

	memcpy(&this_matrix_[0], &temp.this_matrix_[0], sizeof(T) * length_);

	//for (size_t i = 0; i < length_; i++) { this_matrix_[i] = temp.this_matrix_[i]; }

	return *this;
}

template <typename T> KOKKOS_FUNCTION
ViewCArrayMat<T>& ViewCArrayMat<T>::set(const ViewCArrayMat <T>& temp)
{
	if (this != &temp) {
		dim1_ = temp.dim1_;
		dim2_ = temp.dim2_;
		dim3_ = temp.dim3_;
		dim4_ = temp.dim4_;
		dim5_ = temp.dim5_;
		dim6_ = temp.dim6_;

		length_ = temp.length_;

		this_matrix_ = temp.this_matrix_;
	}
	return *this;
}

template <typename T> KOKKOS_FUNCTION
ViewCArrayMat<T>& ViewCArrayMat<T>::operator=  (const T& temp) {

	for (size_t i = 0; i < length_; i++) this_matrix_[i] = temp;

	return *this;
}
template <typename T> KOKKOS_FUNCTION
ViewCArrayMat<T>& ViewCArrayMat<T>::operator=  (T* temp) {

	for (size_t i = 0; i < length_; i++) this_matrix_[i] = temp[i];

	return *this;
}

template <typename T> KOKKOS_FUNCTION
ViewCArrayMat<T>::operator T* () {

	return &this_matrix_[0];
}



template <typename T>
KOKKOS_FUNCTION
ViewCArrayMat<T>::~ViewCArrayMat() {}

////////////////////////////////////////////////////////////////////////////////
// End of ViewCArrayMat
////////////////////////////////////////////////////////////////////////////////




template <typename T>
class ViewCMatrixMat {

private:
    size_t dim1_ = 1;
    size_t dim2_ = 1;
    size_t dim3_ = 1;
    size_t dim4_ = 1;
    size_t dim5_ = 1;
    size_t dim6_ = 1;
    size_t length_ = 1;
	T* this_matrix_;

public:

	KOKKOS_FUNCTION
		ViewCMatrixMat();

	//copy constructor
	KOKKOS_FUNCTION
		ViewCMatrixMat(const ViewCMatrixMat& temp);

	//View Constructors
	KOKKOS_FUNCTION
		ViewCMatrixMat(T* some_matrix, size_t dim1);

	KOKKOS_FUNCTION
		ViewCMatrixMat(T* some_matrix, size_t dim1, size_t dim2);

	KOKKOS_FUNCTION
		ViewCMatrixMat(T* some_matrix, size_t dim1, size_t dim2, size_t dim3);

	KOKKOS_FUNCTION
		ViewCMatrixMat(T* some_matrix, size_t dim1, size_t dim2, size_t dim3,
			size_t dim4);

	KOKKOS_FUNCTION
		ViewCMatrixMat(T* some_matrix, size_t dim1, size_t dim2, size_t dim3,
			size_t dim4, size_t dim5);

	KOKKOS_FUNCTION
		ViewCMatrixMat(T* some_matrix, size_t dim1, size_t dim2, size_t dim3,
			size_t dim4, size_t dim5, size_t dim6);

	KOKKOS_FUNCTION
		ViewCMatrixMat(T* some_matrix, size_t dim1, size_t dim2, size_t dim3,
			size_t dim4, size_t dim5, size_t dim6, size_t dim7);

	//Reverse Constructor for Direct Assignment
	KOKKOS_FUNCTION
		ViewCMatrixMat(size_t length, size_t dim1, size_t dim2, size_t dim3,
			size_t dim4, size_t dim5, size_t dim6, T* some_matrix);

	KOKKOS_FORCEINLINE_FUNCTION
		T& operator()(size_t i) const;

	KOKKOS_FORCEINLINE_FUNCTION
		T& operator()(size_t i, size_t j) const;

	KOKKOS_FORCEINLINE_FUNCTION
		T& operator()(size_t i, size_t j, size_t k) const;

	KOKKOS_FORCEINLINE_FUNCTION
		T& operator()(size_t i, size_t j, size_t k, size_t l) const;

	KOKKOS_FORCEINLINE_FUNCTION
		T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

	KOKKOS_FORCEINLINE_FUNCTION
		T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const;

	KOKKOS_FORCEINLINE_FUNCTION
		T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t o) const;



	//overload = operator
	KOKKOS_FUNCTION
		ViewCMatrixMat& operator=    (const ViewCMatrixMat& temp);
	KOKKOS_FUNCTION
		ViewCMatrixMat& set(const ViewCMatrixMat& temp);
	KOKKOS_FUNCTION
		ViewCMatrixMat& operator=    (const       T temp);
	KOKKOS_FUNCTION
		ViewCMatrixMat& operator=    (T* temp);

	KOKKOS_FUNCTION
		operator T* ();

	KOKKOS_FUNCTION
		operator ViewCArrayMat<T> ();

	KOKKOS_FUNCTION
		size_t size();

	size_t extent();

	KOKKOS_FUNCTION
		~ViewCMatrixMat();

}; // End of ViewCMatrixMat

// Default constructor
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixMat<T>::ViewCMatrixMat() { }


//Copy Constructor
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixMat<T>::ViewCMatrixMat(const ViewCMatrixMat& temp)
{
	dim1_ = temp.dim1_;
	dim2_ = temp.dim2_;
	dim3_ = temp.dim3_;
	dim4_ = temp.dim4_;
	dim5_ = temp.dim5_;
	dim6_ = temp.dim6_;
	length_ = temp.length_;

	this_matrix_ = temp.this_matrix_;
}

// Overloaded 1D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixMat<T>::ViewCMatrixMat(T* some_matrix, size_t dim1) {
	length_ = dim1;
	this_matrix_ = some_matrix;
}

// Overloaded 2D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixMat<T>::ViewCMatrixMat(T* some_matrix, size_t dim1,
	size_t dim2) {
	dim1_ = dim2;
	length_ = (dim1_ * dim1);
	this_matrix_ = some_matrix;
}

// Overloaded 3D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixMat<T>::ViewCMatrixMat(T* some_matrix, size_t dim1, size_t dim2,
	size_t dim3) {

	dim2_ = dim3;
	dim1_ = dim2_ * dim2;
	length_ = (dim1_ * dim1);
	this_matrix_ = some_matrix;
}

// Overloaded 4D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixMat<T>::ViewCMatrixMat(T* some_matrix, size_t dim1, size_t dim2,
	size_t dim3, size_t dim4) {
	
	dim3_ = dim4;
	dim2_ = dim3_ * dim3;
	dim1_ = dim2_ * dim2;
	length_ = (dim1_ * dim1);
	this_matrix_ = some_matrix;
}

// Overloaded 5D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixMat<T>::ViewCMatrixMat(T* some_matrix, size_t dim1, size_t dim2,
	size_t dim3, size_t dim4, size_t dim5) {
	dim4_ = dim5;
	dim3_ = dim4_ * dim4;
	dim2_ = dim3_ * dim3;
	dim1_ = dim2_ * dim2;
	
	length_ = (dim1_ * dim1);
	this_matrix_ = some_matrix;
}

// Overloaded 6D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixMat<T>::ViewCMatrixMat(T* some_matrix, size_t dim1, size_t dim2,
	size_t dim3, size_t dim4, size_t dim5,
	size_t dim6) {
	dim5_ = dim6;
	dim4_ = dim5_ * dim5;
	dim3_ = dim4_ * dim4;
	dim2_ = dim3_ * dim3;
	dim1_ = dim2_ * dim2;
	
	length_ = (dim1_ * dim1);
	this_matrix_ = some_matrix;
}

// Overloaded 7D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixMat<T>::ViewCMatrixMat(T* some_matrix, size_t dim1, size_t dim2,
	size_t dim3, size_t dim4, size_t dim5,
	size_t dim6, size_t dim7) {

	dim6_ = dim7;
	dim5_ = dim6_ * dim6;
	dim4_ = dim5_ * dim5;
	dim3_ = dim4_ * dim4;
	dim2_ = dim3_ * dim3;
	dim1_ = dim2_ * dim2;
	length_ = (dim1_ * dim1);
	this_matrix_ = some_matrix;
}



// Reverse Constructor for direct Assignment
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixMat<T>::ViewCMatrixMat(size_t length, size_t dim1, size_t dim2, size_t dim3, size_t dim4, size_t dim5, size_t dim6, T* some_matrix) {
	length_ = length;
	dim1_ = dim1;
	dim2_ = dim2;
	dim3_ = dim3;
	dim4_ = dim4;
	dim5_ = dim5;
	dim6_ = dim6;
	this_matrix_ = some_matrix;
}


template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
T& ViewCMatrixMat<T>::operator()(size_t i) const {
	//assert(i >= 1 && i <= length_ && "i is out of bounds in ViewCMatrixMat 1D!");
	return this_matrix_[(i - 1)];
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
T& ViewCMatrixMat<T>::operator()(size_t i, size_t j) const {
	//assert(i >= 1 && i <= length_/dim1_ && "i is out of bounds in ViewCMatrixMat 2D!");
	//assert(j >= 1 && j <= dim1_/dim2_ && "j is out of bounds in ViewCMatrixMat 2D!");
	return this_matrix_[(i - 1) * dim1_ + (j - 1)];
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
T& ViewCMatrixMat<T>::operator()(size_t i, size_t j, size_t k) const {
	//assert(i >= 1 && i <= length_/dim1_ && "i is out of bounds in ViewCMatrixMat 3D!");
	//assert(j >= 1 && j <= dim1_/dim2_ && "j is out of bounds in ViewCMatrixMat 3D!");
	//assert(k >= 1 && k <= dim2_/dim3_ && "k is out of bounds in ViewCMatrixMat 3D!");
	return this_matrix_[(i - 1) * dim1_ + (j - 1) * dim2_ + k - 1];
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
T& ViewCMatrixMat<T>::operator()(size_t i, size_t j, size_t k, size_t l) const {
	//assert(i >= 1 && i <= length_/dim1_ && "i is out of bounds in ViewCMatrixMat 4D!");
	//assert(j >= 1 && j <= dim1_/dim2_ && "j is out of bounds in ViewCMatrixMat 4D!");
	//assert(k >= 1 && k <= dim2_/dim3_ && "k is out of bounds in ViewCMatrixMat 4D!");
	//assert(l >= 1 && l <= dim3_/dim4_ && "l is out of bounds in ViewCMatrixMat 4D!");
	return this_matrix_[(i - 1) * dim1_ + (j - 1) * dim2_ + (k - 1) * dim3_ + (l - 1)];
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
T& ViewCMatrixMat<T>::operator()(size_t i, size_t j, size_t k, size_t l,
	size_t m) const {
	//assert(i >= 1 && i <= dim1_ && "i is out of bounds for ViewCMatrixMat 5D!");
	//assert(j >= 1 && j <= dim2_ && "j is out of bounds for ViewCMatrixMat 5D!");
	//assert(k >= 1 && k <= dim3_ && "k is out of bounds for ViewCMatrixMat 5D!");
	//assert(l >= 1 && l <= dim4_ && "l is out of bounds for ViewCMatrixMat 5D!");
	//assert(m >= 1 && m <= dim5_ && "m is out of bounds for ViewCMatrixMat 5D!");
	return this_matrix_[(i - 1) * dim1_ + (j - 1) * dim2_ + (k - 1) * dim3_ + (l - 1) * dim4_ + (m - 1)];
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
T& ViewCMatrixMat<T>::operator()(size_t i, size_t j, size_t k, size_t l,
	size_t m, size_t n) const {
	//assert(i >= 1 && i <= dim1_ && "i is out of bounds for ViewCMatrixMat 6D!");
	//assert(j >= 1 && j <= dim2_ && "j is out of bounds for ViewCMatrixMat 6D!");
	//assert(k >= 1 && k <= dim3_ && "k is out of bounds for ViewCMatrixMat 6D!");
	//assert(l >= 1 && l <= dim4_ && "l is out of bounds for ViewCMatrixMat 6D!");
	//assert(m >= 1 && m <= dim5_ && "m is out of bounds for ViewCMatrixMat 6D!");
	//assert(n >= 1 && n <= dim6_ && "n is out of bounds for ViewCMatrixMat 6D!");
	return this_matrix_[(i - 1) * dim1_ + (j - 1) * dim2_ + (k - 1) * dim3_ + (l - 1) * dim4_ + (m - 1) * dim5_ + (n - 1)];
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
T& ViewCMatrixMat<T>::operator()(size_t i, size_t j, size_t k, size_t l,
	size_t m, size_t n, size_t o) const {
	//assert(i >= 1 && i <= dim1_ && "i is out of bounds for ViewCMatrixMat 7D!");
	//assert(j >= 1 && j <= dim2_ && "j is out of bounds for ViewCMatrixMat 7D!");
	//assert(k >= 1 && k <= dim3_ && "k is out of bounds for ViewCMatrixMat 7D!");
	//assert(l >= 1 && l <= dim4_ && "l is out of bounds for ViewCMatrixMat 7D!");
	//assert(m >= 1 && m <= dim5_ && "m is out of bounds for ViewCMatrixMat 7D!");
	//assert(n >= 1 && n <= dim6_ && "n is out of bounds for ViewCMatrixMat 7D!");
	//assert(o >= 1 && o <= dim7_ && "o is out of bounds for ViewCMatrixMat 7D!");
	return this_matrix_[(i - 1) * dim1_ + (j - 1) * dim2_ + (k - 1) * dim3_ + (l - 1) * dim4_ + (m - 1) * dim5_ + (n - 1) * dim6_ + (o - 1)];
}


template <typename T>
KOKKOS_FUNCTION
size_t ViewCMatrixMat<T>::size() {
	return length_;
}

template <typename T>
size_t ViewCMatrixMat<T>::extent() {
	return length_;
}


template <typename T>
KOKKOS_FUNCTION ViewCMatrixMat<T>& ViewCMatrixMat<T>::operator=  (const ViewCMatrixMat <T>& temp)
{

	//memcpy(&this_matrix_[0], &temp.this_matrix_[0], sizeof(T) * length_);

	for (size_t i = 0; i < length_; i++) { this_matrix_[i] = temp.this_matrix_[i]; }

	return *this;
}

template <typename T> KOKKOS_FUNCTION
ViewCMatrixMat<T>& ViewCMatrixMat<T>::set(const ViewCMatrixMat <T>& temp)
{
	if (this != &temp) {
		dim1_ = temp.dim1_;
		dim2_ = temp.dim2_;
		dim3_ = temp.dim3_;
		dim4_ = temp.dim4_;
		dim5_ = temp.dim5_;
		dim6_ = temp.dim6_;

		length_ = temp.length_;

		this_matrix_ = temp.this_matrix_;
	}
	return *this;
}

template <typename T> KOKKOS_FUNCTION
ViewCMatrixMat<T>& ViewCMatrixMat<T>::operator=  (const T temp) {

	//memset(&this_matrix_[0], temp, (sizeof(T) * length_));

	for (size_t i = 0; i < length_; i++) this_matrix_[i] = temp;

	return *this;
}
template <typename T> KOKKOS_FUNCTION
ViewCMatrixMat<T>& ViewCMatrixMat<T>::operator=  (T* temp) {

	//memcpy(&this_matrix_[0], &temp[0], sizeof(T) * length_);

	for (size_t i = 0; i < length_; i++) this_matrix_[i] = temp[i];

	return *this;
}

template <typename T> KOKKOS_FUNCTION
ViewCMatrixMat<T>::operator T* () {

	return &this_matrix_[0];
}


template <typename T> KOKKOS_FUNCTION
ViewCMatrixMat<T>::operator ViewCArrayMat<T> () {

	ViewCArrayMat<T>  out(length_, dim1_, dim2_, dim3_, dim4_, dim5_, dim6_, &this_matrix_[0]);
	
	return out;
}

template <typename T>
KOKKOS_FUNCTION
ViewCMatrixMat<T>::~ViewCMatrixMat() {}

////////////////////////////////////////////////////////////////////////////////
// End of ViewCMatrixMat
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

    // Data member to access host view
    ViewCArray <T> host;

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
DCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DCArrayKokkos() {}

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
    
    DViewCArrayKokkos(T * inp_array, size_t dim0);

    DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1);

    DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2);

    DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3);

    DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4);

    DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, size_t dim5);

    DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1, size_t dim2,
                 size_t dim3, size_t dim4, size_t dim5,
                 size_t dim6);
    
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
DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCArrayKokkos() {}

// Overloaded 1D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCArrayKokkos(T * inp_array, size_t dim0) {
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
    this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    // Create host ViewCArray. Note: inp_array and this_array_host_.data() are the same pointer 
    host = ViewCArray <T> (inp_array, dim0);
}

// Overloaded 2D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1) {
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
    this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    // Create host ViewCArray
    host = ViewCArray <T> (inp_array, dim0, dim1);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2) {
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
    this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    // Create host ViewCArray
    host = ViewCArray <T> (inp_array, dim0, dim1, dim2);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3) {
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
    this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    // Create host ViewCArray
    host = ViewCArray <T> (inp_array, dim0, dim1, dim2, dim3);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3, 
                              size_t dim4) {

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
    this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    // Create host ViewCArray
    host = ViewCArray <T> (inp_array, dim0, dim1, dim2, dim3, dim4);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3, 
                              size_t dim4, size_t dim5) {
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
    this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
    // Create host ViewCArray
    host = ViewCArray <T> (inp_array, dim0, dim1, dim2, dim3, dim4, dim5);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCArrayKokkos(T * inp_array, size_t dim0, size_t dim1,
                              size_t dim2, size_t dim3,
                              size_t dim4, size_t dim5,
                              size_t dim6) {
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
    this_array_ = create_mirror_view_and_copy(ExecSpace(), this_array_host_);
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

    // Data member to access host view
    ViewCMatrix <T> host;

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
DCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DCMatrixKokkos() {}

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
    
private:
    size_t dims_[7];
    size_t length_;
    size_t order_;  // tensor order (rank)
    TArray1D this_matrix_; 
    TArray1DHost this_matrix_host_; 
    T * temp_inp_matrix_;

public:
    DViewCMatrixKokkos();
    
    DViewCMatrixKokkos(T * inp_matrix, size_t dim1);

    DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2);

    DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3);

    DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4);

    DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, size_t dim5);

    DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, size_t dim5, size_t dim6);

    DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2, size_t dim3,
                 size_t dim4, size_t dim5, size_t dim6,
                 size_t dim7);
    
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
DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCMatrixKokkos() {}

// Overloaded 1D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCMatrixKokkos(T * inp_matrix, size_t dim1) {
    
    dims_[0] = dim1;
    order_ = 1;
    length_ = dim1;
    // Create a 1D host view of the external allocation
    this_matrix_host_ = TArray1DHost(inp_matrix, length_);
    // Assign temp point to inp_matrix pointer that is passed in
    temp_inp_matrix_ = inp_matrix;
    // Create a device copy of that host view
    this_matrix_ = create_mirror_view_and_copy(ExecSpace(), this_matrix_host_);
    // Create host ViewCMatrix. Note: inp_matrix and this_matrix_host_.data() are the same pointer 
    host = ViewCMatrix <T> (inp_matrix, dim1);
}

// Overloaded 2D constructor
template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2) {
    
    dims_[0] = dim1;
    dims_[1] = dim2;
    order_ = 2;
    length_ = (dim1 * dim2);
    // Create a 1D host view of the external allocation
    this_matrix_host_ = TArray1DHost(inp_matrix, length_);
    // Assign temp point to inp_matrix pointer that is passed in
    temp_inp_matrix_ = inp_matrix;
    // Create a device copy of that host view
    this_matrix_ = create_mirror_view_and_copy(ExecSpace(), this_matrix_host_);
    // Create host ViewCMatrix
    host = ViewCMatrix <T> (inp_matrix, dim1, dim2);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3) {
    
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
    this_matrix_ = create_mirror_view_and_copy(ExecSpace(), this_matrix_host_);
    // Create host ViewCMatrix
    host = ViewCMatrix <T> (inp_matrix, dim1, dim2, dim3);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4) {
    
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
    this_matrix_ = create_mirror_view_and_copy(ExecSpace(), this_matrix_host_);
    // Create host ViewCMatrix
    host = ViewCMatrix <T> (inp_matrix, dim1, dim2, dim3, dim4);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4, 
                              size_t dim5) {
    
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
    this_matrix_ = create_mirror_view_and_copy(ExecSpace(), this_matrix_host_);
    // Create host ViewCMatrix
    host = ViewCMatrix <T> (inp_matrix, dim1, dim2, dim3, dim4, dim5);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4, 
                              size_t dim5, size_t dim6) {
    
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
    this_matrix_ = create_mirror_view_and_copy(ExecSpace(), this_matrix_host_);
    // Create host ViewCMatrix
    host = ViewCMatrix <T> (inp_matrix, dim1, dim2, dim3, dim4, dim5, dim6);
}

template <typename T, typename Layout, typename ExecSpace, typename MemoryTraits>
DViewCMatrixKokkos<T,Layout,ExecSpace,MemoryTraits>::DViewCMatrixKokkos(T * inp_matrix, size_t dim1, size_t dim2,
                              size_t dim3, size_t dim4,
                              size_t dim5, size_t dim6,
                              size_t dim7) {
    
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
    this_matrix_ = create_mirror_view_and_copy(ExecSpace(), this_matrix_host_);
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
RaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::RaggedRightArrayKokkos() {}

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
      start_index_((i) = 0;
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
RaggedRightArrayofVectorsKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::RaggedRightArrayofVectorsKokkos() {}

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
      start_index_((i) = 0;
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
    assert(j < vector_dim_ && "k is out of vector_dim bounds in RaggedRightArrayKokkos");  // die if >= vector_dim
    
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
RaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits,ILayout>::RaggedDownArrayKokkos() {}

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
      start_index_((i) = 0;
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
DynamicRaggedRightArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DynamicRaggedRightArrayKokkos () {}

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
      strides_((i) = 0;
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
DynamicRaggedDownArrayKokkos<T,Layout,ExecSpace,MemoryTraits>::DynamicRaggedDownArrayKokkos () {}

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
      strides_((i) = 0;
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


#endif







#endif // MATAR_H
