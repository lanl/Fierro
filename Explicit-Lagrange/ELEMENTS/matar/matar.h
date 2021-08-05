#ifndef MATAR_H
#define MATAR_H

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

#include <assert.h>
//#include "kokkos_alias.h"

#ifdef HAVE_KOKKOS
#include <Kokkos_Core.hpp>

//MACROS to make the code less scary
#define kmalloc(size) ( Kokkos::kokkos_malloc<MemSpace>(size) )
#define kfree(pnt)        (  Kokkos::kokkos_free(pnt) ) 
#define ProfileRegionStart  ( Kokkos::Profiling::pushRegion )
#define ProfileRegionEnd  ( Kokkos::Profiling::popRegion )

using real_t = double;
using u_int  = unsigned int;

#ifdef HAVE_CUDA
//using UVMMemSpace     = Kokkos::CudaUVMSpace;
using MemSpace        = Kokkos::CudaSpace;
using ExecSpace       = Kokkos::Cuda;
using Layout          = Kokkos::LayoutLeft;
#endif

// Won't have both
#if HAVE_OPENMP
using MemSpace        = Kokkos::HostSpace;
using ExecSpace       = Kokkos::OpenMP;
using Layout          = Kokkos::LayoutRight;
#endif

using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
using mdrange_policy2 = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
using mdrange_policy3 = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;

using RMatrix1D    = Kokkos::View<real_t *,Layout,ExecSpace>;
using RMatrix2D    = Kokkos::View<real_t **,Layout,ExecSpace>;
using RMatrix3D    = Kokkos::View<real_t ***,Layout,ExecSpace>;
using RMatrix4D    = Kokkos::View<real_t ****,Layout,ExecSpace>;
using RMatrix5D    = Kokkos::View<real_t *****,Layout,ExecSpace>;
using IMatrix1D    = Kokkos::View<int *,Layout,ExecSpace>;
using IMatrix2D    = Kokkos::View<int **,Layout,ExecSpace>;
using IMatrix3D    = Kokkos::View<int ***,Layout,ExecSpace>;
using IMatrix4D    = Kokkos::View<int ****,Layout,ExecSpace>;
using IMatrix5D    = Kokkos::View<int *****,Layout,ExecSpace>;
using SVar         = Kokkos::View<size_t,Layout,ExecSpace>;
using SArray1D     = Kokkos::View<size_t *,Layout,ExecSpace>;
using SArray2D     = Kokkos::View<size_t **,Layout,ExecSpace>;
using SArray3D     = Kokkos::View<size_t ***,Layout,ExecSpace>;
using SArray4D     = Kokkos::View<size_t ****,Layout,ExecSpace>;
using SArray5D     = Kokkos::View<size_t *****,Layout,ExecSpace>;

using SHArray1D     = Kokkos::View<size_t *,Layout,Kokkos::HostSpace>;

#else

#define KOKKOS_FUNCTION 
#define KOKKOS_INLINE_FUNCTION inline

#endif

//To disable asserts, uncomment the following line
//#define NDEBUG

//---Begin Standard Data Structures---

template <typename T>
class ViewCArrayMat {

private:
    size_t dim1_;
    size_t dim2_;
    size_t dim3_;
    size_t dim4_;
    size_t dim5_;
    size_t dim6_;
    size_t length_;
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
        T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const;

    KOKKOS_INLINE_FUNCTION
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
}

// Overloaded 1D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayMat<T>::ViewCArrayMat(T* some_matrix, size_t dim1) {
    dim6_ = 1;
    dim5_ = 1;
    dim4_ = 1;
    dim3_ = 1;
    dim2_ = 1;
    dim1_ = 1;
    length_ = dim1;
    this_matrix_ = some_matrix;
}

// Overloaded 2D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayMat<T>::ViewCArrayMat(T* some_matrix, size_t dim1,
    size_t dim2) {
    dim6_ = 1;
    dim5_ = 1;
    dim4_ = 1;
    dim3_ = 1;
    dim2_ = 1;
    dim1_ = dim2;
    length_ = (dim1_ * dim1);
    this_matrix_ = some_matrix;
}

// Overloaded 3D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCArrayMat<T>::ViewCArrayMat(T* some_matrix, size_t dim1, size_t dim2,
    size_t dim3) {

    dim6_ = 1;
    dim5_ = 1;
    dim4_ = 1;
    dim3_ = 1;
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

    dim6_ = 1;
    dim5_ = 1;
    dim4_ = 1;
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
    dim6_ = 1;
    dim5_ = 1;
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
    dim6_ = 1;
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
KOKKOS_INLINE_FUNCTION
T& ViewCArrayMat<T>::operator()(size_t i) const {
    //assert(i < length_ && "i is out of bounds in ViewCArrayMat 1D!");
    return this_matrix_[(i)];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCArrayMat<T>::operator()(size_t i, size_t j) const {
    //assert(i < length_/dim1_ && "i is out of bounds in ViewCArrayMat 2D!");
    //assert(j < dim1_/dim2_ && "j is out of bounds in ViewCArrayMat 2D!");
    return this_matrix_[(i)*dim1_ + (j)];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCArrayMat<T>::operator()(size_t i, size_t j, size_t k) const {
    //assert(i < length_ / dim1_ && "i is out of bounds in ViewCArrayMat 3D!");
    //assert(j < dim1_ / dim2_ && "j is out of bounds in ViewCArrayMat 3D!");
    //assert(k < dim2_/dim3_ && "k is out of bounds in ViewCArrayMat 3D!");
    return this_matrix_[(i)*dim1_ + (j)*dim2_ + k];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCArrayMat<T>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    //assert(i < length_/dim1_ && "i is out of bounds in ViewCArrayMat 4D!");
    //assert(j < dim1_ / dim2_ && "j is out of bounds in ViewCArrayMat 4D!");
    //assert(k < dim2_/dim3_ && "k is out of bounds in ViewCArrayMat 4D!");
    //assert(l < dim3_/dim4_ && "l is out of bounds in ViewCArrayMat 4D!");
    return this_matrix_[(i)*dim1_ + (j)*dim2_ + (k)*dim3_ + (l)];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCArrayMat<T>::operator()(size_t i, size_t j, size_t k, size_t l,
    size_t m) const {
    //assert(i >= 1 && i <= dim1_ && "i is out of bounds for ViewCArrayMat 5D!");
    //assert(j >= 1 && j <= dim2_ && "j is out of bounds for ViewCArrayMat 5D!");
    //assert(k >= 1 && k <= dim3_ && "k is out of bounds for ViewCArrayMat 5D!");
    //assert(l >= 1 && l <= dim4_ && "l is out of bounds for ViewCArrayMat 5D!");
    //assert(m >= 1 && m <= dim5_ && "m is out of bounds for ViewCArrayMat 5D!");
    return this_matrix_[(i)*dim1_ + (j)*dim2_ + (k)*dim3_ + (l)*dim4_ + (m)];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCArrayMat<T>::operator()(size_t i, size_t j, size_t k, size_t l,
    size_t m, size_t n) const {
    //assert(i >= 1 && i <= dim1_ && "i is out of bounds for ViewCArrayMat 6D!");
    //assert(j >= 1 && j <= dim2_ && "j is out of bounds for ViewCArrayMat 6D!");
    //assert(k >= 1 && k <= dim3_ && "k is out of bounds for ViewCArrayMat 6D!");
    //assert(l >= 1 && l <= dim4_ && "l is out of bounds for ViewCArrayMat 6D!");
    //assert(m >= 1 && m <= dim5_ && "m is out of bounds for ViewCArrayMat 6D!");
    //assert(n >= 1 && n <= dim6_ && "n is out of bounds for ViewCArrayMat 6D!");
    return this_matrix_[(i)*dim1_ + (j)*dim2_ + (k)*dim3_ + (l)*dim4_ + (m)*dim5_ + (n)];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCArrayMat<T>::operator()(size_t i, size_t j, size_t k, size_t l,
    size_t m, size_t n, size_t o) const {
    //assert(i >= 1 && i <= dim1_ && "i is out of bounds for ViewCArrayMat 7D!");
    //assert(j >= 1 && j <= dim2_ && "j is out of bounds for ViewCArrayMat 7D!");
    //assert(k >= 1 && k <= dim3_ && "k is out of bounds for ViewCArrayMat 7D!");
    //assert(l >= 1 && l <= dim4_ && "l is out of bounds for ViewCArrayMat 7D!");
    //assert(m >= 1 && m <= dim5_ && "m is out of bounds for ViewCArrayMat 7D!");
    //assert(n >= 1 && n <= dim6_ && "n is out of bounds for ViewCArrayMat 7D!");
    //assert(o >= 1 && o <= dim7_ && "o is out of bounds for ViewCArrayMat 7D!");
    return this_matrix_[(i)*dim1_ + (j)*dim2_ + (k)*dim3_ + (l)*dim4_ + (m)*dim5_ + (n)*dim6_ + (o)];
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
    size_t dim1_;
    size_t dim2_;
    size_t dim3_;
    size_t dim4_;
    size_t dim5_;
    size_t dim6_;
    size_t length_;
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
        T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const;

    KOKKOS_INLINE_FUNCTION
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
        operator ViewCArrayMat<T>();

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
    dim6_ = 1;
    dim5_ = 1;
    dim4_ = 1;
    dim3_ = 1;
    dim2_ = 1;
    dim1_ = 1;
    length_ = dim1;
    this_matrix_ = some_matrix;
}

// Overloaded 2D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixMat<T>::ViewCMatrixMat(T* some_matrix, size_t dim1,
    size_t dim2) {
    dim6_ = 1;
    dim5_ = 1;
    dim4_ = 1;
    dim3_ = 1;
    dim2_ = 1;
    dim1_ = dim2;
    length_ = (dim1_ * dim1);
    this_matrix_ = some_matrix;
}

// Overloaded 3D constructor
template <typename T>
KOKKOS_FUNCTION
ViewCMatrixMat<T>::ViewCMatrixMat(T* some_matrix, size_t dim1, size_t dim2,
    size_t dim3) {

    dim6_ = 1;
    dim5_ = 1;
    dim4_ = 1;
    dim3_ = 1;
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

    dim6_ = 1;
    dim5_ = 1;
    dim4_ = 1;
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
    dim6_ = 1;
    dim5_ = 1;
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
    dim6_ = 1;
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
KOKKOS_INLINE_FUNCTION
T& ViewCMatrixMat<T>::operator()(size_t i) const {
    //assert(i >= 1 && i <= length_ && "i is out of bounds in ViewCMatrixMat 1D!");
    return this_matrix_[(i - 1)];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCMatrixMat<T>::operator()(size_t i, size_t j) const {
    //assert(i >= 1 && i <= length_/dim1_ && "i is out of bounds in ViewCMatrixMat 2D!");
    //assert(j >= 1 && j <= dim1_/dim2_ && "j is out of bounds in ViewCMatrixMat 2D!");
    return this_matrix_[(i - 1) * dim1_ + (j - 1)];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCMatrixMat<T>::operator()(size_t i, size_t j, size_t k) const {
    //assert(i >= 1 && i <= length_/dim1_ && "i is out of bounds in ViewCMatrixMat 3D!");
    //assert(j >= 1 && j <= dim1_/dim2_ && "j is out of bounds in ViewCMatrixMat 3D!");
    //assert(k >= 1 && k <= dim2_/dim3_ && "k is out of bounds in ViewCMatrixMat 3D!");
    return this_matrix_[(i - 1) * dim1_ + (j - 1) * dim2_ + k - 1];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
T& ViewCMatrixMat<T>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    //assert(i >= 1 && i <= length_/dim1_ && "i is out of bounds in ViewCMatrixMat 4D!");
    //assert(j >= 1 && j <= dim1_/dim2_ && "j is out of bounds in ViewCMatrixMat 4D!");
    //assert(k >= 1 && k <= dim2_/dim3_ && "k is out of bounds in ViewCMatrixMat 4D!");
    //assert(l >= 1 && l <= dim3_/dim4_ && "l is out of bounds in ViewCMatrixMat 4D!");
    return this_matrix_[(i - 1) * dim1_ + (j - 1) * dim2_ + (k - 1) * dim3_ + (l - 1)];
}

template <typename T>
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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

    memcpy(&this_matrix_[0], &temp.this_matrix_[0], sizeof(T) * length_);

    //for (size_t i = 0; i < length_; i++) { this_matrix_[i] = temp.this_matrix_[i]; }

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
ViewCMatrixMat<T>::operator ViewCArrayMat<T>() {

    ViewCArrayMat<T>  out(length_, dim1_, dim2_, dim3_, dim4_, dim5_, dim6_, &this_matrix_[0]);

    return out;
}

template <typename T>
KOKKOS_FUNCTION
ViewCMatrixMat<T>::~ViewCMatrixMat() {}

////////////////////////////////////////////////////////////////////////////////
// End of ViewCMatrixMat
////////////////////////////////////////////////////////////////////////////////



//1. FArray
template <typename T>
class FArray {
    
private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_, length_;
    T * this_array;
    
public:
    
    // default constructor
   FArray ();
   
    //overload constructors from 1D to 6D
     
   FArray(size_t some_dim1);
   FArray(size_t some_dim1, size_t some_dim2);
   FArray(size_t some_dim1, size_t some_dim2, size_t some_dim3);
   FArray(size_t some_dim1, 
          size_t some_dim2,
          size_t some_dim3,
          size_t some_dim4);
    
   FArray(size_t some_dim1,
          size_t some_dim2,
          size_t some_dim3,
          size_t some_dim4,
          size_t some_dim5);

   FArray(size_t some_dim1,
          size_t some_dim2,
          size_t some_dim3,
          size_t some_dim4,
          size_t some_dim5,
          size_t some_dim6);

    // overload operator() to access data as array(i,....,n);
    T& operator()(size_t i);
    T& operator()(size_t i, size_t j);
    T& operator()(size_t i, size_t j, size_t k);
    T& operator()(size_t i, size_t j, size_t k, size_t l);
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m);
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n);

    //overload = operator
    FArray& operator=(const FArray& temp);

    // deconstructor
    ~FArray ( );
    
}; // end of f_array_t

//---FArray class definnitions----

//constructors
template <typename T>
FArray<T>::FArray(){}

//1D
template <typename T>
FArray<T>::FArray(size_t some_dim1) {
    dim1_ = some_dim1;
    length_ = dim1_;
    this_array = new T[length_];
}

template <typename T>
FArray<T>::FArray(size_t some_dim1, size_t some_dim2) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = dim1_*dim2_;
    this_array = new T[length_];
}

//3D
template <typename T>
FArray<T>::FArray(size_t some_dim1, size_t some_dim2, size_t some_dim3) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = dim1_*dim2_*dim3_;
    this_array = new T[length_];
}

//4D
template <typename T>
FArray<T>::FArray(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = dim1_*dim2_*dim3_*dim4_;
    this_array = new T[length_];
}

//5D
template <typename T>
FArray<T>::FArray(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = dim1_*dim2_*dim3_*dim4_*dim5_;
    this_array = new T[length_];
}

//6D
template <typename T>
FArray<T>::FArray(size_t some_dim1,size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5, size_t some_dim6) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = dim1_*dim2_*dim3_*dim4_*dim5_*dim6_;
    this_array = new T[length_];
}

//overload operator () for 1D to 6D
//indices are from [0:N-1]

//1D
template <typename T>
T& FArray<T>::operator()(size_t i)
{
    assert( i < dim1_ && "i is out of bounds in FArray 1D!");
    return this_array[i];
}

//2D
template <typename T>
T& FArray<T>::operator()(size_t i, size_t j)
{
    assert( i < dim1_ && "i is out of bounds in FArray 2D!");
    assert( j < dim2_ && "j is out of bounds in FArray 2D!");
    return this_array[i + j*dim1_];
}

//3D
template <typename T>
T& FArray<T>::operator()(size_t i, size_t j, size_t k)
{
    assert( i < dim1_ && "i is out of bounds in FArray 3D!");
    assert( j < dim2_ && "j is out of bounds in Farray 3D!");
    assert( k < dim3_ && "k is out of bounds in FArray 3D!");
    return this_array[i + j*dim1_ + k*dim1_*dim2_];
}

//4D
template <typename T>
T& FArray<T>::operator()(size_t i, size_t j, size_t k, size_t l)
{
    assert( i < dim1_ && "i is out of bounds in FArray 4D!");
    assert( j < dim2_ && "j is out of bounds in FArray 4D!");
    assert( k < dim3_ && "k is out of bounds in FArray 4D!");
    assert( l < dim4_ && "l is out of bounds in FArray 4D!");
    return this_array[ i + j*dim1_ + k*dim1_*dim2_ + l*dim1_*dim2_*dim3_];
}

//5D
template <typename T>
T& FArray<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m)
{
    assert( i < dim1_ && "i is out of bounds in FArray 5D!");
    assert( j < dim2_ && "j is out of bounds in FArray 5D!");
    assert( k < dim3_ && "k is out of bounds in FArray 5D!");
    assert( l < dim4_ && "l is out of bounds in FArray 5D!");
    assert( m < dim5_ && "m is out of bounds in FArray 5D!");
    return this_array[i + j*dim1_ + k*dim1_*dim2_ + l*dim1_*dim2_*dim3_ + m*dim1_*dim2_*dim3_*dim4_];
}

//6D
template <typename T>
T& FArray<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n)
{

    assert( i < dim1_ && "i is out of bounds in FArray 5D!");
    assert( j < dim2_ && "j is out of bounds in FArray 5D!");
    assert( k < dim3_ && "k is out of bounds in FArray 5D!");
    assert( l < dim4_ && "l is out of bounds in FArray 5D!");
    assert( m < dim5_ && "m is out of bounds in FArray 5D!");
    assert( n < dim6_ && "n is out of bounds in FArray 6D!");
    return this_array[i + j*dim1_ + k*dim1_*dim2_ + l*dim1_*dim2_*dim3_ + m*dim1_*dim2_*dim3_*dim4_ + n*dim1_*dim2_*dim3_*dim4_*dim5_];
}

// = operator
//THIS = FArray <> TEMP(n,m,...)
template <typename T>
FArray<T>& FArray<T>::operator= (const FArray& temp)
{
	if(this != & temp) {
	  dim1_ = temp.dim1_;
	  dim2_ = temp.dim2_;
	  dim3_ = temp.dim3_;
	  dim4_ = temp.dim4_;
	  dim5_ = temp.dim5_;
	  dim6_ = temp.dim6_;
	  length_ = temp.length_;
	  this_array = new T[length_];
	}
  return *this;
} 

//delete FArray
template <typename T>
FArray<T>::~FArray(){
    delete [] this_array;
}

//---end of FArray class definitions----

//2. ViewFArray
template <typename T>
class ViewFArray {

private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_;
    T * this_array;
    
public:
    
    // default constructor
    ViewFArray ();

    //---1D array---
    ViewFArray(T *some_array, size_t some_dim1);
    T& operator()(size_t i);
    
    //--- 2D array ---
    
    // overloaded constructor
    ViewFArray (T *some_array, size_t some_dim1, size_t some_dim2);
    T& operator()(size_t i, size_t j);
    
    //--- 3D array ---
    ViewFArray (T *some_array, size_t some_dim1, size_t some_dim2, size_t some_dim3);
    T& operator()(size_t i, size_t j, size_t k);
    
    //--- 4D array ---
    // overloaded constructor
    ViewFArray (T *some_array, size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4);
    T& operator()(size_t i, size_t j, size_t k, size_t l);

    //--- 5D array ---
    // overloaded constructor
    ViewFArray (T *some_array, size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5);
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m);
   
    //------6D -----

    ViewFArray (T *some_array,size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5, size_t some_dim6);
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n);

    
}; // end of viewFArray

//class definitions for viewFArray

//~~~~constructors for viewFArray for 1D to 6D~~~~~~~

//no dimension
template <typename T>
ViewFArray<T>::ViewFArray(){}

//1D
template <typename T>
ViewFArray<T>::ViewFArray(T *some_array, size_t some_dim1){
	dim1_ = some_dim1;
	this_array = some_array;
}

//2D
template <typename T>
ViewFArray<T>::ViewFArray(T *some_array, size_t some_dim1, size_t some_dim2) {
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	this_array = some_array;
}

//3D
template <typename T>
ViewFArray<T>::ViewFArray(T *some_array, size_t some_dim1, size_t some_dim2, size_t some_dim3){
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	dim3_ = some_dim3;
	this_array = some_array;
}

//4D
template <typename T>
ViewFArray<T>::ViewFArray(T *some_array, size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4){
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	dim3_ = some_dim3;
	dim4_ = some_dim4;
	this_array = some_array;
}

//5D
template <typename T>
ViewFArray<T>::ViewFArray(T *some_array, size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5){
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	dim3_ = some_dim3;
	dim4_ = some_dim4;
	dim5_ = some_dim5;
	this_array = some_array;
}

//6D
template <typename T>
ViewFArray<T>::ViewFArray(T *some_array, size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5, size_t some_dim6){
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	dim3_ = some_dim3;
	dim4_ = some_dim4;
	dim5_ = some_dim5;
	dim6_ = some_dim6;
	this_array = some_array;
}

//~~~~~~operator () overload 
//for dimensions 1D to 6D
//indices for array are from 0...N-1

//1D
template <typename T>
T& ViewFArray<T>::operator()(size_t i)
{
	assert(i < dim1_ && "i is out of bounds in ViewFArray 1D!");
	return this_array[i];
}

//2D
template <typename T>
T& ViewFArray<T>::operator()(size_t i, size_t j) 
{
	assert( i < dim1_ && "i is out of bounds in ViewFArray 2D!");
	assert( j < dim2_ && "j is out of bounds in ViewFArray 2D!");
	return this_array[i + j*dim1_];
}

//3D
template <typename T>
T& ViewFArray<T>::operator()(size_t i, size_t j,size_t k)
{
	assert(i < dim1_ && "i is out of bounds in ViewFArray 3D!");
	assert(j < dim2_ && "j is out of bounds in ViewFArray 3D!");
	assert(k < dim3_ && "k is out of bounds in ViewFArray 3D!");
	return this_array[i + j*dim1_ + k*dim1_*dim2_];
}

//4D
template <typename T>
T& ViewFArray<T>::operator()(size_t i, size_t j, size_t k, size_t l)
{
	assert(i < dim1_ && "i is out of bounds in ViewFArray 4D!");
	assert(j < dim2_ && "j is out of bounds in ViewFArray 4D!");
	assert(k < dim3_ && "k is out of bounds in ViewFArray 4D!");
	assert(l < dim4_ && "l is out of bounds in ViewFArray 4D!");
	return this_array[i + j*dim1_ + k*dim1_*dim2_ + l*dim1_*dim2_*dim3_];
}

//5D
template <typename T>
T& ViewFArray<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m)
{
	assert(i < dim1_ && "i is out of bounds in ViewFArray 5D!");
	assert(j < dim2_ && "j is out of bounds in ViewFArray 5D!");
	assert(k < dim3_ && "k is out of bounds in ViewFArray 5D!");
	assert(l < dim4_ && "l is out of bounds in ViewFArray 5D!");
	assert(m < dim5_ && "m is out of bounds in ViewFArray 5D!");
	return this_array[i + j*dim1_ + k*dim1_*dim2_ + l*dim1_*dim2_*dim3_ + m*dim1_*dim2_*dim3_*dim4_];
}

//6D
template <typename T>
T& ViewFArray<T>:: operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n)
{
	assert(i < dim1_ && "i is out of bounds in ViewFArray 6D!");
	assert(j < dim2_ && "j is out of bounds in ViewFArray 6D!");
	assert(k < dim3_ && "k is out of bounds in ViewFArray 6D!");
	assert(l < dim4_ && "l is out of bounds in ViewFArray 6D!");
	assert(m < dim5_ && "m is out of bounds in ViewFArray 6D!");
	assert(n < dim6_ && "n is out of bounds in ViewFArray 6D!");
	return this_array[i + j*dim1_ + k*dim1_*dim2_ + l*dim1_*dim2_*dim3_ + m*dim1_*dim2_*dim3_*dim4_ + n*dim1_*dim2_*dim3_*dim4_*dim5_];
}

//---end of ViewFArray class definitions---

//3. FMatrix
template <typename T>
class FMatrix {
private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_;
    size_t length_; // Length of 1D array
    T* this_matrix_;

public:
    // Default constructor
    FMatrix ();

    // --- 1D matrix ---
    
    // Overloaded constructor
    FMatrix (size_t some_dim1);

    // Overload operator() to access data as matrix(i), where i = [1:N]
    T& operator() (size_t i) const;

    // --- 2D matrix ---
    
    // Overloaded constructor
    FMatrix (size_t some_dim1, size_t some_dim2);

    // Overload operator() to access data as matrix(i, j),
    // where i = [1:N], j = [1:N]
    T& operator() (size_t i, size_t j) const;

    // --- 3D matrix ---
    
    // Overloaded constructor
    FMatrix (size_t some_dim1, size_t some_dim2, size_t some_dim3);

    // Overload operator() to access data as matrix(i, j, k),
    // where i = [1:N], j = [1:N], k = [1:N]
    T& operator() (size_t i, size_t j, size_t k) const;

    // --- 4D matrix ---
    
    // Overloaded constructor
    FMatrix (size_t some_dim1, size_t some_dim2, size_t some_dim3, 
               size_t some_dim4);

    // Overload operator() to access data as matrix(i, j, k, l),
    // where i = [1:N], j = [1:N], k = [1:N], l = [1:N]
    T& operator() (size_t i, size_t j, size_t k, size_t l) const;

    // --- 5D matrix ---
    
    // Overloaded constructor
    FMatrix (size_t some_dim1, size_t some_dim2, size_t some_dim3,
               size_t some_dim4, size_t some_dim5);

    // Overload operator() to access data as matrix(i, j, k, l, m),
    // where i = [1:N], j = [1:N], k = [1:N], l = [1:N],
    // m = [1:N]
    T& operator() (size_t i, size_t j, size_t k, size_t l,
                          size_t m) const;

    // --- 6D matrix ---
    
    // Overloaded constructor
    FMatrix (size_t some_dim1, size_t some_dim2, size_t some_dim3,
               size_t some_dim4, size_t some_dim5, size_t some_dim6);

    // Overload operator() to access data as matrix(i, j, k, l, m, n),
    // where i = [1:N], j = [1:N], k = [1:N], l = [1:N], 
    // m = [1:N], n = [1:N]
    T& operator() (size_t i, size_t j, size_t k, size_t l,
                          size_t m, size_t n) const;

    // Overload copy assignment operator
    FMatrix& operator=(const FMatrix& temp);

    size_t size() const;

    // Deconstructor
    ~FMatrix ();

}; // End of FMatrix

//---FMatrix class definitions---

//constructors

//1D
template <typename T>
FMatrix<T>::FMatrix(size_t some_dim1) {
    // assert(some_dim1 > 0);
    dim1_ = some_dim1;
    length_ = dim1_;
    this_matrix_ = new T[length_];
}

//2D
template <typename T>
FMatrix<T>::FMatrix(size_t some_dim1, size_t some_dim2) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = dim1_ * dim2_;
    this_matrix_ = new T[length_];
}

//3D
template <typename T>
FMatrix<T>::FMatrix(size_t some_dim1, size_t some_dim2, size_t some_dim3) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = dim1_ * dim2_ * dim3_;
    this_matrix_ = new T[length_];
}

//4D
template <typename T>
FMatrix<T>::FMatrix(size_t some_dim1, size_t some_dim2, size_t some_dim3, 
                        size_t some_dim4) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    // assert(some_dim4 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = dim1_ * dim2_ * dim3_ * dim4_;
    this_matrix_ = new T[length_];
}

//5D
template <typename T>
FMatrix<T>::FMatrix(size_t some_dim1, size_t some_dim2, size_t some_dim3,
                        size_t some_dim4, size_t some_dim5) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    // assert(some_dim4 > 0):
    // assert(some_dim5 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_;
    this_matrix_ = new T[length_];
}

//6D
template <typename T>
FMatrix<T>::FMatrix(size_t some_dim1, size_t some_dim2, size_t some_dim3,
                        size_t some_dim4, size_t some_dim5, size_t some_dim6) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    // assert(some_dim4 > 0):
    // assert(some_dim5 > 0);
    // assert(some_dim6 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_;
    this_matrix_ = new T[length_];

}

//overload operators

//1D
template <typename T>
inline T& FMatrix<T>::operator() (size_t i) const {
    assert(i >= 1 && i <= dim1_);
    return this_matrix_[i - 1];
}

//2D
template <typename T>
inline T& FMatrix<T>::operator() (size_t i, size_t j) const {
    assert(i >= 1 && i <= dim1_);
    assert(j >= 1 && j <= dim2_);
    return this_matrix_[(i - 1) + ((j - 1) * dim1_)];
}

//3D
template <typename T>
inline T& FMatrix<T>::operator() (size_t i, size_t j, size_t k) const {
    assert(i >= 1 && i <= dim1_);
    assert(j >= 1 && j <= dim2_);
    assert(k >= 1 && k <= dim3_);
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_)];
}

//4D
template <typename T>
inline T& FMatrix<T>::operator() (size_t i, size_t j, size_t k, size_t l) const {
    assert(i >= 1 && i <= dim1_);
    assert(j >= 1 && j <= dim2_);
    assert(k >= 1 && k <= dim3_);
    assert(l >= 1 && l <= dim4_);
    return this_matrix_[(i - 1) + ((j - 1) * dim1_)  
                                + ((k - 1) * dim1_ * dim2_)  
                                + ((l - 1) * dim1_ * dim2_ * dim3_)];
}

//5D
template <typename T>
inline T& FMatrix<T>::operator() (size_t i, size_t j, size_t k, size_t l,
                                    size_t m) const {
    assert(i >= 1 && i <= dim1_);
    assert(j >= 1 && j <= dim2_);
    assert(k >= 1 && k <= dim3_);
    assert(l >= 1 && l <= dim4_);
    assert(m >= 1 && m <= dim5_);
    return this_matrix_[(i - 1) + ((j - 1) * dim1_)  
                                + ((k - 1) * dim1_ * dim2_)  
                                + ((l - 1) * dim1_ * dim2_ * dim3_) 
                                + ((m - 1) * dim1_ * dim2_ * dim3_ * dim4_)];
}

//6D
template <typename T>
inline T& FMatrix<T>::operator() (size_t i, size_t j, size_t k, size_t l,
                                    size_t m, size_t n) const {
    assert(i >= 1 && i <= dim1_);
    assert(j >= 1 && j <= dim2_);
    assert(k >= 1 && k <= dim3_);
    assert(l >= 1 && l <= dim4_);
    assert(m >= 1 && m <= dim5_);
    assert(n >= 1 && n <= dim6_);
    return this_matrix_[(i - 1) + ((j - 1) * dim1_)  
                                + ((k - 1) * dim1_ * dim2_)  
                                + ((l - 1) * dim1_ * dim2_ * dim3_)  
                                + ((m - 1) * dim1_ * dim2_ * dim3_ * dim4_)  
                                + ((n - 1) * dim1_ * dim2_ * dim3_ * dim4_ * dim5_)];
}

template <typename T>
inline FMatrix<T>& FMatrix<T>::operator= (const FMatrix& temp)
{
    // Do nothing if assignment is of the form x = x
    if (this != &temp) {
        dim1_ = temp.dim1_;
        dim2_ = temp.dim2_;
        dim3_ = temp.dim3_;
        dim4_ = temp.dim4_;
        dim5_ = temp.dim5_;
        dim6_ = temp.dim6_;
        length_ = temp.length_;
        this_matrix_ = new T[length_];
    }
    
    return *this;
}

template <typename T>
inline size_t FMatrix<T>::size() const {
    return length_;
}

template <typename T>
FMatrix<T>::~FMatrix() {
    delete[] this_matrix_;
}

//----end of FMatrix class definitions----

//4. ViewFMatrix
template <typename T>
class ViewFMatrix {

private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_;
    size_t length_; // Length of 1D array
    T * this_matrix_;
    
public:
    
    // Default constructor
    ViewFMatrix ();
    
    //--- 1D matrix ---
    
    // Overloaded constructor
    ViewFMatrix(T *some_matrix,
                size_t some_dim1);
    
    // Overload  operator() to access data as matrix(i,j),
    // where i = [1:N], j = [1:N]
    T& operator()(size_t i) const;
    
    //--- 2D matrix ---
    
    // Overloaded constructor
    ViewFMatrix(T *some_matrix,
                size_t some_dim1,
                size_t some_dim2);
    
    // Overload operator() to access data as matrix(i,j),
    //  where i=[1:N], j=[1:N]
    T& operator()(size_t i, size_t j) const;
    
    //--- 3D matrix ---
    
    // Overloaded constructor
    ViewFMatrix(T *some_matrix,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3);
    
    // Overload operator() to access data as matrix(i,j,k),
    // where i = [1:N], j = [1:N], k = [1:N]
    T& operator()(size_t i, size_t j, size_t k) const;
    
    //--- 4D matrix ---
    
    // Overloaded constructor
    ViewFMatrix(T *some_matrix,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4);
        
    // Overload operator() to access data as matrix(i, j, k, l),
    // where i = [0:n-1], j = [1:N], k = [1:N], l = [1:N]
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;
    
    //--- 5D matrix ---
    
    // Overloaded constructor
    ViewFMatrix (T *some_matrix,
                 size_t some_dim1,
                 size_t some_dim2,
                 size_t some_dim3,
                 size_t some_dim4,
                 size_t some_dim5);
        
    // Overload operator() to access data as matrix(i,j,k,l,m),
    // where i = [1:N], j = [1:N], k = [1:N], l = [1:N], m = [1:N]
    T& operator() (size_t i, 
                   size_t j, 
                   size_t k, 
                   size_t l, 
                   size_t m) const;
    
    //--- 6D matrix ---
    
    // Overloaded constructor
    ViewFMatrix (T *some_matrix,
                 size_t some_dim1,
                 size_t some_dim2,
                 size_t some_dim3,
                 size_t some_dim4,
                 size_t some_dim5,
                 size_t some_dim6);
        
    // Overload operator() to access data as matrix(i,j,k,l,m,n),
    // where i = [1:N], j = [1:N], k = [1:N], l = [1:N], m = [1:N]
    T& operator()(size_t i, 
                  size_t j, 
                  size_t k, 
                  size_t l, 
                  size_t m, 
                  size_t n) const;

    size_t size() const;
}; // end of ViewFMatrix

//constructors

//no dimension
template <typename T>
ViewFMatrix<T>::ViewFMatrix() {}

//1D
template <typename T>
ViewFMatrix<T>::ViewFMatrix(T *some_matrix,
                            size_t some_dim1)
{
    dim1_ = some_dim1;
    length_ = dim1_;
    this_matrix_ = some_matrix;
}

//2D
template <typename T>
ViewFMatrix<T>::ViewFMatrix(T *some_matrix,
                            size_t some_dim1,
                            size_t some_dim2)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = dim1_ * dim2_;
    this_matrix_ = some_matrix;
}

//3D
template <typename T>
ViewFMatrix<T>::ViewFMatrix (T *some_matrix,
                             size_t some_dim1,
                             size_t some_dim2,
                             size_t some_dim3)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = dim1_ * dim2_ * dim3_;
    this_matrix_ = some_matrix;
}

//4D
template <typename T>
ViewFMatrix<T>::ViewFMatrix(T *some_matrix,
                            size_t some_dim1,
                            size_t some_dim2,
                            size_t some_dim3,
                            size_t some_dim4)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = dim1_ * dim2_ * dim3_ * dim4_;
    this_matrix_ = some_matrix;
}

//5D
template <typename T>
ViewFMatrix<T>::ViewFMatrix(T *some_matrix,
                            size_t some_dim1,
                            size_t some_dim2,
                            size_t some_dim3,
                            size_t some_dim4,
                            size_t some_dim5)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_;
    this_matrix_ = some_matrix;
}

//6D
template <typename T>
ViewFMatrix<T>::ViewFMatrix(T *some_matrix,
                            size_t some_dim1,
                            size_t some_dim2,
                            size_t some_dim3,
                            size_t some_dim4,
                            size_t some_dim5,
                            size_t some_dim6)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_;
    this_matrix_ = some_matrix;
}

//overload operator ()

//1D
template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i) const
{
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in c_matrix 1D");  // die if >= dim1
        
    return this_matrix_[(i - 1)];
}

//2D
template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i, 
                                     size_t j) const
{
       
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in c_matrix 2D");  // die if >= dim1
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in c_matrix 2D");  // die if >= dim2
        
    return this_matrix_[(i - 1) + ((j - 1) * dim1_)];
}

//3D
template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i, 
                                     size_t j, 
                                     size_t k) const
{
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in c_matrix 3D");  // die if >= dim1
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in c_matrix 3D");  // die if >= dim2
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in c_matrix 3D");  // die if >= dim3
        
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_)];
}

//4D
template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i, 
                                     size_t j, 
                                     size_t k, 
                                     size_t l) const
{
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in c_matrix 4D");  // die if >= dim1
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in c_matrix 4D");  // die if >= dim2
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in c_matrix 4D");  // die if >= dim3
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in c_matrix 4D");  // die if >= dim4
        
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_) 
                                + ((l - 1) * dim1_ * dim2_ * dim3_)];
}

//5D
template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i, 
                                     size_t j, 
                                     size_t k, 
                                     size_t l, 
                                     size_t m) const
{
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in c_matrix 5D");  // die if >= dim1
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in c_matrix 5D");  // die if >= dim2
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in c_matrix 5D");  // die if >= dim3
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in c_matrix 5D");  // die if >= dim4
    assert(m >= 1 && m <= dim5_ && "m is out of bounds in c_matrix 5D");  // die if >= dim5
       
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_) 
                                + ((l - 1) * dim1_ * dim2_ * dim3_)
                                + ((m - 1) * dim1_ * dim2_ * dim3_ * dim4_)];
}

//6D
template <typename T>
inline T& ViewFMatrix<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
{
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in c_matrix 6D");  // die if >= dim1
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in c_matrix 6D");  // die if >= dim2
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in c_matrix 6D");  // die if >= dim3
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in c_matrix 6D");  // die if >= dim4
    assert(m >= 1 && m <= dim5_ && "m is out of bounds in c_matrix 6D");  // die if >= dim5
    assert(n >= 1 && n <= dim6_ && "n is out of bounds in c_matrix 6D");  // die if >= dim6
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_) 
                                + ((l - 1) * dim1_ * dim2_ * dim3_)
                                + ((m - 1) * dim1_ * dim2_ * dim3_ * dim4_)
                                + ((n - 1) * dim1_ * dim2_ * dim3_ * dim4_ * dim5_)];
}

//-----end ViewFMatrix-----

//5. CArray
template <typename T>
class CArray {
private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_;
    size_t length_; // Length of 1D array
    T* this_array_;

public:
    // Default constructor
    CArray ();

    // --- 1D array ---
    
    // Overloaded constructor
    CArray (size_t some_dim1);

    // Overload operator() to access data as array(i), where i = [0:N-1]
    T& operator() (size_t i) const;

    // --- 2D array ---
    
    // Overloaded constructor
    CArray (size_t some_dim1, size_t some_dim2);

    // Overload operator() to access data as array(i, j),
    // where i = [0:N-1], j = [0:N-1]
    T& operator() (size_t i, size_t j) const;

    // --- 3D array ---
    
    // Overloaded constructor
    CArray (size_t some_dim1, size_t some_dim2, size_t some_dim3);

    // Overload operator() to access data as array(i, j, k),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1]
    T& operator() (size_t i, size_t j, size_t k) const;

    // --- 4D array ---
    
    // Overloaded constructor
    CArray (size_t some_dim1, size_t some_dim2, size_t some_dim3, 
               size_t some_dim4);

    // Overload operator() to access data as array(i, j, k, l),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1]
    T& operator() (size_t i, size_t j, size_t k, size_t l) const;

    // --- 5D array ---
    
    // Overloaded constructor
    CArray (size_t some_dim1, size_t some_dim2, size_t some_dim3,
               size_t some_dim4, size_t some_dim5);

    // Overload operator() to access data as array(i, j, k, l, m),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1],
    // m = [0:N-1]
    T& operator() (size_t i, size_t j, size_t k, size_t l,
                          size_t m) const;

    // --- 6D array ---
    
    // Overloaded constructor
    CArray (size_t some_dim1, size_t some_dim2, size_t some_dim3,
               size_t some_dim4, size_t some_dim5, size_t some_dim6);

    // Overload operator() to access data as array(i, j, k, l, m, n),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1], 
    // m = [0:N-1], n = [0:N-1]
    T& operator() (size_t i, size_t j, size_t k, size_t l,
                          size_t m, size_t n) const;

    // Overload copy assignment operator
    CArray& operator= (const CArray& temp); 

    size_t size() const;

    // Deconstructor
    ~CArray ();

}; // End of CArray

//---carray class declarations---

//constructors

//no dim
template <typename T>
CArray<T>::CArray() {}

//1D
template <typename T>
CArray<T>::CArray(size_t some_dim1) {
    // assert(some_dim1 > 0);
    dim1_ = some_dim1;
    length_ = dim1_;
    this_array_ = new T[length_];
}

//2D
template <typename T>
CArray<T>::CArray(size_t some_dim1, size_t some_dim2) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = dim1_ * dim2_;
    this_array_ = new T[length_];
}

//3D
template <typename T>
CArray<T>::CArray(size_t some_dim1, size_t some_dim2, size_t some_dim3) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = dim1_ * dim2_ * dim3_;
    this_array_ = new T[length_];
}

//4D
template <typename T>
CArray<T>::CArray(size_t some_dim1, size_t some_dim2, size_t some_dim3, 
                        size_t some_dim4) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    // assert(some_dim4 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = dim1_ * dim2_ * dim3_ * dim4_;
    this_array_ = new T[length_];
}

//5D
template <typename T>
CArray<T>::CArray(size_t some_dim1, size_t some_dim2, size_t some_dim3,
                        size_t some_dim4, size_t some_dim5) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    // assert(some_dim4 > 0):
    // assert(some_dim5 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_;
    this_array_ = new T[length_];
}

//6D
template <typename T>
CArray<T>::CArray(size_t some_dim1, size_t some_dim2, size_t some_dim3,
                        size_t some_dim4, size_t some_dim5, size_t some_dim6) {
    // assert(some_dim1 > 0);
    // assert(some_dim2 > 0);
    // assert(some_dim3 > 0);
    // assert(some_dim4 > 0):
    // assert(some_dim5 > 0);
    // assert(some_dim6 > 0);
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_;
    this_array_ = new T[length_];

}

//overload () operator

//1D
template <typename T>
inline T& CArray<T>::operator() (size_t i) const {
    assert(i < dim1_);
    return this_array_[i];
}

//2D
template <typename T>
inline T& CArray<T>::operator() (size_t i, size_t j) const {
    assert(i < dim1_);
    assert(j < dim2_);
    return this_array_[j + (i * dim2_)];
}

//3D
template <typename T>
inline T& CArray<T>::operator() (size_t i, size_t j, size_t k) const {
    assert(i < dim1_);
    assert(j < dim2_);
    assert(k < dim3_);
    return this_array_[k + (j * dim3_) + (i * dim3_ * dim2_)];
}

//4D
template <typename T>
inline T& CArray<T>::operator() (size_t i, size_t j, size_t k, size_t l) const {
    assert(i < dim1_);
    assert(j < dim2_);
    assert(k < dim3_);
    assert(l < dim4_);
    return this_array_[l + (k * dim4_) 
                         + (j * dim4_ * dim3_)  
                         + (i * dim4_ * dim3_ * dim2_)];
}

//5D
template <typename T>
inline T& CArray<T>::operator() (size_t i, size_t j, size_t k, size_t l,
                                    size_t m) const {
    assert(i < dim1_);
    assert(j < dim2_);
    assert(k < dim3_);
    assert(l < dim4_);
    assert(m < dim5_);
    return this_array_[m + (l * dim5_) 
                         + (k * dim5_ * dim4_) 
                         + (j * dim5_ * dim4_ * dim3_) 
                         + (i * dim5_ * dim4_ * dim3_ * dim2_)];
}

//6D
template <typename T>
inline T& CArray<T>::operator() (size_t i, size_t j, size_t k, size_t l,
                                    size_t m, size_t n) const {
    assert(i < dim1_);
    assert(j < dim2_);
    assert(k < dim3_);
    assert(l < dim4_);
    assert(m < dim5_);
    assert(n < dim6_);
    return this_array_[n + (m * dim6_) 
                         + (l * dim6_ * dim5_)  
                         + (k * dim6_ * dim5_ * dim4_) 
                         + (j * dim6_ * dim5_ * dim4_ * dim3_)  
                         + (i * dim6_ * dim5_ * dim4_ * dim3_ * dim2_)];

}

//overload = operator
template <typename T>
inline CArray<T>& CArray<T>::operator= (const CArray& temp)
{
    // Do nothing if the assignment is of the form x = x
    if (this != &temp) {
        dim1_ = temp.dim1_;
        dim2_ = temp.dim2_;
        dim3_ = temp.dim3_;
        dim4_ = temp.dim4_;
        dim5_ = temp.dim5_;
        dim6_ = temp.dim6_;
        length_ = temp.length_;
        this_array_ = new T[length_];
    }
    return *this;
}

//return size
template <typename T>
inline size_t CArray<T>::size() const {
    return length_;
}

//destructor
template <typename T>
CArray<T>::~CArray() {
    delete[] this_array_;
}

//----endof carray class definitions----

//6. ViewCArray
template <typename T>
class ViewCArray {

private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_;
    size_t length_; // Length of 1D array
    T * this_array_;
    
public:
    
    // Default constructor
    ViewCArray ();
    
    //--- 1D array ---
    
    // Overloaded constructor
    ViewCArray(T *some_array,
               size_t some_dim1);
    
    // Overload  operator() to access data as array(i,j),
    // where i = [0:N-1], j = [0:N-1]
    T& operator()(size_t i) const;
    
    //--- 2D array ---
    
    // Overloaded constructor
    ViewCArray(T *some_array,
               size_t some_dim1,
               size_t some_dim2);
    
    // Overload operator() to access data as array(i,j),
    //  where i=[0:N-1], j=[0:N-1]
    T& operator()(size_t i, size_t j) const;
    
    //--- 3D array ---
    
    // Overloaded constructor
    ViewCArray(T *some_array,
               size_t some_dim1,
               size_t some_dim2,
               size_t some_dim3);
    
    // Overload operator() to access data as array(i,j,k),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1]
    T& operator()(size_t i, size_t j, size_t k) const;
    
    //--- 4D array ---
    
    // Overloaded constructor
    ViewCArray(T *some_array,
               size_t some_dim1,
               size_t some_dim2,
               size_t some_dim3,
               size_t some_dim4);
        
    // Overload operator() to access data as array(i, j, k, l),
    // where i = [0:n-1], j = [0:N-1], k = [0:N-1], l = [0:N-1]
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;
    
    //--- 5D array ---
    
    // Overloaded constructor
    ViewCArray (T *some_array,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4,
                size_t some_dim5);
        
    // Overload operator() to access data as array(i,j,k,l,m),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1], m = [0:N-1]
    T& operator() (size_t i, size_t j, size_t k, size_t l, size_t m) const;
    
    //--- 6D array ---
    
    // Overloaded constructor
    ViewCArray (T *some_array,
                size_t some_dim1,
                size_t some_dim2,
                size_t some_dim3,
                size_t some_dim4,
                size_t some_dim5,
                size_t some_dim6);
        
    // Overload operator() to access data as array(i,j,k,l,m,n),
    // where i = [0:N-1], j = [0:N-1], k = [0:N-1], l = [0:N-1], m = [0:N-1]
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const;
    
    size_t size() const;
    
}; // end of ViewCArray

//class definitions

//constructors

//no dim
template <typename T>
ViewCArray<T>::ViewCArray() {}

//1D
template <typename T>
ViewCArray<T>::ViewCArray(T *some_array,
                          size_t some_dim1)
{
    dim1_ = some_dim1;
    length_ = dim1_;
    this_array_ = some_array;
}

//2D
template <typename T>
ViewCArray<T>::ViewCArray(T *some_array,
                          size_t some_dim1,
                          size_t some_dim2)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = dim1_ * dim2_;
    this_array_ = some_array;
}

//3D
template <typename T>
ViewCArray<T>::ViewCArray (T *some_array,
                           size_t some_dim1,
                           size_t some_dim2,
                           size_t some_dim3)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = dim1_ * dim2_ * dim3_;
    this_array_ = some_array;
}

//4D
template <typename T>
ViewCArray<T>::ViewCArray(T *some_array,
                          size_t some_dim1,
                          size_t some_dim2,
                          size_t some_dim3,
                          size_t some_dim4)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_);
    this_array_ = some_array;
}

//5D
template <typename T>
ViewCArray<T>::ViewCArray(T *some_array,
                          size_t some_dim1,
                          size_t some_dim2,
                          size_t some_dim3,
                          size_t some_dim4,
                          size_t some_dim5)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_);
    this_array_ = some_array;
}

//6D
template <typename T>
ViewCArray<T>::ViewCArray(T *some_array,
                          size_t some_dim1,
                          size_t some_dim2,
                          size_t some_dim3,
                          size_t some_dim4,
                          size_t some_dim5,
                          size_t some_dim6)
{
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_;
    this_array_ = some_array;
}

//overload () operator

//1D
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i) const
{
    assert(i < dim1_ && "i is out of bounds in c_array 1D");  // die if >= dim1
    
    return this_array_[i];
}

//2D
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i, 
                                    size_t j) const
{
   
    assert(i < dim1_ && "i is out of bounds in c_array 2D");  // die if >= dim1
    assert(j < dim2_ && "j is out of bounds in c_array 2D");  // die if >= dim2
    
    return this_array_[j + (i * dim2_)];
}

//3D
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i, 
                                    size_t j, 
                                    size_t k) const
{
    assert(i < dim1_ && "i is out of bounds in c_array 3D");  // die if >= dim1
    assert(j < dim2_ && "j is out of bounds in c_array 3D");  // die if >= dim2
    assert(k < dim3_ && "k is out of bounds in c_array 3D");  // die if >= dim3
    
    return this_array_[k + (j * dim3_) 
                         + (i * dim3_ * dim2_)];
}

//4D
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i, 
                                    size_t j, 
                                    size_t k, 
                                    size_t l) const
{
    assert(i < dim1_ && "i is out of bounds in c_array 4D");  // die if >= dim1
    assert(j < dim2_ && "j is out of bounds in c_array 4D");  // die if >= dim2
    assert(k < dim3_ && "k is out of bounds in c_array 4D");  // die if >= dim3
    assert(l < dim4_ && "l is out of bounds in c_array 4D");  // die if >= dim4
    
    return this_array_[l + (k * dim4_) 
                         + (j * dim4_ * dim3_) 
                         + (i * dim4_ * dim3_ * dim2_)];
}

//5D
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i, 
                                    size_t j, 
                                    size_t k, 
                                    size_t l, 
                                    size_t m) const
{
    assert(i < dim1_ && "i is out of bounds in c_array 5D");  // die if >= dim1
    assert(j < dim2_ && "j is out of bounds in c_array 5D");  // die if >= dim2
    assert(k < dim3_ && "k is out of bounds in c_array 5D");  // die if >= dim3
    assert(l < dim4_ && "l is out of bounds in c_array 5D");  // die if >= dim4
    assert(m < dim5_ && "m is out of bounds in c_array 5D");  // die if >= dim5
    
    return this_array_[m + (l * dim5_) 
                         + (k * dim5_ * dim4_) 
                         + (j * dim5_ * dim4_ * dim3_)
                         + (i * dim5_ * dim4_ * dim3_ * dim2_)];
}

//6D
template <typename T>
inline T& ViewCArray<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const
{
    assert(i < dim1_ && "i is out of bounds in c_array 6D");  // die if >= dim1
    assert(j < dim2_ && "j is out of bounds in c_array 6D");  // die if >= dim2
    assert(k < dim3_ && "k is out of bounds in c_array 6D");  // die if >= dim3
    assert(l < dim4_ && "l is out of bounds in c_array 6D");  // die if >= dim4
    assert(m < dim5_ && "m is out of bounds in c_array 6D");  // die if >= dim5
    assert(n < dim6_ && "n is out of bounds in c_array 6D");  // die if >= dim6
    
    return this_array_[n + (m * dim6_) 
                         + (l * dim6_ * dim5_) 
                         + (k * dim6_ * dim5_ * dim4_)
                         + (j * dim6_ * dim5_ * dim4_ * dim3_) 
                         + (i * dim6_ * dim5_ * dim4_ * dim3_ * dim2_)];
}

//return size    
template <typename T>
inline size_t ViewCArray<T>::size() const {
    return length_;
}

//---end of ViewCArray class definitions----


//7. CMatrix
template <typename T>
class CMatrix {
        
private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_,length_;
    T * this_matrix;
            
public:
        
       // default constructor
       CMatrix();
       CMatrix(size_t some_dim1);
       CMatrix(size_t some_dim1, size_t some_dim2);
       CMatrix(size_t some_dim1, size_t some_dim2, size_t some_dim3);
       CMatrix(size_t some_dim1,
           size_t some_dim2,
           size_t some_dim3,
           size_t some_dim4);
       CMatrix(size_t some_dim1,
           size_t some_dim2,
           size_t some_dim3,
           size_t some_dim4,
           size_t some_dim5);
       CMatrix (size_t some_dim1,
            size_t some_dim2,
            size_t some_dim3,
            size_t some_dim4,
            size_t some_dim5,
            size_t some_dim6);
           
    //overload operators to access data
       T& operator()(size_t i);
       T& operator()(size_t i, size_t j);
       T& operator()(size_t i, size_t j, size_t k);
       T& operator()(size_t i, size_t j, size_t k, size_t l);
       T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m);
       T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n);

       //overload = operator
	CMatrix& operator= (const CMatrix &temp);
            
       // deconstructor
       ~CMatrix( );
        
}; // end of CMatrix

// CMatrix class definitions

//constructors

//no dim

//1D
template <typename T>
CMatrix<T>::CMatrix() {}

//1D
template <typename T>
CMatrix<T>::CMatrix(size_t some_dim1) {
    dim1_ = some_dim1;
    length_ = dim1_;
    this_matrix = new T[length_];
}

//2D
template <typename T>
CMatrix<T>::CMatrix(size_t some_dim1, size_t some_dim2) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = dim1_*dim2_;
    this_matrix = new T[length_];
}

//3D
template <typename T>
CMatrix<T>::CMatrix(size_t some_dim1, size_t some_dim2, size_t some_dim3) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = dim1_*dim2_*dim3_;
    this_matrix = new T[length_];
}

//4D
template <typename T>
CMatrix<T>::CMatrix(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = dim1_*dim2_*dim3_*dim4_;
    this_matrix= new T[length_];
}   

//5D
template <typename T>
CMatrix<T>::CMatrix(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = dim1_*dim2_*dim3_*dim4_*dim5_;
    this_matrix = new T[length_];
}

//6D
template <typename T>
CMatrix<T>::CMatrix(size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5, size_t some_dim6) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = dim1_*dim2_*dim3_*dim4_*dim5_*dim6_;
    this_matrix = new T[length_];
}


//overload () operator

//1D
template <typename T>
T& CMatrix<T>::operator()(size_t i)
{
    assert( i < dim1_+1 && "i is out of bounds in CMatrix 1D!");
    return this_matrix[i-1];
}

//2D
template <typename T>
T& CMatrix<T>::operator()(size_t i, size_t j)
{
    assert( i < dim1_+1 && "i is out of bounds in CMatrix 2D!");
    assert( j < dim2_+1 && "j is out of bounds in CMatrix 2D!");
    return this_matrix[(j-1) + (i-1)*dim2_];
}

//3D
template <typename T>
T& CMatrix<T>::operator()(size_t i, size_t j, size_t k)
{
    assert( i < dim1_+1 && "i is out of bounds in CMatrix 3D!");
    assert( j < dim2_+1 && "j is out of bounds in CMatrix 3D!");
    assert( k < dim3_+1 && "k is out of bounds in CMatrix 3D!");
    return this_matrix[(k-1) + (j-1)*dim3_ + (i-1)*dim3_*dim2_];
}

//4D
template <typename T>
T& CMatrix<T>::operator()(size_t i, size_t j, size_t k, size_t l)
{
    assert( i < dim1_+1 && "i is out of bounds in CMatrix 4D!");
    assert( j < dim2_+1 && "j is out of bounds in CMatrix 4D!");
    assert( k < dim3_+1 && "k is out of bounds in CMatrix 4D!");
    assert( l < dim4_+1 && "l is out of bounds in CMatrix 4D!");
    return this_matrix[ (l-1) + (k-1)*dim4_ + (j-1)*dim4_*dim3_ + (i-1)*dim4_*dim3_*dim2_];
}

//5D
template <typename T>
T& CMatrix<T>::operator()(size_t i, size_t j, size_t k, size_t l,size_t m)
{
    assert( i < dim1_+1 && "i is out of bounds in CMatrix 5D!");
    assert( j < dim2_+1 && "j is out of bounds in CMatrix 5D!");
    assert( k < dim3_+1 && "k is out of bounds in CMatrix 5D!");
    assert( l < dim4_+1 && "l is out of bounds in CMatrix 5D!");
    assert( m < dim5_+1 && "m is out of bounds in CMatrix 5D!");
    return this_matrix[(m-1) + (l-1)*dim5_ + (k-1)*dim5_*dim4_ + (j-1)*dim5_*dim4_*dim3_ + (i-1)*dim5_*dim4_*dim3_*dim2_];
}

//6D
template <typename T>
T& CMatrix<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n)
{
    assert( i < dim1_+1 && "i is out of bounds in CMatrix 6D!");
    assert( j < dim2_+1 && "j is out of bounds in CMatrix 6D!");
    assert( k < dim3_+1 && "k is out of bounds in CMatrix 6D!");
    assert( l < dim4_+1 && "l is out of bounds in CMatrix 6D!");
    assert( m < dim5_+1 && "m is out of bounds in CMatrix 6D!");
    assert( n < dim6_+1 && "n is out of bounds in CMatrix 6D!");
    return this_matrix[ (n-1) + (m-1)*dim6_ + (l-1)*dim6_*dim5_ + (k-1)*dim6_*dim5_*dim4_ + (j-1)*dim6_*dim5_*dim4_*dim3_ + (i-1)*dim6_*dim5_*dim4_*dim3_*dim2_];
}

//overload = operator
//THIS = CMatrix<> temp
template <typename T>
CMatrix<T> &CMatrix<T>::operator= (const CMatrix &temp) {
	if(this != &temp) {
        dim1_ = temp.dim1_;
        dim2_ = temp.dim2_;
        dim3_ = temp.dim3_;
        dim4_ = temp.dim4_;
        dim5_ = temp.dim5_;
        dim6_ = temp.dim6_;
        length_ = temp.length_;
        this_matrix = new T[length_];
	}
  return *this;
}

// Destructor
template <typename T>
CMatrix<T>::~CMatrix(){
    delete[] this_matrix;
}

//----end of CMatrix class definitions----

//8. ViewCMatrix
template <typename T>
class ViewCMatrix {

private:
    size_t dim1_, dim2_, dim3_, dim4_, dim5_, dim6_;
     T * this_matrix;
		    
public:
		    
    // default constructor
    ViewCMatrix();
		    
		    
    //--- 1D array ---	   	    
    // overloaded constructor
    ViewCMatrix (T *some_matrix,size_t some_dim1);
    T& operator() (size_t i);
		    
    //--- 2D array ---	    
    // overloaded constructor
    ViewCMatrix (T *some_matrix, size_t some_dim1, size_t some_dim2);
		    
    T& operator() (size_t i, size_t j);
		    
    //--- 3D array ---	    
    // overloaded constructor
    ViewCMatrix (T *some_matrix,
		size_t some_dim1,
		size_t some_dim2,
		size_t some_dim3);
    T& operator() (size_t i, size_t j, size_t k);
		    
    //--- 4D array ---
		    
    // overloaded constructor
    ViewCMatrix (T *some_matrix,
		size_t some_dim1,
		size_t some_dim2,
		size_t some_dim3,
		size_t some_dim4);
		    
    T& operator() (size_t i, size_t j, size_t k, size_t l);

		    
    //--- 5D array ---
    // overloaded constructor
    ViewCMatrix (T *some_matrix,
		size_t some_dim1,
		size_t some_dim2,
		size_t some_dim3,
		size_t some_dim4,
		size_t some_dim5);
    T& operator() (size_t i, size_t j, size_t k, size_t l, size_t m);

    //--- 6D array ---		    
    // overloaded constructor
    ViewCMatrix (T *some_matrix,
		   size_t some_dim1,
		   size_t some_dim2,
		   size_t some_dim3,
		   size_t some_dim4,
		   size_t some_dim5,
		   size_t some_dim6);
		    
   T& operator() (size_t i, size_t j, size_t k, size_t l, size_t m, size_t n);

		    
}; // end of ViewCMatrix

//class definitions

//constructors

//no dim
template <typename T>
ViewCMatrix<T>::ViewCMatrix(){}

//1D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *some_matrix,size_t some_dim1) {
	dim1_ = some_dim1;
	this_matrix = some_matrix;
}

//2D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *some_matrix, size_t some_dim1, size_t some_dim2) {
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	this_matrix = some_matrix;
}

//3D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *some_matrix, size_t some_dim1, size_t some_dim2, size_t some_dim3) {
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	dim3_ = some_dim3;
	this_matrix = some_matrix;
}

//4D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *some_matrix, size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4){
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	dim3_ = some_dim3;
	dim4_ = some_dim4;
	this_matrix = some_matrix;
}

//5D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *some_matrix, size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5){
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	dim3_ = some_dim3;
	dim4_ = some_dim4;
	dim5_ = some_dim5;
	this_matrix = some_matrix;
}

//6D
template <typename T>
ViewCMatrix<T>::ViewCMatrix(T *some_matrix, size_t some_dim1, size_t some_dim2, size_t some_dim3, size_t some_dim4, size_t some_dim5, size_t some_dim6) {
	dim1_ = some_dim1;
	dim2_ = some_dim2;
	dim3_ = some_dim3;
	dim4_ = some_dim4;
	dim5_ = some_dim5;
	dim6_ = some_dim6;
	this_matrix = some_matrix;
}

//overload () operator

//1D
template <typename T>
T& ViewCMatrix<T>:: operator() (size_t i)
{
	assert(i < dim1_+1 && "i is out of bounds for ViewCMatrix 1D!");
	return this_matrix[i-1];
}

//2D
template <typename T>
T& ViewCMatrix<T>::operator() (size_t i, size_t j)
{
	assert(i < dim1_+1 && "i is out of bounds for ViewCMatrix 2D!");
	assert(j < dim2_+1 && "j is out of bounds for ViewCMatrix 2D!");
	return this_matrix[(i-1)*dim2_ + (j-1)];
}

//3D
template <typename T>
T& ViewCMatrix<T>::operator () (size_t i, size_t j, size_t k)
{
	assert(i < dim1_+1 && "i is out of bounds for ViewCMatrix 3D!");
	assert(j < dim2_+1 && "j is out of bounds for ViewCMatrix 3D!");
	assert(k < dim3_+1 && "k is out of bounds for ViewCMatrix 3D!");
	return this_matrix[(k-1) + (j-1)*dim3_ + (i-1)*dim3_*dim2_];
}

//4D
template <typename T>
T& ViewCMatrix<T>::operator()(size_t i, size_t j, size_t k, size_t l)
{
	assert(i < dim1_+1 && "i is out of bounds for ViewCMatrix 4D!");
	assert(j < dim2_+1 && "j is out of bounds for ViewCMatrix 4D!");
	assert(k < dim3_+1 && "k is out of bounds for ViewCMatrix 4D!");
	assert(l < dim4_+1 && "l is out of bounds for ViewCMatrix 4D!");
	return this_matrix[(l-1) + (k-1)*dim4_ + (j-1)*dim4_*dim3_ + (i-1)*dim4_*dim3_*dim2_];
}

//5D
template <typename T>
T& ViewCMatrix<T>::operator()(size_t i, size_t j, size_t k,size_t l, size_t m)
{
	assert(i < dim1_+1 && "i is out of bounds for ViewCMatrix 5D!");
	assert(j < dim2_+1 && "j is out of bounds for ViewCMatrix 5D!");
	assert(k < dim3_+1 && "k is out of bounds for ViewCMatrix 5D!");
	assert(l < dim4_+1 && "l is out of bounds for ViewCMatrix 5D!");
	assert(m < dim5_+1 && "m is out of bounds for ViewCMatrix 5D!");
	return this_matrix[(m-1) + (l-1)*dim5_ + (k-1)*dim5_*dim4_ + (j-1)*dim5_*dim4_*dim3_ + (i-1)*dim5_*dim4_*dim3_*dim2_];
}

//6D
template <typename T>
T& ViewCMatrix<T>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n)
{
	assert(i < dim1_+1 && "i is out of bounds for ViewCMatrix 6D!");
	assert(j < dim2_+1 && "j is out of bounds for ViewCMatrix 6D!");
	assert(k < dim3_+1 && "k is out of bounds for ViewCMatrix 6D!");
	assert(l < dim4_+1 && "l is out of bounds for ViewCMatrix 6D!");
	assert(m < dim5_+1 && "m is out of bounds for ViewCMatrix 6D!");
	assert(n < dim6_+1 && "n is out of bounds for ViewCMatrix 6D!");
	return this_matrix[(n-1)+ (m-1)*dim6_ + (l-1)*dim5_*dim6_ + (k-1)*dim6_*dim5_*dim4_ + (j-1)*dim6_*dim5_*dim4_*dim3_ + (i-1)*dim5_*dim6_*dim4_*dim3_*dim2_];
}


//----end of ViewCMatrix class definitions----

//9. RaggedRightArray
template <typename T>
class RaggedRightArray {
private:
    size_t *start_index_;
    T * array_;
    
    size_t dim1_, length_;
    
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
    
    // A method to return the stride size
    size_t stride(size_t i) const;
    
    // Overload operator() to access data as array(i,j)
    // where i=[0:N-1], j=[stride(i)]
    T& operator()(size_t i, size_t j) const;

    RaggedRightArray& operator= (const RaggedRightArray &temp);

    // Destructor
    ~RaggedRightArray ( );
}; // End of RaggedRightArray

// Default constructor
template <typename T>
RaggedRightArray<T>::RaggedRightArray () {}


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
    
    array_ = new T[count];
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
    
    array_ = new T[count];
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
    
    array_ = new T[count];
} // End constructor

// A method to return the stride size
template <typename T>
inline size_t RaggedRightArray<T>::stride(size_t i) const {
    // Ensure that i is within bounds
    assert(i < (dim1_ + 1) && "i is greater than dim1_ in RaggedRightArray");

    return start_index_[(i + 1)] - start_index_[i];
}

// Overload operator() to access data as array(i,j)
// where i=[0:N-1], j=[0:stride(i)]
template <typename T>
inline T& RaggedRightArray<T>::operator()(size_t i, size_t j) const {
    // get the 1D array index
    size_t start = start_index_[i];
    
    // asserts
    assert(i < dim1_ && "i is out of dim1 bounds in RaggedRightArray");  // die if >= dim1
    assert(j < stride(i) && "j is out of stride bounds in RaggedRightArray");  // die if >= stride
    
    return array_[j + start];
} // End operator()

//overload = operator
template <typename T>
RaggedRightArray<T> & RaggedRightArray<T>::operator= (const RaggedRightArray &temp) {

    if( this != &temp) {
        dim1_ = temp.dim1_;
        length_ = temp.length_;
        start_index_ = new size_t[dim1_ + 1];
        for (int j = 0; j < dim1_; j++) {
            start_index_[j] = temp.start_index_[j];  
        }
        array_ = new T[length_];
    }
	
    return *this;
}

// Destructor
template <typename T>
RaggedRightArray<T>::~RaggedRightArray () {
    delete[] array_;
    delete[] start_index_;
}

//----end of RaggedRightArray class definitions----

//10. RaggedDownArray
template <typename T>
class RaggedDownArray { 
private:
    size_t *start_index_;
	T * array_;

	size_t dim2_;
    size_t length_;

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

	//method to return stride size
	size_t stride(size_t j);

	//overload () operator to access data as array (i,j)
	T& operator()(size_t i, size_t j);

    // method to return total size
    size_t size();

    //overload = operator
    RaggedDownArray& operator= (const RaggedDownArray &temp);

    //destructor
    ~RaggedDownArray();

}; //~~~~~end of RaggedDownArray class declarations~~~~~~~~	

//no dims
template <typename T>
RaggedDownArray<T>::RaggedDownArray() {}

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

    array_ = new T[count];

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

// Check the stride size
template <typename T>
size_t RaggedDownArray<T>::stride(size_t j) {
    return start_index_[j+1] - start_index_[j];
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
		
    return array_[i + start];

} // End () operator

//overload = operator
template <typename T>
RaggedDownArray<T> & RaggedDownArray<T>::operator= (const RaggedDownArray &temp) {

    if( this != &temp) {
        dim2_ = temp.dim2_;
        length_ = temp.length_;
        start_index_ = new size_t[dim2_ + 1];
        for (int j = 0; j < dim2_; j++) {
            start_index_[j] = temp.start_index_[j];  
        }
        array_ = new T[length_];
    }
	
    return *this;
}

// Destructor
template <typename T>
RaggedDownArray<T>::~RaggedDownArray() {
    delete[] array_;
    delete[] start_index_;

} // End destructor


//----end of RaggedDownArray----

//11. DynamicRaggedRightArray
/*
template <typename T>
class DynamicRaggedRightArrayKokkos {
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
    
    // A method to increase the stride size
    void push_back(size_t i) const;
    
    // Overload operator() to access data as array(i,j),
    // where i=[0:N-1], j=[stride(i)]
    T& operator()(size_t i, size_t j) const;
    
    // Destructor
    ~DynamicRaggedRightArray ();
}; 

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


// A method to increase the stride size for row i
template <typename T>
void DynamicRaggedRightArray<T>::push_back(size_t i) const {
    stride_[i]++;
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

// Destructor
template <typename T>
DynamicRaggedRightArray<T>::~DynamicRaggedRightArray() {
    delete[] array_;
    delete[] stride_;
}
*/



//----end DynamicRaggedRightArray class definitions----

//12. DynamicRaggedDownArray


//----end of DynamicRaggedDownArray class definitions-----

//13. SparseRowArray
template <typename T>
class SparseRowArray {
private:
    size_t *start_index_;
    size_t *column_index_;
    
    T * array_;
    
    size_t dim1_;
    
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
    size_t column_index(size_t i, size_t j) const;
    
    // A method to access data as array.value(i,j),
    // where i=[0:N-1], j=[stride(i)]
    T& value(size_t i, size_t j) const;
    
    // Destructor
    ~SparseRowArray ();
}; 


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
size_t SparseRowArray<T>::column_index(size_t i, size_t j) const {
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

// Destructor
template <typename T>
SparseRowArray<T>::~SparseRowArray() {
    delete[] array_;
    delete[] start_index_;
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

	size_t dim2_;

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

	//destructor
	~SparseColArray();
};

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

//destructor
template <typename T>
SparseColArray<T>::~SparseColArray() {
	delete [] array_;
	delete [] start_index_;
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
template <typename T>
class FArrayKokkos {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
private:

    /*!
     * Private variable that specifies the length of the first dimension
     * of the FArrayKokkos object.
     */
    size_t dim1_;
    /*!
     * Private variable that specifies the length of the second dimension
     * of the FArrayKokkos object.
     */
    size_t dim2_;
    size_t dim3_;
    size_t dim4_;
    size_t dim5_;
    size_t dim6_;
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

        \param some_dim1 the length of the first dimension
     */
    FArrayKokkos(size_t some_dim1);

    /*!
     * \brief An overloaded constructor used to construct a 2D FArrayKokkos
              object.

        \param some_dim1 the length of the first dimension
        \param some_dim2 the length of the second dimension
     */
    FArrayKokkos(size_t some_dim1, size_t some_dim2);

    /*!
     * \brief An overloaded constructor used to construct a 3D FArrayKokkos
              object.

        \param some_dim1 the length of the first dimension
        \param some_dim2 the length of the second dimension
        \param some_dim3 the length of the third dimension
     */
    FArrayKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3);

    FArrayKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, 
                 size_t some_dim4);

    FArrayKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, 
                 size_t some_dim4, size_t some_dim5); 

    FArrayKokkos(size_t some_dim1, size_t sone_dim2, size_t some_dim3, 
                 size_t some_dim4, size_t some_dim5, size_t some_dim6);

    // Overload operator() to acces data
    // from 1D to 6D
    
    KOKKOS_FUNCTION
    T& operator()(size_t i) const;
    
    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j, size_t k) const;

    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const;

    // Overload = operator
    FArrayKokkos& operator= (const FArrayKokkos &temp);

    KOKKOS_FUNCTION
    size_t size();

    size_t extent();

    T* pointer();

    // Destructor
    KOKKOS_FUNCTION
    ~FArrayKokkos();    

}; //end of FArrayKokkos declarations

// Default constructor
template <typename T>
FArrayKokkos<T>::FArrayKokkos() {}

// Overloaded 1D constructor
template <typename T>
FArrayKokkos<T>::FArrayKokkos(size_t some_dim1){
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    length_ = dim1_;
    this_array_ = TArray1D("this_array_", length_);
}

// Overloaded 2D constructor
template <typename T>
FArrayKokkos<T>::FArrayKokkos(size_t some_dim1, size_t some_dim2) {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = (dim1_ * dim2_);
    this_array_ = TArray1D("this_array_", length_);
}

// Overloaded 3D constructor
template <typename T>
FArrayKokkos<T>::FArrayKokkos(size_t some_dim1, size_t some_dim2, 
                              size_t some_dim3) {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = (dim1_ * dim2_ * dim3_);
    this_array_ = TArray1D("this_array_", length_);
}

// Overloaded 4D constructor
template <typename T>
FArrayKokkos<T>::FArrayKokkos(size_t some_dim1, size_t some_dim2, 
                              size_t some_dim3, size_t some_dim4) {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_);
    this_array_ = TArray1D("this_array_", length_);
}

// Overloaded 5D constructor
template <typename T>
FArrayKokkos<T>::FArrayKokkos(size_t some_dim1, size_t some_dim2, 
                              size_t some_dim3, size_t some_dim4, 
                              size_t some_dim5) {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_);
    this_array_ = TArray1D("this_array_", length_);
}

// Overloaded 6D constructor
template <typename T>
FArrayKokkos<T>::FArrayKokkos(size_t some_dim1, size_t some_dim2, 
                              size_t some_dim3, size_t some_dim4, 
                              size_t some_dim5, size_t some_dim6) {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_);
    this_array_ = TArray1D("this_array_", length_);
}

// Definitions of overload operator()
// for 1D to 6D
// Note: the indices for array all start at 0

// 1D
template<typename T>
KOKKOS_FUNCTION
T& FArrayKokkos<T>::operator()( size_t i) const {
    assert( i < dim1_ && "i is out of bounds in FArrayKokkos 1D!");
    return this_array_(i);
}

// 2D
template <typename T>
KOKKOS_FUNCTION
T& FArrayKokkos<T>::operator()(size_t i, size_t j) const {
    assert( i < dim1_ && "i is out of bounds in FArrayKokkos 2D!");
    assert( j < dim2_ && "j is out of bounds in FArrayKokkos 2D!");
    return this_array_(i + (j * dim1_));
}

// 3D
template <typename T>
KOKKOS_FUNCTION
T& FArrayKokkos<T>::operator()(size_t i, size_t j, size_t k) const {
    assert( i < dim1_ && "i is out of bounds in FArrayKokkos 3D!");
    assert( j < dim2_ && "j is out of bounds in FArrayKokkos 3D!");
    assert( k < dim3_ && "k is out of bounds in FArrayKokkos 3D!");
    return this_array_(i + (j * dim1_) 
                         + (k * dim1_ * dim2_));
}

// 4D
template <typename T>
KOKKOS_FUNCTION
T& FArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    assert( i < dim1_ && "i is out of bounds in FArrayKokkos 4D!");
    assert( j < dim2_ && "j is out of bounds in FArrayKokkos 4D!");
    assert( k < dim3_ && "k is out of bounds in FArrayKokkos 4D!");
    assert( l < dim4_ && "l is out of bounds in FArrayKokkos 4D!");
    return this_array_(i + (j * dim1_) 
                         + (k * dim1_ * dim2_) 
                         + (l * dim1_ * dim2_ * dim3_));
}

// 5D
template <typename T>
KOKKOS_FUNCTION
T& FArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, 
                               size_t m) const {
    assert( i < dim1_ && "i is out of bounds in FArrayKokkos 5D!");
    assert( j < dim2_ && "j is out of bounds in FArrayKokkos 5D!");
    assert( k < dim3_ && "k is out of bounds in FArrayKokkos 5D!");
    assert( l < dim4_ && "l is out of bounds in FArrayKokkos 5D!");
    assert( m < dim5_ && "m is out of bounds in FArrayKokkos 5D!");
    return this_array_(i + (j * dim1_) 
                         + (k * dim1_ * dim2_) 
                         + (l * dim1_ * dim2_ * dim3_) 
                         + (m * dim1_ * dim2_ * dim3_ * dim4_));
}

// 6D
template <typename T>
KOKKOS_FUNCTION
T& FArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, 
                               size_t m, size_t n) const {
    assert( i < dim1_ && "i is out of bounds in FArrayKokkos 6D!");
    assert( j < dim2_ && "j is out of bounds in FArrayKokkos 6D!");
    assert( k < dim3_ && "k is out of bounds in FArrayKokkos 6D!");
    assert( l < dim4_ && "l is out of bounds in FArrayKokkos 6D!");
    assert( m < dim5_ && "m is out of bounds in FArrayKokkos 6D!");
    assert( n < dim6_ && "n is out of bounds in FArrayKokkos 6D!");
    return this_array_(i + (j * dim1_) 
                         + (k * dim1_ * dim2_) 
                         + (l * dim1_ * dim2_ * dim3_) 
                         + (m * dim1_ * dim2_ * dim3_ * dim4_) 
                         + (n * dim1_ * dim2_ * dim3_ * dim4_ * dim5_));
}

// Overload = operator
// for object assingment THIS = FArrayKokkos<> TEMP(n,m,,,,)
template <typename T>
FArrayKokkos<T>& FArrayKokkos<T>::operator= (const FArrayKokkos& temp) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    if (this != &temp) {
      dim1_ = temp.dim1_;
      dim2_ = temp.dim2_;
      dim3_ = temp.dim3_;
      dim4_ = temp.dim4_;
      dim5_ = temp.dim5_;
      dim6_ = temp.dim6_;
      length_ = temp.length_;
      this_array_ = TArray1D("this_array_", length_);
    }
    return *this;
}

template <typename T>
KOKKOS_FUNCTION
size_t FArrayKokkos<T>::size() {
    return length_;
}

template <typename T>
size_t FArrayKokkos<T>::extent() {
    return length_;
}

template <typename T>
T* FArrayKokkos<T>::pointer() {
    return this_array_.data();
}

// Destructor
template <typename T>
KOKKOS_FUNCTION
FArrayKokkos<T>::~FArrayKokkos() {}

////////////////////////////////////////////////////////////////////////////////
// End of FArrayKokkos
////////////////////////////////////////////////////////////////////////////////

/*! \brief Kokkos version of the serial ViewFArray class.
 *
 */
template <typename T>
class ViewFArrayKokkos {

private: 
    size_t dim1_;
    size_t dim2_;
    size_t dim3_;
    size_t dim4_;
    size_t dim5_;
    size_t dim6_;
    size_t length_;
    T* this_array_;

public:
    ViewFArrayKokkos();

    ViewFArrayKokkos(T* some_array, size_t dim1);
    
    ViewFArrayKokkos(T* some_array, size_t dim1, size_t dim2);
    
    ViewFArrayKokkos(T* some_array, size_t dim1, size_t dim2, size_t dim3);
    
    ViewFArrayKokkos(T* some_array, size_t dim1, size_t dim2, size_t dim3, 
                     size_t dim4);
    
    ViewFArrayKokkos(T* some_array, size_t dim1, size_t dim2, size_t dim3, 
                     size_t dim4, size_t dim5);
    
    ViewFArrayKokkos(T* some_array, size_t dim1, size_t dim2, size_t dim3, 
                     size_t dim4, size_t dim5, size_t dim6);

    KOKKOS_FUNCTION
    T& operator()(size_t i) const; 

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const; 

    KOKKOS_FUNCTION
    size_t size();

    size_t extent();

    KOKKOS_FUNCTION
    ~ViewFArrayKokkos();

}; // End of ViewFArrayKokkos declarations

// Default constructor
template <typename T>
ViewFArrayKokkos<T>::ViewFArrayKokkos() {}

// Overloaded 1D constructor
template <typename T>
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim1) {
    dim1_ = dim1;
    length_ = dim1_;
    this_array_ = some_array;
}

// Overloaded 2D constructor
template <typename T>
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim1, size_t dim2) {
    dim1_ = dim1;
    dim2_ = dim2;
    length_ = (dim1_ * dim2_);
    this_array_ = some_array;
}

// Overloaded 3D constructor
template <typename T>
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim1, size_t dim2, 
                                      size_t dim3) {
    dim1_ = dim1;
    dim2_ = dim2;
    dim3_ = dim3;
    length_ = (dim1_ * dim2_ * dim3_);
    this_array_ = some_array;
}

// Overloaded 4D constructor
template <typename T>
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim1, size_t dim2, 
                                      size_t dim3, size_t dim4) {
    dim1_ = dim1;
    dim2_ = dim2;
    dim3_ = dim3;
    dim4_ = dim4;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_);
    this_array_ = some_array;
}

// Overloaded 5D constructor
template <typename T>
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim1, size_t dim2, 
                                      size_t dim3, size_t dim4, size_t dim5) {
    dim1_ = dim1;
    dim2_ = dim2;
    dim3_ = dim3;
    dim4_ = dim4;
    dim5_ = dim5;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_);
    this_array_ = some_array;
}

// Overloaded 6D constructor
template <typename T>
ViewFArrayKokkos<T>::ViewFArrayKokkos(T *some_array, size_t dim1, size_t dim2, 
                                      size_t dim3, size_t dim4, size_t dim5, 
                                      size_t dim6) {
    dim1_ = dim1;
    dim2_ = dim2;
    dim3_ = dim3;
    dim4_ = dim4;
    dim5_ = dim5;
    dim6_ = dim6;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_);
    this_array_ = some_array;
}

// Overloaded operator() for 1D array access
template <typename T>
KOKKOS_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i) const {
    assert( i < dim1_ && "i is out of bounds in ViewFArrayKokkos 1D!");
    return this_array_[i];
}

//2D
template <typename T>
KOKKOS_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i, size_t j) const {
    assert(i < dim1_ && "i is out of bounds in ViewFArrayKokkos 2D!");
    assert(j < dim2_ && "j is out of bounds in ViewFArrayKokkos 2D!");
    return this_array_[i + (j * dim1_)];
}

//3D
template <typename T>
KOKKOS_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i, size_t j, size_t k) const {
    assert(i < dim1_ && "i is out of bounds in ViewFArrayKokkos 3D!");
    assert(j < dim2_ && "j is out of bounds in ViewFArrayKokkos 3D!");
    assert(k < dim3_ && "k is out of bounds in ViewFArrayKokkos 3D!");
    return this_array_[i + (j * dim1_) 
                         + (k * dim1_ * dim2_)];
}

//4D
template <typename T>
KOKKOS_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, 
                                   size_t l) const {
    assert(i < dim1_ && "i is out of bounds in ViewFArrayKokkos 4D!");
    assert(j < dim2_ && "j is out of bounds in ViewFArrayKokkos 4D!");
    assert(k < dim3_ && "k is out of bounds in ViewFArrayKokkos 4D!");
    assert(l < dim4_ && "l is out of bounds in ViewFArrayKokkos 4D!");
    return this_array_[i + (j * dim1_) 
                         + (k * dim1_ * dim2_) 
                         + (l * dim1_ * dim2_ *dim3_)];
}

//5D
template <typename T>
KOKKOS_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, 
                                   size_t l, size_t m) const {
    assert(i < dim1_ && "i is out of bounds in ViewFArrayKokkos 5D!");
    assert(j < dim2_ && "j is out of bounds in ViewFArrayKokkos 5D!");
    assert(k < dim3_ && "k is out of bounds in ViewFArrayKokkos 5D!");
    assert(l < dim4_ && "l is out of bounds in ViewFArrayKokkos 5D!");
    assert(m < dim5_ && "m is out of bounds in ViewFArrayKokkos 5D!");
    return this_array_[i + (j * dim1_) 
                         + (k * dim1_ * dim2_) 
                         + (l * dim1_ * dim2_ * dim3_) 
                         + (m * dim1_ * dim2_ * dim3_ * dim4_)];
}

//6D
template <typename T>
KOKKOS_FUNCTION
T& ViewFArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, 
                                   size_t l, size_t m, size_t n) const {
    assert(i < dim1_ && "i is out of bounds in ViewFArrayKokkos 6D!");
    assert(j < dim2_ && "j is out of bounds in ViewFArrayKokkos 6D!");
    assert(k < dim3_ && "k is out of bounds in ViewFArrayKokkos 6D!");
    assert(l < dim4_ && "l is out of bounds in ViewFArrayKokkos 6D!");
    assert(m < dim5_ && "m is out of bounds in ViewFArrayKokkos 6D!");
    assert(n < dim6_ && "n is out of bounds in ViewFArrayKokkos 6D!");
    return this_array_[i + (j * dim1_) 
                         + (k * dim1_ * dim2_) 
                         + (l * dim1_ * dim2_ * dim3_) 
                         + (m * dim1_ * dim2_ * dim3_ * dim4_)
                         + (n * dim1_ * dim2_ * dim3_ * dim4_ * dim5_)];
}

template <typename T>
KOKKOS_FUNCTION
size_t ViewFArrayKokkos<T>::size() {
    return length_;
}

template <typename T>
size_t ViewFArrayKokkos<T>::extent() {
    return length_;
}

template <typename T>
KOKKOS_FUNCTION
ViewFArrayKokkos<T>::~ViewFArrayKokkos() {}

////////////////////////////////////////////////////////////////////////////////
// End of ViewFArrayKokkos
////////////////////////////////////////////////////////////////////////////////

/*! \brief Kokkos version of the serial FMatrix class.
 *
 */
template <typename T>
class FMatrixKokkos {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
private:

    size_t dim1_;
    size_t dim2_;
    size_t dim3_;
    size_t dim4_;
    size_t dim5_;
    size_t dim6_;
    size_t length_; 
    TArray1D this_matrix_; 

public:
    FMatrixKokkos();

    FMatrixKokkos(size_t some_dim1);

    FMatrixKokkos(size_t some_dim1, size_t some_dim2);

    FMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3);

    FMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, 
                  size_t some_dim4);

    FMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3,
                  size_t some_dim4, size_t some_dim5);

    FMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3,
                  size_t some_dim4, size_t some_dim5, size_t some_dim6);

    KOKKOS_FUNCTION
    T& operator()(size_t i) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, 
                  size_t n) const;

    FMatrixKokkos& operator=(const FMatrixKokkos& temp);

    KOKKOS_FUNCTION
    size_t size();

    size_t extent();

    T* pointer();

    KOKKOS_FUNCTION
    ~FMatrixKokkos();

}; // End of FMatrixKokkos

// Default constructor
template <typename T>
FMatrixKokkos<T>::FMatrixKokkos() {}

// Overloaded 1D constructor
template <typename T>
FMatrixKokkos<T>::FMatrixKokkos(size_t some_dim1) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    length_ = dim1_;
    this_matrix_ = TArray1D("this_matrix_", length_);
}

// Overloaded 2D constructor
template <typename T>
FMatrixKokkos<T>::FMatrixKokkos(size_t some_dim1, size_t some_dim2) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = (dim1_ * dim2_);
    this_matrix_ = TArray1D("this_matrix_", length_);
}

// Overloaded 3D constructor
template <typename T>
FMatrixKokkos<T>::FMatrixKokkos(size_t some_dim1, size_t some_dim2, 
                                size_t some_dim3) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = (dim1_ * dim2_ * dim3_);
    this_matrix_ = TArray1D("this_matrix_", length_);
}

// Overloaded 4D constructor
template <typename T>
FMatrixKokkos<T>::FMatrixKokkos(size_t some_dim1, size_t some_dim2, 
                                size_t some_dim3, size_t some_dim4) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_);
    this_matrix_ = TArray1D("this_matrix_", length_);
}

// Overloaded 5D constructor
template <typename T>
FMatrixKokkos<T>::FMatrixKokkos(size_t some_dim1, size_t some_dim2, 
                                size_t some_dim3, size_t some_dim4, 
                                size_t some_dim5) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_);
    this_matrix_ = TArray1D("this_matrix_", length_);
}

// Overloaded 5D constructor
template <typename T>
FMatrixKokkos<T>::FMatrixKokkos(size_t some_dim1, size_t some_dim2, 
                                size_t some_dim3, size_t some_dim4, 
                                size_t some_dim5, size_t some_dim6) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_);
    this_matrix_ = TArray1D("this_matrix_", length_);
}

template <typename T>
KOKKOS_FUNCTION
T& FMatrixKokkos<T>::operator()(size_t i) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in FMatrixKokkos in 1D!");
    return this_matrix_((i - 1));
}

template <typename T>
KOKKOS_FUNCTION
T& FMatrixKokkos<T>::operator()(size_t i, size_t j) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in FMatrixKokkos in 2D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in FMatrixKokkos in 2D!");
    return this_matrix_((i - 1) + ((j - 1) * dim1_));
}

template <typename T>
KOKKOS_FUNCTION
T& FMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in FMatrixKokkos in 3D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in FMatrixKokkos in 3D!");
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in FMatrixKokkos in 3D!");
    return this_matrix_((i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_));
}

template <typename T>
KOKKOS_FUNCTION
T& FMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in FMatrixKokkos in 4D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in FMatrixKokkos in 4D!");
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in FMatrixKokkos in 4D!");
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in FMatrixKokkos in 4D!");
    return this_matrix_((i - 1) + ((j - 1) * dim1_)  
                                + ((k - 1) * dim1_ * dim2_)  
                                + ((l - 1) * dim1_ * dim2_ * dim3_));
}

template <typename T>
KOKKOS_FUNCTION
T& FMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l,
                                size_t m) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in FMatrixKokkos in 5D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in FMatrixKokkos in 5D!");
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in FMatrixKokkos in 5D!");
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in FMatrixKokkos in 5D!");
    assert(m >= 1 && m <= dim5_ && "m is out of bounds in FMatrixKokkos in 5D!");
    return this_matrix_((i - 1) + ((j - 1) * dim1_)  
                                + ((k - 1) * dim1_ * dim2_)  
                                + ((l - 1) * dim1_ * dim2_ * dim3_) 
                                + ((m - 1) * dim1_ * dim2_ * dim3_ * dim4_));
}

template <typename T>
KOKKOS_FUNCTION
T& FMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l,
                                size_t m, size_t n) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in FMatrixKokkos in 6D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in FMatrixKokkos in 6D!");
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in FMatrixKokkos in 6D!");
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in FMatrixKokkos in 6D!");
    assert(m >= 1 && m <= dim5_ && "m is out of bounds in FMatrixKokkos in 6D!");
    assert(n >= 1 && n <= dim6_ && "n is out of bounds in FMatrixKokkos in 6D!");
    return this_matrix_[(i - 1) + ((j - 1) * dim1_)  
                                + ((k - 1) * dim1_ * dim2_)  
                                + ((l - 1) * dim1_ * dim2_ * dim3_)  
                                + ((m - 1) * dim1_ * dim2_ * dim3_ * dim4_)  
                                + ((n - 1) * dim1_ * dim2_ * dim3_ * dim4_ * dim5_)];
}

template <typename T>
FMatrixKokkos<T>& FMatrixKokkos<T>::operator=(const FMatrixKokkos& temp) {
    // Do nothing if the assignment is of the form x = x
    if (this != &temp) {
        dim1_ = temp.dim1_;
        dim2_ = temp.dim2_;
        dim3_ = temp.dim3_;
        dim4_ = temp.dim4_;
        dim5_ = temp.dim5_;
        dim6_ = temp.dim6_;
        length_ = temp.length_;
        this_matrix_ = TArray1D("this_matrix_", length_);
    }
    return *this;
}

template <typename T>
KOKKOS_FUNCTION
size_t FMatrixKokkos<T>::size() {
    return length_;
}

template <typename T>
size_t FMatrixKokkos<T>::extent() {
    return length_;
}

template <typename T>
T* FMatrixKokkos<T>::pointer() {
    return this_matrix_.data();
}

template <typename T>
KOKKOS_FUNCTION
FMatrixKokkos<T>::~FMatrixKokkos() {}

////////////////////////////////////////////////////////////////////////////////
// End of FMatrixKokkos
////////////////////////////////////////////////////////////////////////////////

/*! \brief Kokkos version of the serial ViewFMatrix class.
 * 
 */
template <typename T>
class ViewFMatrixKokkos {

private:

    size_t dim1_;
    size_t dim2_;
    size_t dim3_;
    size_t dim4_;
    size_t dim5_;
    size_t dim6_;
    size_t length_; 
    T* this_matrix_;
    
public:
    
    ViewFMatrixKokkos();
    
    ViewFMatrixKokkos(T* some_matrix, size_t some_dim1);

    ViewFMatrixKokkos(T* some_matrix, size_t some_dim1, size_t some_dim2);

    ViewFMatrixKokkos(T* some_matrix, size_t some_dim1, size_t some_dim2,
                      size_t some_dim3);

    ViewFMatrixKokkos(T* some_matrix, size_t some_dim1, size_t some_dim2,
                      size_t some_dim3, size_t some_dim4);

    ViewFMatrixKokkos(T* some_matrix, size_t some_dim1, size_t some_dim2,
                      size_t some_dim3, size_t some_dim4, size_t some_dim5);

    ViewFMatrixKokkos(T* some_matrix, size_t some_dim1, size_t some_dim2, 
                      size_t some_dim3, size_t some_dim4, size_t some_dim5,
                      size_t some_dim6);
    
    KOKKOS_FUNCTION
    T& operator()(size_t i) const;
    
    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j) const;
    
    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;
        
    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, 
                  size_t n) const;
    
    KOKKOS_FUNCTION
    size_t size();

    size_t extent();

    KOKKOS_FUNCTION
    ~ViewFMatrixKokkos();
    
}; // end of ViewFMatrixKokkos

// Default constructor
template <typename T>
ViewFMatrixKokkos<T>::ViewFMatrixKokkos() {}

// Overloaded 1D constructor
template <typename T>
ViewFMatrixKokkos<T>::ViewFMatrixKokkos(T* some_matrix, size_t some_dim1) {
    dim1_ = some_dim1;
    length_ = dim1_;
    this_matrix_ = some_matrix;
}

// Overloaded 2D constructor
template <typename T>
ViewFMatrixKokkos<T>::ViewFMatrixKokkos(T* some_matrix, size_t some_dim1,
                                        size_t some_dim2) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = (dim1_ * dim2_);
    this_matrix_ = some_matrix;
}

// Overloaded 3D constructor
template <typename T>
ViewFMatrixKokkos<T>::ViewFMatrixKokkos(T* some_matrix, size_t some_dim1,
                                        size_t some_dim2, size_t some_dim3) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = (dim1_ * dim2_ * dim3_);
    this_matrix_ = some_matrix;
}

// Overloaded 4D constructor
template <typename T>
ViewFMatrixKokkos<T>::ViewFMatrixKokkos(T* some_matrix, size_t some_dim1,
                                        size_t some_dim2, size_t some_dim3,
                                        size_t some_dim4) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_);
    this_matrix_ = some_matrix;
}

// Overloaded 5D constructor
template <typename T>
ViewFMatrixKokkos<T>::ViewFMatrixKokkos(T* some_matrix, size_t some_dim1,
                                        size_t some_dim2, size_t some_dim3,
                                        size_t some_dim4, size_t some_dim5) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_);
    this_matrix_ = some_matrix;
}

// Overloaded 6D constructor
template <typename T>
ViewFMatrixKokkos<T>::ViewFMatrixKokkos(T* some_matrix, size_t some_dim1,
                                        size_t some_dim2, size_t some_dim3,
                                        size_t some_dim4, size_t some_dim5,
                                        size_t some_dim6) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_);
    this_matrix_ = some_matrix;
}

template <typename T>
KOKKOS_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in ViewFMatrixKokkos 1D!"); 
    return this_matrix_[(i - 1)];
}

template <typename T>
KOKKOS_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i, size_t j) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in ViewFMatrixKokkos 2D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in ViewFMatrixKokkos 2D!");  
    return this_matrix_[(i - 1) + ((j - 1) * dim1_)];
}

template <typename T>
KOKKOS_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k) const
{
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in ViewFMatrixKokkos 3D!");  
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in ViewFMatrixKokkos 3D!");  
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in ViewFMatrixKokkos 3D!"); 
    
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_)];
}

template <typename T>
KOKKOS_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, 
                                    size_t l) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in ViewFMatrixKokkos 4D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in ViewFMatrixKokkos 4D!");
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in ViewFMatrixKokkos 4D!");
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in ViewFMatrixKokkos 4D!");
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_)
                                + ((l - 1) * dim1_ * dim2_ * dim3_)];
}

template <typename T>
KOKKOS_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, 
                                    size_t m) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in ViewFMatrixKokkos 5D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in ViewFMatrixKokkos 5D!");
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in ViewFMatrixKokkos 5D!");
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in ViewFMatrixKokkos 5D!");
    assert(m >= 1 && m <= dim5_ && "m is out of bounds in ViewFMatrixKokkos 5D!");
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_) 
                                + ((l - 1) * dim1_ * dim2_ * dim3_)
                                + ((m - 1) * dim1_ * dim2_ * dim3_ * dim4_)];
}

template <typename T>
KOKKOS_FUNCTION
T& ViewFMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, 
                                    size_t m, size_t n) const
{
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in ViewFMatrixKokkos 6D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in ViewFMatrixKokkos 6D!");
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in ViewFMatrixKokkos 6D!");
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in ViewFMatrixKokkos 6D!");
    assert(m >= 1 && m <= dim5_ && "m is out of bounds in ViewFMatrixKokkos 6D!");
    assert(n >= 1 && n <= dim6_ && "n is out of bounds in ViewFMatrixKokkos 6D!");
    return this_matrix_[(i - 1) + ((j - 1) * dim1_) 
                                + ((k - 1) * dim1_ * dim2_) 
                                + ((l - 1) * dim1_ * dim2_ * dim3_)
                                + ((m - 1) * dim1_ * dim2_ * dim3_ * dim4_)
                                + ((n - 1) * dim1_ * dim2_ * dim3_ * dim4_ * dim5_)];
}

template <typename T>
KOKKOS_FUNCTION
size_t ViewFMatrixKokkos<T>::size() {
    return length_;
}

template <typename T>
size_t ViewFMatrixKokkos<T>::extent() {
    return length_;
}

template <typename T>
KOKKOS_FUNCTION
ViewFMatrixKokkos<T>::~ViewFMatrixKokkos() {}

////////////////////////////////////////////////////////////////////////////////
// End of ViewFMatrixKokkos
////////////////////////////////////////////////////////////////////////////////

/*! \brief Kokkos version of the serial CArray class.
 *
 */
template <typename T>
class CArrayKokkos {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
private:
    size_t dim1_;
    size_t dim2_;
    size_t dim3_;
    size_t dim4_;
    size_t dim5_;
    size_t dim6_;
    size_t length_;
    TArray1D this_array_; 

public:
    CArrayKokkos();
    
    CArrayKokkos(size_t some_dim1);

    CArrayKokkos(size_t some_dim1, size_t some_dim2);

    CArrayKokkos (size_t some_dim1, size_t some_dim2, size_t some_dim3);

    CArrayKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, 
                 size_t some_dim4);

    CArrayKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3,
                 size_t some_dim4, size_t some_dim5);

    CArrayKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3,
                 size_t some_dim4, size_t some_dim5, size_t some_dim6);

    KOKKOS_FUNCTION
    T& operator()(size_t i) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, 
                  size_t n) const;

    CArrayKokkos& operator=(const CArrayKokkos& temp);

    // GPU Method
    // Method that returns size
    KOKKOS_FUNCTION
    size_t size();

    // Host Method
    // Method that returns size
    size_t extent();

    // Methods returns the raw pointer (most likely GPU) of the Kokkos View
    T* pointer();

    // Deconstructor
    KOKKOS_FUNCTION
    ~CArrayKokkos ();
}; // End of CArrayKokkos

// Default constructor
template <typename T>
CArrayKokkos<T>::CArrayKokkos() {}

// Overloaded 1D constructor
template <typename T>
CArrayKokkos<T>::CArrayKokkos(size_t some_dim1) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    length_ = dim1_;
    this_array_ = TArray1D("this_array_", length_);
}

// Overloaded 2D constructor
template <typename T>
CArrayKokkos<T>::CArrayKokkos(size_t some_dim1, size_t some_dim2) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = (dim1_ * dim2_);
    this_array_ = TArray1D("this_array_", length_);
}

template <typename T>
CArrayKokkos<T>::CArrayKokkos(size_t some_dim1, size_t some_dim2, 
                              size_t some_dim3) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = (dim1_ * dim2_ * dim3_);
    this_array_ = TArray1D("this_array_", length_);
}

template <typename T>
CArrayKokkos<T>::CArrayKokkos(size_t some_dim1, size_t some_dim2, 
                              size_t some_dim3, size_t some_dim4) {
    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_);
    this_array_ = TArray1D("this_array_", length_);
}

template <typename T>
CArrayKokkos<T>::CArrayKokkos(size_t some_dim1, size_t some_dim2, 
                              size_t some_dim3, size_t some_dim4, 
                              size_t some_dim5) {

    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_);
    this_array_ = TArray1D("this_array_", length_);
}

template <typename T>
CArrayKokkos<T>::CArrayKokkos(size_t some_dim1, size_t some_dim2, 
                              size_t some_dim3, size_t some_dim4, 
                              size_t some_dim5, size_t some_dim6) {
    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_);
    this_array_ = TArray1D("this_array_", length_);
}

template <typename T>
KOKKOS_FUNCTION
T& CArrayKokkos<T>::operator()(size_t i) const {
    assert(i < dim1_ && "i is out of bounds in CArrayKokkos 1D!");
    return this_array_(i);
}

template <typename T>
KOKKOS_FUNCTION
T& CArrayKokkos<T>::operator()(size_t i, size_t j) const {
    assert(i < dim1_ && "i is out of bounds in CArrayKokkos 2D!");
    assert(j < dim2_ && "j is out of bounds in CArrayKokkos 2D!");
    return this_array_(j + (i * dim2_));
}

template <typename T>
KOKKOS_FUNCTION
T& CArrayKokkos<T>::operator()(size_t i, size_t j, size_t k) const {
    assert(i < dim1_ && "i is out of bounds in CArrayKokkos 3D!");
    assert(j < dim2_ && "j is out of bounds in CArrayKokkos 3D!");
    assert(k < dim3_ && "k is out of bounds in CArrayKokkos 3D!");
    return this_array_(k + (j * dim3_) 
                         + (i * dim3_ * dim2_));
}

template <typename T>
KOKKOS_FUNCTION
T& CArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    assert(i < dim1_ && "i is out of bounds in CArrayKokkos 4D!");
    assert(j < dim2_ && "j is out of bounds in CArrayKokkos 4D!");
    assert(k < dim3_ && "k is out of bounds in CArrayKokkos 4D!");
    assert(l < dim4_ && "l is out of bounds in CArrayKokkos 4D!");
    return this_array_(l + (k * dim4_) 
                         + (j * dim4_ * dim3_)  
                         + (i * dim4_ * dim3_ * dim2_));
}

template <typename T>
KOKKOS_FUNCTION
T& CArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m) const {
    assert(i < dim1_ && "i is out of bounds in CArrayKokkos 5D!");
    assert(j < dim2_ && "j is out of bounds in CArrayKokkos 5D!");
    assert(k < dim3_ && "k is out of bounds in CArrayKokkos 5D!");
    assert(l < dim4_ && "l is out of bounds in CArrayKokkos 5D!");
    assert(m < dim5_ && "m is out of bounds in CArrayKokkos 5D!");
    return this_array_(m + (l * dim5_) 
                         + (k * dim5_ * dim4_) 
                         + (j * dim5_ * dim4_ * dim3_) 
                         + (i * dim5_ * dim4_ * dim3_ * dim2_));
}

template <typename T>
KOKKOS_FUNCTION
T& CArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l,
                               size_t m, size_t n) const {
    assert(i < dim1_ && "i is out of bounds in CArrayKokkos 6D!");
    assert(j < dim2_ && "j is out of bounds in CArrayKokkos 6D!");
    assert(k < dim3_ && "k is out of bounds in CArrayKokkos 6D!");
    assert(l < dim4_ && "l is out of bounds in CArrayKokkos 6D!");
    assert(m < dim5_ && "m is out of bounds in CArrayKokkos 6D!");
    assert(n < dim6_ && "n is out of bounds in CArrayKokkos 6D!");
    return this_array_(n + (m * dim6_) 
                         + (l * dim6_ * dim5_)  
                         + (k * dim6_ * dim5_ * dim4_) 
                         + (j * dim6_ * dim5_ * dim4_ * dim3_)  
                         + (i * dim6_ * dim5_ * dim4_ * dim3_ * dim2_));
}

template <typename T>
CArrayKokkos<T>& CArrayKokkos<T>::operator= (const CArrayKokkos& temp) {
    using TArray1D = Kokkos::View<T *,Layout,ExecSpace>;
    
    // Do nothing if the assignment is of the form x = x
    if (this != &temp) {
        dim1_ = temp.dim1_;
        dim2_ = temp.dim2_;
        dim3_ = temp.dim3_;
        dim4_ = temp.dim4_;
        dim5_ = temp.dim5_;
        dim6_ = temp.dim6_;
        length_ = temp.length_;
        this_array_ = TArray1D("this_array_", length_);
    }
    
    return *this;
}

// Return size
template <typename T>
KOKKOS_FUNCTION
size_t CArrayKokkos<T>::size() {
    return length_;
}

template <typename T>
size_t CArrayKokkos<T>::extent() {
    return length_;
}

template <typename T>
T* CArrayKokkos<T>::pointer() {
    return this_array_.data();
}

template <typename T>
KOKKOS_FUNCTION
CArrayKokkos<T>::~CArrayKokkos() {}

////////////////////////////////////////////////////////////////////////////////
// End of CArrayKokkos
////////////////////////////////////////////////////////////////////////////////

/*! \brief Kokkos version of the serial ViewCArray class.
 *
 */
template <typename T>
class ViewCArrayKokkos {

private:
    size_t dim1_;
    size_t dim2_;
    size_t dim3_;
    size_t dim4_;
    size_t dim5_;
    size_t dim6_;
    size_t length_;  // Length of 1D array
    T* this_array_;
    
public:
    ViewCArrayKokkos();

    ViewCArrayKokkos(T* some_array, size_t some_dim1);

    ViewCArrayKokkos(T* some_array, size_t some_dim1, size_t some_dim2);

    ViewCArrayKokkos(T* some_array, size_t some_dim1, size_t some_dim2,
                     size_t some_dim3);

    ViewCArrayKokkos(T* some_array, size_t some_dim1, size_t some_dim2,
                     size_t some_dim3, size_t some_dim4);

    ViewCArrayKokkos(T* some_array, size_t some_dim1, size_t some_dim2,
                     size_t some_dim3, size_t some_dim4, size_t some_dim5);

    ViewCArrayKokkos(T* some_array, size_t some_dim1, size_t some_dim2,
                     size_t some_dim3, size_t some_dim4, size_t some_dim5,
                     size_t some_dim6);
    
    KOKKOS_FUNCTION
    T& operator()(size_t i) const;
    
    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j) const;
    
    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;
        
    KOKKOS_FUNCTION
    T& operator() (size_t i, size_t j, size_t k, size_t l, size_t m) const;
        
    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, 
                  size_t n) const;
    
    KOKKOS_FUNCTION
    size_t size();

    size_t extent();

    KOKKOS_FUNCTION
    ~ViewCArrayKokkos();
    
}; // end of ViewCArrayKokkos

// Default constructor
template <typename T>
ViewCArrayKokkos<T>::ViewCArrayKokkos() {}

// Overloaded 1D constructor
template <typename T>
ViewCArrayKokkos<T>::ViewCArrayKokkos(T* some_array, size_t some_dim1) {
    dim1_ = some_dim1;
    length_ = dim1_;
    this_array_ = some_array;
}

// Overloaded 2D constructor
template <typename T>
ViewCArrayKokkos<T>::ViewCArrayKokkos(T* some_array, size_t some_dim1, 
                                      size_t some_dim2) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = (dim1_ * dim2_);
    this_array_ = some_array;
}

// Overloaded 3D constructor
template <typename T>
ViewCArrayKokkos<T>::ViewCArrayKokkos(T* some_array, size_t some_dim1,
                                      size_t some_dim2, size_t some_dim3) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = (dim1_ * dim2_ * dim3_);
    this_array_ = some_array;
}

// Overloaded 4D constructor
template <typename T>
ViewCArrayKokkos<T>::ViewCArrayKokkos(T* some_array, size_t some_dim1,
                                      size_t some_dim2, size_t some_dim3,
                                      size_t some_dim4) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_);
    this_array_ = some_array;
}

// Overloaded 5D constructor
template <typename T>
ViewCArrayKokkos<T>::ViewCArrayKokkos(T* some_array, size_t some_dim1,
                                      size_t some_dim2, size_t some_dim3,
                                      size_t some_dim4, size_t some_dim5) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_);
    this_array_ = some_array;
}

// Overloaded 6D constructor
template <typename T>
ViewCArrayKokkos<T>::ViewCArrayKokkos(T* some_array, size_t some_dim1,
                                      size_t some_dim2, size_t some_dim3,
                                      size_t some_dim4, size_t some_dim5,
                                      size_t some_dim6) {
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_);
    this_array_ = some_array;
}

template <typename T>
KOKKOS_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i) const {
    assert(i < dim1_ && "i is out of bounds in ViewCArrayKokkos 1D!");
    return this_array_[i];
}

template <typename T>
KOKKOS_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i, size_t j) const {
    assert(i < dim1_ && "i is out of bounds in ViewCArrayKokkos 2D!");
    assert(j < dim2_ && "j is out of bounds in ViewCArrayKokkos 2D!");  
    return this_array_[j + (i * dim2_)];
}

template <typename T>
KOKKOS_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i, size_t j, size_t k) const {
    assert(i < dim1_ && "i is out of bounds in ViewCArrayKokkos 3D!");
    assert(j < dim2_ && "j is out of bounds in ViewCArrayKokkos 3D!");
    assert(k < dim3_ && "k is out of bounds in ViewCArrayKokkos 3D!");
    return this_array_[k + (j * dim3_) 
                         + (i * dim3_ * dim2_)];
}

template <typename T>
KOKKOS_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, 
                                   size_t l) const {
    assert(i < dim1_ && "i is out of bounds in ViewCArrayKokkos 4D!");
    assert(j < dim2_ && "j is out of bounds in ViewCArrayKokkos 4D!");
    assert(k < dim3_ && "k is out of bounds in ViewCArrayKokkos 4D!");
    assert(l < dim4_ && "l is out of bounds in ViewCArrayKokkos 4D!");
    return this_array_[l + (k * dim4_) 
                         + (j * dim4_ * dim3_) 
                         + (i * dim4_ * dim3_ * dim2_)];
}

template <typename T>
KOKKOS_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, 
                                   size_t m) const {
    assert(i < dim1_ && "i is out of bounds in ViewCArrayKokkos 5D!");
    assert(j < dim2_ && "j is out of bounds in ViewCArrayKokkos 5D!");
    assert(k < dim3_ && "k is out of bounds in ViewCArrayKokkos 5D!");
    assert(l < dim4_ && "l is out of bounds in ViewCArrayKokkos 5D!");
    assert(m < dim5_ && "m is out of bounds in ViewCArrayKokkos 5D!");
    return this_array_[m + (l * dim5_) 
                         + (k * dim5_ * dim4_) 
                         + (j * dim5_ * dim4_ * dim3_)
                         + (i * dim5_ * dim4_ * dim3_ * dim2_)];
}

template <typename T>
KOKKOS_FUNCTION
T& ViewCArrayKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, 
                                   size_t m, size_t n) const {
    assert(i < dim1_ && "i is out of bounds in ViewCArrayKokkos 6D!");
    assert(j < dim2_ && "j is out of bounds in ViewCArrayKokkos 6D!");
    assert(k < dim3_ && "k is out of bounds in ViewCArrayKokkos 6D!");
    assert(l < dim4_ && "l is out of bounds in ViewCArrayKokkos 6D!");
    assert(m < dim5_ && "m is out of bounds in ViewCArrayKokkos 6D!");
    assert(n < dim6_ && "n is out of bounds in ViewCArrayKokkos 6D!");
    return this_array_[n + (m * dim6_) 
                         + (l * dim6_ * dim5_) 
                         + (k * dim6_ * dim5_ * dim4_)
                         + (j * dim6_ * dim5_ * dim4_ * dim3_) 
                         + (i * dim6_ * dim5_ * dim4_ * dim3_ * dim2_)];
}

template <typename T>
KOKKOS_FUNCTION
size_t ViewCArrayKokkos<T>::size() {
    return length_;
}

template <typename T>
size_t ViewCArrayKokkos<T>::extent() {
    return length_;
}

template <typename T>
KOKKOS_FUNCTION
ViewCArrayKokkos<T>::~ViewCArrayKokkos() {}

////////////////////////////////////////////////////////////////////////////////
// End of ViewCArrayKokkos
////////////////////////////////////////////////////////////////////////////////

/*! \brief Kokkos version of the serial CMatrix class.
 *
 */
template <typename T>
class CMatrixKokkos {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
private:
    size_t dim1_;
    size_t dim2_;
    size_t dim3_;
    size_t dim4_;
    size_t dim5_;
    size_t dim6_;
    size_t length_;
    TArray1D this_matrix_; 

public:
    CMatrixKokkos();

    CMatrixKokkos(size_t some_dim1);

    CMatrixKokkos(size_t some_dim1, size_t some_dim2);

    CMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3);    

    CMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, 
                  size_t some_dim4);

    CMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, 
                  size_t some_dim4, size_t some_dim5);

    CMatrixKokkos(size_t some_dim1, size_t some_dim2, size_t some_dim3, 
                  size_t some_dim4, size_t some_dim5, size_t some_dim6);

    KOKKOS_FUNCTION
    T& operator()(size_t i) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, 
                  size_t n) const;

    CMatrixKokkos& operator=(const CMatrixKokkos &temp);

    KOKKOS_FUNCTION
    size_t size();

    size_t extent();

    T* pointer();

    KOKKOS_FUNCTION
    ~CMatrixKokkos();

}; // End of CMatrixKokkos

// Default constructor
template <typename T>
CMatrixKokkos<T>::CMatrixKokkos() {}

// Overloaded 1D constructor
template <typename T>
CMatrixKokkos<T>::CMatrixKokkos(size_t some_dim1) { 
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    length_ = dim1_;
    this_matrix_ = TArray1D("this_matrix_", length_);
}

// Overloaded 2D constructor
template <typename T>
CMatrixKokkos<T>::CMatrixKokkos(size_t some_dim1, size_t some_dim2) { 
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    length_ = (dim1_ * dim2_);
    this_matrix_ = TArray1D("this_matrix_", length_);
}

// Overloaded 3D constructor
template <typename T>
CMatrixKokkos<T>::CMatrixKokkos(size_t some_dim1, size_t some_dim2, 
                                size_t some_dim3) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    length_ = (dim1_ * dim2_ * dim3_);
    this_matrix_ = TArray1D("this_matrix_", length_);
}

// Overloaded 4D constructor
template <typename T>
CMatrixKokkos<T>::CMatrixKokkos(size_t some_dim1, size_t some_dim2, 
                                size_t some_dim3, size_t some_dim4) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_);
    this_matrix_ = TArray1D("this_matrix_", length_);
}

// Overloaded 5D constructor
template <typename T>
CMatrixKokkos<T>::CMatrixKokkos(size_t some_dim1, size_t some_dim2, 
                                size_t some_dim3, size_t some_dim4, 
                                size_t some_dim5) {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_);
    this_matrix_ = TArray1D("this_matrix_", length_);
}

// Overloaded 6D constructor
template <typename T>
CMatrixKokkos<T>::CMatrixKokkos(size_t some_dim1, size_t some_dim2, 
                                size_t some_dim3, size_t some_dim4, 
                                size_t some_dim5, size_t some_dim6) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
    dim1_ = some_dim1;
    dim2_ = some_dim2;
    dim3_ = some_dim3;
    dim4_ = some_dim4;
    dim5_ = some_dim5;
    dim6_ = some_dim6;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_);
    this_matrix_ = TArray1D("this_matrix_", length_);
}

template<typename T>
KOKKOS_FUNCTION
T& CMatrixKokkos<T>::operator()(size_t i) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in CMatrixKokkos 1D!");
    return this_matrix_((i - 1));
}

template <typename T>
KOKKOS_FUNCTION
T& CMatrixKokkos<T>::operator()(size_t i, size_t j) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in CMatrixKokkos 2D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in CMatrixKokkos 2D!");
    return this_matrix_((j - 1) + ((i - 1) * dim2_));
}

template<typename T>
KOKKOS_FUNCTION
T& CMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in CMatrixKokkos 3D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in CMatrixKokkos 3D!");
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in CMatrixKokkos 3D!");
    return this_matrix_((k - 1) + ((j - 1) * dim3_) 
                                + ((i - 1) * dim3_ * dim2_));
}

template <typename T>
KOKKOS_FUNCTION
T& CMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in CMatrixKokkos 4D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in CMatrixKokkos 4D!");
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in CMatrixKokkos 4D!");
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in CMatrixKokkos 4D!");
    return this_matrix_((l - 1) + ((k - 1) * dim4_) 
                                + ((j - 1) * dim4_ * dim3_) 
                                + ((i - 1) * dim4_ * dim3_ * dim2_));
}

template <typename T>
KOKKOS_FUNCTION
T& CMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, 
                                size_t m) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in CMatrixKokkos 5D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in CMatrixKokkos 5D!");
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in CMatrixKokkos 5D!");
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in CMatrixKokkos 5D!");
    assert(m >= 1 && m <= dim5_ && "m is out of bounds in CMatrixKokkos 5D!");
    return this_matrix_((m - 1) + ((l - 1) * dim5_) 
                                + ((k - 1) * dim5_ * dim4_) 
                                + ((j - 1) * dim5_ * dim4_ * dim3_) 
                                + ((i - 1) * dim5_ * dim4_ * dim3_ * dim2_));
}

template <typename T>
KOKKOS_FUNCTION
T& CMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, 
                                size_t m, size_t n) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in CMatrixKokkos 6D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in CMatrixKokkos 6D!");
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in CMatrixKokkos 6D!");
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in CMatrixKokkos 6D!");
    assert(m >= 1 && m <= dim5_ && "m is out of bounds in CMatrixKokkos 6D!");
    assert(n >= 1 && n <= dim6_ && "n is out of bounds in CMatrixKokkos 6D!");
    return this_matrix_((n - 1) + ((m - 1) * dim6_) 
                                + ((l - 1) * dim6_ * dim5_) 
                                + ((k - 1) * dim6_ * dim5_ * dim4_) 
                                + ((j - 1) * dim6_ * dim5_ * dim4_ * dim3_) 
                                + ((i - 1) * dim6_ * dim5_ * dim4_ * dim3_ * dim2_));
}

// Overload = operator
// for object assignment THIS = CMatrixKokkos <> temp
template <typename T>
CMatrixKokkos<T> & CMatrixKokkos<T>::operator=(const CMatrixKokkos &temp) {
    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;

    if( this != &temp) {
        dim1_ = temp.dim1_;
        dim2_ = temp.dim2_;
        dim3_ = temp.dim3_;
        dim4_ = temp.dim4_;
        dim5_ = temp.dim5_;
        dim6_ = temp.dim6_;
        length_ = temp.length_;
        this_matrix_ = TArray1D("this_matrix_", length_);
    }
    
    return *this;
}

template <typename T>
KOKKOS_FUNCTION
size_t CMatrixKokkos<T>::size() {
    return length_;
}

template <typename T>
size_t CMatrixKokkos<T>::extent() {
    return length_;
}

template <typename T>
T* CMatrixKokkos<T>::pointer() {
    return this_matrix_.data();
}

// Deconstructor
template <typename T>
KOKKOS_FUNCTION
CMatrixKokkos<T>::~CMatrixKokkos() {}

////////////////////////////////////////////////////////////////////////////////
// End of CMatrixKokkos
////////////////////////////////////////////////////////////////////////////////

/*! \brief Kokkos version of the serial ViewCMatrix class.
 *
 */
template <typename T>
class ViewCMatrixKokkos {

private:
    size_t dim1_;
    size_t dim2_;
    size_t dim3_;
    size_t dim4_;
    size_t dim5_;
    size_t dim6_;
    size_t length_;
    T* this_matrix_;

public:
    ViewCMatrixKokkos();

    ViewCMatrixKokkos(T* some_matrix, size_t dim1);

    ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2);

    ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2, size_t dim3);

    ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2, size_t dim3, 
                      size_t dim4);

    ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2, size_t dim3, 
                      size_t dim4, size_t dim5);

    ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2, size_t dim3,
                      size_t dim4, size_t dim5, size_t dim6);

    KOKKOS_FUNCTION
    T& operator()(size_t i) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j , size_t k) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k , size_t l) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;

    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n) const;

    KOKKOS_FUNCTION
    size_t size();

    size_t extent();

    KOKKOS_FUNCTION
    ~ViewCMatrixKokkos();

}; // End of ViewCMatrixKokkos

// Default constructor
template <typename T>
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(){ }

// Overloaded 1D constructor
template <typename T>
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T* some_matrix, size_t dim1) {
    dim1_ = dim1;
    length_ = dim1_;
    this_matrix_ = some_matrix;
}

// Overloaded 2D constructor
template <typename T>
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T* some_matrix, size_t dim1, 
                                        size_t dim2) {
    dim1_ = dim1;
    dim2_ = dim2;
    length_ = (dim1_ * dim2_);
    this_matrix_ = some_matrix;
}

// Overloaded 3D constructor
template <typename T>
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2,
                                        size_t dim3) {
    dim1_ = dim1;
    dim2_ = dim2;
    dim3_ = dim3;
    length_ = (dim1_ * dim2_ * dim3_);
    this_matrix_ = some_matrix;
}

// Overloaded 4D constructor
template <typename T>
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2,
                                        size_t dim3, size_t dim4) {
    dim1_ = dim1;
    dim2_ = dim2;
    dim3_ = dim3;
    dim4_ = dim4;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_);
    this_matrix_ = some_matrix;
}

// Overloaded 5D constructor
template <typename T>
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2,
                                        size_t dim3, size_t dim4, size_t dim5) {
    dim1_ = dim1;
    dim2_ = dim2;
    dim3_ = dim3;
    dim4_ = dim4;
    dim5_ = dim5;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_);
    this_matrix_ = some_matrix;
}

// Overloaded 6D constructor
template <typename T>
ViewCMatrixKokkos<T>::ViewCMatrixKokkos(T* some_matrix, size_t dim1, size_t dim2,
                                        size_t dim3, size_t dim4, size_t dim5,
                                        size_t dim6) {
    dim1_ = dim1;
    dim2_ = dim2;
    dim3_ = dim3;
    dim4_ = dim4;
    dim5_ = dim5;
    dim6_ = dim6;
    length_ = (dim1_ * dim2_ * dim3_ * dim4_ * dim5_ * dim6_);
    this_matrix_ = some_matrix;
}

template <typename T>
KOKKOS_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in ViewCMatrixKokkos 1D!");
    return this_matrix_[(i - 1)];
}

template <typename T>
KOKKOS_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i, size_t j) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in ViewCMatrixKokkos 2D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in ViewCMatrixKokkos 2D!");
    return this_matrix_[(j - 1) + ((i - 1) * dim2_)];
}

template <typename T>
KOKKOS_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in ViewCMatrixKokkos 3D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in ViewCMatrixKokkos 3D!");
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in ViewCMatrixKokkos 3D!");
    return this_matrix_[(k - 1) + ((j - 1) * dim3_) 
                                + ((i - 1) * dim3_ * dim2_)];
}

template <typename T>
KOKKOS_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i, size_t j , size_t k, size_t l) const { 
    assert(i >= 1 && i <= dim1_ && "i is out of bounds in ViewCMatrixKokkos 4D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds in ViewCMatrixKokkos 4D!");
    assert(k >= 1 && k <= dim3_ && "k is out of bounds in ViewCMatrixKokkos 4D!");
    assert(l >= 1 && l <= dim4_ && "l is out of bounds in ViewCMatrixKokkos 4D!");
    return this_matrix_[(l - 1) + ((k - 1) * dim4_) 
                                + ((j - 1) * dim4_ * dim3_) 
                                + ((i - 1) * dim4_ * dim3_ * dim2_)];
}

template <typename T>
KOKKOS_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, 
                                    size_t m) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds for ViewCMatrixKokkos 5D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds for ViewCMatrixKokkos 5D!");
    assert(k >= 1 && k <= dim3_ && "k is out of bounds for ViewCMatrixKokkos 5D!");
    assert(l >= 1 && l <= dim4_ && "l is out of bounds for ViewCMatrixKokkos 5D!");
    assert(m >= 1 && m <= dim5_ && "m is out of bounds for ViewCMatrixKokkos 5D!");
    return this_matrix_[(m - 1) + ((l - 1) * dim5_)
                                + ((k - 1) * dim5_ * dim4_)
                                + ((j - 1) * dim5_ * dim4_ * dim3_)
                                + ((i - 1) * dim5_ * dim4_ * dim3_ * dim2_)];
}

template <typename T>
KOKKOS_FUNCTION
T& ViewCMatrixKokkos<T>::operator()(size_t i, size_t j, size_t k, size_t l, 
                                    size_t m, size_t n) const {
    assert(i >= 1 && i <= dim1_ && "i is out of bounds for ViewCMatrixKokkos 6D!");
    assert(j >= 1 && j <= dim2_ && "j is out of bounds for ViewCMatrixKokkos 6D!");
    assert(k >= 1 && k <= dim3_ && "k is out of bounds for ViewCMatrixKokkos 6D!");
    assert(l >= 1 && l <= dim4_ && "l is out of bounds for ViewCMatrixKokkos 6D!");
    assert(m >= 1 && m <= dim5_ && "m is out of bounds for ViewCMatrixKokkos 6D!");
    assert(n >= 1 && n <= dim6_ && "n is out of bounds for ViewCMatrixKokkos 6D!");
    return this_matrix_[(n - 1) + ((m - 1) * dim6_)
                                + ((l - 1) * dim6_ * dim5_)
                                + ((k - 1) * dim6_ * dim5_ * dim4_)
                                + ((j - 1) * dim6_ * dim5_ * dim4_ * dim3_)
                                + ((i - 1) * dim6_ * dim5_ * dim4_ * dim3_ * dim2_)];
}

template <typename T>
KOKKOS_FUNCTION
size_t ViewCMatrixKokkos<T>::size() {
    return length_;
}

template <typename T>
size_t ViewCMatrixKokkos<T>::extent() {
    return length_;
}

template <typename T>
KOKKOS_FUNCTION
ViewCMatrixKokkos<T>::~ViewCMatrixKokkos() {}

////////////////////////////////////////////////////////////////////////////////
// End of ViewCMatrixKokkos
////////////////////////////////////////////////////////////////////////////////

/*! \brief Kokkos version of the serial RaggedRightArray class.
 *
 */
template <typename T>
class RaggedRightArrayKokkos {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
private:
    SArray1D start_index_;
    TArray1D array_; 
    
    size_t dim1_;
    size_t length_;

    // THIS WILL BE A GPU POINTER!
    size_t* mystrides_;
    
public:
    // Default constructor
    RaggedRightArrayKokkos();
    
    //--- 2D array access of a ragged right array ---
    
    // Overload constructor for a CArray
    RaggedRightArrayKokkos(CArrayKokkos<size_t> &strides_array);
    
    // Overload constructor for a ViewCArray
    RaggedRightArrayKokkos(ViewCArray<size_t> &strides_array);
    
    // Overloaded constructor for a traditional array
    RaggedRightArrayKokkos(size_t* strides_array, size_t some_dim1);

    // A method to return the stride size
    KOKKOS_FUNCTION
    size_t stride(size_t i) const;
    
    // Overload operator() to access data as array(i,j)
    // where i=[0:N-1], j=[stride(i)]
    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j) const;

    T& pointer();

    RaggedRightArrayKokkos& operator= (const RaggedRightArrayKokkos &temp);

    // Destructor
    KOKKOS_FUNCTION
    ~RaggedRightArrayKokkos ( );
}; // End of RaggedRightArray

template <typename T>
RaggedRightArrayKokkos<T>::RaggedRightArrayKokkos() {}

// Overloaded constructor
template <typename T>
RaggedRightArrayKokkos<T>::RaggedRightArrayKokkos(CArrayKokkos<size_t> &strides_array) {
    mystrides_ = strides_array.pointer();
    dim1_ = strides_array.extent();
} // End constructor

// Overloaded constructor
template <typename T>
RaggedRightArrayKokkos<T>::RaggedRightArrayKokkos(ViewCArray<size_t> &strides_array) {
} // End constructor

// Overloaded constructor
template <typename T>
RaggedRightArrayKokkos<T>::RaggedRightArrayKokkos(size_t* strides_array, 
                                                  size_t some_dim1) {
    mystrides_ = strides_array;
    dim1_ = some_dim1;
} // End constructor

// A method to return the stride size
template <typename T>
KOKKOS_FUNCTION
size_t RaggedRightArrayKokkos<T>::stride(size_t i) const {
    // Ensure that i is within bounds
    assert(i < (dim1_ + 1) && "i is greater than dim1_ in RaggedRightArray");

    return start_index_((i + 1)) - start_index_(i);
}

// Overload operator() to access data as array(i,j)
// where i=[0:N-1], j=[0:stride(i)]
template <typename T>
KOKKOS_FUNCTION
T& RaggedRightArrayKokkos<T>::operator()(size_t i, size_t j) const {
    // Get the 1D array index
    size_t start = start_index_(i);
    
    // asserts
    assert(i < dim1_ && "i is out of dim1 bounds in RaggedRightArrayKokkos");  // die if >= dim1
    assert(j < stride(i) && "j is out of stride bounds in RaggedRightArrayKokkos");  // die if >= stride
    
    return array_(j + start);
} // End operator()

template <typename T>
RaggedRightArrayKokkos<T> & RaggedRightArrayKokkos<T>::operator= (const RaggedRightArrayKokkos &temp) {

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
    start_index_ = SArray1D("start_index_", dim1_ + 1);
    //start_index_(0) = 0; // the 1D array starts at 0
    Kokkos::parallel_for("StartFirst", 1, KOKKOS_CLASS_LAMBDA(const int&) {
            start_index_(0) = 0;
        });
    Kokkos::fence();
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    Kokkos::parallel_scan("StartValues", dim1_, KOKKOS_CLASS_LAMBDA(const int i, double& update, const bool final) {
            // Load old value in case we update it before accumulating
            const size_t count = temp.mystrides_[i];
            update += count;
            if (final) {
                start_index_((i+1)) = update;
            }       

        });
    Kokkos::fence();

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

    SArray1D templen = SArray1D("templen", 1);
    auto h_templen = Kokkos::create_mirror_view(templen);
    Kokkos::parallel_for("ArrayLength", 1, KOKKOS_CLASS_LAMBDA(const int&) {
            templen(0) = start_index_(dim1_);
            //length_ = start_index_(dim1_);
        });
    Kokkos::fence();
    deep_copy(h_templen, templen);
    length_ = h_templen(0);

    //printf("Length %ld\n", length_);

    //Kokkos::parallel_for("StartCheck", dim1_+1, KOKKOS_CLASS_LAMBDA(const int i) {
    //        printf("%d) Start %ld\n", i, start_index_(i));
    //    });
    //Kokkos::fence();
    
    array_ = TArray1D("array_", length_);

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

// Destructor
template <typename T>
KOKKOS_FUNCTION
RaggedRightArrayKokkos<T>::~RaggedRightArrayKokkos() { }

////////////////////////////////////////////////////////////////////////////////
// End of RaggedRightArrayKokkos
////////////////////////////////////////////////////////////////////////////////

/*! \brief Kokkos version of the serial RaggedDownArray class.
 *
 */
template <typename T>
class RaggedDownArrayKokkos {

    using TArray1D = Kokkos::View<T*, Layout, ExecSpace>;
    
private:
    SArray1D start_index_;
    TArray1D array_; 
    
    size_t dim2_;
    size_t length_;

    // THIS WILL BE A GPU POINTER!
    size_t* mystrides_;
    
public:
    // Default constructor
    RaggedDownArrayKokkos();
    
    //--- 2D array access of a ragged right array ---
    
    // Overload constructor for a CArray
    RaggedDownArrayKokkos(CArrayKokkos<size_t> &strides_array);
    
    // Overload constructor for a ViewCArray
    RaggedDownArrayKokkos(ViewCArray<size_t> &strides_array);
    
    // Overloaded constructor for a traditional array
    RaggedDownArrayKokkos(size_t* strides_array, size_t some_dim2);

    // A method to return the stride size
    KOKKOS_FUNCTION
    size_t stride(size_t j) const;
    
    // Overload operator() to access data as array(i,j)
    // where i=[0:N-1], j=[stride(i)]
    KOKKOS_FUNCTION
    T& operator()(size_t i, size_t j) const;

    T& pointer();

    RaggedDownArrayKokkos& operator= (const RaggedDownArrayKokkos &temp);

    // Destructor
    KOKKOS_FUNCTION
    ~RaggedDownArrayKokkos ( );
}; // End of RaggedDownArray

template <typename T>
RaggedDownArrayKokkos<T>::RaggedDownArrayKokkos() {}

// Overloaded constructor
template <typename T>
RaggedDownArrayKokkos<T>::RaggedDownArrayKokkos(CArrayKokkos<size_t> &strides_array) {
    mystrides_ = strides_array.pointer();
    dim2_ = strides_array.extent();
} // End constructor

// Overloaded constructor
template <typename T>
RaggedDownArrayKokkos<T>::RaggedDownArrayKokkos(ViewCArray<size_t> &strides_array) {
} // End constructor

// Overloaded constructor
template <typename T>
RaggedDownArrayKokkos<T>::RaggedDownArrayKokkos(size_t* strides_array, 
                                                  size_t some_dim2) {
    mystrides_ = strides_array;
    dim2_ = some_dim2;
} // End constructor

// A method to return the stride size
template <typename T>
KOKKOS_FUNCTION
size_t RaggedDownArrayKokkos<T>::stride(size_t j) const {
    // Ensure that j is within bounds
    assert(j < (dim2_ + 1) && "j is greater than dim1_ in RaggedDownArray");

    return start_index_((j + 1)) - start_index_(j);
}

// Overload operator() to access data as array(i,j)
// where i=[0:N-1], j=[0:stride(i)]
template <typename T>
KOKKOS_FUNCTION
T& RaggedDownArrayKokkos<T>::operator()(size_t i, size_t j) const {
    // Get the 1D array index
    size_t start = start_index_(j);
    
    // asserts
    assert(i < stride(j) && "i is out of stride bounds in RaggedDownArrayKokkos");  // die if >= stride
    assert(j < dim2_ && "j is out of dim1 bounds in RaggedDownArrayKokkos");  // die if >= dim1
    
    return array_(i + start);
} // End operator()

template <typename T>
RaggedDownArrayKokkos<T> & RaggedDownArrayKokkos<T>::operator= (const RaggedDownArrayKokkos &temp) {

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
    start_index_ = SArray1D("start_index_", dim2_ + 1);
    //start_index_(0) = 0; // the 1D array starts at 0
    Kokkos::parallel_for("StartFirst", 1, KOKKOS_CLASS_LAMBDA(const int&) {
            start_index_(0) = 0;
        });
    Kokkos::fence();
    
    // Loop over to find the total length of the 1D array to
    // represent the ragged-right array and set the starting 1D index
    Kokkos::parallel_scan("StartValues", dim2_, KOKKOS_CLASS_LAMBDA(const int j, double& update, const bool final) {
            // Load old value in case we update it before accumulating
            const size_t count = temp.mystrides_[j];
            update += count;
            if (final) {
                start_index_((j+1)) = update;
            }       

        });
    Kokkos::fence();

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

    SArray1D templen = SArray1D("templen", 1);
    auto h_templen = Kokkos::create_mirror_view(templen);
    Kokkos::parallel_for("ArrayLength", 1, KOKKOS_CLASS_LAMBDA(const int&) {
            templen(0) = start_index_(dim2_);
            //length_ = start_index_(dim2_);
        });
    Kokkos::fence();
    deep_copy(h_templen, templen);
    length_ = h_templen(0);

    printf("Length %ld\n", length_);

    Kokkos::parallel_for("StartCheck", dim2_+1, KOKKOS_CLASS_LAMBDA(const int j) {
            printf("%d) Start %ld\n", j, start_index_(j));
        });
    Kokkos::fence();
    
    array_ = TArray1D("array_", length_);

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

// Destructor
template <typename T>
KOKKOS_FUNCTION
RaggedDownArrayKokkos<T>::~RaggedDownArrayKokkos() { }

////////////////////////////////////////////////////////////////////////////////
// End of RaggedDownArrayKokkos
////////////////////////////////////////////////////////////////////////////////

//////////////////////////
// Inherited Class Array
//////////////////////////

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

    KOKKOS_FUNCTION
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
    KOKKOS_FUNCTION
    size_t size();

    // Host Method
    // Method that returns size
    size_t extent();

    // Methods returns the raw pointer (most likely GPU) of the Kokkos View
    T* pointer();

    // Deconstructor
    KOKKOS_FUNCTION
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
KOKKOS_FUNCTION
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
KOKKOS_FUNCTION
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
KOKKOS_FUNCTION
InheritedArray2L<T>::~InheritedArray2L() {}

////////////////////////////////////////////////////////////////////////////////
// End of InheritedArray2L
////////////////////////////////////////////////////////////////////////////////

#endif







#endif // MATAR_H
