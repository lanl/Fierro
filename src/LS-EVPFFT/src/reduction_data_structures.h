#pragma once

#include "matar.h"



//========================================================================================
// ArrayType : used for performing reduction into an array in parallel region.
template <typename T, int N>
struct ArrayType {
  int size=N;
  T array[N];

  KOKKOS_INLINE_FUNCTION    // Default constructor
  ArrayType();

  KOKKOS_INLINE_FUNCTION    // Copy constructor
  ArrayType(const ArrayType & rhs);

  KOKKOS_INLINE_FUNCTION    // add operator
  ArrayType& operator+=(const ArrayType & rhs);

  KOKKOS_INLINE_FUNCTION    // volatile add operator
  void operator+=(const volatile ArrayType & rhs) volatile;

  //KOKKOS_INLINE_FUNCTION    // accessor
  //T& operator[](int i) const;

};

template <typename T, int N>
KOKKOS_INLINE_FUNCTION
ArrayType<T,N>::ArrayType()
{
  for (int i = 0; i < N; i++) {
    this->array[i] = T(0);
  }
}

template <typename T, int N>
KOKKOS_INLINE_FUNCTION
ArrayType<T,N>::ArrayType(const ArrayType & rhs)
{
  for (int i = 0; i < N; i++) {
    this->array[i] = rhs.array[i];
  }
}

template <typename T, int N>
KOKKOS_INLINE_FUNCTION
ArrayType<T,N>& ArrayType<T,N>::operator+=(const ArrayType & rhs)
{
  for (int i = 0; i < N; i++) {
    this->array[i] += rhs.array[i];
  }
  return *this;
}

template <typename T, int N>
KOKKOS_INLINE_FUNCTION
void ArrayType<T,N>::operator+=(const volatile ArrayType & rhs) volatile
{
  for (int i = 0; i < N; i++) {
    this->array[i] += rhs.array[i];
  }
}

//template <typename T, int N>
//KOKKOS_INLINE_FUNCTION
//T& ArrayType<T,N>::operator[](int i) const
//{
//  return this->array[i];
//}


//========================================================================================
// ArrayOfArrayType : used for performing reduction into an array of ArrayType in parallel region.
template <int M, typename T, int N>
struct ArrayOfArrayType {
  int size=M;
  ArrayType <T,N> array[M];

  KOKKOS_INLINE_FUNCTION    // Default constructor
  ArrayOfArrayType();

  KOKKOS_INLINE_FUNCTION    // Copy constructor
  ArrayOfArrayType(const ArrayOfArrayType & rhs);

  KOKKOS_INLINE_FUNCTION    // add operator
  ArrayOfArrayType& operator+=(const ArrayOfArrayType & rhs);

  KOKKOS_INLINE_FUNCTION    // volatile add operator
  void operator+=(const volatile ArrayOfArrayType & rhs) volatile;

  //KOKKOS_INLINE_FUNCTION    // accessor
  //ArrayType<T,N>& operator[](int i) const;

};

template <int M, typename T, int N>
KOKKOS_INLINE_FUNCTION
ArrayOfArrayType<M,T,N>::ArrayOfArrayType()
{}

template <int M, typename T, int N>
KOKKOS_INLINE_FUNCTION
ArrayOfArrayType<M,T,N>::ArrayOfArrayType(const ArrayOfArrayType & rhs)
{
  for (int j = 0; j < M; j++) {
    for (int i = 0; i < N; i++) {
      this->array[j].array[i] = rhs.array[j].array[i];
    }
  }
}

template <int M, typename T, int N>
KOKKOS_INLINE_FUNCTION
ArrayOfArrayType<M,T,N>& ArrayOfArrayType<M,T,N>::operator+=(const ArrayOfArrayType & rhs)
{
  for (int j = 0; j < M; j++) {
    this->array[j] += rhs.array[j];
  }
  return *this;
}

template <int M, typename T, int N>
KOKKOS_INLINE_FUNCTION
void ArrayOfArrayType<M,T,N>::operator+=(const volatile ArrayOfArrayType & rhs) volatile
{
  for (int j = 0; j < M; j++) {
    this->array[j] += rhs.array[j];
  }
}

//template <int M, typename T, int N>
//KOKKOS_INLINE_FUNCTION
//ArrayType<T,N>& ArrayOfArrayType<M,T,N>::operator[](int i) const
//{
//  return this->array[i];
//}




