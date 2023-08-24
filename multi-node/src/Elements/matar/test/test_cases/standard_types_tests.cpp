#include "matar.h"
#include "gtest/gtest.h"
#include <stdio.h>

using namespace mtr; // matar namespace

template <typename T>
T return_dense_type()
{
  size_t size = 5;

  T dense_type(size);

  for (size_t i = 0; i < size; i++) {
    dense_type.pointer()[i] = i;
  }

  return dense_type;
}

template <typename T>
void modify_dense_type(T &dense_type)
{
  for (size_t i = 0; i < dense_type.size(); i++) {
    dense_type.pointer()[i] += 1;
  }
}


template <typename T>
T return_ragged_type(size_t *strides, size_t dim)
{
  // if T is RaggedRightArray
  // create a ragged-right array of integers
  //
  // [0, 1, 2]
  // [3, 4]
  // [5]
  // [6, 7, 8, 9]
 
  // if T is RaggedDownArray
  // create a ragged-down array of integers
  // _ _ _ _
  // 0 3 5 6
  // 1 4   7
  // 2     8
  //       9
 
  T ragged_type(strides, dim);

  size_t size;
  for (size_t i = 0; i < dim; i++) {
    size += strides[i];
  }

  for (size_t i = 0; i < size; i++) {
    ragged_type.pointer()[i] = i;
  }

  return ragged_type;

}


template <typename T>
void modify_ragged_type(T &ragged_type)
{
  for (size_t i = 0; i < ragged_type.size(); i++) {
    ragged_type.pointer()[i] += 1;
  }
}

template <typename T>
void empty_function(T x_type)
{}


TEST(StandaredTypesTests, FunctionReturnDenseTypes)
{
  // FArray
  FArray<int> farray = return_dense_type <FArray<int>>();
  for (size_t i = 0; i < farray.size(); i++) {
    EXPECT_EQ(i, farray(i));
  }

  // FMatrix
  FMatrix<int> fmatrix = return_dense_type <FMatrix<int>>();
  for (size_t i = 1; i <= fmatrix.size(); i++) {
    EXPECT_EQ(i-1, fmatrix(i));
  }

  // CArray
  CArray<int> carray = return_dense_type <CArray<int>>();
  for (size_t i = 0; i < carray.size(); i++) {
    EXPECT_EQ(i, carray(i));
  }

  // CMatrix
  CMatrix<int> cmatrix = return_dense_type <CMatrix<int>>();
  for (size_t i = 1; i <= cmatrix.size(); i++) {
    EXPECT_EQ(i-1, cmatrix(i));
  }
}


TEST(StandaredTypesTests, FunctionReturnRaggedTypes)
{
  const size_t dim = 4;
  size_t strides[dim] = {3,2,1,4};

  // RaggedRightArray 
  RaggedRightArray <int> ragged_right = return_ragged_type <RaggedRightArray<int>> (strides, dim);
  int value = 0;
  for (size_t i = 0; i < dim; i++) {
    for (size_t j = 0; j < ragged_right.stride(i); j++) {
      EXPECT_EQ(value, ragged_right(i,j));
      value++;
    }
  }

  // RaggedDownArray
  RaggedDownArray <int> ragged_down = return_ragged_type <RaggedDownArray<int>> (strides, dim);
  value = 0;
  for (size_t i = 0; i < dim; i++) {
    for (size_t j = 0; j < ragged_down.stride(i); j++) {
      EXPECT_EQ(value, ragged_down(j,i));
      value++;
    }
  }
}


TEST(StandaredTypesTests, PassDenseTypeToEmptyFunctionByValue)
{
  // FArray
  FArray<int> farray = return_dense_type <FArray<int>>();
  empty_function(farray);
  for (size_t i = 0; i < farray.size(); i++) {
    EXPECT_EQ(i, farray(i));
  }

  // FMatrix
  FMatrix<int> fmatrix = return_dense_type <FMatrix<int>>();
  empty_function(fmatrix);
  for (size_t i = 1; i <= fmatrix.size(); i++) {
    EXPECT_EQ(i-1, fmatrix(i));
  }

  // CArray
  CArray<int> carray = return_dense_type <CArray<int>>();
  empty_function(carray);
  for (size_t i = 0; i < carray.size(); i++) {
    EXPECT_EQ(i, carray(i));
  }

  // CMatrix
  CMatrix<int> cmatrix = return_dense_type <CMatrix<int>>();
  empty_function(cmatrix);
  for (size_t i = 1; i <= cmatrix.size(); i++) {
    EXPECT_EQ(i-1, cmatrix(i));
  }
}


TEST(StandaredTypesTests, PassRaggedTypeToEmptyFunctionByValue)
{
  const size_t dim = 4;
  size_t strides[dim] = {3,2,1,4};
  
  // RaggedRightArray
  RaggedRightArray <int> ragged_right = return_ragged_type <RaggedRightArray<int>> (strides, dim);
  empty_function(ragged_right);
  int value = 0;
  for (size_t i = 0; i < dim; i++) {
    for (size_t j = 0; j < ragged_right.stride(i); j++) {
      EXPECT_EQ(value, ragged_right(i,j));
      value++;
    }
  }
  
  // RaggedDownArray
  RaggedDownArray <int> ragged_down = return_ragged_type <RaggedDownArray<int>> (strides, dim);
  empty_function(ragged_down);
  value = 0;
  for (size_t i = 0; i < dim; i++) {
    for (size_t j = 0; j < ragged_down.stride(i); j++) {
      EXPECT_EQ(value, ragged_down(j,i));
      value++;
    }
  }
}


TEST(StandaredTypesTests, FunctionModifyDenseTypes)
{
  // FArray
  FArray<int> farray = return_dense_type <FArray<int>>();
  modify_dense_type(farray); // add 1 to all element
  for (size_t i = 0; i < farray.size(); i++) {
    EXPECT_EQ(i+1, farray(i));
  }

  // FMatrix
  FMatrix<int> fmatrix = return_dense_type <FMatrix<int>>();
  modify_dense_type(fmatrix); // add 1 to all element
  for (size_t i = 1; i <= fmatrix.size(); i++) {
    EXPECT_EQ(i, fmatrix(i));
  }

  // CArray
  CArray<int> carray = return_dense_type <CArray<int>>();
  modify_dense_type(carray); // add 1 to all element
  for (size_t i = 0; i < carray.size(); i++) {
    EXPECT_EQ(i+1, carray(i));
  }

  // CMatrix
  CMatrix<int> cmatrix = return_dense_type <CMatrix<int>>();
  modify_dense_type(cmatrix); // add 1 to all element
  for (size_t i = 1; i <= cmatrix.size(); i++) {
    EXPECT_EQ(i, cmatrix(i));
  }
}


TEST(StandaredTypesTests, FunctionModifyRaggedTypes)
{
  const size_t dim = 4;
  size_t strides[dim] = {3,2,1,4};

  // RaggedRightArray 
  RaggedRightArray <int> ragged_right = return_ragged_type <RaggedRightArray<int>> (strides, dim);
  modify_ragged_type(ragged_right); // add 1 to all element
  int value = 1;
  for (size_t i = 0; i < dim; i++) {
    for (size_t j = 0; j < ragged_right.stride(i); j++) {
      EXPECT_EQ(value, ragged_right(i,j));
      value++;
    }
  }

  // RaggedDownArray
  RaggedDownArray <int> ragged_down = return_ragged_type <RaggedDownArray<int>> (strides, dim);
  modify_ragged_type(ragged_down); // add 1 to all element
  value = 1;
  for (size_t i = 0; i < dim; i++) {
    for (size_t j = 0; j < ragged_down.stride(i); j++) {
      EXPECT_EQ(value, ragged_down(j,i));
      value++;
    }
  }
}

int main(int argc, char* argv[])
{

  int result = 0;

  testing::InitGoogleTest(&argc, argv);

  result = RUN_ALL_TESTS();

  return result;

}
