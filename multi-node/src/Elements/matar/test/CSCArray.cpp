#include <stdio.h>
#include <math.h>
#include <matar.h>
#include <limits.h>
#include <gtest/gtest.h>

using namespace mtr; // matar namespace

void buildTestArray(CArray<int> &data, CArray<size_t>& cols, CArray<size_t>& rows){
    int i = 0;
    for(i = 0; i < 8; i++){
        data(i) = i+1;
    }
    rows(0) = 0;
    rows(1) = 0;
    rows(2) = 1;
    rows(3) = 2;
    rows(4) = 1;
    rows(5) = 2;
    rows(6) = 2;
    rows(7) = 3;
    cols(0) = 0;
    cols(1) = 1;
    cols(2) = 3;
    cols(3) = 4;
    cols(4) = 6;
    cols(5) = 7;
    cols(6) = 8;
}


// The array here looks like 
//  | 1 4 7 | 
//  | 2 5 8 |
//  | 3 6 0 | 
TEST(CSCArray, AccessingData){
    CArray<int> data(8);
    CArray<size_t> cols(4);
    CArray<size_t> rows(8);
    int i, j;
    for(i = 0; i < 8; i++){
        data(i) = i+1;
    }
    cols(0) = 0;
    cols(1) = 3;
    cols(2) = 6;
    cols(3) = 8;
    rows(0) = 0;
    rows(1) = 1;
    rows(2) = 2;
    rows(3) = 0;
    rows(4) = 1;
    rows(5) = 2;
    rows(6) = 0;
    rows(7) = 1;

    CSCArray<int> A = CSCArray<int> (data, rows, cols, 3,3);
    
    
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            if( (i != 2) || (j != 2)){
                EXPECT_EQ(A(i,j), (j*3)+i + 1) << "Element access is different than expected at " << i << " " << j;
            } else {
                EXPECT_EQ(A(i,j), 0) << "Element was expected to be 0 at " << i << " " << j; 
            }
        }
    }
    
}

// The array here and below looks like 
// | 1 2 0 0 0 0 |
// | 0 3 0 4 0 0 |
// | 0 0 5 6 7 0 |
// | 0 0 0 0 0 8 |
// Initially but we then modify it during the test
TEST(CSCArray,  ModifyValue){
    CArray<int> data(8);
    CArray<size_t> cols(7);
    CArray<size_t> rows(8);
    buildTestArray(data, cols, rows);
    CSCArray<int> B  = CSCArray<int>(data, rows, cols, 4, 6 );
    //B.setVal(3,5,99);
    B(3,5) = 99;
    EXPECT_EQ(B(3,5), 99) << "Element was expected to be 99 at " << 3 << " " << 5 ;
    // This test below should not change anything as the code currently exists
    B(0,5) = 66;
    B(1,0) = 199;
    EXPECT_EQ(B(0,5), 0) << "This operation should not change anything but we got value " << B(0,5);    
}

TEST(CSCArray, NonZeroRow){
    CArray<int> data(8);
    CArray<size_t> cols(7);
    CArray<size_t> rows(8);
    buildTestArray(data, cols, rows);
    CSCArray<int> B  = CSCArray<int>(data, rows, cols, 4, 6 );
  
    EXPECT_EQ(B.nnz(0), 1) << "expected there to be " << 1 << " nonzero elements in col " << 0; 
    EXPECT_EQ(B.nnz(1), 2) << "expected there to be " << 2 << " nonzero elements in col " << 1; 
    EXPECT_EQ(B.nnz(2), 1) << "expected there to be " << 2 << " nonzero elements in col " << 2; 
    EXPECT_EQ(B.nnz(3), 2) << "expected there to be " << 1 << " nonzero elements in col " << 3; 

}

TEST(CSCArray, NonZeroTotal){
    CArray<int> data(8);
    CArray<size_t> cols(7);
    CArray<size_t> rows(8);
    buildTestArray(data, cols, rows);
    CSCArray<int> B  = CSCArray<int>(data, rows, cols, 4, 6 );
    EXPECT_EQ(B.nnz(), 8) << "expected there to be " << 8 << " nonzero elements in total "; 
}

// Flat index take i j and returns the position it is in the underlying 1d array 
// if A(i,j) is unfilled then it isn't represented in the 1d array so we opt to 
// return -1
TEST(CSCArray, FlatIndex){
    CArray<int> data(8);
    CArray<size_t> cols(7);
    CArray<size_t> rows(8);
    buildTestArray(data, cols, rows);
    CSCArray<int> B  = CSCArray<int>(data, rows, cols, 4, 6 );
    EXPECT_EQ(B.flat_index(3,5), 7) << "Expected the location of 3,5 to be at " << 7 ;
    EXPECT_EQ(B.flat_index(3,0), -1) << "There shouldn't be a proper value here, we are supposed to return -1 in this case";

}

// .begin(i) and .end(i) work as expected if you have seen how the iterators from std datastructures like vector 
// work. Specifically both return a pointer to the data array
TEST(CSCArray, ColAccess){
    CArray<int> data(8);
    CArray<size_t> cols(7);
    CArray<size_t> rows(8);
    buildTestArray(data, cols, rows);
    CSCArray<int> B  = CSCArray<int>(data, rows, cols, 4, 6 );  
    
    size_t elements_row = B.nnz(1);    
    int i = 0; 
    for(auto start = B.begin(1); start != B.end(1); ++start){    
        if(i == 0){
            EXPECT_EQ(2, *start) << "Wrong value of " << *start << " in iteration 0";
            i++; 
        } else if(i == 1) {
            EXPECT_EQ(3, *start) << "Wrong value of " << *start << " in iterion 1";
        } else {
            EXPECT_EQ(0,1) << "We've looped too far :( ";
         }   
    }
}

// This code goes through the flat arrays. 
// The important difference here is we return the index to the start and end position instead of a pointer
// One on hand this is useful because we can access both the column information and the actual data if we 
// return the integer
TEST(CSCArray, RowAccessFlat){
    CArray<int> data(8);
    CArray<size_t> cols(7);
    CArray<size_t> rows(8);
    buildTestArray(data, cols, rows);
    CSCArray<int> B  = CSCArray<int>(data, rows, cols, 4, 6 );
     
    size_t elements_row = B.nnz(1);    
    int i = 0; 
    for(auto start = B.begin_index(1); start != B.end_index(1); ++start){    
        if(i == 0){
            EXPECT_EQ(2, B.get_val_flat(start)) << "Wrong value of " << start << " in iteration 0";
            i++; 
        } else if(i == 1) {
            EXPECT_EQ(3, B.get_val_flat(start)) << "Wrong value of " << start << " in iterion 1";
        } else {
            EXPECT_EQ(0,1) << "We've looped too far :( ";
         }   
    }
}

TEST(CSCArray, toSparseColArray){
    CArray<int> data(8);
    CArray<size_t> cols(7);
    CArray<size_t> rows(8);
    buildTestArray(data, cols, rows);
    CSCArray<int> B  = CSCArray<int>(data, rows, cols, 4, 6 );
    size_t elements_row = B.nnz(1);
    size_t nnz = B.nnz();
    CArray<int> out_data(nnz);
    CArray<size_t> out_cols(nnz);
    CArray<size_t> out_rows(nnz); // change this
    B.toCSR(out_data, out_cols, out_rows);
    int expected_data[] = {1, 2, 3, 5, 4, 6, 7, 8};
    for(int i = 0; i < 8; i++){
        EXPECT_EQ(expected_data[i], out_data(i)) << "Wrong value at index: " << i << " got " << out_data(i) << " but wanted: " << expected_data[i];
    }
}

int main(int argc, char* argv[]){
    int result = 0;
        
    testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
    return result; 
    
} 
