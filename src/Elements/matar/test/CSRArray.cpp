#include <stdio.h>
#include <math.h>
#include <matar.h>
#include <limits.h>
#include <gtest/gtest.h> 

using namespace mtr; // matar namespace

void voidTestArray(CArray<int> &data, CArray<size_t> &cols, CArray<size_t> &rows){
    int i;
    for(i = 0; i < 8; i++){
        data(i) = i+1;
        cols(i) = (i % 3);
    }
    rows(0) = 0;
    rows(1) = 2;
    rows(2) = 4;
    rows(3) = 7;
    rows(4) = 8;
    cols(0) = 0;
    cols(1) = 1;
    cols(2) = 1;
    cols(3) = 3;
    cols(4) = 2;
    cols(5) = 3;
    cols(6) = 4;
    cols(7) = 5;
    CSRArray<int> B  = CSRArray<int>(data, cols, rows, 4, 6);
}
// The array here looks like 
//  | 1 2 3 | 
//  | 4 5 6 |
//  | 7 8 0 | 
TEST(CSRArray, AccessinData){
    CArray<int> data(8);
    CArray<size_t> cols(8);
    CArray<size_t> rows(4);
    int i, j;
    for(i = 0; i < 8; i++){
        data(i) = i+1;
        cols(i) = (i % 3);
    }
    rows(0) = 0;
    rows(1) = 3;
    rows(2) = 6;
    rows(3) = 8;


    CSRArray<int> A = CSRArray<int> (data, cols, rows, 3,3);

    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            if( (i != 2) || (j != 2)){
                EXPECT_EQ(A(i,j), (i*3)+j + 1) << "Element access is different than expected at " << i << " " << j;
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
TEST(CSRArray,  ModifyValue){
    CArray<int> data(8);
    CArray<size_t> cols(8);
    CArray<size_t> rows(5);
    voidTestArray(data, cols, rows);
    CSRArray<int> B(data, cols, rows, 4,6);   
    B(3,5) = 99;
    EXPECT_EQ(B(3,5), 99) << "Element was expected to be 99 at " << 3 << " " << 5 ;
    // This test below should not change anything as the code currently exists
    //B.setVal(0,5, 66);
    B(0,5) = 66;
    B(1,0) = 199;
    EXPECT_EQ(B(0,5), 0) << "This operation should not change anything but we got value " << B(0,5);    
}

TEST(CSRArray, NonZeroRow){
    CArray<int> data(8);
    CArray<size_t> cols(8);
    CArray<size_t> rows(5);
    voidTestArray(data, cols, rows);
    CSRArray<int> B(data, cols, rows, 4,6);   

    EXPECT_EQ(B.nnz(0), 2) << "expected there to be " << 2 << " nonzero elements in row " << 0; 
    EXPECT_EQ(B.nnz(1), 2) << "expected there to be " << 2 << " nonzero elements in row " << 1; 
    EXPECT_EQ(B.nnz(2), 3) << "expected there to be " << 3 << " nonzero elements in row " << 2; 
    EXPECT_EQ(B.nnz(3), 1) << "expected there to be " << 1 << " nonzero elements in row " << 3; 

}

TEST(CSRArray, NonZeroTotal){
    CArray<int> data(8);
    CArray<size_t> cols(8);
    CArray<size_t> rows(5);
    voidTestArray(data, cols, rows);
    CSRArray<int> B(data, cols, rows, 4,6);   
    EXPECT_EQ(B.nnz(), 8) << "expected there to be " << 8 << " nonzero elements"; 

}

// Flat index take i j and returns the position it is in the underlying 1d array 
// if A(i,j) is unfilled then it isn't represented in the 1d array so we opt to 
// return -1
TEST(CSRArray, FlatIndex){
    CArray<int> data(8);
    CArray<size_t> cols(8);
    CArray<size_t> rows(5);
    voidTestArray(data, cols, rows);
    CSRArray<int> B(data, cols, rows, 4,6);   
    EXPECT_EQ(B.flat_index(3,5), 7) << "Expected the location of 3,5 to be at " << 7 ;
    EXPECT_EQ(B.flat_index(3,0), 8) << "There shouldn't be a proper value here, we are supposed to return -1 in this case";

}

// .begin(i) and .end(i) work as expected if you have seen how the iterators from std datastructures like vector 
// work. Specifically both return a pointer to the data array
TEST(CSRArray, RowAccess){
    CArray<int> data(8);
    CArray<size_t> cols(8);
    CArray<size_t> rows(5);
    voidTestArray(data, cols, rows);
    CSRArray<int> B(data, cols, rows, 4,6);   

    size_t elements_row = B.nnz(1);    
    int i = 0; 
    for(auto start = B.begin(1); start != B.end(1); ++start){    
        if(i == 0){
            EXPECT_EQ(3, *start) << "Wrong value of " << *start << " in iteration 0";
            i++; 
        } else if(i == 1) {
            EXPECT_EQ(4, *start) << "Wrong value of " << *start << " in iterion 1";
        } else {
            EXPECT_EQ(0,1) << "We've looped too far :( ";
         }   
    }
}

// This code goes through the flat arrays. 
// The important difference here is we return the index to the start and end position instead of a pointer
// One on hand this is useful because we can access both the column information and the actual data if we 
// return the integer
TEST(CSRArray, RowAccessFlat){
    CArray<int> data(8);
    CArray<size_t> cols(8);
    CArray<size_t> rows(5);
    voidTestArray(data, cols, rows);
    CSRArray<int> B(data, cols, rows, 4,6);   

     
    size_t elements_row = B.nnz(1);    
    int i = 0; 
    for(auto start = B.begin_index(1); start != B.end_index(1); ++start){    
        if(i == 0){
            EXPECT_EQ(3, B.get_val_flat(start)) << "Wrong value of " << start << " in iteration 0";
            i++; 
        } else if(i == 1) {
            EXPECT_EQ(4, B.get_val_flat(start)) << "Wrong value of " << start << " in iterion 1";
        } else {
            EXPECT_EQ(0,1) << "We've looped too far :( ";
         }   
    }
}

TEST(CSRArray, toCSCArray){
    CArray<int> data(8);
    CArray<size_t> cols(8);
    CArray<size_t> rows(5);
    voidTestArray(data, cols, rows);
    CSRArray<int> B(data, cols, rows, 4,6);   

     
    int i;
    
    size_t elements_row = B.nnz(1);
    size_t nnz = B.nnz();
    CArray<int> out_data(nnz);
    CArray<size_t> out_cols(6);
    CArray<size_t> out_rows(nnz); // change this
    B.toCSC(out_data, out_cols, out_rows);
    int expected_data[] = {1, 2, 3, 5, 4, 6, 7, 8}; 
    for(i = 0; i < 8; i++){
        EXPECT_EQ(expected_data[i], out_data(i)) << "Wrong value at index: " << i << " got " << out_data(i) << " but wanted: " << expected_data[i];
    }
}


int main(int argc, char* argv[]){
    int result = 0;
        
    testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
    return result; 
    
} 
