#include <stdio.h>
#include <math.h>
#include <matar.h>
#include <limits.h>

using namespace mtr; // matar namespace
    

int main(int argc, char* argv[]){
    
    Kokkos::initialize(); { 
    size_t nnz = 6 ;
    size_t dim1 = 3;
    size_t dim2 = 10;
    CArrayKokkos<size_t> starts(dim2 + 1);
    CArrayKokkos<size_t> rows(nnz);
    CArrayKokkos<int> array(nnz + 1); 
    RUN ({
        starts(1) = 1;
        starts(2) = 2;
        starts(3) = 3;
        starts(4) = 4;
        starts(5) = 5; 
        starts(6) = 6; 
        starts(7) = 6;
        starts(8) = 6;
        starts(9) = 6;

        rows(0) = 0;
        rows(1) = 0;
        rows(2) = 1;
        rows(3) = 1;
        rows(4) = 2;
        rows(5) = 2;

        array(0) = 1;
        array(1) = 2;
        array(2) = 3;
        array(3) = 4;
        array(4) = 5;
        array(5) = 6;
        array(6) = 0;
                    
    });

    /*
    |1 2 2 0 0 0 0 0 0 0|
    |0 0 3 4 0 0 0 0 0 0|
    |0 0 0 0 5 6 0 0 0 0|
    */
    
    const std::string s = "hello";   
    // Testing = op 
    auto pre_A = CSCArrayKokkos<int>(array, starts,rows, dim1, dim2, s);
    auto A = pre_A;
    int* values = A.pointer();
    auto a_start = A.get_starts();
    int total = 0;
  
    RUN ({
         printf("This matix is %ld x %ld \n" , A.dim1(), A.dim2());
    });

    RUN ({
         printf("nnz : %ld \n", A.nnz());
    });


    int loc_total = 0;  
    loc_total += 0; // Get rid of warning
    REDUCE_SUM(i, 0, nnz,
                loc_total, {
                    loc_total += values[i];
                    }, total);    
    printf("Sum of nnz from pointer method %d\n", total);
    total = 0;
    REDUCE_SUM(i, 0, nnz,
                loc_total, {
                    loc_total += a_start[i];
                    }, total);    
    printf("Sum of start indices form .get_starts() %d\n", total); 
    total = 0;
    
    REDUCE_SUM(i, 0, dim1,
               j, 0, dim2-1,
                loc_total, {
                    loc_total += A(i,j);
                    }, total);    
    printf("Sum of nnz in array notation %d\n", total);
    
    } Kokkos::finalize();
    return 0; 

    
}
