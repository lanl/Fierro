#include <stdio.h>
#include <math.h>
#include <matar.h>
#include <limits.h>
#include <chrono> 
#include <time.h>

#define EXPORT true 

using namespace mtr; // matar namespace

void matVecSparse(CSRArrayKokkos<double> &A, CArrayKokkos<double> &v, CArrayKokkos<double> &b){
        size_t m = A.dim2();
        size_t n = A.dim1();
        FOR_ALL(i, 0, n, { 
                size_t col;
                for(auto j = A.begin_index(i); j <  A.end_index(i); j++){
                    col = A.get_col_flat(j);    
                    b(i) += A(i,col) * v(col);
                }
            }
        );
        Kokkos::fence();
}
 

int main(int argc, char** argv){
    Kokkos::initialize(); {
    int nrows = 55;
    int ncols = 55;
    size_t n;
    if(argc != 2){
            printf("Usage is .powerTest <MatrixSize> using default of 5000\n");
            n = 5000;
        } else{
             n = (size_t) atoi(argv[1]);
        } 
    nrows = n;
    ncols = n;
    CArrayKokkos<double> data(3*nrows);
    CArrayKokkos<size_t> starts(nrows+1);
    CArrayKokkos<size_t> cols(3*nrows);
    CArrayKokkos<double> v1(ncols);
    CArrayKokkos<double> v2(ncols);
    CArrayKokkos<double> b1(nrows);
    CArrayKokkos<double> b2(nrows);
    
    int i;
    i = 0;
    FOR_ALL(i, 0, ncols,{
                    v1(i) = 1;
                    v2(i) = 1;
                    b1(i) = 0;
                    b2(i) = 0;
                    });
    FOR_ALL(i, 0, nrows,{
                    
                    if(i == nrows -2){
                        data(3*i) = i;
                        data(3*i+1) = i;
                        data(3*i+2) = i;
                        cols(3*i) = i-1;
                        cols(3*i+1) = i;
                        cols(3*i+2) = i+1;
                        b1(i) = 0;
                        b2(i) = 0;
                        starts(i) = 3*i;
                    }
                    else if(i == nrows -1){
                        data(3*i) = i;
                        data(3*i+1) = i;
                        data(3*i+2) = i;
                        cols(3*i) = i-2;
                        cols(3*i+1) = i-1;
                        cols(3*i+2) = i;
                        b1(i) = 0;
                        b2(i) = 0;
                        starts(i) = 3*i;
                    }
                    else {
                        data(3*i) = i ;
                        data(3*i+1) = i;
                        data(3*i+2) = i;
                        cols(3*i) = i;
                        cols(3*i+1) = i+1;
                        cols(3*i+2) = i+2;
                        b1(i) = 0;
                        b2(i) = 0;
                        starts(i) = 3*i;
                    }
    }); 
    RUN({ 
        starts(0) = 0;
        starts(nrows) = 3*nrows;
    });
    CSRArrayKokkos<double> B (data, starts, cols, nrows, ncols);
    auto start = std::chrono::high_resolution_clock::now();

    matVecSparse(B,v2,b2);
    Kokkos::fence(); 
    auto lap1 = std::chrono::high_resolution_clock::now();
    auto time1 =  std::chrono::duration_cast<std::chrono::nanoseconds>(lap1 - start); 

    if(!EXPORT){
            RUN({ printf("Size: %ld,  Sparse: %.2e, %f, %f \n", n, time1.count() * 1e-9 , b1(57980), b2(57980) ); });
        } else {
                RUN({
                    for(int i = 0; i < n; i++){
                        if(abs(b1(i) - b2(i) > 1e-7)){
                            printf("b1(%d) - b2(%d) = %.2e\n", i, i, b1(i)-b2(i));
                        }
                    }
                    printf("%ld, %.2e, %.2e, %f, %f, %f \n", n, time1.count() * 1e-9, b1(25), b2(25) );
                }); 
       }
    }Kokkos::finalize();
}
