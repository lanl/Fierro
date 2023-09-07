
#include <stdio.h>
#include <iostream>
#include <matar.h>

#include <algorithm>  // std::max, std::min, etc.

using namespace mtr; // matar namespace

// main
int main(){

    // A view example
    int A[10];
    ViewCArray <int> arr(A,10);
    FOR_ALL (i, 0, 10, {
        arr(i) = 314;
    });


    // A 2D array example
    CArray <int> arr_2D(10,10);
    FOR_ALL (i, 0, 10,
             j, 0, 10,{
        arr_2D(i,j) = 314;
    });
    

    // A 3D array example
    CArray <int> arr_3D(10,10,10);
    FOR_ALL (i, 0, 10,
             j, 0, 10,
             k, 0, 10,{
        arr_3D(i,j,k) = 314;
    });


    int loc_sum = 0;
    int result = 0;
    REDUCE_SUM(i, 0, 10,
               loc_sum, {
                   loc_sum += arr(i)*arr(i);
               }, result);
    
    // testing
    loc_sum = 0;
    for (int i=0; i<10; i++){
        loc_sum += 314*314;
    }
    std::cout << "1D reduce : " << result << " vs. " << loc_sum << " \n";
    
    
    loc_sum = 0;
    result = 0;
    REDUCE_SUM(i, 0, 10,
               j, 0, 10,
               loc_sum, {
                   loc_sum += arr_2D(i,j)*arr_2D(i,j);
               }, result);
    
    // testing
    loc_sum = 0;
    for (int i=0; i<10; i++){
        for (int j=0; j<10; j++){
            loc_sum += 314*314;
        }
    }
    std::cout << "2D reduce : " << result << " vs. " << loc_sum << " \n";
    
    
    loc_sum = 0;
    result = 0;
    REDUCE_SUM(i, 0, 10,
               j, 0, 10,
               k, 0, 10,
               loc_sum, {
                   loc_sum += arr_3D(i,j,k)*arr_3D(i,j,k);
               }, result);
    
    // testing
    loc_sum = 0;
    for (int i=0; i<10; i++){
        for (int j=0; j<10; j++){
            for (int k=0; k<10; k++){
                loc_sum += 314*314;
            }
        }
    }
    std::cout << "3D reduce : " << result << " vs. " << loc_sum << " \n";
    
    
    int loc_max;
    REDUCE_MAX(i, 0, 10,
               j, 0, 10,
               k, 0, 10,
               loc_max, {
                   loc_max = std::max<int>(arr_3D(i,j,k), loc_max);
               }, result);
    
    std::cout << "3D reduce MAX : " << result << " \n";
    
    
    
    int loc_min;
    REDUCE_MIN(i, 0, 10,
               j, 0, 10,
               k, 0, 10,
               loc_min, {
                   loc_min = std::min<int>(arr_3D(i,j,k), loc_min);
               }, result);
    
    std::cout << "3D reduce MIN : " << result << " \n";
    
    
    REDUCE_MIN_CLASS(i, 0, 10,
               j, 0, 10,
               k, 0, 10,
               loc_min, {
                   loc_min = std::min<int>(arr_3D(i,j,k), loc_min);
               }, result);
    
    std::cout << "3D reduce MIN CLASS : " << result << " \n";
    
    
    std::cout << "done" << std::endl;

    return 0;
}


