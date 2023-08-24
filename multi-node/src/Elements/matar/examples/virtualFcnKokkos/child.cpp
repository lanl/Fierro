//
//  child.h
//
//
//

#include <stdio.h>
#include "child.hpp"



//----------------------------
    KOKKOS_FUNCTION
    child_variables::child_variables() {};

//----------------------------
    KOKKOS_FUNCTION
    child_models::child_models() {};

    KOKKOS_FUNCTION
    baby1::baby1(double glitter, double food){
        this_glitter = glitter;
        this_food = food;
    }
    
    KOKKOS_FUNCTION
    double baby1::math(double jump, double bounce){
        return this_glitter + this_food;
    }

    KOKKOS_FUNCTION
    baby2::baby2() {
        this_glitter = 0.0;
        this_food = 0.0;
    }
    
    KOKKOS_FUNCTION
    double baby2::math(double jump, double bounce){
        double sum = jump * bounce;
        return sum;
    }










