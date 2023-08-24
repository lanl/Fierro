//
//  parent.h
//  
//
//
#ifndef PARENT_H
#define PARENT_H
#include "child.hpp"


class parent_variables{
    
    public:
        //child_variables *child_var;
        int num_pnts;
        int type;
    
        // child variables
        double *child_p; // var


        // variables
        double *var_fake1;


    // ...
    
        // default constructor
        KOKKOS_FUNCTION
        parent_variables() {};


    
        // deconstructor
        KOKKOS_FUNCTION
        ~parent_variables(){};
    
};

class parent_models{
    
    public:
        child_models *child;

    // ...
    
    
    // deconstructor
        KOKKOS_FUNCTION
        ~parent_models(){};
    
    
}; // end of parent


#endif
