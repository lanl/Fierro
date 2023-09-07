//
//  child.h
//
//
//
#ifndef CHILD_H
#define CHILD_H

#include <stdio.h>
#include <Kokkos_Core.hpp>

enum child_type
{
    BABY1,
    BABY2
};

class child_variables {
    protected:

    public:
        int num_pnts;
        int type;
        double *p;

        KOKKOS_FUNCTION
        child_variables();

        KOKKOS_FUNCTION
        ~child_variables(){};
};


//----------------------------

//----------------------------
class child_models{

    protected:
        // local parameters
        double this_glitter;
        double this_food;
    
    public:
        KOKKOS_FUNCTION
        child_models();

        // child functions
        KOKKOS_FUNCTION
        virtual double math(double jump, double bounce) {return 0.0;};
        
        KOKKOS_FUNCTION
        virtual ~child_models(){};
};


class baby1: public child_models{
    
    public:
    
        // constructor
        KOKKOS_FUNCTION
        baby1(double glitter, double food);
    
        KOKKOS_FUNCTION
        double math(double jump, double bounce);
};

class baby2: public child_models{

    public:
    
        // constructor
        KOKKOS_FUNCTION
        baby2();
    
        KOKKOS_FUNCTION
        double math(double jump, double bounce);
};


#endif
