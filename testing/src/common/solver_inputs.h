

#ifndef FIERRO_SOLVER_INPUT_OPTIONS_H
#define FIERRO_SOLVER_INPUT_OPTIONS_H
#include <stdio.h>
#include "matar.h"

namespace solver_input
{

    // solver method
    enum method
    {
        SGH = 0,        
    };
} // end of namespace

static std::map <std::string, solver_input::method> solver_map
{
    {"SGH", solver_input::SGH}
};


// solver input parameters
struct solver_input_t{

    solver_input::method method;

}; // solver_input_t

// ----------------------------------
// valid inputs for solver options
// ----------------------------------
static std::vector <std::string> str_solver_inps
{
    "method"
};



#endif // end Header Guard