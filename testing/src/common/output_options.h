

#ifndef FIERRO_OUTPUT_OPTIONS_H
#define FIERRO_OUTPUT_OPTIONS_H
#include <stdio.h>
#include "matar.h"

namespace output_options
{

    // output file options
    enum format
    {
        vtk = 0, 
        ensight = 1,
    };


    // timer output level
    enum timer_output_level
    {
        thorough = 0,   
    };

} // end of namespace

static std::map <std::string, output_options::format> output_format_map
{
    {"vtk", output_options::vtk},
    {"ensight", output_options::ensight}
};

static std::map <std::string, output_options::timer_output_level> timer_output_level_map
{
    {"thorough", output_options::thorough}
};

// mmeshing input parameters
struct output_options_t {

    output_options::format format;
    output_options::timer_output_level timer_level;

    real_t graphics_time_step = 1.0;
    int graphics_iteration_step = 2000000;
    


}; // output_options_t

// ----------------------------------
// valid inputs for output options
// ----------------------------------
static std::vector <std::string> str_output_options_inps
{
    "timer_output_level",
    "output_file_format",
    "graphics_time_step",
    "graphics_iteration_step"
};



#endif // end Header Guard