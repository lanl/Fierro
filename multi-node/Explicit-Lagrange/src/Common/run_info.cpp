// Prints cycle info

#include <stdio.h>
#include <iostream>  // std::cout etc.


#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

void run_info(int cycle){

// std::cout<<"(TIME-graphics_time)"<<(TIME-graphics_time)<<std::endl;
	if ((cycle%10) == 0) {
        real_t percent_comp = (TIME/TFINAL)*100.;
        printf("Percent complete = %.2f \n", percent_comp); 
    }

    
    if ((cycle%50) == 0) {
        printf("Step = %d   dt = %e  time = %e  \n ", cycle, dt, TIME);
    }

    // std::cout<<"DT = "<<dt<<std::endl;
    
    // output a graphics dump
    if ((cycle%graphics_cyc_ival) == 0 || (TIME-graphics_time) >= 0.0) {
        printf("****************Graphics write*************** \n");
        ensight();
        
        // set next time to dump a graphics output
        graphics_time += graphics_dt_ival;
    }

}
