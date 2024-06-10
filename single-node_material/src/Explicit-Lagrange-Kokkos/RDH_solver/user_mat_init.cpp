// -----------------------------------------------------------------------------
// This code contains the initialization of state vars for supplied models
//------------------------------------------------------------------------------
#include <string.h>
#include <sys/stat.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm>



#include "mesh.h"
#include "state.h"




// -----------------------------------------------------------------------------
// The function to read in the state vars for a user supplied model
//------------------------------------------------------------------------------
void user_model_init(const DCArrayKokkos <double> &file_state_vars,
                     const size_t num_state_vars,
                     const size_t mat_id,
                     const size_t num_elems) {

    // initialize to zero
    for (size_t elem_gid = 0; elem_gid<num_elems; elem_gid++) {
        for(size_t var=0; var<num_state_vars; var++){
            file_state_vars.host(mat_id,elem_gid,var) = 0.0;
        }
    }
    

	return;
}


