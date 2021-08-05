
#include <string.h>
#include <sys/stat.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>


#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"
#include "user_mat.h"


using namespace utils;

void setup_material() {

    material_t mat;

    size_t mesh_mat_size = mesh.num_gauss_pts();
    size_t mat_sz = user_material_number_vars();

    printf("NSTATV = %i\nTSTATV = %i\n", mat_sz, mat_sz * mesh_mat_size);


    //Kokkos::View<double**, Layout, ExecSpace, void> state_vars("state_vars", mesh_mat_size, mat_sz);
    state_vars = (real_t*)malloc((size_t)(mat_sz * mesh_mat_size * sizeof(real_t)));

    for (size_t i = 0; i < (mat_sz * mesh_mat_size); i++) {
        state_vars[i] = 0.0;
    }

    //auto state_vars_host = Kokkos::create_mirror(state_vars);

    mat.init(&state_vars[0]);


    // WARNING; INITIALIZING ALL MATERIAL POINTS WITH SAME PROPERTIES
    for (size_t i = 1; i < mesh_mat_size; i++) {
        for (size_t j = 0; j < mat_sz; j++) {
            state_vars[i * mat_sz + j] = state_vars[j];
        }
    }



    //Kokkos::deep_copy(state_vars, state_vars_host);


	return;
}