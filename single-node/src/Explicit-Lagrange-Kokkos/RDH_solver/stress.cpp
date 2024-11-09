#include "matar.h"
#include "mesh.h"
#include "rdh.h"
#include "state.h"
#include "ref_elem.h"
#include <cmath>

#define PI 3.141592653589793

using namespace mtr;

void get_stress(const mesh_t &mesh,
                const mat_pt_t &mat_pt,
                const DViewCArrayKokkos <double> &pressure,
                DViewCArrayKokkos <double> &sigma,
                const int stage){
    
    FOR_ALL(gauss_gid, 0, mat_pt.num_leg_pts, {
        
        for (int i = 0; i < mesh.num_dims; i++){

            for (int j = 0; j < mesh.num_dims; j++){

                sigma(stage, gauss_gid, i, j) = 0.0;

            }

        }

    });// FOR_ALL
    
    FOR_ALL(gauss_gid, 0, mat_pt.num_leg_pts, {
        
        for (int i = 0; i < mesh.num_dims; i++){

                sigma(stage, gauss_gid, i, i) = -1.0*pressure(gauss_gid);

        }
        
    });// FOR_ALL
    
}// end get_stress

void get_SigmaJacInv(const mesh_t &mesh,
                     const mat_pt_t &mat_pt,
                     const DViewCArrayKokkos <double> &sigma,
                     const DViewCArrayKokkos <double> &JacInv,
                     DViewCArrayKokkos <double> &SigmaJacInv,
                     const int stage){

    FOR_ALL( gauss_gid, 0, mat_pt.num_leg_pts, {
        for (int j = 0; j < mesh.num_dims; j++){
            for (int k = 0; k < mesh.num_dims; k++){
                SigmaJacInv(stage, gauss_gid, j, k) = 0.0;
            }
        }
    });

    FOR_ALL( gauss_gid, 0, mat_pt.num_leg_pts, {
        //  printf("-------------\n");
        for (int k = 0; k < mesh.num_dims; k++){
            for (int m = 0; m < mesh.num_dims; m++){
                
                SigmaJacInv(stage, gauss_gid, k, m) = sigma(stage, gauss_gid, k, 0)*JacInv(gauss_gid, m, 0)
                                                      + sigma(stage, gauss_gid, k, 1)*JacInv(gauss_gid, m, 1)
                                                      + sigma(stage, gauss_gid, k, 2)*JacInv(gauss_gid, m, 2);

                // printf(" %f ", SigmaJacInv(stage, gauss_gid, k, m));
            }// m
            // printf("\n");
        }// k
        // printf("-------------\n");
    });

}
