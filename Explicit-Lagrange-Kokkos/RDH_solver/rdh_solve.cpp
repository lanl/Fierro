// Call cycle loop for the SGH solver


#include "state.h"
#include "mesh.h"
#include <chrono>

void rdh_solve(CArrayKokkos <material_t> &material,
               CArrayKokkos <boundary_t> &boundary,
               mesh_t &mesh,
               DViewCArrayKokkos <double> &node_coords,
               DViewCArrayKokkos <double> &node_vel,
               DViewCArrayKokkos <double> &node_mass,
               DViewCArrayKokkos <double> &elem_den,
               DViewCArrayKokkos <double> &elem_pres,
               DViewCArrayKokkos <double> &elem_stress,
               DViewCArrayKokkos <double> &elem_sspd,
               DViewCArrayKokkos <double> &elem_sie,
               DViewCArrayKokkos <double> &elem_vol,
               DViewCArrayKokkos <double> &elem_div,
               DViewCArrayKokkos <double> &elem_mass,
               DViewCArrayKokkos <size_t> &elem_mat_id,
               DViewCArrayKokkos <double> &elem_statev,
               DViewCArrayKokkos <double> &corner_force,
               DViewCArrayKokkos <double> &corner_mass,
               double &time_value,
               const double time_final,
               const double dt_max,
               const double dt_min,
               const double dt_cfl,
               double &graphics_time,
               size_t graphics_cyc_ival,
               double graphics_dt_ival,
               const size_t cycle_stop,
               const size_t rk_num_stages,
               double dt,
               const double fuzz,
               const double tiny,
               const double small,
               CArray <double> &graphics_times,
               size_t &graphics_id){
    
    
    printf("Writing outputs to file at %f \n", time_value);
    write_outputs(mesh,
                  node_coords,
                  node_vel,
                  node_mass,
                  elem_den,
                  elem_pres,
                  elem_stress,
                  elem_sspd,
                  elem_sie,
                  elem_vol,
                  elem_mass,
                  elem_mat_id,
                  graphics_times,
                  graphics_id,
                  time_value);
    
        VTKHexN(mesh,
                node_coords,
                node_vel,
                node_mass,
                elem_den,
                elem_pres,
                elem_stress,
                elem_sspd,
                elem_sie,
                elem_vol,
                elem_mass,
                elem_mat_id,
                graphics_times,
                graphics_id,
                time_value);
    
    return;
    
} // end of RDH solve

