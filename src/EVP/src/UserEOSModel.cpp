#include "UserEOSModel.h"

KOKKOS_FUNCTION
UserEOSModel::UserEOSModel(
  const DCArrayKokkos <material_t> &material,
  const DCArrayKokkos <double> &elem_state_vars,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t mat_id,
  const size_t elem_gid)
{
}

KOKKOS_FUNCTION
UserEOSModel::~UserEOSModel()
{
}

KOKKOS_FUNCTION
int UserEOSModel::calc_sound_speed(
  const DViewCArrayKokkos <double> &elem_pres,
  const DViewCArrayKokkos <double> &elem_stress,
  const size_t elem_gid,
  const size_t mat_id,
  const DCArrayKokkos <double> &elem_state_vars,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const DViewCArrayKokkos <double> &elem_sspd,
  const double den,
  const double sie)
{
    elem_sspd(elem_gid) = 0.357;//1.0e-8;  // sound speed
    
    return 0;
}

KOKKOS_FUNCTION
int UserEOSModel::calc_pressure(
  const DViewCArrayKokkos <double> &elem_pres,
  const DViewCArrayKokkos <double> &elem_stress,
  const size_t elem_gid,
  const size_t mat_id,
  const DCArrayKokkos <double> &elem_state_vars,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const DViewCArrayKokkos <double> &elem_sspd,
  const double den,
  const double sie)
{
    const int num_dims = 3;
    
    elem_pres(elem_gid) = 0.0;  // pressure
    
    // pressure = 1/3tr(stress)
    for (int i=0; i<num_dims; i++){
        elem_pres(elem_gid) -= elem_stress(1,elem_gid,i,i);
    }
    elem_pres(elem_gid) *= 1.0/3.0;
 
    return 0;
}


