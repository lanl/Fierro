// -----------------------------------------------------------------------------
// This code contains the equation of state information (constitutive relations)
//------------------------------------------------------------------------------
#include "state.h"




// -----------------------------------------------------------------------------
// This is the gamma-law gas eos
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
void ideal_gas(const DViewCArrayKokkos <double> &elem_pres,
               const size_t elem_gid,
               const DViewCArrayKokkos <size_t> &mat_id,
               const DViewCArrayKokkos <double> &elem_state_vars,
               const DViewCArrayKokkos <double> &elem_sspd,
               const DViewCArrayKokkos <double> &elem_den,
               const DViewCArrayKokkos <double> &elem_sie){
    
    // statev(0) = specific heat
    // statev(1) = ref temperature
    // statev(2) = ref density
    // statev(3) = ref specific internal energy
    // statev(4) = gamma
    // statev(5) = minimum sound speed
    
    double gamma = elem_state_vars(elem_gid,4);
    double csmin = elem_state_vars(elem_gid,5);
    
    
    // pressure
    elem_pres(elem_gid) = (gamma - 1.0)*
                           elem_sie(1,elem_gid)*elem_den(elem_gid);
    
    // sound speed
    elem_sspd(elem_gid) = gamma*(gamma - 1.0)*elem_sie(1,elem_gid);
    elem_sspd(elem_gid) = sqrt(elem_sspd(elem_gid));
    
    // ensure soundspeed is great than min specified
    if (elem_sspd(elem_gid) < csmin){
        elem_sspd(elem_gid) = csmin;
    } // end if

} // end of ideal_gas

