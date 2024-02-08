// -----------------------------------------------------------------------------
// This code contains the equation of state information (constitutive relations)
// ------------------------------------------------------------------------------
#include "state.h"

// -----------------------------------------------------------------------------
// This is the gamma-law gas eos
// ------------------------------------------------------------------------------
KOKKOS_FUNCTION
void ideal_gas(const DViewCArrayKokkos<double>& elem_pres,
               const DViewCArrayKokkos<double>& elem_stress,
               const size_t                     elem_gid,
               const size_t                     mat_id,
               const DViewCArrayKokkos<double>& elem_state_vars,
               const DViewCArrayKokkos<double>& elem_sspd,
               const double                     den,
               const double                     sie)
{
    // statev(0) = gamma
    // statev(1) = minimum sound speed
    // statev(2) = specific heat c_v
    // statev(3) = ref temperature
    // statev(4) = ref density
    // statev(5) = ref specific internal energy

    double gamma = elem_state_vars(elem_gid, 0);
    double csmin = elem_state_vars(elem_gid, 1);

    // pressure
    elem_pres(elem_gid) = (gamma - 1.0) * sie * den;

    // sound speed
    elem_sspd(elem_gid) = sqrt(gamma * (gamma - 1.0) * sie);

    // ensure soundspeed is great than min specified
    if (elem_sspd(elem_gid) < csmin)
    {
        elem_sspd(elem_gid) = csmin;
    } // end if

    return;
} // end of ideal_gas
