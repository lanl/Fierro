#include "matar.h"
#include "mesh.h"
#include "rdh.h"
#include "state.h"
#include "ref_elem.h"
#include "ref_surf_elem.h"
#include <cmath>

#define PI 3.141592653589793

using namespace mtr;


// These functions are called inside get_*_rhs, where *  is either momentum or energy

KOKKOS_FUNCTION
void get_viscosity_coefficient(const mesh_t &mesh,
                               const fe_ref_elem_t &ref_elem,
                               DViewCarrayKokkos <double> &alpha,
                               const DViewCArrayKokkos <double> &sspd){
                             
    // get maximum wave speed in each element

}

KOKKOS_FUNCTION
void get_viscosity(const mesh_t &mesh,
                    const fe_ref_elem_t &ref_elem,
                    DViewCArrayKokkos <double> &F){
    
}