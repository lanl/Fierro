// -----------------------------------------------------------------------------
// This code contains the equation of state information (constitutive relations)
//------------------------------------------------------------------------------
#include "state.h"
#include "rdh.h"





// -----------------------------------------------------------------------------
// This is the gamma-law gas eos
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
void ideal_gas(const DViewCArrayKokkos <double> &elem_pres,
               const DViewCArrayKokkos <double> &elem_stress,
               const size_t elem_gid,
               const size_t legendre_gid,
               const DViewCArrayKokkos <double> &elem_state_vars,
               const DViewCArrayKokkos <double> &elem_sspd,
               const double den,
               const double sie){
    

    // statev(0) = gamma
    // statev(1) = minimum sound speed
    // statev(2) = specific heat c_v
    // statev(3) = ref temperature
    // statev(4) = ref density
    // statev(5) = ref specific internal energy
    
    double gamma = elem_state_vars(elem_gid,0);
    double csmin = elem_state_vars(elem_gid,1);
    // printf("gamma : %f \n", gamma);
    // printf("csmin : %f \n", csmin);
    
    // pressure
    elem_pres(legendre_gid) = (gamma - 1.0)*sie*den;
    // printf("p : %f \n", elem_pres(legendre_gid));
    
    // sound speed
    elem_sspd(legendre_gid) = sqrt(gamma*(gamma - 1.0)*sie);
    // printf("c_s : %f \n", elem_sspd(legendre_gid));
    
    // ensure soundspeed is greater than min specified
    if (elem_sspd(legendre_gid) < csmin){
        elem_sspd(legendre_gid) = csmin;
    } // end if

    return;
} // end of ideal_gas


void update_thermo(const mesh_t &mesh, 
                   const fe_ref_elem_t &ref_elem,
                   const DViewCArrayKokkos <double> sie,
                   const DViewCArrayKokkos <double> den,
                   DViewCArrayKokkos <double> pres,
                   DViewCArrayKokkos <double> stress,
                   DViewCArrayKokkos <double> sspd,
                   const int stage,
                   const DViewCArrayKokkos <double> &elem_state_vars){
    FOR_ALL(elem_gid, 0, mesh.num_elems, {

        for (int gauss_lid = 0; gauss_lid < ref_elem.num_gauss_leg_in_elem; gauss_lid++){
            int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);

            double sie_val = 0.0;
			double val = 0.0;
            eval_sie(sie, elem_gid, gauss_lid, mesh, ref_elem, val, stage);
            //if (val < 0.0){
			//	printf("sie_val : %.25f \n", val);
			//}
			sie_val = fmax(0.0, val);
            double den_val = den(gauss_gid);
            
            ideal_gas(pres, stress, elem_gid, gauss_gid,
                      elem_state_vars, sspd, den_val, sie_val);

        }
    });// FOR_ALL
    Kokkos::fence();
}// end update thermo
