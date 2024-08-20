// -----------------------------------------------------------------------------
// This code contains the constitutive relation for a user supplied model
//------------------------------------------------------------------------------
#include "state.h"
#include "mesh.h"
#include "VPSC7.h"

// -----------------------------------------------------------------------------
// This is the user material model function for the equation of state
// An eos function must be supplied or the code will fail to run.
// The pressure and sound speed can be calculated from an analytic eos.
// The pressure can also be calculated using p = -1/3 Trace(Stress)
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
void user_eos_model(const DViewCArrayKokkos <double> &elem_pres,
                    const DViewCArrayKokkos <double> &elem_stress,
                    const size_t elem_gid,
                    const DViewCArrayKokkos <double> &elem_sspd){
    
    const size_t num_dims = 3;
    
    // -----------------------------------------------------------------------------
    // Required variables are here
    //------------------------------------------------------------------------------
    elem_pres(elem_gid) = 0.0;  // pressure
    elem_sspd(elem_gid) = 2400.0;  // sound speed
    
    // pressure = 1/3tr(stress)
    for (size_t i=0; i<num_dims; i++){
        elem_pres(elem_gid) -= elem_stress(1,elem_gid,i,i);
    }
    elem_pres(elem_gid) *= 1.0/3.0;
    
    
    return;
    
} // end for user_eos_model





// -----------------------------------------------------------------------------
// This is the user material model function for the stress tensor
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
void user_strength_model(const DViewCArrayKokkos <double> &elem_pres,
                         const DViewCArrayKokkos <double> &elem_stress,
                         const size_t elem_gid,
                         const size_t mat_id,
                         const DViewCArrayKokkos <double> &elem_state_vars,
                         const DViewCArrayKokkos <double> &elem_sspd,
                         const double den,
                         const double sie,
                         const ViewCArrayKokkos <double> &vel_grad,
                         const ViewCArrayKokkos <size_t> &elem_node_gids,
                         const DViewCArrayKokkos <double> &node_coords,
                         const DViewCArrayKokkos <double> &node_vel,
                         const double vol,
                         const double dt,
                         const double rk_alpha){
    

    // statev(0) = var_1
    //   :
    //   :
    //   :
    // statev(N) = var_N

    int num_dims = 3;

    
    // -----------------------------------------------------------------------------
    // The user must coding goes here
    //------------------------------------------------------------------------------

    
    return;
    
} // end of user mat


// -----------------------------------------------------------------------------
// This is the user material model function
//------------------------------------------------------------------------------
KOKKOS_FUNCTION
void user_strength_model_vpsc(const DViewCArrayKokkos <double> &elem_pres,
                              const DViewCArrayKokkos <double> &elem_stress,
                              const size_t elem_gid,
                              const size_t mat_id,
                              const DViewCArrayKokkos <double> &elem_state_vars,
                              const DViewCArrayKokkos <double> &elem_sspd,
                              const double den,
                              const double sie,
                              const ViewCArrayKokkos <double> &vel_grad,
                              const ViewCArrayKokkos <size_t> &elem_node_gids,
                              const DViewCArrayKokkos <double> &node_coords,
                              const DViewCArrayKokkos <double> &node_vel,
                              const double vol,
                              const double dt,
                              const double rk_alpha){
    

    // statev(0) = var_1
    //   :
    //   :
    //   :
    // statev(N) = var_N

    const int num_dims = 3;
    
    // -----------------------------------------------------------------------------
    // Required variables are here
    //------------------------------------------------------------------------------
    //elem_pres(elem_gid) = 1.0e-15;  // pressure
    //elem_sspd(elem_gid) = 1.0e-15;  // sound speed

    
    // -----------------------------------------------------------------------------
    // The user must coding goes here
    //------------------------------------------------------------------------------
    
    // For hypo-elastic models
	double ddsdde_values[36];
	
	double stress_scale;
	double Lnorm;
	
	stress_scale = 1.0;
	
    ViewCMatrixKokkos <double> ddsdde(&ddsdde_values[0], 6, 6);
    
	Lnorm = 0.0;
	
    // calculate strain and rotation increments
    for (size_t i = 1; i <= 3; i++) {
        for (size_t j = 1; j <= 3; j++) {
			Lnorm = Lnorm +  vel_grad(i-1, j-1)*vel_grad(i-1, j-1);
        }
    } // end for
    
	
	
    //interpret the state variable array
	if (sqrt(Lnorm) > 1.0e-16) {
        size_t cnt = 0;
		
		PropType props;
        StateVType statv;
        IntVType intv;


        props.init(&elem_state_vars(elem_gid, 0), cnt);
        statv.init(&elem_state_vars(elem_gid, 0), cnt);
        intv.init(&elem_state_vars(elem_gid, 0), cnt);


        //supply statv with defined vars
        // update the temperature
        update_temperature(statv.svm, statv.evmp, statv.evmpold, statv.temp, statv.dtemp, statv.wrplastic, statv.wplastic, props.b_int);

        // temp = statv.temp;
        evol_vpsc(statv, props, intv, vel_grad, &elem_stress(1, elem_gid, 0, 0), dt, elem_gid, ddsdde);

        statv.tempold = statv.temp;

        // write the scalars to the state variable array
        props.write(&elem_state_vars(elem_gid, 0));
        statv.write(&elem_state_vars(elem_gid, 0));
		
		
		
		if (elem_stress(1, elem_gid, 0, 0) != elem_stress(1, elem_gid, 0, 0)) {
            printf("L(1, 1) : %16.4e; stress(1, 1) : %16.8f;\n", vel_grad(0, 0), elem_stress(1, elem_gid, 0, 0));
            exit(1);
        }
    }
	
	//printf("deps(3, 3) : %16.4e; stress(3, 3) : %16.8f;\n", deps(3, 3), elem_stress(1, elem_gid, 2, 2));
    return;
    
} // end of ideal_gas

