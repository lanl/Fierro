// -----------------------------------------------------------------------------
// This code contains the initialization of state vars for supplied models
//------------------------------------------------------------------------------
#include <string.h>
#include <sys/stat.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>

#include "user_mat.h"
#include "Simulation_Parameters/Material.h"

#define CLEAR_LINE(ifstream) ( ifstream.ignore(std::numeric_limits<std::streamsize>::max(), '\n') );
// for string delimiter parsing
std::vector<std::string> split (std::string s, std::string delimiter);

	//This goes in input.cpp
//	material(0).num_state_vars = 163;  // actual num_state_vars 
//	material(0).eos_model = user_eos_model; // EOS model is required
//	material(0).strength_type = model::hypo;
//	material(0).strength_setup = model_init::input;
//	material(0).strength_model = user_strength_model;

// -----------------------------------------------------------------------------
// The function to read in the state vars for a user supplied model
//------------------------------------------------------------------------------


void init_state_vars(
  const DCArrayKokkos <material_t> &material,
  const DViewCArrayKokkos <size_t> &elem_mat_id,
  const DCArrayKokkos <double> &state_vars,
  const DCArrayKokkos <double> &global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t num_elems)
{
    const real_t SQR2_ = 1.41421356237309;
    const real_t RSQ2_ = 0.70710678118654744;
    const real_t RSQ3_ = 0.57735026918962584;
    const real_t RSQ6_ = 0.40824829046386304;

	CMatrix <double> cc66v(6,6);
	CMatrix <double> cdim(3);
	CMatrix <int> isn(24,3);  //assumed BCC with 2 modes here
	CMatrix <int> isb(24,3);  //assumed BCC with 2 modes here
	CMatrix <double> sn(3);
	CMatrix <double> sb(3);
	CMatrix <double> aux5(5);
	CMatrix <double> aux33(3,3);
	// Accessed within state_vars for loop
	DFMatrixKokkos <double> dnca(3,24);  //assumed BCC with 2 modes here
	DFMatrixKokkos <double> dbca(3,24);  //assumed BCC with 2 modes here
	DCMatrixKokkos <double> cc66(6,6);
	DCMatrixKokkos <double> B_basis(3,3,6);
	// FMatrix <double> schca(5,12*2);  //assumed BCC with 2 modes here
	std::ifstream ur1;
	int iso,nmodesx,nmodes,modex,nsmx,nrsx,nsysx;
	real_t gamd0x,tau0xf,tau0xb,tau1x,thet0x,thet1x; //real_t may be something else in fierro
	real_t hselfx,hlatex,snor,qnor,prod,ph,th,om;
	std::string junk,icryst;
	 
	real_t stress_conv = 1.0e-5;

    for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
    		for (int k = 1; k <= 6; k++) {
    			B_basis.host(i,j,k) = 0.0;
  			}
		}
		}
    B_basis.host(1,1,2) = -RSQ6_;
    B_basis.host(2,2,2) = -RSQ6_;
    B_basis.host(3,3,2) =  real_t(2.0)*RSQ6_;

    B_basis.host(1,1,1) = -RSQ2_;
    B_basis.host(2,2,1) =  RSQ2_;

    B_basis.host(2,3,3) =  RSQ2_;
    B_basis.host(3,2,3) =  RSQ2_;

    B_basis.host(1,3,4) =  RSQ2_;
    B_basis.host(3,1,4) =  RSQ2_;

    B_basis.host(1,2,5) =  RSQ2_;
    B_basis.host(2,1,5) =  RSQ2_;

    B_basis.host(1,1,6) =  RSQ3_;
    B_basis.host(2,2,6) =  RSQ3_;
    B_basis.host(3,3,6) =  RSQ3_;
    B_basis.update_device();

	ur1.open("./ta_el.sx");//("./cuel1.txt");
	// read el txt file here
	ur1 >> iso; CLEAR_LINE(ur1);

	if (iso == 0) {
		for (int i = 1; i <= 6; i++) {
	    	for (int j = 1; j <= 6; j++) {
	      		ur1 >> cc66v(i,j);
	    	}
	    	CLEAR_LINE(ur1);
	  	}
	}

	// CMatrix <double> dde(3,3);
	// CMatrix <double> xid4(3,3,3,3);
	// CMatrix <double> cc3333(3,3,3,3);

	// for (int i = 1; i <= 3; i++) {
	//     for (int j = 1; j <= 3; j++) {
	//         dde(i,j) = 0.0;
	// 	    if (i == j) dde(i,j) = 1.0;
	//     }
	// }

	// for (int i = 1; i <= 3; i++) {
	//   for (int j = 1; j <= 3; j++) {
	//     for (int k = 1; k <= 3; k++) {
	//       for (int l = 1; l <= 3; l++) {
	//         xid4(i,j,k,l) = (dde(i,k)*dde(j,l) + dde(i,l)*dde(j,k)) / real_t(2.0);
	//       }
	//     }
	//   }
	// }
	// real_t tmu = 119800.0 / (real_t(2.0)*(real_t(1.0) + 0.35));
	// real_t tla = real_t(2.0)*tmu*0.35 / (real_t(1.0) - real_t(2.0)*0.35);

	// printf("%g %g\n",tmu,tla );

 //    for (int i = 1; i <= 3; i++) {
 //      for (int j = 1; j <= 3; j++) {
 //        for (int k = 1; k <= 3; k++) {
 //          for (int l = 1; l <= 3; l++) {
 //            cc3333(i,j,k,l) = tla*dde(i,j)*dde(k,l) + real_t(2.0)*tmu*xid4(i,j,k,l);
 //            printf("%d %d %d %d %g \n",i,j,k,l,cc3333(i,j,k,l) );
 //          }
 //        }
 //      }
 //    }
 //    chg_basis_4(cc66.pointer(), cc3333.pointer(), 4, 6, B_basis.pointer());

    int i1,i2,j1,j2;
	int ijv_[6*2] = {1,2,3,2,1,1,1,2,3,3,3,2};
	CMatrix <double> cc3333(3,3,3,3);
	ViewFMatrixKokkos <int> ijv(ijv_,6,2);
	for (int i = 1; i <= 6; i++) {
        i1 = ijv(i,1);
        i2 = ijv(i,2);
        for (int j = 1; j <= 6; j++) {
        	j1 = ijv(j,1);
        	j2 = ijv(j,2);
        	cc3333(i1,i2,j1,j2) = cc66v(i,j);
	        cc3333(i2,i1,j1,j2) = cc66v(i,j);
	        cc3333(i1,i2,j2,j1) = cc66v(i,j);
	        cc3333(i2,i1,j2,j1) = cc66v(i,j);
	    }
    }

    for (int i = 1; i <= 3; i++) {
    	for (int j = 1; j <= 3; j++) {
    		for (int k = 1; k <= 3; k++) {
    			for (int l = 1; l <= 3; l++) {
      				printf("%d %d %d %d %g \n", i,j,k,l,cc3333(i,j,k,l));
      			}
      		}
      	}
    }

    chg_basis_4(cc66.host.pointer(), cc3333.pointer(), 4, 6, B_basis.host.pointer());

    cc66.update_device();

    // for (int i = 1; i <= 6; i++) {
    //     printf("cc=%.8f %.8f %.8f %.8f %.8f %.8f\n", cc66(i,1),cc66(i,2),cc66(i,3),cc66(i,4),cc66(i,5),cc66(i,6));
    // }


	ur1.close();

	ur1.open("./ta_pl.sx");//("./cupl3.txt");
	// read pl txt file here, up to slip systems
	ur1 >> junk; CLEAR_LINE(ur1);
	ur1 >> icryst; CLEAR_LINE(ur1);
	for (int i = 1; i <= 3; i++) {
		ur1 >> cdim(i);
	}
	CLEAR_LINE(ur1);

	ur1 >> nmodesx; CLEAR_LINE(ur1);
	ur1 >> nmodes; CLEAR_LINE(ur1);
	// temporary check until general crystal file read in is implemented
	// if (nmodes > 1) {
	// 	printf("nmodes IN IS > 1\n");
	// 	printf("ONLY 1 SLIP MODE ALLOWED FOR NOW\n");
	//   	exit(1);
	// }
	ur1 >> junk;
	CLEAR_LINE(ur1);

	ur1 >> junk; CLEAR_LINE(ur1);
	ur1 >> modex >> nsmx >> nrsx >> gamd0x; CLEAR_LINE(ur1);
	ur1 >> tau0xf >> tau0xb >> tau1x >> thet0x >> thet1x; CLEAR_LINE(ur1);
	ur1 >> hselfx >> hlatex; CLEAR_LINE(ur1);
	ur1 >> junk; CLEAR_LINE(ur1);

	if (thet0x < thet1x) {
	  	printf("INITIAL HARDENING LOWER THAN FINAL HARDENING\n");
	  	exit(1);
	}
	//
	//  CASE TAU1=0 CORRESPONDS TO LINEAR HARDENING AND IS INDEPENDENT OF TAU0.
	//  AVOID DIVISION BY ZERO
	if (tau1x <= 1.0e-6) {
	  	tau1x = 1.0e-6;
	  	thet0x = thet1x;
	}

	for (int j = 1; j <= nsmx; j++) {
	  	for (int k = 1; k <= 3; k++) ur1 >> isn(j,k);
	  	for (int k = 1; k <= 3; k++) ur1 >> isb(j,k);
	  	CLEAR_LINE(ur1);
	}

	for (int js = 1; js <= nsmx; js++) {
	  	nsysx = js;
	  	//
	  	//   DEFINES RATE SENSITIVITY AND CRSS FOR EACH SYSTEM IN THE MODE
	  	//
	  	for (int m = 1; m <= 3; m++) {
	    	sn(m) = isn(js,m) / cdim(m);
	    	sb(m) = isb(js,m) * cdim(m);
	  	}
	    //
	    // *** NORMALIZES SYSTEM VECTORS AND CHECKS NORMALITY
	    //
	    snor = sqrt( sn(1)*sn(1) + sn(2)*sn(2) + sn(3)*sn(3) );
	    qnor = sqrt( sb(1)*sb(1) + sb(2)*sb(2) + sb(3)*sb(3) );
	    prod = 0.0;
	    for (int j = 1; j <= 3; j++) {
		dnca.host(j,nsysx) = sn(j) / snor;
	        dbca.host(j,nsysx) = sb(j) / qnor;
	        if (abs(dnca.host(j,nsysx)) < 1.e-03) dnca.host(j,nsysx) = 0.0;
	        if (abs(dbca.host(j,nsysx)) < 1.e-03) dbca.host(j,nsysx) = 0.0;
	        prod += dnca.host(j,nsysx) * dbca.host(j,nsysx);
	    }

	    if (prod >= 1.0e-3) {
			printf("SYSTEM %d IS NOT ORTHOGONAL !!\n", js);
			printf("%d %d %d,", isb(js,1), isb(js,2), isb(js,3));
			printf("%d %d %d", isn(js,1), isn(js,2), isn(js,3));
			exit(1);
		}      
		//
		//   DEFINE SCHMID VECTOR IN CRYSTAL AXES FOR EACH SYSTEM
		//
		//for (int i = 1; i <= 3; i++) {
		//	for (int j = 1; j <= 3; j++) {
		//		aux33(i,j) = (dnca(i,nsysx)*dbca(j,nsysx) + dnca(j,nsysx)*dbca(i,nsysx)) / 2.0;
		//	}
		//}
	 
	    //cb.chg_basis_2(aux5.pointer(), aux33.pointer(), 2, 5, cb.B_basis_host_pointer()); //heavy modification needed here

		//for (int i = 1; i <= 5; i++) {
	    //    schca(i,nsysx) = aux5(i);
	    //}
	} // end for js

	// second set of slip systems for Ta BCC
	ur1 >> junk; CLEAR_LINE(ur1);
	ur1 >> junk; CLEAR_LINE(ur1);
	ur1 >> junk; CLEAR_LINE(ur1);
	ur1 >> junk; CLEAR_LINE(ur1);
	ur1 >> junk; CLEAR_LINE(ur1);

	for (int j = nsmx+1; j <= nsmx*2; j++) {
	  	for (int k = 1; k <= 3; k++) ur1 >> isn(j,k);
	  	for (int k = 1; k <= 3; k++) ur1 >> isb(j,k);
	  	CLEAR_LINE(ur1);
	}
	for (int js = nsmx+1; js <= nsmx*2; js++) {
	  	nsysx = js;
	  	//
	  	//   DEFINES RATE SENSITIVITY AND CRSS FOR EACH SYSTEM IN THE MODE
	  	//
	  	for (int m = 1; m <= 3; m++) {
	    	sn(m) = isn(js,m) / cdim(m);
	    	sb(m) = isb(js,m) * cdim(m);
	  	}
	    //
	    // *** NORMALIZES SYSTEM VECTORS AND CHECKS NORMALITY
	    //
	    snor = sqrt( sn(1)*sn(1) + sn(2)*sn(2) + sn(3)*sn(3) );
	    qnor = sqrt( sb(1)*sb(1) + sb(2)*sb(2) + sb(3)*sb(3) );
	    prod = 0.0;
	    for (int j = 1; j <= 3; j++) {
		dnca.host(j,nsysx) = sn(j) / snor;
	        dbca.host(j,nsysx) = sb(j) / qnor;
	        if (abs(dnca.host(j,nsysx)) < 1.e-03) dnca.host(j,nsysx) = 0.0;
	        if (abs(dbca.host(j,nsysx)) < 1.e-03) dbca.host(j,nsysx) = 0.0;
	        prod += dnca.host(j,nsysx) * dbca.host(j,nsysx);
	    }

	    if (prod >= 1.0e-3) {
			printf("SYSTEM %d IS NOT ORTHOGONAL !!\n", js);
			printf("%d %d %d,", isb(js,1), isb(js,2), isb(js,3));
			printf("%d %d %d", isn(js,1), isn(js,2), isn(js,3));
			exit(1);
		}
	}
        dnca.update_device();
        dbca.update_device();

	nsmx = nsmx*2;
	ur1.close();

	//=====================================================================

    for (size_t elem_gid = 0; elem_gid<num_elems; elem_gid++) {
    // Should be able to parallelize this
    //FOR_ALL(elem_gid, 0, num_elems, {
    	//----Determine grain id of element----
    	const size_t rk_level = 1;

        int grain_id = 0;

    	// 0-35 Crystal elastic stiffness tensor
		int c = 0;
		for (int i = 1; i <= 6; i++) {
		  for (int j = 1; j <= 6; j++) {
		    state_vars(elem_gid,c) = cc66(i,j)*stress_conv;  // C_ij
		    c += 1;
		  }
		}

		// 36-45 Crystal plasticity parameters
		state_vars(elem_gid,36) = nsmx;  // nsm # slip systems
		state_vars(elem_gid,37) = nrsx;  // nrs rate sensitivity exponent
		state_vars(elem_gid,38) = gamd0x*1.e-6;  // gamd0 reference slip shear rate
		state_vars(elem_gid,39) = tau0xf*stress_conv;  // tau0f Voce parameter
		state_vars(elem_gid,40) = tau0xb*stress_conv;  // tau0b Voce parameter
		state_vars(elem_gid,41) = tau1x*stress_conv;  // tau1 Voce parameter
		state_vars(elem_gid,42) = thet0x*stress_conv;  // thet0 Voce parameter
		state_vars(elem_gid,43) = thet1x*stress_conv;  // thet1 Voce parameter
		state_vars(elem_gid,44) = hselfx;  // hself Self Hardening
		state_vars(elem_gid,45) = hlatex;  // hlate Latent Hardening

		// 46-117 Crystal slip normal vectors (BCC - 24 systems)
		c = 46;
		for (int i = 1; i <= nsmx; i++) {
			for (int j = 1; j <= 3; j++) {
				state_vars(elem_gid,c) = dnca(j,i); // slip normal component j for slip system i
				c += 1;
			}
		}
		// 118-189 Crystal slip direction vectors (BCC - 24 systems)
		for (int i = 1; i <= nsmx; i++) {
			for (int j = 1; j <= 3; j++) {
				state_vars(elem_gid,c) = dbca(j,i); // slip direction component j for slip system i
				c += 1;
			}
		}
		// 46-105 Crystal Schmid vectors [5 length] (FCC - 12 systems)
		//for (int i = 0; i < nsmx; i++) {
		//    for (int j = 0; j < 5; j++) {
		//        file_state_vars.host(mat_id,elem_gid,c) = schca(j,i);  // Schmid tensor component j for slip system i
		//        c += 1;
		//    }
		//}

		// 190-237 Variables to remember per slip system per element
		c = 190;
		for (int i = 0; i < nsmx; i++) {
			state_vars(elem_gid,c) = tau0xf*stress_conv;  // crss Critical resolved shear stress (forward) on slip system i
			c += 1;
		}
		for (int i = 0; i < nsmx; i++) {
			state_vars(elem_gid,c) = tau0xb*stress_conv;  // crss Critical resolved shear stress (backward) on slip system i
			c += 1;
		}
		// 238-255 Variables to remember per element
		for (int i = 238; i <= 246; i++) {
			state_vars(elem_gid,i) = 0.0;  // etot Total strain tensor 238-246 [9 terms]
		//} lets get fancy
		//for (int i = 151; i <= 159, i++) {
			state_vars(elem_gid,i+9) = 0.0;  // ept Plastic strain tensor 247-255 [9 terms]
		}
		state_vars(elem_gid,256) = 0.0;      // gacumgr Accumulated shear across all slip systems and time 

		// 161-163 Euler angles for each element. Only initial term for now that isn't flat
		//ur1 >> ph >> th >> om; CLEAR_LINE(ur1);
		// single crystal test
		ph = 0.0; th = 0.0; om = 0.0; //001
    	state_vars(elem_gid,257) = ph;
    	state_vars(elem_gid,258) = th;
    	state_vars(elem_gid,259) = om;

		// 260-313 Basis change array values
	    c = 260;
	    for (int i = 1; i <= 3; i++) {
    		for (int j = 1; j <= 3; j++) {
        		for (int k = 1; k <= 6; k++) {
        			state_vars(elem_gid,c) = B_basis(i,j,k);
        			c += 1;
      			}
    		}
  		}
  		state_vars(elem_gid,314) = grain_id;
  //});
	}
  Kokkos::fence();
  //state_vars.update_host();

    printf("user_mat_init completed\n");
    
	return;
}
