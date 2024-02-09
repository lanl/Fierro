#include "evpfft.h"
#include "voigt.h"
#include "definitions.h"
#include "utilities.h"

using namespace utils;


void EVPFFT::data_crystal_elast(int iph, const std::string & filecrysel)
{
  std::ifstream ur1;
  ur1.open(filecrysel);
  check_that_file_is_open(ur1, filecrysel.c_str());

  MatrixTypeRealHost dde(3,3);
  MatrixTypeRealHost xid4(3,3,3,3);
  MatrixTypeRealHost cc66v(6,6);
  MatrixTypeRealHost ccaux(3,3,3,3);
  MatrixTypeRealHost aux6(6);
  MatrixTypeRealHost aux33(3,3);
  MatrixTypeRealHost aux66(6,6);
  MatrixTypeRealHost aux3333(3,3,3,3); 
  int iso;
  real_t young,tnu,tmu,tla;
  
  //       UNITARY TENSORS
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      dde(i,j) = 0.0;
	    if (i == j) dde(i,j) = 1.0;
    }
  }

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      for (int k = 1; k <= 3; k++) {
        for (int l = 1; l <= 3; l++) {
	        xid4(i,j,k,l) = (dde(i,k)*dde(j,l) + dde(i,l)*dde(j,k)) / real_t(2.0);
        }
      }
    }
  }

  ur1 >> iso; CLEAR_LINE(ur1);

	if (iso == 0) {
    for (int i = 1; i <= 6; i++) {
      for (int j = 1; j <= 6; j++) {
        ur1 >> cc66v(i,j);
      }
      CLEAR_LINE(ur1);
    }

    voigt(cc66v.pointer(),ccaux.pointer(),1);

    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
        for (int k = 1; k <= 3; k++) {
          for (int l = 1; l <= 3; l++) {
            cc(i,j,k,l,iph) = ccaux(i,j,k,l);
          }
        }
      }
    }

  } else {

  	ur1 >> young >> tnu; CLEAR_LINE(ur1);
    tmu = young / (real_t(2.0)*(real_t(1.0) + tnu));
	  tla = real_t(2.0)*tmu*tnu / (real_t(1.0) - real_t(2.0)*tnu);

    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
        for (int k = 1; k <= 3; k++) {
          for (int l = 1; l <= 3; l++) {
            cc(i,j,k,l,iph) = tla*dde(i,j)*dde(k,l) + real_t(2.0)*tmu*xid4(i,j,k,l);
          }
        }
      }
    }

  } // end if (iso == 0)


  // scale stiffness tensor
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      for (int k = 1; k <= 3; k++) {
        for (int l = 1; l <= 3; l++) {
          cc(i,j,k,l,iph) *= stress_scale;
        }
      }
    }
  } 

  ur1.close();
}
