#include <iostream>

#include "evpfft.h"
#include "definitions.h"
#include "utilities.h"

using namespace utils;
      
void EVPFFT::data_crystal(int iph, const std::string & filecryspl)
{

  std::ifstream ur1;
  ur1.open(filecryspl);
  check_that_file_is_open(ur1, filecryspl.c_str());
  std::string prosa;

  MatrixTypeIntHost isn(NSYSMX,4);
  MatrixTypeIntHost isb(NSYSMX,4);
  MatrixTypeIntHost mode(10);
  MatrixTypeRealHost sn(3);
  MatrixTypeRealHost sb(3);
  MatrixTypeRealHost cdim(3);
  MatrixTypeRealHost aux5(5);
  MatrixTypeRealHost aux33(3,3);
  MatrixTypeRealHost aux55(5,5);
  MatrixTypeRealHost aux3333(3,3,3,3);
  MatrixTypeRealHost hselfx(10);
  MatrixTypeRealHost hlatex(10,10);
  int nmodesx,kount,nm,modex,nsmx,nrsx;
  int isectwx,nsysx,nrs_j2x;
  real_t covera,gamd0x,twshx,tau0xf,tau0xb,tau1x,thet0x,thet1x;
  real_t snor,qnor,prod;
  real_t edotp0_j2x,sigma0x,sigma1x,thet0_j2x,thet1_j2x;

  ur1 >> prosa; CLEAR_LINE(ur1);
  ur1 >> icryst(iph); CLEAR_LINE(ur1);
  if (icryst(iph).substr(0,2) == "J2" || icryst(iph).substr(0,2) == "j2") {

    iJ2.host(iph) = 1;

    ur1 >> nrs_j2x >> edotp0_j2x; CLEAR_LINE(ur1);
    ur1 >> sigma0x >> sigma1x >> thet0_j2x >> thet1_j2x; CLEAR_LINE(ur1);

    nrs_j2.host(iph) = nrs_j2x;
    edotp0_j2.host(iph) = edotp0_j2x;
    sigma0.host(iph) = sigma0x;
    sigma1.host(iph) = sigma1x;
    thet0_j2.host(iph) = thet0_j2x;
    thet1_j2.host(iph) = thet1_j2x;

  } else { // else for if (icryst(iph).substr(0,2) == "J2" || icryst(iph).substr(0,2) == "j2") {

    iJ2.host(iph) = 0;

    for (int i = 1; i <= 3; i++) {
      ur1 >> cdim(i);
    }
    CLEAR_LINE(ur1);
    covera = cdim(3)/cdim(1);
  
    ur1 >> nmodesx; CLEAR_LINE(ur1);
    ur1 >> nmodes(iph); CLEAR_LINE(ur1);
    for (int i = 1; i <= mode.size(); i++) mode(i) = 0;
    for (int i = 1; i <= nmodes(iph); i++) {
      ur1 >> mode(i);
    }
    CLEAR_LINE(ur1);
  
    ntwmod(iph) = 0;
    nsyst.host(iph)  = 0;
    ntwsys(iph) = 0;
    kount = 1;
  
    //
    //     START READING DEFORMATION MODES FROM FILECRYS
    //
    //nm = 0;
    //label_100: {
    //  do {
    //    nm++;
    for (int nm = 1; nm <= nmodesx; nm++) {
  
        ur1 >> prosa; CLEAR_LINE(ur1);
        ur1 >> modex >> nsmx >> nrsx >> gamd0x >> twshx >> isectwx; CLEAR_LINE(ur1);
        ur1 >> tau0xf >> tau0xb >> tau1x >> thet0x >> thet1x; CLEAR_LINE(ur1);
        ur1 >> hselfx(nm);
        for (int jm = 1; jm <= nmodesx; jm++) {
          ur1 >> hlatex(nm,jm);
        }
        CLEAR_LINE(ur1);
  
        // scale variables
        gamd0x /= time_scale;
        tau0xf  *= stress_scale;
        tau0xb *= stress_scale;
        tau1x *= stress_scale;
        thet0x *= stress_scale;
        thet1x *= stress_scale;
  
        //     SKIPS nsmx LINES IF THE MODE IS NOT IN THE LIST.
        if (modex != mode(kount)) {
          for (int iz = 1; iz <= nsmx; iz++) {
            CLEAR_LINE(ur1);
          } 
          //goto label_100;
          continue;
        }
  
        if (thet0x < thet1x) {
          printf("INITIAL HARDENING LOWER THAN FINAL HARDENING FOR MODE %d IN PHASE %d\n", kount, iph);
          exit(1);
        }
  
        //
        //  CASE TAU1=0 CORRESPONDS TO LINEAR HARDENING AND IS INDEPENDENT OF TAU0.
        //  AVOID DIVISION BY ZERO
        if (tau1x <= 1.0e-6) {
          tau1x = 1.0e-6;
          thet0x = thet1x;
        }
        
        // REORDER HARDENING COEFFICIENTS TO ACCOUNT ONLY FOR ACTIVE MODES
        hselfx(kount) = hselfx(nm);
        for (int i = 1; i <= nmodes(iph); i++) {
          hlatex(kount,i) = hlatex(nm,mode(i));
        }
  
        // SYSTEMS GIVEN IN FOUR INDEX NOTATION: HEXAGONALS AND TRIGONALS
        // SYSTEMS GIVEN IN THREE INDEX NOTATION: CUBIC AND ORTHORHOMBIC
        if (icryst(iph).substr(0,3) == "HEX" || icryst(iph).substr(0,3) == "TRI" ||
            icryst(iph).substr(0,3) == "hex" || icryst(iph).substr(0,3) == "tri") {
          for (int j = 1; j <= nsmx; j++) {
            for (int k = 1; k <= 4; k++) ur1 >> isn(j,k);
            for (int k = 1; k <= 4; k++) ur1 >> isb(j,k);
            CLEAR_LINE(ur1);
          }
        } else if (icryst(iph).substr(0,3) == "CUB" || icryst(iph).substr(0,3) == "ORT" ||
                   icryst(iph).substr(0,3) == "cub" || icryst(iph).substr(0,3) == "ort") {
          for (int j = 1; j <= nsmx; j++) {
            for (int k = 1; k <= 3; k++) ur1 >> isn(j,k);
            for (int k = 1; k <= 3; k++) ur1 >> isb(j,k);
            CLEAR_LINE(ur1);
          }
        } else {
          printf(" cannot identify the crystal symmetry of phase  %d\n", iph);
          exit(1);
        }
  
        nsm(kount,iph) = nsmx;
        if (twshx != 0) ntwmod(iph) = ntwmod(iph) + 1;
  
#ifdef NON_SCHMID_EFFECTS
      //===============non_schmid_effects ===================================    
      // NOTE: the reading sequence here might be wrong.
      // But was duplicated form the original fortran code verbatim. 
      for (int kmo = 1; kmo <= nmodes(iph); kmo++) {
        CLEAR_LINE(ur1);
        ur1 >> cns(1,kmo,iph) >> cns(2,kmo,iph) 
               >> cns(3,kmo,iph) >> cns(4,kmo,iph);
        // FIXME: the commented out line below leads to use of uninitialised cns(_,_,_)
        // in the later part of this function. But this is the way it appears
        // in the original fortran code. It was modified as below. 
        //cns(5,1,1) = -( cns(3,kmo,iph) + cns(4,kmo,iph) );
        cns(5,kmo,iph) = -( cns(3,kmo,iph) + cns(4,kmo,iph) );
      }
      //=====================================================================
#endif
  
        for (int js = 1; js <= nsm(kount,iph); js++) {
          nsyst.host(iph) = nsyst.host(iph) + 1;
          nsysx = nsyst.host(iph);
          if (twshx != 0) ntwsys(iph) = ntwsys(iph) + 1;
  
          //
          //   DEFINES RATE SENSITIVITY AND CRSS FOR EACH SYSTEM IN THE MODE
          //
          gamd0.host(nsysx,iph)  = gamd0x;
          nrs.host(nsysx,iph)    = nrsx;
          twsh(nsysx,iph)   = twshx;
          tau.host(nsysx,1,iph)  = tau0xf;
          tau.host(nsysx,2,iph)  = tau0xb;
          tau.host(nsysx,3,iph)  = tau1x;
          thet.host(nsysx,1,iph) = thet0x;
          thet.host(nsysx,2,iph) = thet1x;
  
          //
          isectw(nsysx,iph) = isectwx;
          //
  
          if (icryst(iph).substr(0,3) == "HEX" || icryst(iph).substr(0,3) == "TRI" ||
              icryst(iph).substr(0,3) == "hex" || icryst(iph).substr(0,3) == "tri") {
            sn(1) = isn(js,1);
            sn(2) = (isn(js,1)+2.0*isn(js,2))/sqrt(3.0);
            sn(3) = isn(js,4)/covera;
            sb(1) = 3.0/2.0*isb(js,1);
            sb(2) = (isb(js,1)/2.0+isb(js,2))*sqrt(3.0);
            sb(3) = isb(js,4)*covera;
          } else if (icryst(iph).substr(0,3) == "CUB" || icryst(iph).substr(0,3) == "ORT" ||
                     icryst(iph).substr(0,3) == "cub" || icryst(iph).substr(0,3) == "ort") {
            for (int m = 1; m <= 3; m++) {
              sn(m) = isn(js,m) / cdim(m);
              sb(m) = isb(js,m) * cdim(m);
            }
          }
  
          //
          // *** NORMALIZES SYSTEM VECTORS AND CHECKS NORMALITY
          //
          snor = sqrt( sn(1)*sn(1) + sn(2)*sn(2) + sn(3)*sn(3) );
          qnor = sqrt( sb(1)*sb(1) + sb(2)*sb(2) + sb(3)*sb(3) );
          prod = 0.0;
          for (int j = 1; j <= 3; j++) {
            dnca.host(j,nsysx,iph) = sn(j) / snor;
            dbca.host(j,nsysx,iph) = sb(j) / qnor;
            if (abs(dnca.host(j,nsysx,iph)) < 1.e-03) dnca.host(j,nsysx,iph) = 0.0;
            if (abs(dbca.host(j,nsysx,iph)) < 1.e-03) dbca.host(j,nsysx,iph) = 0.0;
            prod += dnca.host(j,nsysx,iph) * dbca.host(j,nsysx,iph);
          }
  
          if (prod >= 1.0e-3) {
            printf("SYSTEM %d IN MODE %d IN PHASE %d IS NOT ORTHOGONAL !!\n", js, nm, iph);
            printf("%d %d %d", isb(js,1), isb(js,2), isb(js,3));
            printf("%d %d %d", isn(js,1), isn(js,2), isn(js,3));
            exit(1);
          }
  
#ifdef NON_SCHMID_EFFECTS
          dtca(1,nsysx,iph) =  dnca.host(2,nsysx,iph)*dbca.host(3,nsysx,iph) - dnca.host(3,nsysx,iph)*dbca.host(2,nsysx,iph);
          dtca(2,nsysx,iph) = -dnca.host(1,nsysx,iph)*dbca.host(3,nsysx,iph) + dnca.host(3,nsysx,iph)*dbca.host(1,nsysx,iph);
          dtca(3,nsysx,iph) =  dnca.host(1,nsysx,iph)*dbca.host(2,nsysx,iph) - dnca.host(2,nsysx,iph)*dbca.host(1,nsysx,iph);
#endif
             
          //
          //   DEFINE SCHMID VECTOR IN CRYSTAL AXES FOR EACH SYSTEM
          //
  
          for (int i = 1; i <= 3; i++) {
            for (int j = 1; j <= 3; j++) {
              aux33(i,j) = (dnca.host(i,nsysx,iph)*dbca.host(j,nsysx,iph) + dnca.host(j,nsysx,iph)*dbca.host(i,nsysx,iph)) / 2.0;
            }
          }
   
          cb.chg_basis_2(aux5.pointer(), aux33.pointer(), 2, 5, cb.B_basis_host_pointer());
  
          for (int i = 1; i <= 5; i++) {
            schca.host(i,nsysx,iph) = aux5(i);
          }
  
#ifdef NON_SCHMID_EFFECTS
          for (int i = 1; i <= 3; i++) {
            for (int j = 1; j <= 3; j++) {
              aux33(i,j) = cns(1,nm,iph)*(dtca(i,nsysx,iph)*dbca.host(j,nsysx,iph) + 
                          dtca(j,nsysx,iph)*dbca.host(i,nsysx,iph))/2.0 +
                          cns(2,nm,iph)*(dtca(i,nsysx,iph)*dnca.host(j,nsysx,iph) +
                          dtca(j,nsysx,iph)*dnca.host(i,nsysx,iph))/2.0 +
                          cns(3,nm,iph)*(dnca.host(i,nsysx,iph)*dnca.host(j,nsysx,iph) +
                          dnca.host(j,nsysx,iph)*dnca.host(i,nsysx,iph))/2.0 +
                          cns(4,nm,iph)*(dtca(i,nsysx,iph)*dtca(j,nsysx,iph) +
                          dtca(j,nsysx,iph)*dtca(i,nsysx,iph))/2.0 +
                          cns(5,nm,iph)*(dbca.host(i,nsysx,iph)*dbca.host(j,nsysx,iph) +
                          dbca.host(j,nsysx,iph)*dbca.host(i,nsysx,iph))/2.0 +
                          (dnca.host(i,nsysx,iph)*dbca.host(j,nsysx,iph) + 
                          dnca.host(j,nsysx,iph)*dbca.host(i,nsysx,iph))/2.0;
            }
          }
  
          cb.chg_basis_2(aux5.pointer(), aux33.pointer(), 2, 5, cb.B_basis_host_pointer());
  
          for (int i = 1; i <= 5; i++) {
            schcnon.host(i,nsysx,iph) = aux5(i);
          }
#endif
  
        } // end for js
  
        kount++;
      //} while (nm < nmodesx); // end of do/while loop
    //} // end of label_100 scope
    } // end for (nm)
  
    //     INITIALIZE SELF & LATENT HARDENING COEFS FOR EACH SYSTEM OF THE PHASE.
    //     ABSOLUTE UNITS ARE ACCOUNTED FOR BY MODULATING FACTOR IN HARDENING LAW.
  
    int i = 0;
    int j;
    for (int im = 1; im <= nmodes(iph); im++) {
      for (int is = 1; is <= nsm(im,iph); is++) {
        i += 1;
        j = 0;
        for (int jm = 1; jm <= nmodes(iph); jm++) {
          for (int js = 1; js <= nsm(jm,iph); js++) {
            j += 1;
            hard.host(i,j,iph) = hlatex(im,jm);
          }
        }
        hard.host(i,i,iph) = hselfx(im);
      }
    }
  
    //     VERIFICATION OF TWINNING DATA TO BE SURE PROGRAM WILL RUN PROPERLY
    if (nmodes(iph) > 1) {
      for (int i = 2; i <= nsyst.host(iph); i++) {
        if (twsh(i,iph) == 0.0 && twsh(i-1,iph) != 0.0) {
          printf(" WARNING! THE TWINNING MODES MUST FOLLOW THE \n");
          printf(" SLIP MODES   -->   REORDER CRYSTAL FILE ");
          exit(1);
        }
      }
    }
  } // end for if (icryst(iph).substr(0,2) == "J2" || icryst(iph).substr(0,2) == "j2") {

  ur1.close();  
}
