// *****************************************************************************
//     Harden_Routine    Harden_Routine  Harden_Routine  Harden_Routine    
//     SUBROUTINE UPDATE_CRSS_VOCE     --->      VERSION OF 23/APR/00
//
//     A VOCE LAW FUNCTION OF THE ACCUMULATED SHEAR IN EACH GRAIN IS ADDED
//     AS A MULTIPLICATIVE FACTOR THAT MODULATES THE ORIGINAL LINEAR HARDENING.
//     THE UNITS (AND THE STRENGTH) OF THE HARDENING ARE CARRIED BY THE
//     MULTIPLICATIVE FACTOR 'VOCE' (THE SAME FOR EVERY MODE).
//     THE SELF & LATENT COUPLING COEFFICIENTS 'HARD' ARE DIMENSIONLESS
//     CONSTANTS RELATIVE TO THE FACTOR 'VOCE'.
//*******************************************************************************


#include "evpfft.h"
#include "utilities.h"
#include "Profiler.h"
#include "vm.h"

void EVPFFT::harden(int imicro)
{
  Profiler profiler(__FUNCTION__);
  update_crss_voce();
}

void EVPFFT::update_crss_voce()
{
  Kokkos::parallel_for(
    Kokkos::MDRangePolicy<Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER>>({1,1,1}, {npts3+1,npts2+1,npts1+1}),
    KOKKOS_CLASS_LAMBDA(const int kk, const int jj, const int ii) {
    
    int iph;
    real_t gamtotx;
    real_t deltgam;
    real_t dtau;
    real_t tau0;
    real_t tau1;
    real_t thet0;
    real_t thet1;
    real_t tiny;
    real_t voce;
    real_t fact;
    real_t expini;
    real_t expdel;   
    
    real_t aux_[3*3];
    ViewMatrixTypeReal aux(aux_,3,3);

    iph = jphase(ii,jj,kk);

    if (igas(iph) == 0) {

      if (iJ2(iph) == 0) { 

        gamtotx = gacumgr(ii,jj,kk);
        deltgam = 0.0;
        for (int is = 1; is <= nsyst(iph); is++) {
          deltgam += ABS(gamdot(is,ii,jj,kk)) * tdot;
        }
  
        for (int is = 1; is <= nsyst(iph); is++) {
          dtau = 0.0;
          for (int js = 1; js <= nsyst(iph); js++) {
            dtau += hard(is,js,iph) * ABS(gamdot(js,ii,jj,kk)) * tdot;
          }
          tau0  = tau(is,1,iph);
          tau1  = tau(is,3,iph);
          thet0 = thet(is,1,iph);
          thet1 = thet(is,2,iph);        
          tiny = 1.0e-4*tau0;
  
          voce = 0.0;
          if (ABS(thet0) > tiny) {
            voce = thet1 * deltgam;
            if (ABS(tau1) > tiny) {
              fact = ABS(thet0/tau1);
              expini = EXP(-gamtotx*fact);
              expdel = EXP(-deltgam*fact);
              voce += -(fact*tau1-thet1)/fact*expini*(expdel-1.0) - 
                        thet1/fact*expini*(expdel*((gamtotx+deltgam)*fact+1.0)-(gamtotx*fact+1.0));
            } // end if (ABS(tau1) > tiny)
          } // end if (ABS(thet0) > tiny)
  
          if (ABS(deltgam) > 0.0) {
            crss(is,1,ii,jj,kk) += dtau*voce/deltgam;
            crss(is,2,ii,jj,kk) += dtau*voce/deltgam;
          }
  
          //trialtau(is,1,ii,jj,kk) = crss(is,1,ii,jj,kk);
          //trialtau(is,2,ii,jj,kk) = crss(is,2,ii,jj,kk);
        } // end for is
  
        gacumgr(ii,jj,kk) = gamtotx + deltgam;

      } else { // else for if (iJ2(jph) == 0)

        for (int i = 1; i <= 3; i++) {
          for (int j = 1; j <= 3; j++) {
            aux(i,j) = edotp(i,j,ii,jj,kk);
          }
        }

        gacumgr(ii,jj,kk) = gacumgr(ii,jj,kk) + vm(aux.pointer())*tdot;

        sigma0gr(ii,jj,kk) = sigma0(iph) + (sigma1(iph) + thet1_j2(iph)*gacumgr(ii,jj,kk))*
          (1.0 - exp(- gacumgr(ii,jj,kk)*ABS(thet0_j2(iph)/sigma1(iph))));

      } // end for if (iJ2(jph) == 0)

    } // end if (igas(iph) == 0)

  });

}

