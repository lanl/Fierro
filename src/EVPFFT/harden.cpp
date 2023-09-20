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


real_t specific_heat_tantalum(real_t T);
real_t specific_heat_copper(real_t T);

void EVPFFT::harden(int imicro)
{
  update_crss_voce();
  if (itemphard == 1) {
    update_temperature();
    update_crss_temp(); // apply softerning due to temperature
  }
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

    iph = jphase(ii,jj,kk);

    if (igas(iph) == 0) {

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

    } // end if (igas(iph) == 0)

  });

}

void EVPFFT::update_crss_temp()
{
  Kokkos::parallel_for(
    Kokkos::MDRangePolicy<Kokkos::Rank<3,LOOP_ORDER,LOOP_ORDER>>({1,1,1}, {npts3+1,npts2+1,npts1+1}),
    KOKKOS_CLASS_LAMBDA(const int kk, const int jj, const int ii) {

    int iph;
    real_t tau0;
    iph = jphase(ii,jj,kk);

    if (igas(iph) == 0) {
      for (int is = 1; is <= nsyst(iph); is++) {
        tau0 = tau0_mode(is,1,iph) *
               (exp(-(temp - tau0_mode(is,3,iph)) / tau0_mode(is,2,iph)) -
               exp(-(tempold - tau0_mode(is,3,iph)) / tau0_mode(is,2,iph)));
        crss(is,1,ii,jj,kk) += tau0;
        crss(is,2,ii,jj,kk) += tau0;
      } // end for is
    } // end if (igas(iph) == 0)
  }); 
}

void EVPFFT::update_temperature()
{
    real_t rho, eta, cp;
    real_t devmp, dwplastic, dtemp;

//#if 0
    // tantalum
    eta = 0.9; // 1 default for poly tantalum, adjusted based on number of active slips;
    rho = 16640; // kg/m^3;
    cp = specific_heat_tantalum(temp);
//#endif
#if 0
    // copper
    eta = 0.9;
    rho = 8960; // kg/m^3
    cp  = specific_heat_copper(temp);
#endif

    devmp = evmp - evmpold;
    evmpold = evmp;

    dwplastic = svm * devmp;
    wplastic +=  dwplastic;

    dtemp = eta * dwplastic * 1e6 / (rho * cp);

    tempold = temp;
    temp += dtemp;
}

real_t specific_heat_tantalum(real_t T)
{
  /*
  * Chen, S., Gray, G., Bingert, S.R. 1996. 
  * Mechanical properties and constitutive relations for tantalum and tantalum alloys under high-rate deformation. 
  * Los Alamos National Laboratory Report LA-UR-96-0602. doi: http://dx.doi.org/10.2172/226058.
  */

  real_t cpa0 = 145.5; // j/kgk;
  real_t cpa1 = 0.009544; // j/kgk^2;
  real_t cpa2 = -68900; // jk/kg;
 
  real_t cp = cpa0 + cpa1 * T + cpa2 / (T * T);
  return cp;
}

real_t specific_heat_copper(real_t T)
{
  /*
  * Banerjee, B. (2005).
  * An evaluation of plastic flow stress models for the simulation of high-temperature and high-strain-rate deformation of metals.
  * arXiv preprint cond-mat/0512466.
  */

  real_t cp;
  if (T < 270) {
    cp = 0.0000416*pow(T,3) - 0.027*pow(T,2) + 6.21*T - 142.6;
  }
  else {
    cp = 0.1009*T + 358.4;
  }
  return cp;
}
