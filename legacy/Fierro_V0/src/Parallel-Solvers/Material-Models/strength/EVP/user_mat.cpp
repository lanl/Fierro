#include "user_mat.h"

// using namespace to avoid collision of functions that are the same in other models like EVPFFT or LS-EVPFFT
namespace EVP
{

// euler <-> rotation matrix transformation function
KOKKOS_FUNCTION
void euler(int iopt, real_t &ph, real_t &th, real_t &tm, real_t *a_){
    ViewCMatrixKokkos <double> a(a_,3,3);
    real_t sth,sph,cph,cth,stm,ctm; 
    const real_t pi = 4.*atan(1.0);

    //  CALCULATE THE EULER ANGLES ASSOCIATED WITH THE TRANSFORMATION
    //  MATRIX A(I,J) IF IOPT=1 AND VICEVERSA IF IOPT=2
    //  A(i,j) TRANSFORMS FROM SYSTEM sa TO SYSTEM ca.
    //  ph,th,om ARE THE EULER ANGLES OF ca REFERRED TO sa.

    if (iopt == 1) {
      th = acos(a(3,3));
      if(fabs(a(3,3)) >= 0.9999) {
        tm = 0.0;
        ph = atan2(a(1,2),a(1,1));
      } else {
        sth = sin(th);
        tm  = atan2(a(1,3)/sth,a(2,3)/sth);
        ph  = atan2(a(3,1)/sth,-a(3,2)/sth);
      }
      // th = th*180.0/pi;
      // ph = ph*180.0/pi;
      // tm = tm*180.0/pi;
    } else if (iopt == 2) {
      sph = sin(ph);
      cph = cos(ph);
      sth = sin(th);
      cth = cos(th);
      stm = sin(tm);
      ctm = cos(tm);
      a(1,1) =  ctm*cph-sph*stm*cth;
      a(2,1) = -stm*cph-sph*ctm*cth;
      a(3,1) =  sph*sth;
      a(1,2) =  ctm*sph+cph*stm*cth;
      a(2,2) = -sph*stm+cph*ctm*cth;
      a(3,2) = -sth*cph;
      a(1,3) =  sth*stm;
      a(2,3) =  ctm*sth;
      a(3,3) =  cth;
    } // end if (iopt == 1)
} // end of euler
// update a 3x3 strain tensor
KOKKOS_FUNCTION
void update_strain(real_t *strain_, real_t *strain_n_, real_t *rate_, const real_t dt){
    ViewCMatrixKokkos <double> strain_n(strain_n_,3,3);
    ViewCMatrixKokkos <double> strain(strain_,3,3);
    ViewCMatrixKokkos <double> rate(rate_,3,3);
    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        strain(ii,jj) = strain_n(ii,jj) + rate(ii,jj) * dt;
      } // end for ii
    } // end for jj
} // end update strain
// rotate crystal elastic stiffness tensor with local rotation
KOKKOS_FUNCTION
void update_stiff(real_t *cg66_, real_t *cc66_, real_t *aa_, real_t *B_basis_){
    real_t caux66_[6*6];
    real_t caux3333_[3*3*3*3];

    ViewCMatrixKokkos <double> cg66(cg66_,6,6);
    ViewCMatrixKokkos <double> cc66(cc66_,6,6);
    ViewCMatrixKokkos <double> aa(aa_,3,3);
    ViewCMatrixKokkos <double> B_basis(B_basis_,3,3,6);
    ViewCMatrixKokkos <double> caux66(caux66_,6,6);
    ViewCMatrixKokkos <double> cc3333(caux3333_,3,3,3,3);
    ViewCMatrixKokkos <double> caux3333(caux3333_,3,3,3,3);

    for (int i = 1; i <= 6; i++) {
      for (int j = 1; j <= 6; j++) {
        caux66(i,j) = cc66(i,j);
      }
    }

    chg_basis_3(caux66.pointer(), cc3333.pointer(), 3, 6, B_basis.pointer());

    for (int i1 = 1; i1 <= 3; i1++) {
      for (int j1 = 1; j1 <= 3; j1++) {
        for (int k1 = 1; k1 <= 3; k1++) {
          for (int l1 = 1; l1 <= 3; l1++) {      
            double dum = 0.0;
            for (int i2 = 1; i2 <= 3; i2++) {
              for (int j2 = 1; j2 <= 3; j2++) {
                for (int k2 = 1; k2 <= 3; k2++) {
                  for (int l2 = 1; l2 <= 3; l2++) {        
                    dum += aa(i2,i1) *
                           aa(j2,j1) *
                           aa(k2,k1) *
                           aa(l2,l1) *
                           cc3333(i2,j2,k2,l2);         
                  }
                }
              }
            }              
            caux3333(i1,j1,k1,l1) = dum;
          }
        }
      }
    }

    chg_basis_4(cg66.pointer(), caux3333.pointer(), 4, 6, B_basis.pointer());

} // end update stiffness
// rotate crystal schmid tensor with local rotation
KOKKOS_FUNCTION
void update_schmid(real_t *sch_, real_t *schca_, real_t *aa_, const int nsm, real_t *B_basis_){

    real_t aux5_[5];
    real_t aux33_[3*3]; 
    real_t aux33r_[3*3];
    ViewFMatrixKokkos <double> sch(sch_,5,nsm);
    ViewFMatrixKokkos <double> schca(schca_,5,nsm);
    ViewCMatrixKokkos <double> aa(aa_,3,3);
    ViewCMatrixKokkos <double> B_basis(B_basis_,3,3,6);
    ViewCMatrixKokkos <double> aux5    (aux5_,5);
    ViewCMatrixKokkos <double> aux33   (aux33_,3,3);
    ViewCMatrixKokkos <double> aux33r  (aux33r_,3,3);
    
    // EL: Maybe don't need this loop in this version
    for (int is = 1; is <= nsm; is++) {
      for (int j = 1; j <= 5; j++) {
        sch(j,is) = 0.0;
      } // end for j
    } // end for is

    for (int is = 1; is <= nsm; is++) {
      for (int j = 1; j <= 5; j++) {
        aux5(j) = schca(j,is);
      } // end for j

      chg_basis_1(aux5.pointer(), aux33.pointer(), 1, 5, B_basis.pointer());

      for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
          aux33r(i,j) = 0.0;
          for (int i1 = 1; i1 <= 3; i1++) {
            for (int j1 = 1; j1 <= 3; j1++) {
              aux33r(i,j) += aa(i,i1)*aa(j,j1)*aux33(i1,j1);
            } // end for j1
          } // end for i1
        } // end for j
      } // end for i

      chg_basis_2(aux5.pointer(), aux33r.pointer(), 2, 5, B_basis.pointer());

      for (int j = 1; j <= 5; j++) {
        sch(j,is) = aux5(j);
      } // end for j

    } // end for is 
} // end update schmid
// update crystal orientation with calculated body and plastic rotations
KOKKOS_FUNCTION
void update_orient(real_t *ag_, real_t *vel_grad_, real_t *gamdot_, real_t *dnca_, real_t *dbca_, const int nsm, const real_t dt){
    // thread private arrays
    real_t aa_[3*3];
    real_t distor_[3*3];
    real_t dnsa_[3];
    real_t dbsa_[3];
    real_t rotslip_[3*3];
    real_t rotloc_[3*3];
    real_t rot_[3*3];

    ViewCMatrixKokkos <double> ag(ag_,3,3);
    ViewCMatrixKokkos <double> vel_grad(vel_grad_,3,3);
    ViewCMatrixKokkos <double> gamdot(gamdot_,nsm);
    ViewFMatrixKokkos <double> dnca(dnca_,3,nsm);
    ViewFMatrixKokkos <double> dbca(dbca_,3,nsm);
    // create views of thread private arrays
    ViewCMatrixKokkos <double> aa(aa_,3,3);
    ViewCMatrixKokkos <double> distor(distor_,3,3);
    ViewCMatrixKokkos <double> dnsa(dnsa_,3);
    ViewCMatrixKokkos <double> dbsa(dbsa_,3);
    ViewCMatrixKokkos <double> rotslip(rotslip_,3,3);
    ViewCMatrixKokkos <double> rotloc(rotloc_,3,3);
    ViewCMatrixKokkos <double> rot(rot_,3,3);

    //
    //      LOCAL ROTATION RATE: ANTISYM(VELGRAD)
    //
    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        rotloc(ii,jj) = (vel_grad(ii,jj) - vel_grad(jj,ii)) / 2.0;
      }
    }
    //
    //      SLIP ROTATION RATE
    //
    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        aa(ii,jj) = ag(ii,jj);
        distor(ii,jj) = 0.0;
      }
    }

    for (int is = 1; is <= nsm; is++) {
      for (int ii = 1; ii <= 3; ii++) {
        dnsa(ii) = 0.0;
        dbsa(ii) = 0.0;
        for (int jj = 1; jj <= 3; jj++) {
          dnsa(ii) += aa(ii,jj) * dnca(jj,is);
          dbsa(ii) += aa(ii,jj) * dbca(jj,is);
        }
      }

      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          distor(ii,jj) += dbsa(ii) * dnsa(jj) * gamdot(is);
        }
      }
    } // end for is

    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        rotslip(ii,jj) = (distor(ii,jj)-distor(jj,ii)) / 2.0;
      }
    }
    //
    //      TOTAL ROTATION
    //
    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        //rot(ii,jj) = (fd.tomtot(ii,jj)+rotloc(ii,jj)-rotslip(ii,jj)) * dt; // No total body rotation tomtot (yet?)
        rot(ii,jj) = (rotloc(ii,jj)-rotslip(ii,jj)) * dt;
      }
    }

    //
    //      REORIENTATION
    //
    orient(aa.pointer(), rot.pointer());

    //
    //      UPDATE ORIENTATION MATRIX
    //
    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        ag(ii,jj) = aa(ii,jj);
      }
    }
} // end update orientation
// Update orientation with rotation
KOKKOS_FUNCTION
void orient(real_t *a_, real_t *c_){
    real_t snorm;
    real_t snorm1;
    // thread private arrays
    real_t th2_[3*3];
    real_t v_[3];
    real_t vbar_[3];
    real_t th_[3*3];
    real_t rot_[3*3];
    real_t anew_[3*3];

    ViewCMatrixKokkos <double> a(a_,3,3);
    ViewCMatrixKokkos <double> c(c_,3,3);
    // create views of thread private arrays
    ViewCMatrixKokkos <double> th2(th2_,3,3);
    ViewCMatrixKokkos <double> v(v_,3);
    ViewCMatrixKokkos <double> vbar(vbar_,3);
    ViewCMatrixKokkos <double> th(th_,3,3);
    ViewCMatrixKokkos <double> rot(rot_,3,3);
    ViewCMatrixKokkos <double> anew(anew_,3,3);

    // BUILD ROTATION TENSOR BASED ON RODRIGUES FORMULA
    v(1) = c(3,2);
    v(2) = c(1,3);
    v(3) = c(2,1);

    snorm = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3));
    snorm1 = tan(snorm/2.0);
    if (snorm <= 1.0E-06) {
      snorm = 1.0;
    }

    for (int i = 1; i <= 3; i++) {
      vbar(i) = snorm1*v(i)/snorm;
    }

    snorm = vbar(1)*vbar(1)+vbar(2)*vbar(2)+vbar(3)*vbar(3);
    th(3,2) =  vbar(1);
    th(1,3) =  vbar(2);
    th(2,1) =  vbar(3);
    th(2,3) = -vbar(1);
    th(3,1) = -vbar(2);
    th(1,2) = -vbar(3);

    for (int i = 1; i <= 3; i++) {
      th(i,i) = 0.0;
    }

    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
        th2(i,j) = 0.0;
        for (int k = 1; k <= 3; k++) {
          th2(i,j) += th(i,k)*th(k,j);
        }
      }
    }

    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
        rot(i,j) = (i/j)*(j/i)+2.0*(th(i,j)+th2(i,j))/(1.0+snorm);
      }
    }

    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
        anew(i,j) = 0.0;
        for (int k = 1; k <= 3; k++) {
          anew(i,j) += rot(i,k)*a(k,j);
        }
      }
    }

    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
        a(i,j) = anew(i,j);
      }
    }
} // end orient
// Update critical resolved shear stress on each slip system
KOKKOS_FUNCTION
void harden(real_t *crss_, real_t &gacumgr, real_t *gamdot_, real_t *voceParam_, real_t *hard_, const int nsm, const real_t dt){//(int imicro)
    real_t gamtotx;
    real_t deltgam;
    real_t tau0;
    real_t tau1;
    real_t thet0;
    real_t thet1;
    real_t dtau;
    real_t tiny;
    real_t voce;
    real_t fact;
    real_t expini;
    real_t expdel;   

    ViewCMatrixKokkos <double> crss(crss_,2,nsm);
    ViewCMatrixKokkos <double> gamdot(gamdot_,nsm);
    ViewCMatrixKokkos <double> voceParam(voceParam_,5);
    ViewCMatrixKokkos <double> hard(hard_,nsm,nsm);

    gamtotx = gacumgr;
    deltgam = 0.0;
    for (int is = 1; is <= nsm; is++) {
      deltgam += fabs(gamdot(is)) * dt;
    }

    for (int is = 1; is <= nsm; is++) {
      dtau = 0.0;
      for (int js = 1; js <= nsm; js++) {
        dtau += hard(is,js) * fabs(gamdot(js)) * dt;
      }

      tau0  = voceParam(1);
      tau1  = voceParam(3);
      thet0 = voceParam(4);
      thet1 = voceParam(5); 
      tiny = 1.0e-4*tau0;

      voce = 0.0;
      if (fabs(thet0) > tiny) {
        voce = thet1 * deltgam;
        if (fabs(tau1) > tiny) {
          fact = fabs(thet0/tau1);
          expini = exp(-gamtotx*fact);
          expdel = exp(-deltgam*fact);
          voce += -(fact*tau1-thet1)/fact*expini*(expdel-1.0) - 
                    thet1/fact*expini*(expdel*((gamtotx+deltgam)*fact+1.0)-(gamtotx*fact+1.0));
        } // end if (fabs(tau1) > tiny)
      } // end if (fabs(thet0) > tiny)

      crss(1,is) += dtau*voce/deltgam;
      crss(2,is) += dtau*voce/deltgam;
    } // end for is

    gacumgr = gamtotx + deltgam;

    return;
} // end harden

} // end namespace EVP
