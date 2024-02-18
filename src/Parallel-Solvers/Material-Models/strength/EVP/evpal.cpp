#include "user_mat.h"

// evaluate elasto-viscoplastic crystal plasticity model using Newton-Raphson solve
KOKKOS_FUNCTION
void evpal(real_t *stress_, real_t *edotp_, real_t *gamdot_, real_t *stress_n_, real_t *strain_,
		real_t* ept_n_, real_t *cg66_, real_t *sc_, real_t *crss_, const real_t gamd0, const real_t nrs, 
		const int nsmx, const real_t dt, const size_t cycle, real_t *B_basis_, const size_t elnum){ //const real_t gacumgr, real_t *voceParam, real_t *hard, 

		//ChgBasis cb;

    int itmaxal;
    int iter1;
    int isign;

    real_t sgnorm;
    real_t dgnorm;
    real_t erroral;
    real_t erral;
    real_t dsgnorm1;
    real_t ddgnorm;
    real_t errald;
    real_t resn;
    real_t resnold;

    const int nsm=12;

    // thread private arrays
    real_t sgaux_[3*3];
    real_t sg6_[6];
    real_t sg6old_[6];
    real_t eptaux_[3*3];
    real_t ept6_[6];
    real_t edotpaux_[3*3];
    real_t edotp6_[6];
    real_t dedotp66_[6*6];
    real_t strainaux_[3*3];
    real_t strain6_[6];
    real_t sg66_[6*6];
    real_t strainceq_[3*3];
    real_t strainceq6_[6];
    real_t res_[6];
    real_t xjacobinv_[6*6];
    real_t rss_[nsm];
    real_t rss1_[nsm];
    real_t rss2_[nsm];
    real_t taux_[nsm*2];

    ViewCMatrixKokkos <double> stress(stress_,3,3);
    ViewCMatrixKokkos <double> edotp(edotp_,3,3);
    ViewCMatrixKokkos <double> gamdot(gamdot_,nsm);
    ViewCMatrixKokkos <double> stress_n(stress_n_,3,3);
    ViewCMatrixKokkos <double> strain(strain_,3,3);
    ViewCMatrixKokkos <double> ept_n(ept_n_,3,3);
    ViewCMatrixKokkos <double> cg66(cg66_,6,6);
    //ViewCMatrixKokkos <double> sg66(sg66_,6,6);
    ViewFMatrixKokkos <double> sc(sc_,5,nsm);
    ViewCMatrixKokkos <double> crss(crss_,2,nsm);
    ViewCMatrixKokkos <double> B_basis(B_basis_,3,3,6);
    // create views of thread private arrays
    ViewCMatrixKokkos <double> sgaux(sgaux_,3,3);
    ViewCMatrixKokkos <double> sg6(sg6_,6);
    ViewCMatrixKokkos <double> sg6old(sg6old_,6);
    ViewCMatrixKokkos <double> eptaux(eptaux_,3,3);
    ViewCMatrixKokkos <double> ept6(ept6_,6);
    ViewCMatrixKokkos <double> edotpaux(edotpaux_,3,3);
    ViewCMatrixKokkos <double> edotp6(edotp6_,6);
    ViewCMatrixKokkos <double> dedotp66(dedotp66_,6,6);
    ViewCMatrixKokkos <double> strainaux(strainaux_,3,3);
    ViewCMatrixKokkos <double> strain6(strain6_,6);
    ViewCMatrixKokkos <double> sg66(sg66_,6,6);
    ViewCMatrixKokkos <double> strainceq(strainceq_,3,3); //is this ever used?
    ViewCMatrixKokkos <double> strainceq6(strainceq6_,6);
    ViewCMatrixKokkos <double> res(res_,6);
    ViewCMatrixKokkos <double> xjacobinv(xjacobinv_,6,6);
    ViewCMatrixKokkos <double> rss(rss_,nsm);
    ViewCMatrixKokkos <double> rss1(rss1_,nsm);
    ViewCMatrixKokkos <double> rss2(rss2_,nsm);
    ViewCMatrixKokkos <double> taux(taux_,2,nsm);

    // TODO: optimize indexing of this loop
    for (int ii = 1; ii <= 6; ii++) {
      for (int jj = 1; jj <= 6; jj++) {
        sg66(ii,jj) = cg66(ii,jj);
      }
    }
    //lu_inverse(sg66.pointer(), 6);
    inverse_gj(sg66.pointer(), 6);

    // TODO: optimize indexing of this loop
    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        sgaux(ii,jj)     = stress_n(ii,jj);
        eptaux(ii,jj)    = ept_n(ii,jj);
        strainaux(ii,jj) = strain(ii,jj);
      }
    }

    chg_basis_2(sg6.pointer(), sgaux.pointer(), 2, 6, B_basis.pointer());
    chg_basis_2(ept6.pointer(), eptaux.pointer(), 2, 6, B_basis.pointer());
    chg_basis_2(strain6.pointer(), strainaux.pointer(), 2, 6, B_basis.pointer());
    
    sgnorm = 0.0;
    dgnorm = 0.0;
    // TODO: optimize indexing of this loop
    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        sgnorm += stress_n(ii,jj)*stress_n(ii,jj);
        dgnorm += strainaux(ii,jj)*strainaux(ii,jj);
      }
    }
    sgnorm = sqrt(sgnorm);
    dgnorm = sqrt(dgnorm);

    //making prediction for starting from 0 stress
    if (sgnorm == 0.0) {
      for (int ii = 1; ii <= 6; ii++) {
          sg6(ii) = 0.0;
        for (int jj = 1; jj <= 6; jj++) {
          sg6(ii) += cg66(ii,jj)*strain6(jj);
        } // end for jj
      } // end for ii

    }

    erroral = 0.0000001;
    erral   = 2.0*erroral;
    errald  = erral;
    itmaxal = 500;
    iter1 = 0;
    resnold = 1.0;

    for (int is = 1; is <= nsm; is++) {
      taux(1,is) = crss(1,is);
      taux(2,is) = crss(2,is);
    } // end for is

	  while (iter1 < itmaxal && erral > erroral) {
	    iter1 += 1;
      int icut = 0;
      resn = resnold*2.0;

	    //if (cycle != 0) {
      while (resn > resnold and icut < 10) {
	      //     GET RESOLVED SHEAR STRESSES 'rss' AND SHEAR RATES 'gamdot'.
	      //     SIGN(GAMDOT)=SIGN(RSS).
	      //     NRS CAN BE EVEN OR ODD.
	      //     RSS1 IS ALWAYS > 0 AND IS USED TO CALCULATE VISCOUS COMPLIANCE.

	      for (int is = 1; is <= nsm; is++) {
	        rss(is) = sc(1,is)*sg6(1) +
	                  sc(2,is)*sg6(2) + 
	                  sc(3,is)*sg6(3) + 
	                  sc(4,is)*sg6(4) + 
	                  sc(5,is)*sg6(5);
	        isign = 1;
	        if ( rss(is) < 0.0 ) {
	          isign = 2;
	        }

	        rss(is) = rss(is)/taux(isign,is);

	        rss1(is) = gamd0 * nrs * abs(pow(rss(is),(nrs-1))) / taux(isign,is); //in plastic strain rate derivative
	        rss2(is) = gamd0 * abs(pow(rss(is),nrs)) * copysign(1.0,rss(is)); //in plastic strain rate

	        gamdot(is) = rss2(is);
	      } // end for is

	      for (int ii = 1; ii <= 5; ii++) {
	        edotp6(ii) = 0.0;
	        for (int k1 = 1; k1 <= nsm; k1++) {
	          edotp6(ii) += sc(ii,k1) * rss2(k1);
	        }
	      }
	      edotp6(6) = 0.0;

	      for (int jj = 1; jj <= 5; jj++) {
	        for (int ii = 1; ii <= 5; ii++) {

	          dedotp66(ii,jj) = 0.0;

	          for (int k1 = 1; k1 <= nsm; k1++) {       
	            dedotp66(ii,jj) += sc(ii,k1)*sc(jj,k1)*rss1(k1);
	          } // end for k1

	        } // end for ii
	      } // end for jj

	      for (int ii = 1; ii <= 6; ii++) {
	        dedotp66(ii,6) = 0.0;
	        dedotp66(6,ii) = 0.0;
	      }

	    // } else { // else for if (fd.ithermo != 1 || imicro != 1)
	    //   for (int ii = 1; ii <= 6; ii++) {
	    //     edotp6(ii) = 0.0;
	    //   } // end for ii

	    // } // end for else (fd.ithermo != 1 || imicro != 1)


  	    for (int ii = 1; ii <= 6; ii++) {
  	      //if (cycle == 0) {//(fd.ithermo == 1 && imicro == 1) {
  	      //   strainceq6(ii) = ept6(ii);
  	      // } else {
  	        strainceq6(ii) = ept6(ii)+edotp6(ii)*dt;
  	      // }
  	      for (int jj = 1; jj <= 6; jj++) {
  	        strainceq6(ii) += sg66(ii,jj)*sg6(jj);
  	      } // end for jj
  	    } // end for ii

  	    chg_basis_1(strainceq6.pointer(), strainceq.pointer(), 1, 6, B_basis.pointer());
        resn=0.0;
  	    for (int ii = 1; ii <= 6; ii++) {
  	      res(ii) = strainceq6(ii) - strain6(ii);
          resn += res(ii)*res(ii);
  	    }
        resn = sqrt(resn);

        if (iter1==1){
          resnold = resn;
        }

        chg_basis_1(sg6.pointer(), sgaux.pointer(), 1, 6, B_basis.pointer());

        if (resn > resnold){
          for (int ii = 1; ii <= 6; ii++) {
            sg6(ii) = sg6old(ii) + (sg6(ii)-sg6old(ii))*0.5;
          }
        }

        icut += 1;
        chg_basis_1(edotp6.pointer(), edotpaux.pointer(), 1, 6, B_basis.pointer());
      }

      resnold = resn;
      for (int ii = 1; ii <= 6; ii++) {
        sg6old(ii) = sg6(ii);
      }

	    for (int ii = 1; ii <= 6; ii++) {
	      for (int jj = 1; jj <= 6; jj++) {

	        xjacobinv(ii,jj) = 0.0;
	        //if (cycle == 0) {//(fd.ithermo == 1 && imicro == 1) {
	        //   xjacobinv(ii,jj) = sg66(ii,jj);
	        // } else {
	          xjacobinv(ii,jj) = sg66(ii,jj)+dedotp66(ii,jj)*dt;
	        // }

	      } // end for jj
	    } // end for ii

	    //lu_inverse(xjacobinv.pointer(), 6);
	    inverse_gj(xjacobinv.pointer(), 6);

	    // TODO: optimize indexing of this loop
	    for (int ii = 1; ii <= 6; ii++) {
	      for (int jj = 1; jj <= 6; jj++) {
	        sg6(ii) += -xjacobinv(ii,jj)*res(jj);
	      } // end for jj
	    } // end for ii

	    dsgnorm1 = 0.0;
	    ddgnorm  = 0.0;

	    for (int ii = 1; ii <= 6; ii++) {
	//#ifndef NR_correction_cut      
	      dsgnorm1 += (sg6(ii)-sg6old(ii))*(sg6(ii)-sg6old(ii));
	//#else
	//            dsgnorm1 += ((sg6(ii)-sg6old(ii))/NR_CUT_fac)*((sg6(ii)-sg6old(ii))/NR_CUT_fac);
	//#endif
	      ddgnorm  += (strainceq6(ii)-strain6(ii))*(strainceq6(ii)-strain6(ii));
	    } // end for ii

      dsgnorm1 = sqrt(dsgnorm1);
      ddgnorm = sqrt(ddgnorm);

	    erral  = dsgnorm1 / sgnorm;
	    errald = ddgnorm / dgnorm;
      //if (sgnorm==0.0) {
        erral = errald;
      //}
	  } // end if (iter1 < itmaxal && erral > erroral)

    chg_basis_1(sg6.pointer(), sgaux.pointer(), 1, 6, B_basis.pointer());
    chg_basis_1(edotp6.pointer(), edotpaux.pointer(), 1, 6, B_basis.pointer());

    // TODO: optimize indexing of this loop
    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        stress(ii,jj) = sgaux(ii,jj);
        edotp(ii,jj)  = edotpaux(ii,jj);
      } // end for jj
    } // end for ii
} // end evpal

//
//  inverse_gj
//
KOKKOS_FUNCTION
void inverse_gj(real_t *a_, int n)
{
  ViewCMatrixKokkos <double> a(a_,n,n);
  real_t             tmp;
  real_t             fac;

  for (int i = 1; i <= n; i++) { // loop 1
    fac = 1.0 / a(i,i); // This would become inverse if done on blocks
    a(i,i) = 1.0;
    for (int j = 1; j <= n; j++) {
      a(j,i) = a(j,i)*fac; 
    }

    // NOTE:These loops don't shrink as i increases so load should
    // be balanced?
    for (int k = 1; k <= i-1; k++) {
      tmp = a(i,k);
      a(i,k) = 0.0;
      for (int j = 1; j <= n; j++) {
        a(j,k) += -a(j,i)*tmp;
      }
    }
    for (int k = i+1; k <= n; k++) {
      tmp = a(i,k);
      a(i,k) = 0.0;
      for (int j = 1; j <= n; j++) {
        a(j,k) += -a(j,i)*tmp;
      }
    } 
  } // end loop 1
}
//
//  lu_inverse
//
KOKKOS_FUNCTION
void lu_inverse(real_t *a_, int n)
{
  real_t              *y_ = new real_t[n*n];
  int                 *indx_ = new int[n];
  ViewCMatrixKokkos <double> y(y_,n,n);
  ViewCMatrixKokkos <int> indx(indx_,n);
  ViewCMatrixKokkos <double> a(a_,n,n);
  real_t              d = 0.0;
  int                 isingular = 0;

  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      y(i,j) = 0.0;
    }
    y(i,i) = 1.0;
  }

  ludcmp(a_,n,n,indx_,d,isingular);
  if (isingular) {
    printf("singular matrix in ludcmp\n");
    exit(1);
    return;
  }

  for (int j = 1; j <= n; j++) {
    lubksb(a_,n,n,indx_,&y(1,j));
  }

  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      a(i,j) = y(i,j);
    }
  }
  delete [] y_;
  delete [] indx_;
}
//
//  ludcmp
//
KOKKOS_FUNCTION
void ludcmp(real_t *a_, int n, int np, int *indx_, real_t d, int &isingular)
{
  ViewCMatrixKokkos <double> a(a_,n,n);
  ViewCMatrixKokkos <int> indx(indx_,n);
  int                 imax;
  const int           nmax = 500;
  real_t              vv_[nmax];
  ViewCMatrixKokkos <double> vv(vv_,nmax);
  real_t              aamax;
  real_t              dum;
  real_t              sum;

  d = 1.0;
  for (int i = 1; i <= n; i++) { // loop 12
    aamax = 0.0;
    for (int j = 1; j <= n; j++) { // loop 11
      if (abs(a(i,j)) > aamax) aamax = abs(a(i,j));
    } // end loop 11
    //if (aamax == 0.0) {
    //  printf("singular matrix in ludcmp\n");
    //  PAUSE;
    //} // end if (aamax == 0.0)
    if (aamax == 0.0) {
      //printf("singular matrix in ludcmp\n");
      isingular = 1;
      return;
    }
    vv(i) = 1.0 / aamax;
  } // end loop 12

  for (int j = 1; j <= n; j++) { // loop 19

    for (int i = 1; i <= j-1; i++) { // loop 14
      sum = a(i,j);
      for (int k = 1; k <= i-1; k++) { // loop 13
        sum += - a(i,k)*a(k,j);
      } // end loop 13
      a(i,j) = sum;
    } // end loop 14
    aamax = 0.0;

    for (int i = j; i <= n; i++) { // loop 16
      sum = a(i,j);
      for (int k = 1; k <= j-1; k++) { // loop 15
        sum += - a(i,k)*a(k,j);
      } // end loop 15
      a(i,j) = sum;
      dum = vv(i)*abs(sum);
      if (dum >= aamax) {
        imax = i;
        aamax = dum;
      } // end if (dum >= aamax)
    } // end loop 16

    if (j != imax) {
      for (int k = 1; k <= n; k++) { // loop 17
        dum = a(imax,k);
        a(imax,k) = a(j,k);
        a(j,k) = dum;
      } // end loop 17
      d = -d;
      vv(imax) = vv(j);
    } // end if (j != imax)
    indx(j) = imax;

    //if (a(j,j) == 0.0) a(j,j) = TINY;
    if (aamax == 0.0) {
    //if (a(j,j) == 0.0) {
      //printf("singular matrix in ludcmp\n");
      isingular = 1;
      return;
    }

    if (j != n) {
      dum = 1.0 / a(j,j);
      for (int i = j+1; i <= n; i++) { // loop 18
        a(i,j) = a(i,j)*dum;
      } // end loop 18
    } // end if (j != n)

  } // end loop 19

  isingular = 0;
}
//
//  lubksb
//
KOKKOS_FUNCTION
void lubksb(real_t *a_, int n, int np, int *indx_, real_t *b_)
{
  ViewCMatrixKokkos <double> a(a_,n,n);
  ViewCMatrixKokkos <int> indx(indx_,n);
  ViewCMatrixKokkos <double> b(b_,n);
  real_t              sum;
  int                 ii;
  int                 ll;

  ii = 0;
  for (int i = 1; i <= n; i++) { // loop 12
    ll = indx(i);
    sum = b(ll);
    b(ll) = b(i);
    if (ii != 0) {
      for (int j  = ii; j <= i-1; j++) { // loop 11
        sum += - a(i,j)*b(j);
      } // end loop 11
    } else if (sum != 0.0) {
      ii = i;
    } // end else if (sum != 0.0)
    b(i) = sum;
  } // end loop 12
  for (int i = n; i >= 1; i--) { // loop 14
    sum = b(i);
    for (int j = i+1; j <= n; j++) { // loop 13
      sum += - a(i,j)*b(j);
    } // end loop 13
    b(i) = sum / a(i,i);
  } // end loop 14

}
