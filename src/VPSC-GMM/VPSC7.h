#ifndef VPSC7_H
#define VPSC7_H

//

#include <assert.h>
#include <math.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "matar.h"

#define MATARIX ViewCMatrixMat
#define MATARRY ViewCArrayMat
#define NGR_		30
#define NPH_		1
#define NSYS_		24
#define NMOD_		2
#define NTWS_		0
#define NTWM_		0
#define NELTEX_		10
#define NGAUMX_		32
#define NGAUMX2_	1792
#define NSTATV_		NELTEX_ + 4 * NGAUMX2_ + 105 * NGR_ + 264 * NPH_ + 7 * NGR_ * NSYS_ + NGR_ * NTWM_ + NGR_ * NTWS_ + NMOD_ * NPH_ + 21 * NPH_ * NSYS_ + 5 * NPH_ * NTWM_ + NPH_ * NSYS_ * NSYS_ + 679


#define INTERACTION_    1
#define IBCINV_			1


//1. PropType
// Properties of the material.  I am unsure why these are separated from the state variables in the code.
class PropType {

private:

	size_t scalar_start_;

public:

	double cv, errs, errd, errm, errso, stranvm_pl, dstranvmth_elpl, dstranvmth_max, dstranvmth_min, b_int;

	size_t nsymop, nph, iphbot, iphtop, itmaxext, itmaxint, itmaxso, irecover, isave, icauchy, interaction, nunique, nrsmin,
		ihardlaw, iratesens, iupdori, iupdshp, iupdhar, itemphard, ibcinv, irsvar, ipr, ipru, kinctex, neltex, strain_control,
		grbot, grtop;

	MATARIX <double>  dnca;        // (nph, nsys, 3)
	MATARIX <double>  dbca;        // (nph, nsys, 3)
	MATARIX <double>  schca;       // (nph, nsys, 5)
	MATARIX <double>  hard;        // (nph, nsys, nsys)
	MATARIX <double>  tau;         // (nph, nsys, 2)
	MATARIX <double>  tau0_mode;   // (nph, nsys, 3)
	MATARIX <double>  thet;        // (nph, nsys, 2)
	MATARIX <double>  twthres;     // (nph, ntwm, 2)
	MATARIX <double>  twsh;        // (nph, ntwm)
	MATARIX <double>  c2ca;        // (nph, 6, 6)
	MATARIX <double>  alpha;       // (ngaumx2, 3)
	MATARIX <double>  ww;          // (ngaumx2)
	MATARIX <double>  h;           // (24, 3, 3)
	MATARIX <double>  cvec;        // (3, 3)


	MATARIX <double>  ngr;         // (nph + 1)
	MATARIX <double>  nrs;         // (nph, nsys)
	MATARIX <double>  isense;      // (nph, nsys)
	MATARIX <double>  nsm;         // (nph, nmod)
	MATARIX <double>  isectw;      // (nph, nsys)
	MATARIX <double>  nmodes;      // (nph)
	MATARIX <double>  nsyst;       // (nph)
	MATARIX <double>  ntwmod;      // (nph)
	MATARIX <double>  ntwsys;      // (nph)
	MATARIX <double>  ishape;      // (nph + 1)
	MATARIX <double>  idsim;       // (6)
	MATARIX <double>  iscau;       // (6)
	MATARIX <double>  iudot;       // (3, 3)
	MATARIX <double>  ngaussph;    // (3)
	MATARIX <double>  ngaussth;    // (3)
	MATARIX <double>  eltex;       // (neltex)


	// default constructor
	KOKKOS_FUNCTION  PropType();

	//View Constructors
	KOKKOS_FUNCTION  PropType& init(double* arr_in, size_t& cnt);

	//View Constructors
	KOKKOS_FUNCTION  PropType& write(double* arr_in);

	// deconstructor
	KOKKOS_FUNCTION  ~PropType();

}; // end of PropType


//End of PropType Definition


//1. StateVType
// State Vars
class StateVType {

private:

	size_t scalar_start_;
public:

	//scalars
	double kinc, evmt, evmp, svm, dvm, dstranaccvm, dtimeacc, dstranvmth, temp_ini, temp_fact, temp_amp, temp, tempold, dtemp, evmpold, wrplastic, wplastic;

	size_t jran;

	// arrays
	MATARIX <double> dbart;         // (5)
	MATARIX <double> rot_acum;      // (3, 3)
	MATARIX <double> xmtg;          // (5, 5)
	MATARIX <double> xltg;          // (5, 5)
	MATARIX <double> dzero;         // (5)
	MATARIX <double> eulerph;       // (nph + 1, 3)
	MATARIX <double> axisph;        // (nph + 1, 4, 3)
	MATARIX <double> fijph;         // (nph + 1, 3, 3)
	MATARIX <double> ag;            // (ngr, 3, 3)
	MATARIX <double> crss;          // (ngr, nsys)
	MATARIX <double> sg;            // (ngr, 5)
	MATARIX <double> sghyd;         // (ngr)
	MATARIX <double> gtotgr;        // (ngr)
	MATARIX <double> xmctg;         // (ngr, 5, 5)
	MATARIX <double> dczero;        // (ngr, 5)
	MATARIX <double> eftwfr;        // (nph, ntwm)
	MATARIX <double> twfrph;        // (nph, ntwm)
	MATARIX <double> twfrsy;        // (ngr, ntws)
	MATARIX <double> twfgr;			// (ntwm, ngr)
	MATARIX <double> wgt;           // (ngr)
	MATARIX <double> wph;           // (nph)
	MATARIX <double> gamd0g;        // (ngr)
	MATARIX <double> gamdot;        // (ngr, nsys)
	MATARIX <double> sch;           // (ngr, 5, nsys)
	MATARIX <double> udot;          // (3, 3)
	MATARIX <double> dav;           // (5)
	MATARIX <double> asph;          // (nph, 3, 3, 3, 3)
	MATARIX <double> dg;            // (ngr, 5)
	MATARIX <double> dstranacc;     // (3, 3)
	MATARIX <double> drotacc;       // (3, 3)
	MATARIX <double> xmtgdot;       // (5, 5)
	MATARIX <double> dzerodot;      // (5)
	MATARIX <double> cscold;        // (6, 6)
	MATARIX <double> stressold;     // (3, 3)


	MATARIX <double> ktwsmx;        // (ngr)
	MATARIX <double> ntwevents;     // (ngr)
	MATARIX <double> giph;          // (ngr)
	MATARIX <double> iflat;         // (nph + 1)
	MATARIX <double> imark;         // (ngr)
	MATARIX <double> iptsgr;        // (ngr)
	MATARIX <double> iptmgr;        // (ngr)

	// default constructor
	KOKKOS_FUNCTION  StateVType();

	//View Constructors
	KOKKOS_FUNCTION  StateVType& init(double* arr_in, size_t& cnt);

	//View Constructors
	KOKKOS_FUNCTION  StateVType& write(double* arr_in);


	// deconstructor
	KOKKOS_FUNCTION  ~StateVType();

}; // end of StateVType


//End of PropType Definition


//3. IntVType
// Variables associated with each step, allocated and de-allocated each increment
class IntVType {
private:
	size_t scalar_start_;

public:

	//scalars
	double svm, dvm, epsacu, tincr;

	size_t jjxrs, kinc, k_cell;

	// arrays
	MATARIX <double> sastav;         // (5)
	MATARIX <double> sastbar;        // (5)
	MATARIX <double> dast;           // (5)
	MATARIX <double> sav;            // (5)
	MATARIX <double> dav;            // (5)
	MATARIX <double> sbar;           // (5)
	MATARIX <double> dbar;           // (5)
	MATARIX <double> xlijph;         // (nph, 3, 3)
	MATARIX <double> asph;           // (nph, 3, 3, 3, 3)
	MATARIX <double> xmastph;        // (nph, 5, 5)
	MATARIX <double> dg;             // (ngr, 5)
	MATARIX <double> stry;           // (ngr, 5)
	MATARIX <double> scauchy;        // (3, 3)
	MATARIX <double> udot;           // (3, 3)
	MATARIX <double> dsim;           // (3, 3)
	MATARIX <double> rotbar;         // (3, 3)
	MATARIX <double> sdeviat;        // (3, 3)
	MATARIX <double> csc;            // (6, 6)
	MATARIX <double> ssc;            // (6, 6)
	MATARIX <double> dsimtot;        // (3, 3)
	MATARIX <double> sbarold;        // (5)
	MATARIX <double> cgr;            // (ngr, 6, 6)

	//View Constructors
	KOKKOS_FUNCTION  IntVType();

	//View Constructors
	KOKKOS_FUNCTION  IntVType& init(double* this_array_, size_t& cnt);

	// deconstructor
	KOKKOS_FUNCTION  ~IntVType();
};


//End of IntVType Definition

KOKKOS_FUNCTION double sumProd(const size_t sz, 
    double* a1, 
    double* a2, 
    double* a3);


KOKKOS_FUNCTION double sumProd(const size_t sz,
    double* a1,
    double* a2);


KOKKOS_FUNCTION double sumProd(const size_t nk, const size_t nc, double* a, double* b) ;


KOKKOS_FUNCTION double kahanSumProd(const size_t sz,
    double* a1,
    double* a2); 

KOKKOS_FUNCTION double tnorm(MATARIX <double>& mat); 

KOKKOS_FUNCTION double max(MATARIX <double>& mat); 


KOKKOS_FUNCTION double xid(size_t i, size_t j); 



KOKKOS_FUNCTION void operator+= (MATARIX <double>& a, double* temp) ;

KOKKOS_FUNCTION void operator-= (MATARIX <double>& a, double* temp); 

KOKKOS_FUNCTION void operator*= (MATARIX <double>& a, double* temp); 

KOKKOS_FUNCTION void operator/= (MATARIX <double>& a, double* temp); 


KOKKOS_FUNCTION void operator+= (MATARIX <double>& a, const double temp); 

KOKKOS_FUNCTION void operator-= (MATARIX <double>& a, const double temp); 

KOKKOS_FUNCTION void operator*= (MATARIX <double>& a, const double temp); 

KOKKOS_FUNCTION void operator/= (MATARIX <double>& a, const double temp); 


KOKKOS_FUNCTION void operator+= (MATARRY <double>& a, const double temp); 

KOKKOS_FUNCTION void slice(const size_t N, double* a, double* some_matrix); 

KOKKOS_FUNCTION void peye(MATARIX <double>& a, size_t N, double I); 

KOKKOS_FUNCTION void eye(double* y, size_t n, double I) ;

KOKKOS_FUNCTION void transpose(MATARIX <double>& a, 
    const size_t N); 


KOKKOS_FUNCTION double det(MATARIX <double>& a);



KOKKOS_FUNCTION double tmismatch(const size_t sz, 
    MATARIX <double> v1,
    MATARIX <double> v2);


KOKKOS_FUNCTION void fullout(MATARIX <double>& x) ;


KOKKOS_FUNCTION void fullout(const size_t sz, double* x); 

KOKKOS_FUNCTION size_t random2(size_t seed) ;


KOKKOS_FUNCTION void drot2spin(MATARIX <double>& drot, double dtime, MATARIX <double>& w_app);


KOKKOS_FUNCTION void rodrigues(MATARIX <double>& c, MATARIX <double>& arot) ;

KOKKOS_FUNCTION void voigt(MATARIX <double>& t1, 
    MATARIX <double>& t2, 
    size_t iopt);

KOKKOS_FUNCTION void euler(size_t iopt, double& ph, double& th, double& tm, MATARIX <double>& a);



KOKKOS_FUNCTION void eshelby(MATARIX <double>& axis,
    MATARIX <double>& c4,
    double keff,
    MATARIX <double>& esim,
    MATARIX <double>& escr,
    MATARIX <double>& desh,
    MATARIX <double>& pesh,
    double& pdil,
    MATARIX <double>& dldm,
    MATARIX <double>& dsddm,
    PropType& props,
    size_t ioption);


KOKKOS_FUNCTION void elsc(size_t ioption,
    StateVType& statv,
    PropType& props,
    IntVType& intv);


KOKKOS_FUNCTION void newton_raphson(size_t irc,
    size_t kmax,
    double eps,
    double taulim,
    const size_t nsystx, 
    MATARIX <double>& x,
    MATARIX <double>& db,
    MATARIX <double>& xmastx,
    MATARIX <double>& scx,
    MATARIX <double>& itaux,
    MATARIX <double>& gamd0x,
    MATARIX <double>& nrsx,
    MATARIX <double>& isensex,
    size_t& ierror);


KOKKOS_FUNCTION void  scale_3(size_t istep, StateVType& statv, PropType& props, IntVType& intv);



KOKKOS_FUNCTION void twin_orientation(MATARIX <double>& bur, MATARIX <double>& atwin);



KOKKOS_FUNCTION void update_crss_voce(size_t ioption, StateVType& statv, PropType& props, IntVType& intv);



KOKKOS_FUNCTION void update_crss_temp(StateVType& statv, PropType& props, IntVType& intv);



KOKKOS_FUNCTION void update_fij(size_t iph, StateVType& statv, IntVType& intv);



KOKKOS_FUNCTION void update_shape(size_t iph, StateVType& statv, size_t interaction);




KOKKOS_FUNCTION void update_orientation(StateVType& statv, PropType& props, IntVType& intv);



KOKKOS_FUNCTION void update_schmid(size_t ngr, MATARIX <double> nsyst, MATARIX <double> schca, MATARIX <double> ag, MATARIX <double> sch, MATARIX <double> giph);


KOKKOS_FUNCTION void update_twinning(size_t iph,
    StateVType& statv,
    PropType& props,
    IntVType& intv);


KOKKOS_FUNCTION void grain_stress_moduli(size_t interx, size_t jsc, size_t kcalmod, bool calcStress, StateVType& statv, PropType& props, IntVType& intv);






KOKKOS_FUNCTION void  elas_bc(double dt,
    MATARIX <double> xltg,
    MATARIX <double> dzero,
    MATARIX <double> dsimtot,
    MATARIX <double> csc,
    MATARIX <double> dbar,
    MATARIX <double> sbar,
    MATARIX <double> sbarold);



KOKKOS_FUNCTION void vpsc_affine(size_t istep, StateVType& statv, PropType& props, IntVType& intv) ;

//
//



KOKKOS_FUNCTION void vpsc(size_t istep,
    StateVType& statv,
    PropType& props,
    IntVType& intv);





KOKKOS_FUNCTION void update_temperature(double svm, double evmp, double& evmpold, double& temp, double& dtemp, double wrplastic, double& wplastic, double b_int);


KOKKOS_FUNCTION void initial_state_guess(StateVType& statv, PropType& props, IntVType& intv);

KOKKOS_FUNCTION void evol_vpsc(StateVType& statv, PropType& props, IntVType& intv, const ViewCArrayKokkos <double>& vel_grad, double* stress_in, double dtime, size_t k_cell, ViewCMatrixKokkos <double>& ddsdde); 



template <typename T> KOKKOS_FUNCTION bool ismember(T a, MATARIX <T>& b);


KOKKOS_FUNCTION void nextline(std::fstream& in, size_t nline) ;

KOKKOS_FUNCTION void gauss_legendre(double x1, double x2, MATARIX <double>& x, MATARIX <double>& w, size_t n);



KOKKOS_FUNCTION void eshelby_init(PropType& props) ;


KOKKOS_FUNCTION void init_crss_voce(size_t ioption, StateVType& statv, PropType& props);



KOKKOS_FUNCTION void init_crss_temp(StateVType& statv, PropType& props);


KOKKOS_FUNCTION void write_tex(size_t iph, std::string filename, PropType& props, StateVType& statv) ;


KOKKOS_FUNCTION void data_texture(size_t iph, std::string filename, PropType& props, StateVType& statv);


KOKKOS_FUNCTION void crystal_symmetry(size_t ioption, std::fstream& file, PropType& props, size_t& icrysym, MATARIX <double>& sn, MATARIX <double>& sneq, MATARIX <double>& sb, size_t& npol) ;

KOKKOS_FUNCTION void data_crystal(size_t iph, std::string filename, PropType& props, StateVType& statv);



KOKKOS_FUNCTION void user_mat_init(double* statev1D);



KOKKOS_FUNCTION size_t user_mat_number_vars();

#endif 