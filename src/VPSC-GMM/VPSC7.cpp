
// vpsccpp.cpp : this file contains the 'main' function. program execution begins and ends there.
//
#include <assert.h>
#include <math.h>
#include <fstream>

#include "VPSC7.h"
//#include "vpsc_math.cpp"

//fixed sized matrix of dimension 1
#define fmat1(type, name, N1) type name##_1D[N1]; MATARIX <type> name(&name##_1D[0],N1)
//for (size_t ii = 0; ii < (N1); ii++){name##_1D[ii] = 0.0;};
//fixed sized matrix of dimension 2
#define fmat2(type, name, N1, N2) type name##_1D[(N1)*(N2)]; MATARIX <type> name(&name##_1D[0],N1,N2)
//for (size_t ii = 0; ii < (N1)*(N2); ii++){name##_1D[ii] = 0.0;};
//fixed sized matrix of dimension 3
#define fmat3(type, name, N1, N2, N3) type name##_1D[(N1)*(N2)*(N3)]; MATARIX <type> name(&name##_1D[0],N1,N2,N3)
//for (size_t ii = 0; ii < (N1)*(N2)*(N3); ii++){name##_1D[ii] = 0.0;};
//fixed sized matrix of dimension 4
#define fmat4(type, name, N1, N2, N3, N4) type name##_1D[(N1)*(N2)*(N3)*(N4)]; MATARIX <type> name(&name##_1D[0],N1,N2,N3,N4)
//for (size_t ii = 0; ii < (N1)*(N2)*(N3)*(N4); ii++){name##_1D[ii] = 0.0;};

//vector with a max size
#define xmat1(type, name, maxsz, N1) type name##_1D[maxsz]; MATARIX <type> name(&name##_1D[0],N1)
//for (size_t ii = 0; ii < (maxsz); ii++){name##_1D[ii] = 0.0;};
//matrix with a max size of dimension 2
#define xmat2(type, name, maxsz, N1, N2) type name##_1D[maxsz]; MATARIX <type> name(&name##_1D[0],N1,N2)
//for (size_t ii = 0; ii < (maxsz); ii++){name##_1D[ii] = 0.0;};

//allocated matrix of dimension 1
#define vmat1(type, name, N1) type* name##_1D; name##_1D = new type[(N1)]; MATARIX <type> name(&name##_1D[0],N1)
//allocated matrix of dimension 2
#define vmat2(type, name, N1, N2) type* name##_1D; name##_1D = new type[(N1)*(N2)]; MATARIX <type> name(&name##_1D[0],N1,N2)


#define PI   (3.1415926535897932384626433832795)

//KOKKOS_FUNCTION bool chk_sym(MATARIX <double>& a);

//KOKKOS_FUNCTION double tnorm(MATARIX <double>& mat);
//KOKKOS_FUNCTION double   max(MATARIX <double>& mat);
//
//KOKKOS_FUNCTION double xid(size_t i, size_t j);
//KOKKOS_FUNCTION double det(MATARIX <double>& a);
//KOKKOS_FUNCTION double tmismatch(const size_t sz, MATARIX <double> v1, MATARIX <double> v2);
//KOKKOS_FUNCTION double sumProd(const size_t sz, double* a1, double* a2, double* a3);
//KOKKOS_FUNCTION double sumProd(const size_t sz, double* a1, double* a2);
////sumprod for matrix multiplication
//KOKKOS_FUNCTION double sumProd(const size_t nk, const size_t nc, double* a, double* b);
//
//KOKKOS_FUNCTION double kahanSumProd(const size_t sz, double* a1, double* a2);


KOKKOS_FUNCTION  PropType::PropType() {
}

KOKKOS_FUNCTION  PropType& PropType::init(double* this_array_, size_t& cnt)
{
	scalar_start_ = cnt;

	cv = this_array_[cnt];
	cnt++;
	errs = this_array_[cnt];
	cnt++;
	errd = this_array_[cnt];
	cnt++;
	errm = this_array_[cnt];
	cnt++;
	errso = this_array_[cnt];
	cnt++;
	stranvm_pl = this_array_[cnt];
	cnt++;
	dstranvmth_elpl = this_array_[cnt];
	cnt++;
	dstranvmth_max = this_array_[cnt];
	cnt++;
	dstranvmth_min = this_array_[cnt];
	cnt++;
	b_int = this_array_[cnt];
	cnt++;
	nsymop = this_array_[cnt];
	cnt++;
	nph = this_array_[cnt];
	cnt++;
	iphbot = this_array_[cnt];
	cnt++;
	iphtop = this_array_[cnt];
	cnt++;
	itmaxext = this_array_[cnt];
	cnt++;
	itmaxint = this_array_[cnt];
	cnt++;
	itmaxso = this_array_[cnt];
	cnt++;
	irecover = this_array_[cnt];
	cnt++;
	isave = this_array_[cnt];
	cnt++;
	icauchy = this_array_[cnt];
	cnt++;
	interaction = this_array_[cnt];
	cnt++;
	nunique = this_array_[cnt];
	cnt++;
	nrsmin = this_array_[cnt];
	cnt++;
	ihardlaw = this_array_[cnt];
	cnt++;
	iratesens = this_array_[cnt];
	cnt++;
	iupdori = this_array_[cnt];
	cnt++;
	iupdshp = this_array_[cnt];
	cnt++;
	iupdhar = this_array_[cnt];
	cnt++;
	itemphard = this_array_[cnt];
	cnt++;
	ibcinv = this_array_[cnt];
	cnt++;
	irsvar = this_array_[cnt];
	cnt++;
	ipr = this_array_[cnt];
	cnt++;
	ipru = this_array_[cnt];
	cnt++;
	kinctex = this_array_[cnt];
	cnt++;
	neltex = this_array_[cnt];
	cnt++;
	strain_control = this_array_[cnt];
	cnt++;

	grbot = this_array_[cnt];
	cnt++;
	grtop = this_array_[cnt];
	cnt++;
	dnca.set(MATARIX <double>(&this_array_[cnt], NPH_, NSYS_, 3));
	cnt += dnca.size();

	dbca.set(MATARIX <double>(&this_array_[cnt], NPH_, NSYS_, 3));
	cnt += dbca.size();

	schca.set(MATARIX <double>(&this_array_[cnt], NPH_, NSYS_, 5));
	cnt += schca.size();

	hard.set(MATARIX <double>(&this_array_[cnt], NPH_, NSYS_, NSYS_));
	cnt += hard.size();

	tau.set(MATARIX <double>(&this_array_[cnt], NPH_, NSYS_, 2));
	cnt += tau.size();

	tau0_mode.set(MATARIX <double>(&this_array_[cnt], NPH_, NSYS_, 3));
	cnt += tau0_mode.size();

	thet.set(MATARIX <double>(&this_array_[cnt], NPH_, NSYS_, 2));
	cnt += thet.size();

	twthres.set(MATARIX <double>(&this_array_[cnt], NPH_, NTWM_, 2));
	cnt += twthres.size();

	twsh.set(MATARIX <double>(&this_array_[cnt], NPH_, NTWM_));
	cnt += twsh.size();

	c2ca.set(MATARIX <double>(&this_array_[cnt], NPH_, 6, 6));
	cnt += c2ca.size();

	alpha.set(MATARIX <double>(&this_array_[cnt], NGAUMX2_, 3));
	cnt += alpha.size();

	ww.set(MATARIX <double>(&this_array_[cnt], NGAUMX2_));
	cnt += ww.size();

	h.set(MATARIX <double>(&this_array_[cnt], 24, 3, 3));
	cnt += h.size();

	cvec.set(MATARIX <double>(&this_array_[cnt], 3, 3));
	cnt += cvec.size();

	ngr.set(MATARIX <double>(&this_array_[cnt], NPH_ + 1));
	cnt += ngr.size();

	nrs.set(MATARIX <double>(&this_array_[cnt], NPH_, NSYS_));
	cnt += nrs.size();

	isense.set(MATARIX <double>(&this_array_[cnt], NPH_, NSYS_));
	cnt += isense.size();

	nsm.set(MATARIX <double>(&this_array_[cnt], NPH_, NMOD_));
	cnt += nsm.size();

	isectw.set(MATARIX <double>(&this_array_[cnt], NPH_, NSYS_));
	cnt += isectw.size();

	nmodes.set(MATARIX <double>(&this_array_[cnt], NPH_));
	cnt += nmodes.size();

	nsyst.set(MATARIX <double>(&this_array_[cnt], NPH_));
	cnt += nsyst.size();

	ntwmod.set(MATARIX <double>(&this_array_[cnt], NPH_));
	cnt += ntwmod.size();

	ntwsys.set(MATARIX <double>(&this_array_[cnt], NPH_));
	cnt += ntwsys.size();

	ishape.set(MATARIX <double>(&this_array_[cnt], NPH_ + 1));
	cnt += ishape.size();

	idsim.set(MATARIX <double>(&this_array_[cnt], 6));
	cnt += idsim.size();

	iscau.set(MATARIX <double>(&this_array_[cnt], 6));
	cnt += iscau.size();

	iudot.set(MATARIX <double>(&this_array_[cnt], 3, 3));
	cnt += iudot.size();

	ngaussph.set(MATARIX <double>(&this_array_[cnt], 3));
	cnt += ngaussph.size();

	ngaussth.set(MATARIX <double>(&this_array_[cnt], 3));
	cnt += ngaussth.size();

	eltex.set(MATARIX <double>(&this_array_[cnt], NELTEX_));
	cnt += eltex.size();


	return *this;
}

KOKKOS_FUNCTION PropType& PropType::write(double* this_array_) {

	size_t cnt = scalar_start_;

	this_array_[cnt] = cv;
	cnt++;
	this_array_[cnt] = errs;
	cnt++;
	this_array_[cnt] = errd;
	cnt++;
	this_array_[cnt] = errm;
	cnt++;
	this_array_[cnt] = errso;
	cnt++;
	this_array_[cnt] = stranvm_pl;
	cnt++;
	this_array_[cnt] = dstranvmth_elpl;
	cnt++;
	this_array_[cnt] = dstranvmth_max;
	cnt++;
	this_array_[cnt] = dstranvmth_min;
	cnt++;
	this_array_[cnt] = b_int;
	cnt++;
	this_array_[cnt] = nsymop;
	cnt++;
	this_array_[cnt] = nph;
	cnt++;
	this_array_[cnt] = iphbot;
	cnt++;
	this_array_[cnt] = iphtop;
	cnt++;
	this_array_[cnt] = itmaxext;
	cnt++;
	this_array_[cnt] = itmaxint;
	cnt++;
	this_array_[cnt] = itmaxso;
	cnt++;
	this_array_[cnt] = irecover;
	cnt++;
	this_array_[cnt] = isave;
	cnt++;
	this_array_[cnt] = icauchy;
	cnt++;
	this_array_[cnt] = interaction;
	cnt++;
	this_array_[cnt] = nunique;
	cnt++;
	this_array_[cnt] = nrsmin;
	cnt++;
	this_array_[cnt] = ihardlaw;
	cnt++;
	this_array_[cnt] = iratesens;
	cnt++;
	this_array_[cnt] = iupdori;
	cnt++;
	this_array_[cnt] = iupdshp;
	cnt++;
	this_array_[cnt] = iupdhar;
	cnt++;
	this_array_[cnt] = itemphard;
	cnt++;
	this_array_[cnt] = ibcinv;
	cnt++;
	this_array_[cnt] = irsvar;
	cnt++;
	this_array_[cnt] = ipr;
	cnt++;
	this_array_[cnt] = ipru;
	cnt++;
	this_array_[cnt] = kinctex;
	cnt++;
	this_array_[cnt] = neltex;
	cnt++;
	this_array_[cnt] = strain_control;
	cnt++;
	this_array_[cnt] = grbot;
	cnt++;
	this_array_[cnt] = grtop;
	cnt++;
	return *this;
}

KOKKOS_FUNCTION PropType::~PropType() {
}


KOKKOS_FUNCTION  StateVType::StateVType() {
}

KOKKOS_FUNCTION  StateVType& StateVType::init(double* this_array_, size_t& cnt)
{
	scalar_start_ = cnt;

	kinc = this_array_[cnt];
	cnt++;
	evmt = this_array_[cnt];
	cnt++;
	evmp = this_array_[cnt];
	cnt++;
	svm = this_array_[cnt];
	cnt++;
	dvm = this_array_[cnt];
	cnt++;
	dstranaccvm = this_array_[cnt];
	cnt++;
	dtimeacc = this_array_[cnt];
	cnt++;
	dstranvmth = this_array_[cnt];
	cnt++;
	temp_ini = this_array_[cnt];
	cnt++;
	temp_fact = this_array_[cnt];
	cnt++;
	temp_amp = this_array_[cnt];
	cnt++;
	temp = this_array_[cnt];
	cnt++;
	tempold = this_array_[cnt];
	cnt++;
	dtemp = this_array_[cnt];
	cnt++;
	evmpold = this_array_[cnt];
	cnt++;
	wrplastic = this_array_[cnt];
	cnt++;
	wplastic = this_array_[cnt];
	cnt++;
	jran = this_array_[cnt];
	cnt++;
	dbart.set(MATARIX <double>(&this_array_[cnt], 5));
	cnt += dbart.size();

	rot_acum.set(MATARIX <double>(&this_array_[cnt], 3, 3));
	cnt += rot_acum.size();

	xmtg.set(MATARIX <double>(&this_array_[cnt], 5, 5));
	cnt += xmtg.size();

	xltg.set(MATARIX <double>(&this_array_[cnt], 5, 5));
	cnt += xltg.size();

	dzero.set(MATARIX <double>(&this_array_[cnt], 5));
	cnt += dzero.size();

	eulerph.set(MATARIX <double>(&this_array_[cnt], NPH_ + 1, 3));
	cnt += eulerph.size();

	axisph.set(MATARIX <double>(&this_array_[cnt], NPH_ + 1, 4, 3));
	cnt += axisph.size();

	fijph.set(MATARIX <double>(&this_array_[cnt], NPH_ + 1, 3, 3));
	cnt += fijph.size();

	ag.set(MATARIX <double>(&this_array_[cnt], NGR_, 3, 3));
	cnt += ag.size();

	crss.set(MATARIX <double>(&this_array_[cnt], NGR_, NSYS_));
	cnt += crss.size();

	sg.set(MATARIX <double>(&this_array_[cnt], NGR_, 5));
	cnt += sg.size();

	sghyd.set(MATARIX <double>(&this_array_[cnt], NGR_));
	cnt += sghyd.size();

	gtotgr.set(MATARIX <double>(&this_array_[cnt], NGR_));
	cnt += gtotgr.size();

	xmctg.set(MATARIX <double>(&this_array_[cnt], NGR_, 5, 5));
	cnt += xmctg.size();

	dczero.set(MATARIX <double>(&this_array_[cnt], NGR_, 5));
	cnt += dczero.size();

	eftwfr.set(MATARIX <double>(&this_array_[cnt], NPH_, NTWM_));
	cnt += eftwfr.size();

	twfrph.set(MATARIX <double>(&this_array_[cnt], NPH_, NTWM_));
	cnt += twfrph.size();

	twfrsy.set(MATARIX <double>(&this_array_[cnt], NGR_, NTWS_));
	cnt += twfrsy.size();

	twfgr.set(MATARIX <double>(&this_array_[cnt], NTWM_, NGR_));
	cnt += twfgr.size();

	wgt.set(MATARIX <double>(&this_array_[cnt], NGR_));
	cnt += wgt.size();

	wph.set(MATARIX <double>(&this_array_[cnt], NPH_));
	cnt += wph.size();

	gamd0g.set(MATARIX <double>(&this_array_[cnt], NGR_));
	cnt += gamd0g.size();

	gamdot.set(MATARIX <double>(&this_array_[cnt], NGR_, NSYS_));
	cnt += gamdot.size();

	sch.set(MATARIX <double>(&this_array_[cnt], NGR_, 5, NSYS_));
	cnt += sch.size();

	udot.set(MATARIX <double>(&this_array_[cnt], 3, 3));
	cnt += udot.size();

	dav.set(MATARIX <double>(&this_array_[cnt], 5));
	cnt += dav.size();

	asph.set(MATARIX <double>(&this_array_[cnt], NPH_, 3, 3, 3, 3));
	cnt += asph.size();

	dg.set(MATARIX <double>(&this_array_[cnt], NGR_, 5));
	cnt += dg.size();

	dstranacc.set(MATARIX <double>(&this_array_[cnt], 3, 3));
	cnt += dstranacc.size();

	drotacc.set(MATARIX <double>(&this_array_[cnt], 3, 3));
	cnt += drotacc.size();

	xmtgdot.set(MATARIX <double>(&this_array_[cnt], 5, 5));
	cnt += xmtgdot.size();

	dzerodot.set(MATARIX <double>(&this_array_[cnt], 5));
	cnt += dzerodot.size();

	cscold.set(MATARIX <double>(&this_array_[cnt], 6, 6));
	cnt += cscold.size();

	stressold.set(MATARIX <double>(&this_array_[cnt], 3, 3));
	cnt += stressold.size();

	ktwsmx.set(MATARIX <double>(&this_array_[cnt], NGR_));
	cnt += ktwsmx.size();

	ntwevents.set(MATARIX <double>(&this_array_[cnt], NGR_));
	cnt += ntwevents.size();

	giph.set(MATARIX <double>(&this_array_[cnt], NGR_));
	cnt += giph.size();

	iflat.set(MATARIX <double>(&this_array_[cnt], NPH_ + 1));
	cnt += iflat.size();

	imark.set(MATARIX <double>(&this_array_[cnt], NGR_));
	cnt += imark.size();

	iptsgr.set(MATARIX <double>(&this_array_[cnt], NGR_));
	cnt += iptsgr.size();

	iptmgr.set(MATARIX <double>(&this_array_[cnt], NGR_));
	cnt += iptmgr.size();

	return *this;
}

KOKKOS_FUNCTION StateVType& StateVType::write(double* this_array_) {

	size_t cnt = scalar_start_;

	this_array_[cnt] = kinc;
	cnt++;
	this_array_[cnt] = evmt;
	cnt++;
	this_array_[cnt] = evmp;
	cnt++;
	this_array_[cnt] = svm;
	cnt++;
	this_array_[cnt] = dvm;
	cnt++;
	this_array_[cnt] = dstranaccvm;
	cnt++;
	this_array_[cnt] = dtimeacc;
	cnt++;
	this_array_[cnt] = dstranvmth;
	cnt++;
	this_array_[cnt] = temp_ini;
	cnt++;
	this_array_[cnt] = temp_fact;
	cnt++;
	this_array_[cnt] = temp_amp;
	cnt++;
	this_array_[cnt] = temp;
	cnt++;
	this_array_[cnt] = tempold;
	cnt++;
	this_array_[cnt] = dtemp;
	cnt++;
	this_array_[cnt] = evmpold;
	cnt++;
	this_array_[cnt] = wrplastic;
	cnt++;
	this_array_[cnt] = wplastic;
	cnt++;
	this_array_[cnt] = jran;
	cnt++;

	return *this;
}

KOKKOS_FUNCTION StateVType::~StateVType() {
}

KOKKOS_FUNCTION  IntVType::IntVType() {}

KOKKOS_FUNCTION  IntVType& IntVType::init(double* this_array_, size_t& cnt)
{
	sastav.set(MATARIX <double>(&this_array_[cnt], 5));
	cnt += sastav.size();

	sastbar.set(MATARIX <double>(&this_array_[cnt], 5));
	cnt += sastbar.size();

	dast.set(MATARIX <double>(&this_array_[cnt], 5));
	cnt += dast.size();

	sav.set(MATARIX <double>(&this_array_[cnt], 5));
	cnt += sav.size();

	dav.set(MATARIX <double>(&this_array_[cnt], 5));
	cnt += dav.size();

	sbar.set(MATARIX <double>(&this_array_[cnt], 5));
	cnt += sbar.size();

	dbar.set(MATARIX <double>(&this_array_[cnt], 5));
	cnt += dbar.size();

	xlijph.set(MATARIX <double>(&this_array_[cnt], NPH_, 3, 3));
	cnt += xlijph.size();

	asph.set(MATARIX <double>(&this_array_[cnt], NPH_, 3, 3, 3, 3));
	cnt += asph.size();

	xmastph.set(MATARIX <double>(&this_array_[cnt], NPH_, 5, 5));
	cnt += xmastph.size();

	dg.set(MATARIX <double>(&this_array_[cnt], NGR_, 5));
	cnt += dg.size();

	stry.set(MATARIX <double>(&this_array_[cnt], NGR_, 5));
	cnt += stry.size();

	scauchy.set(MATARIX <double>(&this_array_[cnt], 3, 3));
	cnt += scauchy.size();

	udot.set(MATARIX <double>(&this_array_[cnt], 3, 3));
	cnt += udot.size();

	dsim.set(MATARIX <double>(&this_array_[cnt], 3, 3));
	cnt += dsim.size();

	rotbar.set(MATARIX <double>(&this_array_[cnt], 3, 3));
	cnt += rotbar.size();

	sdeviat.set(MATARIX <double>(&this_array_[cnt], 3, 3));
	cnt += sdeviat.size();

	csc.set(MATARIX <double>(&this_array_[cnt], 6, 6));
	cnt += csc.size();

	ssc.set(MATARIX <double>(&this_array_[cnt], 6, 6));
	cnt += ssc.size();

	dsimtot.set(MATARIX <double>(&this_array_[cnt], 3, 3));
	cnt += dsimtot.size();

	sbarold.set(MATARIX <double>(&this_array_[cnt], 5));
	cnt += sbarold.size();

	cgr.set(MATARIX <double>(&this_array_[cnt], NGR_, 6, 6));
	cnt += cgr.size();


	svm = 0.0;
	dvm = 0.0;
	epsacu = 0.0;
	tincr = 0.0;
	jjxrs = 0.0;
	kinc = 0.0;
	k_cell = 0.0;

	return *this;
}

KOKKOS_FUNCTION  IntVType::~IntVType() {}


//Sum of the product of a1, a2, and a3, sz is their length
KOKKOS_FUNCTION double sumProd(const size_t sz, // length of the matrices
    double* a1, // matrix number 1
    double* a2, // matrix number 2
    double* a3) // martrix number 3
{
    double dum = 0.0;

    for (size_t i = 0; i < sz; i++) {
        dum += a1[i] * a2[i] * a3[i];
    }

    return dum;
}

//Sum of the product of a1 and a2, sz is their length
KOKKOS_FUNCTION double sumProd(const size_t sz,
    double* a1,
    double* a2) {
    double dum = 0.0;

    for (size_t i = 0; i < sz; i++) {
        dum += a1[i] * a2[i];
    }

    return dum;
}

//Sum of the product of a1 and a2, sz is their length, nc is the number of columns in b
KOKKOS_FUNCTION double sumProd(const size_t nk, const size_t nc, double* a, double* b) {
    double dum = 0.0;

    for (size_t i = 0; i < nk; i++) {
        dum += a[i] * b[i * nc];
    }

    return dum;
}


//Sum of the product of a1 and a2, sz is their length
KOKKOS_FUNCTION double kahanSumProd(const size_t sz,
    double* a1,
    double* a2) {
    double sum = 0.0, c = 0.0;

    for (size_t i = 0; i < sz; i++) {

        double y = (a1[i] * a2[i]) - c;
        double t = sum + y;

        c = (t - y);
        c -= (y);

        sum = t;
    }

    return sum;
}

KOKKOS_FUNCTION double tnorm(MATARIX <double>& mat) {
    double norm = 0;
    const size_t sz = mat.size();

    for (size_t i = 1; i <= sz; i++) {
        norm += mat(i) * mat(i);
    }
    norm = sqrt(norm);

    return norm;
}

KOKKOS_FUNCTION double max(MATARIX <double>& mat) {

    size_t len = mat.size();
    if (len == 0) return 0;

    double maxval = mat(1);

    for (size_t i = 2; i <= len; i++) {
        if (mat(i) > maxval) maxval = mat(i);
    }

    return maxval;
}




KOKKOS_FUNCTION double xid(size_t i, size_t j) {
    return (double)(i == j);
}


//KOKKOS_FUNCTION double sgn(double val) {
//	double out = (0.0 < val) - (val < 0.0);
//	if (out == 0.0) out = 1.0;
//	return out;
//}


KOKKOS_FUNCTION void operator+= (MATARIX <double>& a, double* temp) {

    double* aa = &a(1);
    const size_t N = a.size();

    for (size_t i = 0; i < N; i++) aa[i] += temp[i];

    return;
}

KOKKOS_FUNCTION void operator-= (MATARIX <double>& a, double* temp) {

    double* aa = &a(1);
    const size_t N = a.size();

    for (size_t i = 0; i < N; i++) aa[i] -= temp[i];

    return;
}

KOKKOS_FUNCTION void operator*= (MATARIX <double>& a, double* temp) {

    double* aa = &a(1);
    const size_t N = a.size();

    for (size_t i = 0; i < N; i++) aa[i] *= temp[i];

    return;
}

KOKKOS_FUNCTION void operator/= (MATARIX <double>& a, double* temp) {

    double* aa = &a(1);
    const size_t N = a.size();

    for (size_t i = 0; i < N; i++) aa[i] /= temp[i];

    return;
}


KOKKOS_FUNCTION void operator+= (MATARIX <double>& a, const double temp) {

    double* aa = &a(1);
    const size_t N = a.size();

    for (size_t i = 0; i < N; i++) aa[i] += temp;

    return;
}

KOKKOS_FUNCTION void operator-= (MATARIX <double>& a, const double temp) {

    double* aa = &a(1);
    const size_t N = a.size();

    for (size_t i = 0; i < N; i++) aa[i] -= temp;

    return;
}

KOKKOS_FUNCTION void operator*= (MATARIX <double>& a, const double temp) {

    double* aa = &a(1);
    const size_t N = a.size();

    for (size_t i = 0; i < N; i++) aa[i] *= temp;

    return;
}

KOKKOS_FUNCTION void operator/= (MATARIX <double>& a, const double temp) {

    double* aa = &a(1);
    const size_t N = a.size();
    const double tt = (1.0 / temp);

    for (size_t i = 0; i < N; i++) aa[i] *= tt;

    return;
}


KOKKOS_FUNCTION void operator+= (MATARRY <double>& a, const double temp) {

    double* aa = &a(0);
    const size_t N = a.size();

    for (size_t i = 0; i < N; i++) aa[i] += temp;

    return;
}
//*/
KOKKOS_FUNCTION void slice(const size_t N, double* a, double* some_matrix) {
    for (size_t i = 0; i < N; i++) a[i] = some_matrix[i];

    //memcpy(&a(1), &some_matrix[0], (8 * a.size()));

    return;
}

KOKKOS_FUNCTION void peye(MATARIX <double>& a, size_t N, double I) {
    const size_t N1 = N + 1;
    const size_t  L = a.size();
    double* aa = &a(1);

    for (size_t i = 0; i < L; i += N1) {
        aa[i] += I;
    }

    return;
}

KOKKOS_FUNCTION void eye(double* y, size_t n, double I) {
    size_t cnt = 0;
    y[cnt] = I;

    for (size_t i = 1; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            cnt++;
            y[cnt] = 0.0;
        }
        cnt++;
        y[cnt] = I;
    }
}

KOKKOS_FUNCTION void transpose(MATARIX <double>& a, // square matrix thatis being transposed
    const size_t N) // size of matrix that is being transposed
{
    //vmat2(double, b, sz1, sz2);
    double swap = 0;
    for (size_t i = 1; i <= N; i++) {
        for (size_t j = i + 1; j <= N; j++) {
            swap = a(i, j);
            a(i, j) = a(j, i);
            a(j, i) = swap;
        }
    }

    return;
}

KOKKOS_FUNCTION double det(MATARIX <double>& a)
// finds determinant of matrix a
{
    return (a(1) * a(5) * a(9) - a(1) * a(6) * a(8) - a(2) * a(4) * a(9) + a(2) * a(6) * a(7) + a(3) * a(4) * a(8) - a(3) * a(5) * a(7));
}


KOKKOS_FUNCTION double sgn(double val) {
	double out = 1.0;
	if (val < 0.0) out = -1.0;
	return out;
}


KOKKOS_FUNCTION double max(double* mat, const size_t len) {

	double maxval = 0.0;
	double maxsgn = 0.0;

	for (size_t i = 0; i < len; i++) {
		if (mat[i] > maxval) {
			maxval = mat[i];
			maxsgn = 1.0;
		}
		else if (mat[i] < -maxval) {
			maxval = -mat[i];
			maxsgn = -1.0;
		}
	}

	maxval *= maxsgn;

	return maxval;
}

KOKKOS_FUNCTION void rotate_mat(MATARRY <double> a, //3x3 matrix that is rotated by rot
	MATARRY <double> rot) {
	//rotates a 3x3 matrix by rotation rot, outputs to c2

	constexpr size_t N = 3;
	double c1_1D[9]; MATARRY <double> c1(&c1_1D[0], 3, 3);
	//fmat2(double, c1, N, N);

	for (size_t j = 0; j < N; j++) {
		for (size_t i = 0; i < N; i++) {
			double dum = 0.0;
			for (size_t k = 0; k < N; k++) {
				dum += a(j, k) * rot(i, k);
			}
			c1(j, i) = dum;
		}
	}

	for (size_t j = 0; j < N; j++) {
		for (size_t i = 0; i < N; i++) {
			double dum = 0.0;
			for (size_t k = 0; k < N; k++) {
				dum += c1(k, i) * rot(j, k);
			}
			a(j, i) = dum;
		}
	}

}

// Made by Matlab Coder
KOKKOS_FUNCTION void mat_eqsystem(double* A_data, double* b_data, const int N, size_t& ierr)
{
	double b_A_data[25];
	int ipiv_data[5];

	double smax;

	double amax;
	double bmax;

	int i;
	int jp1j;
	int k;
	int n;
	int u0;
	int yk;

	int len = N * N;

	ierr = 0;

	amax = 1.0 / max(A_data, len);
	bmax = 1.0 / max(b_data, N);
	
	if ((amax != amax)){
		return;
		
	}
	for (i = 0; i < len; i++) { A_data[i] *= amax; }
	for (i = 0; i < N; i++) { b_data[i] *= bmax; }
    

	for (i = 0; i < N; i++) {
		for (jp1j = 0; jp1j < N; jp1j++) {
			b_A_data[jp1j + N * i] = A_data[jp1j + N * i];//jp1j + N * i
		}
	}

	u0 = N;
	for (k = 0; k < N; k++) {
		ipiv_data[k] = k;
	}

	int Np1;
	Np1 = N;
	u0 = N - 1;
	n = N;

	for (int j = 0; j < u0; j++) {
		int b;
		int jA;
		int jj;
		int mmj_tmp;
		mmj_tmp = N - j;
		b = j * (N + 1);
		jj = j * (Np1 + 1);
		jp1j = b + 1;
		if (mmj_tmp < 1) {
			n = -1;
		}
		else {
			n = 0;
			if (mmj_tmp > 1) {
				smax = fabs(b_A_data[jj]);
				for (k = 1; k < mmj_tmp; k++) {
					double s;
					s = fabs(b_A_data[(b + k)]);
					if (s > smax) {
						n = k;
						smax = s;
					}
				}
			}
		}
		if (b_A_data[jj + n] != 0.0) {
			if (n != 0) {
				yk = j + n;
				ipiv_data[j] = yk;
				for (k = 0; k < N; k++) {
					n = k * N;
					jA = j + n;
					smax = b_A_data[jA];
					i = yk + n;
					b_A_data[jA] = b_A_data[i];
					b_A_data[i] = smax;
				}
			}
			i = jj + mmj_tmp;
			for (yk = jp1j; yk < i; yk++) {
				b_A_data[yk] /= b_A_data[jj];
			}
		}
		yk = b + N;
		jA = jj + Np1;
		for (b = 0; b <= mmj_tmp - 2; b++) {
			n = yk + b * N;
			smax = b_A_data[n];
			if (b_A_data[n] != 0.0) {
				i = jA + 1;
				jp1j = mmj_tmp + jA;
				for (n = i; n < jp1j; n++) {
					b_A_data[n] += b_A_data[((jj + n) - jA)] * -smax;
				}
			}
			jA += N;
		}
	}

	for (yk = 0; yk < N - 1; yk++) {
		i = ipiv_data[yk];
		if (i != yk) {
			smax = b_data[yk];
			b_data[yk] = b_data[i];
			b_data[i] = smax;
		}
	}

	for (k = 0; k < N; k++) {
		n = N * k;
		if (b_data[k] != 0.0) {
			i = k + 1;
			for (yk = i; yk < N; yk++) {
				b_data[yk] -= b_data[k] * b_A_data[(yk + n)];
			}
		}
	}

	for (k = N - 1; k >= 0; k--) {
		n = N * k;
		smax = b_data[k];
		if (smax != 0.0 && b_A_data[(k + n)] != 0.0) {
			b_data[k] = smax / b_A_data[(k + n)];
			for (yk = 0; yk < k; yk++) {
				b_data[yk] -= b_data[k] * b_A_data[yk + n];
			}
		}
	}

	amax = (amax / bmax);

	for (i = 0; i < N; i++) { b_data[i] *= amax; }


	 
}

// Made by Matlab Coder
KOKKOS_FUNCTION void mat_inverse(double* a_data, const int N, size_t& ierr)
{
	double x_data[36];
	double y_data[36];
	int ipiv_data[6];
	int p_data[6];

	int i;
	int j;
	int jA;
	int jp1j;
	int jy;
	int k;
	int n;
	int yk;

	int len = N * N;
	double amax = 0.0;


	ierr = 0.0;

	amax = max(a_data, len);

	if (amax == 0.0) {
		ierr = 1;

		return;
	}

	//for (i = 0; i < len; i++) {
	//	if (a_data[i] > amax) amax = a_data[i];
	//	else if (a_data[i] < -amax) amax = -a_data[i];
	//}

	amax = 1.0 / amax;

	for (i = 0; i < len; i++) a_data[i] *= amax;


	yk = N;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			//x(i, j) = a(j, i)
			x_data[i + N * j] = a_data[j + N * i];
		}
	}

	yk = len;

	for (i = 0; i < len; i++) {
		y_data[i] = 0.0;
	}

	n = N;

	for (k = 0; k < N; k++) {
		ipiv_data[k] = k + 1;
	}

	int u0 = N - 1;

	for (j = 0; j < u0; j++) {
		double smax;
		int b;
		int jj;
		int mmj_tmp;
		mmj_tmp = N - j;
		b = j * (N + 1);
		jj = j * (N + 1);
		jp1j = b + 2;
		if (mmj_tmp < 1) {
			yk = -1;
		}
		else {
			yk = 0;
			if (mmj_tmp > 1) {
				smax = fabs(x_data[jj]);
				for (k = 2; k <= mmj_tmp; k++) {
					double s;
					s = fabs(x_data[(b + k) - 1]);
					if (s > smax) {
						yk = k - 1;
						smax = s;
					}
				}
			}
		}
		if (x_data[jj + yk] != 0.0) {
			if (yk != 0) {
				jy = j + yk;
				ipiv_data[j] = jy + 1;
				for (k = 0; k < N; k++) {
					yk = k * N;
					jA = j + yk;
					smax = x_data[jA];
					i = jy + yk;
					x_data[jA] = x_data[i];
					x_data[i] = smax;
				}
			}
			i = jj + mmj_tmp;
			for (jA = jp1j; jA <= i; jA++) {
				x_data[jA - 1] /= x_data[jj];
			}
		}
		jy = b + N;
		jA = jj + N;
		for (b = 0; b <= mmj_tmp - 2; b++) {
			yk = jy + b * N;
			smax = x_data[yk];
			if (x_data[yk] != 0.0) {
				i = jA + 2;
				jp1j = mmj_tmp + jA;
				for (yk = i; yk <= jp1j; yk++) {
					x_data[yk - 1] += x_data[((jj + yk) - jA) - 1] * -smax;
				}
			}
			jA += N;
		}
	}
	jy = N;

	if (jy > 0) {
		p_data[0] = 1;
		yk = 1;
		for (k = 2; k <= jy; k++) {
			yk++;
			p_data[k - 1] = yk;
		}
	}
	for (k = 0; k < n; k++) {
		i = ipiv_data[k];
		if (i > k + 1) {
			yk = p_data[i - 1];
			p_data[i - 1] = p_data[k];
			p_data[k] = yk;
		}
	}
	for (k = 0; k < N; k++) {
		yk = N * (p_data[k] - 1);
		y_data[k + yk] = 1.0;
		for (j = k + 1; j <= N; j++) {
			i = (j + yk) - 1;
			if (y_data[i] != 0.0) {
				jp1j = j + 1;
				for (jA = jp1j; jA <= N; jA++) {
					jy = (jA + yk) - 1;
					y_data[jy] -= y_data[i] * x_data[(jA + N * (j - 1)) - 1];
				}
			}
		}
	}

	for (j = 0; j < N; j++) {
		yk = N * j - 1;
		for (k = N; k >= 1; k--) {
			jy = N * (k - 1) - 1;
			i = k + yk;

			if (y_data[i] != 0.0 && x_data[k + jy] != 0.0) {
				y_data[i] /= x_data[k + jy];
				for (jA = 0; jA <= k - 2; jA++) {
					jp1j = (jA + yk) + 1;
					y_data[jp1j] -= y_data[i] * x_data[(jA + jy) + 1];
				}
			}
		}
	}

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			a_data[j + N * i] = y_data[i + N * j];
		}
	}



	for (i = 0; i < len; i++) a_data[i] *= amax;
	 
}

KOKKOS_FUNCTION void matinv3(MATARIX <double>& cinv, MATARIX <double>& c)
{

	//real(kind=dp) c(3,3), cinv(3,3)

	//const double i3 = 1.0 / (c(1, 1) * c(2, 2) * c(3, 3) + c(1, 2) * c(2, 3) * c(3, 1) + c(1, 3) * c(2, 1) * c(3, 2) - c(1, 3) * c(2, 2) * c(3, 1) - c(2, 3) * c(3, 2) * c(1, 1) - c(3, 3) * c(1, 2) * c(2, 1));
	//
	//cinv(1, 1) = i3 * (c(2, 2) * c(3, 3) - c(2, 3) * c(3, 2));
	//cinv(2, 2) = i3 * (c(1, 1) * c(3, 3) - c(1, 3) * c(3, 1));
	//cinv(3, 3) = i3 * (c(1, 1) * c(2, 2) - c(1, 2) * c(2, 1));
	//
	//cinv(1, 2) = -i3 * (c(1, 2) * c(3, 3) - c(3, 2) * c(1, 3));
	//cinv(2, 1) = -i3 * (c(2, 1) * c(3, 3) - c(2, 3) * c(3, 1));
	//
	//cinv(1, 3) = i3 * (c(1, 2) * c(2, 3) - c(1, 3) * c(2, 2));
	//cinv(3, 1) = i3 * (c(2, 1) * c(3, 2) - c(2, 2) * c(3, 1));
	//
	//cinv(2, 3) = -i3 * (c(1, 1) * c(2, 3) - c(1, 3) * c(2, 1));
	//cinv(3, 2) = -i3 * (c(1, 1) * c(3, 2) - c(1, 2) * c(3, 1));

	double i3 = 1.0 / (c(1) * c(5) * c(9) + c(2) * c(6) * c(7) + c(3) * c(4) * c(8) - c(3) * c(5) * c(7) - c(6) * c(8) * c(1) - c(9) * c(2) * c(4));

	cinv(1) = i3 * (c(5) * c(9) - c(6) * c(8));
	cinv(2) = -i3 * (c(2) * c(9) - c(8) * c(3));
	cinv(3) = i3 * (c(2) * c(6) - c(3) * c(5));
	cinv(4) = -i3 * (c(4) * c(9) - c(6) * c(7));
	cinv(5) = i3 * (c(1) * c(9) - c(3) * c(7));
	cinv(6) = -i3 * (c(1) * c(6) - c(3) * c(4));
	cinv(7) = i3 * (c(4) * c(8) - c(5) * c(7));
	cinv(8) = -i3 * (c(1) * c(8) - c(2) * c(7));
	cinv(9) = i3 * (c(1) * c(5) - c(2) * c(4));
	 
}

KOKKOS_FUNCTION void chg_basis15(double* ce, double* c) {

	c[0] = -0.70710678118654752440084436210485 * ce[0] - 0.40824829046386301636621401245098 * ce[1];
	c[1] = 0.70710678118654752440084436210485 * ce[4];
	c[2] = 0.70710678118654752440084436210485 * ce[3];
	c[3] = 0.70710678118654752440084436210485 * ce[4];
	c[4] = 0.70710678118654752440084436210485 * ce[0] - 0.40824829046386301636621401245098 * ce[1];
	c[5] = 0.70710678118654752440084436210485 * ce[2];
	c[6] = 0.70710678118654752440084436210485 * ce[3];
	c[7] = 0.70710678118654752440084436210485 * ce[2];
	c[8] = 0.81649658092772603273242802490196 * ce[1];
	 
}

KOKKOS_FUNCTION void chg_basis25(double* ce, double* c) {

	ce[0] = 0.70710678118654752440084436210485 * c[4] - 0.70710678118654752440084436210485 * c[0];
	ce[1] = 0.81649658092772603273242802490196 * c[8] - 0.40824829046386301636621401245098 * c[4] - 0.40824829046386301636621401245098 * c[0];
	ce[2] = 0.70710678118654752440084436210485 * c[5] + 0.70710678118654752440084436210485 * c[7];
	ce[3] = 0.70710678118654752440084436210485 * c[2] + 0.70710678118654752440084436210485 * c[6];
	ce[4] = 0.70710678118654752440084436210485 * c[1] + 0.70710678118654752440084436210485 * c[3];
	 
}

KOKKOS_FUNCTION void chg_basis35(double* ce, double* c) {
	c[0] = 0.5 * ce[0] + 0.28867513459481288225457439025098 * ce[1] + 0.28867513459481288225457439025098 * ce[5] + 0.16666666666666666666666666666667 * ce[6];
	c[1] = -0.5 * ce[4] - 0.28867513459481288225457439025098 * ce[9];
	c[2] = -0.5 * ce[3] - 0.28867513459481288225457439025098 * ce[8];
	c[3] = -0.5 * ce[4] - 0.28867513459481288225457439025098 * ce[9];
	c[4] = 0.28867513459481288225457439025098 * ce[1] - 0.5 * ce[0] - 0.28867513459481288225457439025098 * ce[5] + 0.16666666666666666666666666666667 * ce[6];
	c[5] = -0.5 * ce[2] - 0.28867513459481288225457439025098 * ce[7];
	c[6] = -0.5 * ce[3] - 0.28867513459481288225457439025098 * ce[8];
	c[7] = -0.5 * ce[2] - 0.28867513459481288225457439025098 * ce[7];
	c[8] = -0.57735026918962576450914878050196 * ce[1] - 0.33333333333333333333333333333333 * ce[6];
	c[9] = -0.5 * ce[20] - 0.28867513459481288225457439025098 * ce[21];
	c[10] = 0.5 * ce[24];
	c[11] = 0.5 * ce[23];
	c[12] = 0.5 * ce[24];
	c[13] = 0.5 * ce[20] - 0.28867513459481288225457439025098 * ce[21];
	c[14] = 0.5 * ce[22];
	c[15] = 0.5 * ce[23];
	c[16] = 0.5 * ce[22];
	c[17] = 0.57735026918962576450914878050196 * ce[21];
	c[18] = -0.5 * ce[15] - 0.28867513459481288225457439025098 * ce[16];
	c[19] = 0.5 * ce[19];
	c[20] = 0.5 * ce[18];
	c[21] = 0.5 * ce[19];
	c[22] = 0.5 * ce[15] - 0.28867513459481288225457439025098 * ce[16];
	c[23] = 0.5 * ce[17];
	c[24] = 0.5 * ce[18];
	c[25] = 0.5 * ce[17];
	c[26] = 0.57735026918962576450914878050196 * ce[16];
	c[27] = -0.5 * ce[20] - 0.28867513459481288225457439025098 * ce[21];
	c[28] = 0.5 * ce[24];
	c[29] = 0.5 * ce[23];
	c[30] = 0.5 * ce[24];
	c[31] = 0.5 * ce[20] - 0.28867513459481288225457439025098 * ce[21];
	c[32] = 0.5 * ce[22];
	c[33] = 0.5 * ce[23];
	c[34] = 0.5 * ce[22];
	c[35] = 0.57735026918962576450914878050196 * ce[21];
	c[36] = 0.28867513459481288225457439025098 * ce[5] - 0.28867513459481288225457439025098 * ce[1] - 0.5 * ce[0] + 0.16666666666666666666666666666667 * ce[6];
	c[37] = 0.5 * ce[4] - 0.28867513459481288225457439025098 * ce[9];
	c[38] = 0.5 * ce[3] - 0.28867513459481288225457439025098 * ce[8];
	c[39] = 0.5 * ce[4] - 0.28867513459481288225457439025098 * ce[9];
	c[40] = 0.5 * ce[0] - 0.28867513459481288225457439025098 * ce[1] - 0.28867513459481288225457439025098 * ce[5] + 0.16666666666666666666666666666667 * ce[6];
	c[41] = 0.5 * ce[2] - 0.28867513459481288225457439025098 * ce[7];
	c[42] = 0.5 * ce[3] - 0.28867513459481288225457439025098 * ce[8];
	c[43] = 0.5 * ce[2] - 0.28867513459481288225457439025098 * ce[7];
	c[44] = 0.57735026918962576450914878050196 * ce[1] - 0.33333333333333333333333333333333 * ce[6];
	c[45] = -0.5 * ce[10] - 0.28867513459481288225457439025098 * ce[11];
	c[46] = 0.5 * ce[14];
	c[47] = 0.5 * ce[13];
	c[48] = 0.5 * ce[14];
	c[49] = 0.5 * ce[10] - 0.28867513459481288225457439025098 * ce[11];
	c[50] = 0.5 * ce[12];
	c[51] = 0.5 * ce[13];
	c[52] = 0.5 * ce[12];
	c[53] = 0.57735026918962576450914878050196 * ce[11];
	c[54] = -0.5 * ce[15] - 0.28867513459481288225457439025098 * ce[16];
	c[55] = 0.5 * ce[19];
	c[56] = 0.5 * ce[18];
	c[57] = 0.5 * ce[19];
	c[58] = 0.5 * ce[15] - 0.28867513459481288225457439025098 * ce[16];
	c[59] = 0.5 * ce[17];
	c[60] = 0.5 * ce[18];
	c[61] = 0.5 * ce[17];
	c[62] = 0.57735026918962576450914878050196 * ce[16];
	c[63] = -0.5 * ce[10] - 0.28867513459481288225457439025098 * ce[11];
	c[64] = 0.5 * ce[14];
	c[65] = 0.5 * ce[13];
	c[66] = 0.5 * ce[14];
	c[67] = 0.5 * ce[10] - 0.28867513459481288225457439025098 * ce[11];
	c[68] = 0.5 * ce[12];
	c[69] = 0.5 * ce[13];
	c[70] = 0.5 * ce[12];
	c[71] = 0.57735026918962576450914878050196 * ce[11];
	c[72] = -0.57735026918962576450914878050196 * ce[5] - 0.33333333333333333333333333333333 * ce[6];
	c[73] = 0.57735026918962576450914878050196 * ce[9];
	c[74] = 0.57735026918962576450914878050196 * ce[8];
	c[75] = 0.57735026918962576450914878050196 * ce[9];
	c[76] = 0.57735026918962576450914878050196 * ce[5] - 0.33333333333333333333333333333333 * ce[6];
	c[77] = 0.57735026918962576450914878050196 * ce[7];
	c[78] = 0.57735026918962576450914878050196 * ce[8];
	c[79] = 0.57735026918962576450914878050196 * ce[7];
	c[80] = 0.66666666666666666666666666666667 * ce[6];

	 
}

KOKKOS_FUNCTION void chg_basis45(double* ce, double* c) {
	ce[0] = 0.5 * c[0] - 0.5 * c[4] - 0.5 * c[36] + 0.5 * c[40];
	ce[1] = 0.28867513459481288225457439025098 * c[0] + 0.28867513459481288225457439025098 * c[4] - 0.57735026918962576450914878050196 * c[8] - 0.28867513459481288225457439025098 * c[36] - 0.28867513459481288225457439025098 * c[40] + 0.57735026918962576450914878050196 * c[44];
	ce[2] = 0.5 * c[41] - 0.5 * c[7] - 0.5 * c[5] + 0.5 * c[43];
	ce[3] = 0.5 * c[38] - 0.5 * c[6] - 0.5 * c[2] + 0.5 * c[42];
	ce[4] = 0.5 * c[37] - 0.5 * c[3] - 0.5 * c[1] + 0.5 * c[39];
	ce[5] = 0.28867513459481288225457439025098 * c[0] - 0.28867513459481288225457439025098 * c[4] + 0.28867513459481288225457439025098 * c[36] - 0.28867513459481288225457439025098 * c[40] - 0.57735026918962576450914878050196 * c[72] + 0.57735026918962576450914878050196 * c[76];
	ce[6] = 0.16666666666666666666666666666667 * c[0] + 0.16666666666666666666666666666667 * c[4] - 0.33333333333333333333333333333333 * c[8] + 0.16666666666666666666666666666667 * c[36] + 0.16666666666666666666666666666667 * c[40] - 0.33333333333333333333333333333333 * c[44] - 0.33333333333333333333333333333333 * c[72] - 0.33333333333333333333333333333333 * c[76] + 0.66666666666666666666666666666667 * c[80];
	ce[7] = 0.57735026918962576450914878050196 * c[77] - 0.28867513459481288225457439025098 * c[7] - 0.28867513459481288225457439025098 * c[41] - 0.28867513459481288225457439025098 * c[43] - 0.28867513459481288225457439025098 * c[5] + 0.57735026918962576450914878050196 * c[79];
	ce[8] = 0.57735026918962576450914878050196 * c[74] - 0.28867513459481288225457439025098 * c[6] - 0.28867513459481288225457439025098 * c[38] - 0.28867513459481288225457439025098 * c[42] - 0.28867513459481288225457439025098 * c[2] + 0.57735026918962576450914878050196 * c[78];
	ce[9] = 0.57735026918962576450914878050196 * c[73] - 0.28867513459481288225457439025098 * c[3] - 0.28867513459481288225457439025098 * c[37] - 0.28867513459481288225457439025098 * c[39] - 0.28867513459481288225457439025098 * c[1] + 0.57735026918962576450914878050196 * c[75];
	ce[10] = 0.5 * c[49] - 0.5 * c[45] - 0.5 * c[63] + 0.5 * c[67];
	ce[11] = 0.57735026918962576450914878050196 * c[53] - 0.28867513459481288225457439025098 * c[49] - 0.28867513459481288225457439025098 * c[45] - 0.28867513459481288225457439025098 * c[63] - 0.28867513459481288225457439025098 * c[67] + 0.57735026918962576450914878050196 * c[71];
	ce[12] = 0.5 * c[50] + 0.5 * c[52] + 0.5 * c[68] + 0.5 * c[70];
	ce[13] = 0.5 * c[47] + 0.5 * c[51] + 0.5 * c[65] + 0.5 * c[69];
	ce[14] = 0.5 * c[46] + 0.5 * c[48] + 0.5 * c[64] + 0.5 * c[66];
	ce[15] = 0.5 * c[22] - 0.5 * c[18] - 0.5 * c[54] + 0.5 * c[58];
	ce[16] = 0.57735026918962576450914878050196 * c[26] - 0.28867513459481288225457439025098 * c[22] - 0.28867513459481288225457439025098 * c[18] - 0.28867513459481288225457439025098 * c[54] - 0.28867513459481288225457439025098 * c[58] + 0.57735026918962576450914878050196 * c[62];
	ce[17] = 0.5 * c[23] + 0.5 * c[25] + 0.5 * c[59] + 0.5 * c[61];
	ce[18] = 0.5 * c[20] + 0.5 * c[24] + 0.5 * c[56] + 0.5 * c[60];
	ce[19] = 0.5 * c[19] + 0.5 * c[21] + 0.5 * c[55] + 0.5 * c[57];
	ce[20] = 0.5 * c[13] - 0.5 * c[9] - 0.5 * c[27] + 0.5 * c[31];
	ce[21] = 0.57735026918962576450914878050196 * c[17] - 0.28867513459481288225457439025098 * c[13] - 0.28867513459481288225457439025098 * c[9] - 0.28867513459481288225457439025098 * c[27] - 0.28867513459481288225457439025098 * c[31] + 0.57735026918962576450914878050196 * c[35];
	ce[22] = 0.5 * c[14] + 0.5 * c[16] + 0.5 * c[32] + 0.5 * c[34];
	ce[23] = 0.5 * c[11] + 0.5 * c[15] + 0.5 * c[29] + 0.5 * c[33];
	ce[24] = 0.5 * c[10] + 0.5 * c[12] + 0.5 * c[28] + 0.5 * c[30];

	 
}

KOKKOS_FUNCTION void chg_basis16(double* ce, double* c) {
	c[0] = -0.70710678118654752440084436210485 * ce[0] - 0.40824829046386301636621401245098 * ce[1] + 0.57735026918962576450914878050196 * ce[5];
	c[1] = 0.70710678118654752440084436210485 * ce[4];
	c[2] = 0.70710678118654752440084436210485 * ce[3];
	c[3] = 0.70710678118654752440084436210485 * ce[4];
	c[4] = 0.70710678118654752440084436210485 * ce[0] - 0.40824829046386301636621401245098 * ce[1] + 0.57735026918962576450914878050196 * ce[5];
	c[5] = 0.70710678118654752440084436210485 * ce[2];
	c[6] = 0.70710678118654752440084436210485 * ce[3];
	c[7] = 0.70710678118654752440084436210485 * ce[2];
	c[8] = 0.81649658092772603273242802490196 * ce[1] + 0.57735026918962576450914878050196 * ce[5];

	 
}

KOKKOS_FUNCTION void chg_basis26(double* ce, double* c) {
	ce[0] = -0.70710678118654752440084436210485 * c[0] + 0.70710678118654752440084436210485 * c[4];
	ce[1] = -0.40824829046386301636621401245098 * c[0] - 0.40824829046386301636621401245098 * c[4] + 0.81649658092772603273242802490196 * c[8];
	ce[2] = 0.70710678118654752440084436210485 * c[5] + 0.70710678118654752440084436210485 * c[7];
	ce[3] = 0.70710678118654752440084436210485 * c[2] + 0.70710678118654752440084436210485 * c[6];
	ce[4] = 0.70710678118654752440084436210485 * c[1] + 0.70710678118654752440084436210485 * c[3];
	ce[5] = 0.57735026918962576450914878050196 * c[0] + 0.57735026918962576450914878050196 * c[4] + 0.57735026918962576450914878050196 * c[8];

	 
}

KOKKOS_FUNCTION void chg_basis36(double* ce, double* c) {
	c[0] = 0.5 * ce[0] + 0.28867513459481288225457439025098 * ce[1] - 0.40824829046386301636621401245098 * ce[5] + 0.28867513459481288225457439025098 * ce[6] + 0.16666666666666666666666666666667 * ce[7] - 0.23570226039551584146694812070162 * ce[11] - 0.40824829046386301636621401245098 * ce[30] - 0.23570226039551584146694812070162 * ce[31] + 0.33333333333333333333333333333333 * ce[35];
	c[1] = 0.40824829046386301636621401245098 * ce[34] - 0.28867513459481288225457439025098 * ce[10] - 0.5 * ce[4];
	c[2] = 0.40824829046386301636621401245098 * ce[33] - 0.28867513459481288225457439025098 * ce[9] - 0.5 * ce[3];
	c[3] = 0.40824829046386301636621401245098 * ce[34] - 0.28867513459481288225457439025098 * ce[10] - 0.5 * ce[4];
	c[4] = 0.28867513459481288225457439025098 * ce[1] - 0.5 * ce[0] - 0.40824829046386301636621401245098 * ce[5] - 0.28867513459481288225457439025098 * ce[6] + 0.16666666666666666666666666666667 * ce[7] - 0.23570226039551584146694812070162 * ce[11] + 0.40824829046386301636621401245098 * ce[30] - 0.23570226039551584146694812070162 * ce[31] + 0.33333333333333333333333333333333 * ce[35];
	c[5] = 0.40824829046386301636621401245098 * ce[32] - 0.28867513459481288225457439025098 * ce[8] - 0.5 * ce[2];
	c[6] = 0.40824829046386301636621401245098 * ce[33] - 0.28867513459481288225457439025098 * ce[9] - 0.5 * ce[3];
	c[7] = 0.40824829046386301636621401245098 * ce[32] - 0.28867513459481288225457439025098 * ce[8] - 0.5 * ce[2];
	c[8] = 0.47140452079103168293389624140323 * ce[31] - 0.40824829046386301636621401245098 * ce[5] - 0.33333333333333333333333333333333 * ce[7] - 0.23570226039551584146694812070162 * ce[11] - 0.57735026918962576450914878050196 * ce[1] + 0.33333333333333333333333333333333 * ce[35];
	c[9] = 0.40824829046386301636621401245098 * ce[29] - 0.28867513459481288225457439025098 * ce[25] - 0.5 * ce[24];
	c[10] = 0.5 * ce[28];
	c[11] = 0.5 * ce[27];
	c[12] = 0.5 * ce[28];
	c[13] = 0.5 * ce[24] - 0.28867513459481288225457439025098 * ce[25] + 0.40824829046386301636621401245098 * ce[29];
	c[14] = 0.5 * ce[26];
	c[15] = 0.5 * ce[27];
	c[16] = 0.5 * ce[26];
	c[17] = 0.57735026918962576450914878050196 * ce[25] + 0.40824829046386301636621401245098 * ce[29];
	c[18] = 0.40824829046386301636621401245098 * ce[23] - 0.28867513459481288225457439025098 * ce[19] - 0.5 * ce[18];
	c[19] = 0.5 * ce[22];
	c[20] = 0.5 * ce[21];
	c[21] = 0.5 * ce[22];
	c[22] = 0.5 * ce[18] - 0.28867513459481288225457439025098 * ce[19] + 0.40824829046386301636621401245098 * ce[23];
	c[23] = 0.5 * ce[20];
	c[24] = 0.5 * ce[21];
	c[25] = 0.5 * ce[20];
	c[26] = 0.57735026918962576450914878050196 * ce[19] + 0.40824829046386301636621401245098 * ce[23];
	c[27] = 0.40824829046386301636621401245098 * ce[29] - 0.28867513459481288225457439025098 * ce[25] - 0.5 * ce[24];
	c[28] = 0.5 * ce[28];
	c[29] = 0.5 * ce[27];
	c[30] = 0.5 * ce[28];
	c[31] = 0.5 * ce[24] - 0.28867513459481288225457439025098 * ce[25] + 0.40824829046386301636621401245098 * ce[29];
	c[32] = 0.5 * ce[26];
	c[33] = 0.5 * ce[27];
	c[34] = 0.5 * ce[26];
	c[35] = 0.57735026918962576450914878050196 * ce[25] + 0.40824829046386301636621401245098 * ce[29];
	c[36] = 0.40824829046386301636621401245098 * ce[5] - 0.28867513459481288225457439025098 * ce[1] - 0.5 * ce[0] + 0.28867513459481288225457439025098 * ce[6] + 0.16666666666666666666666666666667 * ce[7] - 0.23570226039551584146694812070162 * ce[11] - 0.40824829046386301636621401245098 * ce[30] - 0.23570226039551584146694812070162 * ce[31] + 0.33333333333333333333333333333333 * ce[35];
	c[37] = 0.5 * ce[4] - 0.28867513459481288225457439025098 * ce[10] + 0.40824829046386301636621401245098 * ce[34];
	c[38] = 0.5 * ce[3] - 0.28867513459481288225457439025098 * ce[9] + 0.40824829046386301636621401245098 * ce[33];
	c[39] = 0.5 * ce[4] - 0.28867513459481288225457439025098 * ce[10] + 0.40824829046386301636621401245098 * ce[34];
	c[40] = 0.5 * ce[0] - 0.28867513459481288225457439025098 * ce[1] + 0.40824829046386301636621401245098 * ce[5] - 0.28867513459481288225457439025098 * ce[6] + 0.16666666666666666666666666666667 * ce[7] - 0.23570226039551584146694812070162 * ce[11] + 0.40824829046386301636621401245098 * ce[30] - 0.23570226039551584146694812070162 * ce[31] + 0.33333333333333333333333333333333 * ce[35];
	c[41] = 0.5 * ce[2] - 0.28867513459481288225457439025098 * ce[8] + 0.40824829046386301636621401245098 * ce[32];
	c[42] = 0.5 * ce[3] - 0.28867513459481288225457439025098 * ce[9] + 0.40824829046386301636621401245098 * ce[33];
	c[43] = 0.5 * ce[2] - 0.28867513459481288225457439025098 * ce[8] + 0.40824829046386301636621401245098 * ce[32];
	c[44] = 0.57735026918962576450914878050196 * ce[1] + 0.40824829046386301636621401245098 * ce[5] - 0.33333333333333333333333333333333 * ce[7] - 0.23570226039551584146694812070162 * ce[11] + 0.47140452079103168293389624140323 * ce[31] + 0.33333333333333333333333333333333 * ce[35];
	c[45] = 0.40824829046386301636621401245098 * ce[17] - 0.28867513459481288225457439025098 * ce[13] - 0.5 * ce[12];
	c[46] = 0.5 * ce[16];
	c[47] = 0.5 * ce[15];
	c[48] = 0.5 * ce[16];
	c[49] = 0.5 * ce[12] - 0.28867513459481288225457439025098 * ce[13] + 0.40824829046386301636621401245098 * ce[17];
	c[50] = 0.5 * ce[14];
	c[51] = 0.5 * ce[15];
	c[52] = 0.5 * ce[14];
	c[53] = 0.57735026918962576450914878050196 * ce[13] + 0.40824829046386301636621401245098 * ce[17];
	c[54] = 0.40824829046386301636621401245098 * ce[23] - 0.28867513459481288225457439025098 * ce[19] - 0.5 * ce[18];
	c[55] = 0.5 * ce[22];
	c[56] = 0.5 * ce[21];
	c[57] = 0.5 * ce[22];
	c[58] = 0.5 * ce[18] - 0.28867513459481288225457439025098 * ce[19] + 0.40824829046386301636621401245098 * ce[23];
	c[59] = 0.5 * ce[20];
	c[60] = 0.5 * ce[21];
	c[61] = 0.5 * ce[20];
	c[62] = 0.57735026918962576450914878050196 * ce[19] + 0.40824829046386301636621401245098 * ce[23];
	c[63] = 0.40824829046386301636621401245098 * ce[17] - 0.28867513459481288225457439025098 * ce[13] - 0.5 * ce[12];
	c[64] = 0.5 * ce[16];
	c[65] = 0.5 * ce[15];
	c[66] = 0.5 * ce[16];
	c[67] = 0.5 * ce[12] - 0.28867513459481288225457439025098 * ce[13] + 0.40824829046386301636621401245098 * ce[17];
	c[68] = 0.5 * ce[14];
	c[69] = 0.5 * ce[15];
	c[70] = 0.5 * ce[14];
	c[71] = 0.57735026918962576450914878050196 * ce[13] + 0.40824829046386301636621401245098 * ce[17];
	c[72] = 0.47140452079103168293389624140323 * ce[11] - 0.33333333333333333333333333333333 * ce[7] - 0.57735026918962576450914878050196 * ce[6] - 0.40824829046386301636621401245098 * ce[30] - 0.23570226039551584146694812070162 * ce[31] + 0.33333333333333333333333333333333 * ce[35];
	c[73] = 0.57735026918962576450914878050196 * ce[10] + 0.40824829046386301636621401245098 * ce[34];
	c[74] = 0.57735026918962576450914878050196 * ce[9] + 0.40824829046386301636621401245098 * ce[33];
	c[75] = 0.57735026918962576450914878050196 * ce[10] + 0.40824829046386301636621401245098 * ce[34];
	c[76] = 0.57735026918962576450914878050196 * ce[6] - 0.33333333333333333333333333333333 * ce[7] + 0.47140452079103168293389624140323 * ce[11] + 0.40824829046386301636621401245098 * ce[30] - 0.23570226039551584146694812070162 * ce[31] + 0.33333333333333333333333333333333 * ce[35];
	c[77] = 0.57735026918962576450914878050196 * ce[8] + 0.40824829046386301636621401245098 * ce[32];
	c[78] = 0.57735026918962576450914878050196 * ce[9] + 0.40824829046386301636621401245098 * ce[33];
	c[79] = 0.57735026918962576450914878050196 * ce[8] + 0.40824829046386301636621401245098 * ce[32];
	c[80] = 0.66666666666666666666666666666667 * ce[7] + 0.47140452079103168293389624140323 * ce[11] + 0.47140452079103168293389624140323 * ce[31] + 0.33333333333333333333333333333333 * ce[35];

	 
}

KOKKOS_FUNCTION void chg_basis46(double* ce, double* c) {
	ce[0] = 0.5 * c[0] - 0.5 * c[4] - 0.5 * c[36] + 0.5 * c[40];
	ce[1] = 0.28867513459481288225457439025098 * c[0] + 0.28867513459481288225457439025098 * c[4] - 0.57735026918962576450914878050196 * c[8] - 0.28867513459481288225457439025098 * c[36] - 0.28867513459481288225457439025098 * c[40] + 0.57735026918962576450914878050196 * c[44];
	ce[2] = 0.5 * c[41] - 0.5 * c[7] - 0.5 * c[5] + 0.5 * c[43];
	ce[3] = 0.5 * c[38] - 0.5 * c[6] - 0.5 * c[2] + 0.5 * c[42];
	ce[4] = 0.5 * c[37] - 0.5 * c[3] - 0.5 * c[1] + 0.5 * c[39];
	ce[5] = 0.40824829046386301636621401245098 * c[36] - 0.40824829046386301636621401245098 * c[4] - 0.40824829046386301636621401245098 * c[8] - 0.40824829046386301636621401245098 * c[0] + 0.40824829046386301636621401245098 * c[40] + 0.40824829046386301636621401245098 * c[44];
	ce[6] = 0.28867513459481288225457439025098 * c[0] - 0.28867513459481288225457439025098 * c[4] + 0.28867513459481288225457439025098 * c[36] - 0.28867513459481288225457439025098 * c[40] - 0.57735026918962576450914878050196 * c[72] + 0.57735026918962576450914878050196 * c[76];
	ce[7] = 0.16666666666666666666666666666667 * c[0] + 0.16666666666666666666666666666667 * c[4] - 0.33333333333333333333333333333333 * c[8] + 0.16666666666666666666666666666667 * c[36] + 0.16666666666666666666666666666667 * c[40] - 0.33333333333333333333333333333333 * c[44] - 0.33333333333333333333333333333333 * c[72] - 0.33333333333333333333333333333333 * c[76] + 0.66666666666666666666666666666667 * c[80];
	ce[8] = 0.57735026918962576450914878050196 * c[77] - 0.28867513459481288225457439025098 * c[7] - 0.28867513459481288225457439025098 * c[41] - 0.28867513459481288225457439025098 * c[43] - 0.28867513459481288225457439025098 * c[5] + 0.57735026918962576450914878050196 * c[79];
	ce[9] = 0.57735026918962576450914878050196 * c[74] - 0.28867513459481288225457439025098 * c[6] - 0.28867513459481288225457439025098 * c[38] - 0.28867513459481288225457439025098 * c[42] - 0.28867513459481288225457439025098 * c[2] + 0.57735026918962576450914878050196 * c[78];
	ce[10] = 0.57735026918962576450914878050196 * c[73] - 0.28867513459481288225457439025098 * c[3] - 0.28867513459481288225457439025098 * c[37] - 0.28867513459481288225457439025098 * c[39] - 0.28867513459481288225457439025098 * c[1] + 0.57735026918962576450914878050196 * c[75];
	ce[11] = 0.47140452079103168293389624140323 * c[72] - 0.23570226039551584146694812070162 * c[4] - 0.23570226039551584146694812070162 * c[8] - 0.23570226039551584146694812070162 * c[36] - 0.23570226039551584146694812070162 * c[40] - 0.23570226039551584146694812070162 * c[44] - 0.23570226039551584146694812070162 * c[0] + 0.47140452079103168293389624140323 * c[76] + 0.47140452079103168293389624140323 * c[80];
	ce[12] = 0.5 * c[49] - 0.5 * c[45] - 0.5 * c[63] + 0.5 * c[67];
	ce[13] = 0.57735026918962576450914878050196 * c[53] - 0.28867513459481288225457439025098 * c[49] - 0.28867513459481288225457439025098 * c[45] - 0.28867513459481288225457439025098 * c[63] - 0.28867513459481288225457439025098 * c[67] + 0.57735026918962576450914878050196 * c[71];
	ce[14] = 0.5 * c[50] + 0.5 * c[52] + 0.5 * c[68] + 0.5 * c[70];
	ce[15] = 0.5 * c[47] + 0.5 * c[51] + 0.5 * c[65] + 0.5 * c[69];
	ce[16] = 0.5 * c[46] + 0.5 * c[48] + 0.5 * c[64] + 0.5 * c[66];
	ce[17] = 0.40824829046386301636621401245098 * c[45] + 0.40824829046386301636621401245098 * c[49] + 0.40824829046386301636621401245098 * c[53] + 0.40824829046386301636621401245098 * c[63] + 0.40824829046386301636621401245098 * c[67] + 0.40824829046386301636621401245098 * c[71];
	ce[18] = 0.5 * c[22] - 0.5 * c[18] - 0.5 * c[54] + 0.5 * c[58];
	ce[19] = 0.57735026918962576450914878050196 * c[26] - 0.28867513459481288225457439025098 * c[22] - 0.28867513459481288225457439025098 * c[18] - 0.28867513459481288225457439025098 * c[54] - 0.28867513459481288225457439025098 * c[58] + 0.57735026918962576450914878050196 * c[62];
	ce[20] = 0.5 * c[23] + 0.5 * c[25] + 0.5 * c[59] + 0.5 * c[61];
	ce[21] = 0.5 * c[20] + 0.5 * c[24] + 0.5 * c[56] + 0.5 * c[60];
	ce[22] = 0.5 * c[19] + 0.5 * c[21] + 0.5 * c[55] + 0.5 * c[57];
	ce[23] = 0.40824829046386301636621401245098 * c[18] + 0.40824829046386301636621401245098 * c[22] + 0.40824829046386301636621401245098 * c[26] + 0.40824829046386301636621401245098 * c[54] + 0.40824829046386301636621401245098 * c[58] + 0.40824829046386301636621401245098 * c[62];
	ce[24] = 0.5 * c[13] - 0.5 * c[9] - 0.5 * c[27] + 0.5 * c[31];
	ce[25] = 0.57735026918962576450914878050196 * c[17] - 0.28867513459481288225457439025098 * c[13] - 0.28867513459481288225457439025098 * c[9] - 0.28867513459481288225457439025098 * c[27] - 0.28867513459481288225457439025098 * c[31] + 0.57735026918962576450914878050196 * c[35];
	ce[26] = 0.5 * c[14] + 0.5 * c[16] + 0.5 * c[32] + 0.5 * c[34];
	ce[27] = 0.5 * c[11] + 0.5 * c[15] + 0.5 * c[29] + 0.5 * c[33];
	ce[28] = 0.5 * c[10] + 0.5 * c[12] + 0.5 * c[28] + 0.5 * c[30];
	ce[29] = 0.40824829046386301636621401245098 * c[9] + 0.40824829046386301636621401245098 * c[13] + 0.40824829046386301636621401245098 * c[17] + 0.40824829046386301636621401245098 * c[27] + 0.40824829046386301636621401245098 * c[31] + 0.40824829046386301636621401245098 * c[35];
	ce[30] = 0.40824829046386301636621401245098 * c[4] - 0.40824829046386301636621401245098 * c[0] - 0.40824829046386301636621401245098 * c[36] + 0.40824829046386301636621401245098 * c[40] - 0.40824829046386301636621401245098 * c[72] + 0.40824829046386301636621401245098 * c[76];
	ce[31] = 0.47140452079103168293389624140323 * c[8] - 0.23570226039551584146694812070162 * c[4] - 0.23570226039551584146694812070162 * c[0] - 0.23570226039551584146694812070162 * c[36] - 0.23570226039551584146694812070162 * c[40] + 0.47140452079103168293389624140323 * c[44] - 0.23570226039551584146694812070162 * c[72] - 0.23570226039551584146694812070162 * c[76] + 0.47140452079103168293389624140323 * c[80];
	ce[32] = 0.40824829046386301636621401245098 * c[5] + 0.40824829046386301636621401245098 * c[7] + 0.40824829046386301636621401245098 * c[41] + 0.40824829046386301636621401245098 * c[43] + 0.40824829046386301636621401245098 * c[77] + 0.40824829046386301636621401245098 * c[79];
	ce[33] = 0.40824829046386301636621401245098 * c[2] + 0.40824829046386301636621401245098 * c[6] + 0.40824829046386301636621401245098 * c[38] + 0.40824829046386301636621401245098 * c[42] + 0.40824829046386301636621401245098 * c[74] + 0.40824829046386301636621401245098 * c[78];
	ce[34] = 0.40824829046386301636621401245098 * c[1] + 0.40824829046386301636621401245098 * c[3] + 0.40824829046386301636621401245098 * c[37] + 0.40824829046386301636621401245098 * c[39] + 0.40824829046386301636621401245098 * c[73] + 0.40824829046386301636621401245098 * c[75];
	ce[35] = 0.33333333333333333333333333333333 * c[0] + 0.33333333333333333333333333333333 * c[4] + 0.33333333333333333333333333333333 * c[8] + 0.33333333333333333333333333333333 * c[36] + 0.33333333333333333333333333333333 * c[40] + 0.33333333333333333333333333333333 * c[44] + 0.33333333333333333333333333333333 * c[72] + 0.33333333333333333333333333333333 * c[76] + 0.33333333333333333333333333333333 * c[80];

	 
}

KOKKOS_FUNCTION void chg_basis(double* ce, double* c, const size_t iopt, const size_t kdim)
{
	// purpose: change of basis
	// calls: nothing else.


	//real (kind=dp) ce2(kdim),c2(3,3),ce4(kdim,kdim),c4(3,3,3,3)

	//      data b /rsq6,0,   0,   0,   rsq6,0,   0,   0,  -2*rsq6,
	//     1        rsq2,0,   0,   0,  -rsq2,0,   0,   0,   0,
	//     1        0,   0,   0,   0,   0,   rsq2,0,   rsq2,0,
	//     1        0,   0,   rsq2,0,   0,   0,   rsq2,0,   0,
	//     1        0,   rsq2,0,   rsq2,0,   0,   0,   0,   0,
	//     1        rsq3,0,   0,   0,   rsq3,0,   0,   0,   rsq3/
	//
	//      common/basis/ b(3,3,6)
	//
	//      if(iopt.eq.0) then
	// *** calculates basis tensors b(n)

	//static double b[54] = { -0.70710678118654746, -0.40824829046386302, 0, 0, 0, 0.57735026918962584, 0, 0, 0, 0, 0.70710678118654746, 0, 0, 0, 0, 0.70710678118654746, 0, 0, 0, 0, 0, 0, 0.70710678118654746, 0, 0.70710678118654746, -0.40824829046386302, 0, 0, 0, 0.57735026918962584, 0, 0, 0.70710678118654746, 0, 0, 0, 0, 0, 0, 0.70710678118654746, 0, 0, 0, 0, 0.70710678118654746, 0, 0, 0, 0, 0.81649658092772603, 0, 0, 0, 0.57735026918962584 };

	//double sqr2 = 1.41421356237309;
	//double rsq2 = 0.70710678118654744;
	//double rsq3 = 0.57735026918962584;
	//double rsq6 = 0.40824829046386304;
	//
	//b = 0.0;
	//
	//b(1, 1, 2) = -rsq6;
	//b(2, 2, 2) = -rsq6;
	//b(3, 3, 2) = 2.0 * rsq6;
	//
	//b(1, 1, 1) = -rsq2;
	//b(2, 2, 1) = rsq2;
	//
	//b(2, 3, 3) = rsq2;
	//b(3, 2, 3) = rsq2;
	//
	//b(1, 3, 4) = rsq2;
	//b(3, 1, 4) = rsq2;
	//
	//b(1, 2, 5) = rsq2;
	//b(2, 1, 5) = rsq2;
	//
	//b(1, 1, 6) = rsq3;
	//b(2, 2, 6) = rsq3;
	//b(3, 3, 6) = rsq3;

	switch (iopt + 4 * (kdim - 5)) {
	case 1:

		// *** calculates cartesian second order tensor from b-components vector.
		chg_basis15(&ce[0], &c[0]);
		break;


	case 2:

		// *** calculates kdimx1 b-components vector from second order tensor.
		chg_basis25(&ce[0], &c[0]);
		break;

	case 3:

		// *** calculates fourth order tensor from b-components matrix.
		chg_basis35(&ce[0], &c[0]);

		break;

	case 4:
		// *** calculates kdimxkdim b-components matrix from fourth order tensor.
		chg_basis45(&ce[0], &c[0]);
		break;

	case 5:

		// *** calculates cartesian second order tensor from b-components vector.
		chg_basis16(&ce[0], &c[0]);
		break;


	case 6:

		// *** calculates kdimx1 b-components vector from second order tensor.
		chg_basis26(&ce[0], &c[0]);
		break;

	case 7:

		// *** calculates fourth order tensor from b-components matrix.
		chg_basis36(&ce[0], &c[0]);

		break;

	case 8:
		// *** calculates kdimxkdim b-components matrix from fourth order tensor.
		chg_basis46(&ce[0], &c[0]);
		break;
	}

	// *** calculates cartesian second order tensor from b-components vector.
	//switch (iopt) {
	//case 1:
	//	for (short i = 0; i < 3; i++) {
	//		for (short j = 0; j < 3; j++) {
	//			c[i * 3 + j] = 0.0;
	//			for (short n = 0; n < kdim; n++) {
	//				c[i * 3 + j] += ce[n] * b[18 * i + 6 * j + n];
	//			}
	//		}
	//	}
	//	break;
	//
	//	// *** calculates kdimx1 b-components vector from second order tensor.
	//case 2:
	//	for (short n = 0; n < kdim; n++) {
	//		ce[n] = 0.0;
	//		for (short i = 0; i < 3; i++) {
	//			for (short j = 0; j < 3; j++) {
	//				ce[n] += c[i * 3 + j] * b[18 * i + 6 * j + n];
	//			}
	//		}
	//	}
	//	break;
	//
	//	// *** calculates fourth order tensor from b-components matrix.
	//case 3:
	//	for (short i = 0; i < 3; i++) {
	//		for (short j = 0; j < 3; j++) {
	//			for (short k = 0; k < 3; k++) {
	//				for (short l = 0; l < 3; l++) {
	//					double dum = 0.0;
	//					for (short n = 0; n < kdim; n++) {
	//						for (short m = 0; m < kdim; m++) {
	//							dum += ce[n * kdim + m] * b[18 * i + 6 * j + n] * b[18 * k + 6 * l + m];
	//						}
	//					}
	//
	//					c[27 * i + 9 * j + 3 * k + l] = dum;
	//				}
	//			}
	//		}
	//	}
	//
	//	break;
	//	// *** calculates kdimxkdim b-components matrix from fourth order tensor.
	//case 4:
	//
	//	for (short n = 0; n < kdim; n++) {
	//		for (short m = 0; m < kdim; m++) {
	//			double dum = 0.0;
	//			for (short i = 0; i < 3; i++) {
	//				for (short j = 0; j < 3; j++) {
	//					for (short k = 0; k < 3; k++) {
	//						for (short l = 0; l < 3; l++) {
	//							dum += c[27 * i + 9 * j + 3 * k + l] * b[18 * i + 6 * j + n] * b[18 * k + 6 * l + m];
	//						}
	//					}
	//				}
	//			}
	//			ce[n * kdim + m] = dum;
	//		}
	//	}
	//	break;
	//}

	 
}

//
// **********************************************************************
//     subroutine inv3_voigt   --->   version 23/jul/01
//
//     inverts the 3x3 symmetric matrix 'a' using explicit voigt notation
//     11->1, 22->2, 33->3, 23=32->4, 31=13->5, 12=21->6
// **********************************************************************
KOKKOS_FUNCTION void inv3_voigt(MATARIX <double>& a, MATARIX <double>& ainv)
{

	//real (kind=dp) a(10),ainv(10)

	double idet;
	idet= 1.0 / (a(1) * a(2) * a(3) + 2 * a(4) * a(5) * a(6) - a(1) * a(4) * a(4)
		- a(2) * a(5) * a(5) - a(3) * a(6) * a(6));

	ainv(1) = (a(2) * a(3) - a(4) * a(4)) * idet;
	ainv(2) = (a(1) * a(3) - a(5) * a(5)) * idet;
	ainv(3) = (a(1) * a(2) - a(6) * a(6)) * idet;
	ainv(4) = (-a(1) * a(4) + a(5) * a(6)) * idet;
	ainv(5) = (a(4) * a(6) - a(2) * a(5)) * idet;
	ainv(6) = (-a(3) * a(6) + a(4) * a(5)) * idet;

	 
}

KOKKOS_FUNCTION void inv4_voigt(MATARIX <double> a, MATARIX <double> ainv)
//
// **********************************************************************
//     subroutine inv4_voigt   --->   version 20/jul/01
//
//     inverts the 4*4 symmetric matrix 'a' using explicit voigt notation
//     11-->1, 22-->2, 33-->3, 23=32-->4, 31=13-->5, 12=21-->6
//     14-->7, 24-->8, 34-->9, 44-->10.
// **********************************************************************
{
	double det;

	ainv(1) = a(2) * a(3) * a(10) + 2 * a(4) * a(8) * a(9) - a(2) * a(9) * a(9) - a(3) * a(8) * a(8) - a(4) * a(4) * a(10);

	ainv(2) = a(1) * a(3) * a(10) + 2 * a(5) * a(7) * a(9) - a(1) * a(9) * a(9) - a(3) * a(7) * a(7) - a(5) * a(5) * a(10);

	ainv(3) = a(1) * a(2) * a(10) + 2 * a(6) * a(7) * a(8) - a(1) * a(8) * a(8) - a(2) * a(7) * a(7) - a(6) * a(6) * a(10);

	ainv(4) = a(1) * a(4) * a(10) + a(5) * a(7) * a(8) + a(6) * a(7) * a(9) -
		a(1) * a(8) * a(9) - a(4) * a(7) * a(7) - a(5) * a(6) * a(10);
	ainv(4) = -ainv(4);

	ainv(5) = a(4) * a(6) * a(10) + a(2) * a(7) * a(9) + a(5) * a(8) * a(8) -
		a(4) * a(7) * a(8) - a(6) * a(8) * a(9) - a(2) * a(5) * a(10);

	ainv(6) = a(3) * a(6) * a(10) + a(5) * a(8) * a(9) + a(4) * a(7) * a(9) -
		a(3) * a(7) * a(8) - a(6) * a(9) * a(9) - a(4) * a(5) * a(10);
	ainv(6) = -ainv(6);

	ainv(7) = a(4) * a(6) * a(9) + a(4) * a(5) * a(8) + a(2) * a(3) * a(7) -
		a(4) * a(4) * a(7) - a(2) * a(5) * a(9) - a(3) * a(6) * a(8);
	ainv(7) = -ainv(7);

	ainv(8) = a(1) * a(4) * a(9) + a(5) * a(5) * a(8) + a(3) * a(6) * a(7) -
		a(4) * a(5) * a(7) - a(5) * a(6) * a(9) - a(1) * a(3) * a(8);

	ainv(9) = a(1) * a(2) * a(9) + a(5) * a(6) * a(8) + a(4) * a(6) * a(7) -
		a(2) * a(5) * a(7) - a(6) * a(6) * a(9) - a(1) * a(4) * a(8);
	ainv(9) = -ainv(9);

	ainv(10) = a(1) * a(2) * a(3) + 2 * a(4) * a(5) * a(6) -
		a(1) * a(4) * a(4) - a(2) * a(5) * a(5) - a(3) * a(6) * a(6);

	det =  (a(1) * ainv(1) + a(2) * ainv(2) + a(3) * ainv(3) +
		2. * a(4) * ainv(4) + 2. * a(5) * ainv(5) + 2. * a(6) * ainv(6) +
		2. * a(7) * ainv(7) + 2. * a(8) * ainv(8) + 2. * a(9) * ainv(9) +
		a(10) * ainv(10));

	//if (det != 0.0) {
	det = 4.0/det;
	for (short i = 1; i <= 10; i++) {
		ainv(i) *= det;
	}
	//}
	//ainv /= det;

	 
}

KOKKOS_FUNCTION void mult_voigt(MATARIX <double> b, MATARIX <double> c, MATARIX <double> a)
{
	//     performs the multiplication:
	//        a(i,k)=b(i,j,k,l)*c(j,l) using voigt's notation
	//        b is a 6*6 symmetric matrix
	//        c is a 3*3 symmetric tensor
	//        a will be a 3*3 symmetric tensor

	//real (kind=dp) b(6,6),c(6),a(9)

  //  a(1) = b(1, 1) * c(1) + b(6, 6) * c(2) + b(5, 5) * c(3)
  //      + 2 * (b(5, 6) * c(4) + b(1, 5) * c(5) + b(1, 6) * c(6));
  //
  //  a(2) = b(6, 6) * c(1) + b(2, 2) * c(2) + b(4, 4) * c(3)
  //      + 2 * (b(2, 4) * c(4) + b(4, 6) * c(5) + b(2, 6) * c(6));
  //
  //  a(3) = b(5, 5) * c(1) + b(4, 4) * c(2) + b(3, 3) * c(3)
  //      + 2 * (b(3, 4) * c(4) + b(3, 5) * c(5) + b(4, 5) * c(6));
  //
  //  a(4) = b(5, 6) * c(1) + b(2, 4) * c(2) + b(3, 4) * c(3)
  //      + (b(2, 3) + b(4, 4)) * c(4)
  //      + (b(3, 6) + b(4, 5)) * c(5)
  //      + (b(4, 6) + b(2, 5)) * c(6);
  //
  //  a(5) = b(1, 5) * c(1) + b(4, 6) * c(2) + b(3, 5) * c(3)
  //      + (b(3, 6) + b(4, 5)) * c(4)
  //      + (b(1, 3) + b(5, 5)) * c(5)
  //      + (b(1, 4) + b(5, 6)) * c(6);
  //
  //  a(6) = b(1, 6) * c(1) + b(2, 6) * c(2) + b(4, 5) * c(3)
  //      + (b(4, 6) + b(2, 5)) * c(4)
  //      + (b(1, 4) + b(5, 6)) * c(5)
  //      + (b(1, 2) + b(6, 6)) * c(6);



	a(1) = b(1) * c(1) + b(36) * c(2) + b(29) * c(3)
		+ 2 * (b(30) * c(4) + b(5) * c(5) + b(6) * c(6));

	a(2) = b(36) * c(1) + b(8) * c(2) + b(22) * c(3)
		+ 2 * (b(10) * c(4) + b(24) * c(5) + b(12) * c(6));

	a(3) = b(29) * c(1) + b(22) * c(2) + b(15) * c(3)
		+ 2 * (b(16) * c(4) + b(17) * c(5) + b(23) * c(6));

	a(4) = b(30) * c(1) + b(10) * c(2) + b(16) * c(3)
		+ (b(9) + b(22)) * c(4)
		+ (b(18) + b(23)) * c(5)
		+ (b(24) + b(11)) * c(6);

	a(5) = b(5) * c(1) + b(24) * c(2) + b(17) * c(3)
		+ (b(18) + b(23)) * c(4)
		+ (b(3) + b(29)) * c(5)
		+ (b(4) + b(30)) * c(6);

	a(6) = b(6) * c(1) + b(12) * c(2) + b(23) * c(3)
		+ (b(24) + b(11)) * c(4)
		+ (b(4) + b(30)) * c(5)
		+ (b(2) + b(36)) * c(6);
	 
}

//----------------------------------------------------------------------
//ROTATE_TENS4
//Rotates a column major rank 4 tensor 
// 
//Nathaniel Morgan - nmorgan@lanl.gov - 6-20-2021
//----------------------------------------------------------------------
KOKKOS_FUNCTION void rotate_tens4(MATARRY <double> ag, MATARRY <double> in, MATARRY <double> out) {

	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < 3; j++) {
			for (size_t m = 0; m < 3; m++) {
				for (size_t n = 0; n < 3; n++) {
					double dum = 0.0, agi, agj, agm;
					for (size_t i1 = 0; i1 < 3; i1++) {
						agi = ag(i, i1);
						for (size_t j1 = 0; j1 < 3; j1++) {
							agj = agi * ag(j, j1);
							for (size_t m1 = 0; m1 < 3; m1++) {
								agm = agj * ag(m, m1);
								for (size_t n1 = 0; n1 < 3; n1++) {
									dum += agm * ag(n, n1) * in(i1, j1, m1, n1); // i1 * 27 + j1 * 9 + m1 * 3 + n1];
								}
							}
						}
					}
					out(i, j, m, n) = dum;
				}
			}
		}
	}

	 
}


template <size_t N> KOKKOS_FUNCTION void matMul(MATARRY <double> a, MATARRY <double> b, MATARRY <double> c)
{
	for (size_t j = 0; j < N; j++) {
		for (size_t i = 0; i < N; i++) {
			// call sumprod written for 2nd index
			double dum = 0.0;

			for (size_t k = 0; k < N; k++) {
				dum += a(j, k) * b(k, i);
			}

			c(j, i) += dum;
		}
	}
	 
}


template <size_t N> KOKKOS_FUNCTION void vecMul(MATARRY <double> a, MATARRY <double> b, MATARRY <double> c)
{
	for (size_t j = 0; j < N; j++) {

		double dum;
		dum = 0.0;

		for (size_t k = 0; k < N; k++) {
			dum += a(j, k) * b(k);
		}

		c(j) += dum;
	}
	 
}


KOKKOS_FUNCTION double powm(double x, size_t y) {
	double x2, x4, xo;
	switch (y) {
	case 0:
		xo = 1.0;
		break;
	case 1:
		break;
	case 2:
		xo = (x * x);
		break;
	case 3:
		xo = (x * x * x);
		break;
	case 4:

		x = x * x;
		xo = (x * x);

		break;

	case 5:

		x2 = x * x;
		xo = (x2 * x2 * x);

		break;
	case 6:

		x = x * x;
		xo = (x * x * x);

		break;
	case 7:

		x2 = x * x;
		xo = (x2 * x2 * x2 * x);

		break;
	case 8:

		x = x * x;
		x = x * x;
		xo = (x * x);

		break;
	case 9:

		x2 = x * x;
		x4 = x2 * x2;
		xo = (x4 * x4 * x);

		break;
	case 10:

		x2 = x * x;
		x4 = x2 * x2;
		xo = (x4 * x4 * x2);

		break;
	case 11:

		x2 = x * x;
		x4 = x2 * x2;
		xo = (x4 * x4 * x2 * x);

		break;
	case 12:

		x2 = x * x;
		x4 = x2 * x2;
		xo = (x4 * x4 * x4);

		break;
	case 13:

		x2 = x * x;
		x4 = x2 * x2;
		xo = (x4 * x4 * x4 * x);

		break;
	case 14:

		x2 = x * x;
		x4 = x2 * x2;
		xo = (x4 * x4 * x4 * x2);

		break;
	case 15:

		x2 = x * x;
		x4 = x2 * x2;
		xo = (x4 * x4 * x4 * x2 * x);

		break;
	case 16:

		x = x * x;
		x = x * x;
		x = x * x;
		xo = (x * x);

		break;

	default:
		xo = pow(x, y);
	}

	return xo;
}


KOKKOS_FUNCTION void eigsort(MATARIX <double>& d,
	MATARIX <double>& v,
	size_t n,
	size_t np)
{
	//real(kind=dp) d(np),v(np,np)
	size_t k;
	double p;
	for (size_t i = 1; i <= n - 1; i++) {
		k = i;
		p = d(i);
		for (size_t j = i + 1; j <= n; j++) {
			if (d(j) >= p) {
				k = j;
				p = d(j);
			}
		}
		if (k != i) {
			d(k) = d(i);
			d(i) = p;
			for (size_t j = 1; j <= n; j++) {
				p = v(j, i);
				v(j, i) = v(j, k);
				v(j, k) = p;
			}
		}
	}
	 
}

//
// **********************************************************************
//
KOKKOS_FUNCTION void eigen(MATARIX <double>& a,
	const size_t n,
	const size_t np,
	MATARIX <double>& d,
	MATARIX <double>& v,
	size_t& nrot,
	size_t& ier) {

	const size_t nmax = n;
	double b1D[10], z1D[10];
	MATARIX <double> b(b1D, nmax), z(z1D, nmax);


	//xmat1(double, b, 10, nmax); xmat1(double, z, 10, nmax);
	double c, g, h, s, sm, t, tau, theta, tresh;

	const size_t len = n * n;
	
	//v      integer i,ip,iq,j
	if (max(a, len) == 0.0) {
		//printf("zeros in eigen\n");
		
		return;
	}


	for (size_t ip = 1; ip <= n; ip++) {
		for (size_t iq = 1; iq <= n; iq++) {
			v(ip, iq) = 0.0;
		}
		v(ip, ip) = 1.0;
	}

	for (size_t ip = 1; ip <= n; ip++) {
		b(ip) = a(ip, ip);
		d(ip) = b(ip);
		z(ip) = 0.0;
	}
	nrot = 0;
	for (size_t i = 1; i <= 50; i++) {
		sm = 0.;
		for (size_t ip = 1; ip <= n - 1; ip++) {
			for (size_t iq = ip + 1; iq <= n; iq++) {
				sm += fabs(a(ip, iq));

			}
		}
		//
		if (sm == 0.) {
			ier = 0;
			return;
		}
		//
		if (i < 4) {
			tresh = 0.2 * sm / n * n;
		}
		else tresh = 0.;

		for (size_t ip = 1; ip <= n - 1; ip++) {
			for (size_t iq = ip + 1; iq <= n; iq++) {
				g = 100. * fabs(a(ip, iq));
				if ((i > 4) && (fabs(d(ip)) + g == fabs(d(ip))) && (fabs(d(iq)) + g == fabs(d(iq)))) {
					a(ip, iq) = 0.;
				}
				else if (fabs(a(ip, iq)) > tresh) {
					h = d(iq) - d(ip);
					if (fabs(h) + g == fabs(h)) {
						t = a(ip, iq) / h;

					}
					else {
						theta = 0.5 * h / a(ip, iq);
						t = 1. / (fabs(theta) + sqrt(1. + theta * theta));
						if (theta < 0.)t = -t;
					}
					c = 1. / sqrt(1 + t * t);
					s = t * c;
					tau = s / (1. + c);
					h = t * a(ip, iq);
					z(ip) -= h;
					z(iq) += h;
					d(ip) -= h;
					d(iq) += h;
					a(ip, iq) = 0.;
					for (size_t j = 1; j <= ip - 1; j++) {
						g = a(j, ip);
						h = a(j, iq);

						a(j, ip) = g - s * (h + g * tau);
						a(j, iq) = h + s * (g - h * tau);
					}
					for (size_t j = ip + 1; j <= iq - 1; j++) {
						g = a(ip, j);
						h = a(j, iq);
						a(ip, j) = g - s * (h + g * tau);
						a(j, iq) = h + s * (g - h * tau);
					}
					for (size_t j = iq + 1; j <= n; j++) {
						g = a(ip, j);
						h = a(iq, j);
						a(ip, j) = g - s * (h + g * tau);
						a(iq, j) = h + s * (g - h * tau);
					}
					for (size_t j = 1; j <= n; j++) {
						g = v(j, ip);
						h = v(j, iq);

						v(j, ip) = g - s * (h + g * tau);
						v(j, iq) = h + s * (g - h * tau);
					}
					nrot = nrot + 1;
				}
			}
		}
		for (size_t ip = 1; ip <= n; ip++) {
			b(ip) = b(ip) + z(ip);
			d(ip) = b(ip);
			z(ip) = 0.;
		}
	}
	//      pause 'too many iterations in eigen'
	//
	ier = 1;
	//
	 
}

KOKKOS_FUNCTION double tmismatch(const size_t sz, //total size of matrices
    MATARIX <double> v1,//first matrix, compared to second matrix
    MATARIX <double> v2)//second matrix
{
    double v_dif, v_ave, out;
v_dif = 0.0;
v_ave = 0.0;

    for (size_t i = 1; i <= sz; i++) {
        v_dif += (v1(i) - v2(i)) * (v1(i) - v2(i));
        v_ave += 0.25 * (v1(i) + v2(i)) * (v1(i) + v2(i));
    }
    //if (v_ave != 0.0) 
    out = (sqrt(v_dif) / sqrt(v_ave));
    //else out = sqrt(v_dif);
    return out;
}


KOKKOS_FUNCTION void fullout(MATARIX <double>& x) {
    const size_t sz = x.size();

    printf("x : ");
    for (size_t i = 1; i <= sz; i++) {
        printf("%26.16E,  ", x(i));
    }

    printf("\n");

}


KOKKOS_FUNCTION void fullout(const size_t sz, double* x) {
    printf("x : ");
    for (size_t i = 0; i < sz; i++) {
        printf("%26.16E,  ", x[i]);
    }

    printf("\n");

}

KOKKOS_FUNCTION size_t random2(size_t seed) {
    return seed;
}

// ********************************************************************
//     subroutine drot2spin
// ********************************************************************
KOKKOS_FUNCTION void drot2spin(MATARIX <double>& drot, double dtime, MATARIX <double>& w_app)
{
    //
    //   calculates spin based on incremental rotation.
    //   input:
    //     drot - incremental rotation (3x3)
    //     dtime - time increment
    //   output
    //     w_app - applied spin
    //
    //real (kind=dp) drot(3,3),dtime,w_app(3,3);
    //MATARIX <double> n_dual(3), n(3, 3);
    fmat1(double, n_dual, 3); fmat2(double, n, 3, 3);

    double n_dual_norm, theta_dot;

    n_dual(1) = drot(3, 2) - drot(2, 3);//dual vector;
    n_dual(2) = -drot(3, 1) + drot(1, 3);
    n_dual(3) = drot(2, 1) - drot(1, 2);
    n_dual_norm = tnorm(n_dual);

    if (n_dual_norm == 0.0) {
        n_dual = 0.0;
    }
    else {
        n_dual *= (1.0 / n_dual_norm);
    }

    // small rotation safe guard
    if (fabs(0.5 * (drot(1, 1) + drot(2, 2) + drot(3, 3) - 1.0)) >= 1.0) {
        theta_dot = 0.0;
    }
    else {
        if (acos(0.5 * (drot(1, 1) + drot(2, 2) + drot(3, 3) - 1.0)) > 1.0e-10) {
            theta_dot = acos(0.5 * (drot(1, 1) + drot(2, 2) + drot(3, 3) - 1.0)) / dtime; //angular velocity;
        }
        else {
            theta_dot = 0.0;
        }
    }

    n = 0.0; //skew tensor;
    n(1, 2) = -n_dual(3);
    n(1, 3) = n_dual(2);
    n(2, 3) = -n_dual(1);

    for (size_t i = 1; i <= 3; i++) {
        for (size_t j = 1; j <= 3; j++) {
            w_app(i, j) = (n(i, j) - n(j, i)) * theta_dot;
        }
    }

    //w_app = n;
    //transpose(n, 3);
    //w_app -= n;
    //w_app *= theta_dot;
    //spin;
    return;
}


KOKKOS_FUNCTION void rodrigues(MATARIX <double>& c, MATARIX <double>& arot) {
    // **********************************************************************
    //     rodrigues --> version 7-2021
    //
    //     builds incremental rotation matrix 'arot' based on rodrigues formu
    //     'c' is the incremental lattice spin.
    //     'arot' transforms from initial to final orientation.
    // **********************************************************************

    fmat2(double, th, 3, 3);
    fmat2(double, th2, 3, 3);
    fmat1(double, v, 3);

    double snorm, snorm1;

    //v(1) = c(3, 2);
    //v(2) = c(1, 3);
    //v(3) = c(2, 1);

    v(1) = c(8);
    v(2) = c(3);
    v(3) = c(4);

    snorm = sqrt(v(1) * v(1) + v(2) * v(2) + v(3) * v(3));

    snorm1 = tan(snorm * 0.5);

    if (snorm < 1.e-06) snorm = 1.0;

    for (size_t i = 1; i <= 3; i++) {
        v(i) = snorm1 * v(i) / snorm;
    }

    snorm = v(1) * v(1) + v(2) * v(2) + v(3) * v(3);


    //th(1, 1) = 0.0;
    //th(1, 2) = -v(3);
    //th(1, 3) = v(2);
    //th(2, 1) = v(3);
    //th(2, 2) = 0.0;
    //th(2, 3) = -v(1);
    //th(3, 1) = -v(2);
    //th(3, 2) = v(1);
    //th(3, 3) = 0.0;


    // Assign theta
    th(1) = 0.0;
    th(2) = -v(3);
    th(3) = v(2);
    th(4) = v(3);
    th(5) = 0.0;
    th(6) = -v(1);
    th(7) = -v(2);
    th(8) = v(1);
    th(9) = 0.0;

    th2 = 0.0;
    matMul<3>(th, th, th2);

    //for (size_t i = 1; i <= 3; i++) {
    //    for (size_t j = 1; j <= 3; j++) {
    //        arot(i, j) = 2. * (th(i, j) + th2(i, j)) / (1.0 + snorm);
    //    }
    //    //arot(i, i) += 1;
    //}

    //no longer snorm
    snorm = (2.0 / (1.0 + snorm));

    for (size_t i = 1; i <= 9; i++) {
        arot(i) = (th(i) + th2(i)) * snorm;
    }

    arot(1) += 1.0;
    arot(5) += 1.0;
    arot(9) += 1.0;
    return;
}

KOKKOS_FUNCTION void voigt(MATARIX <double>& t1, //either (3x3) or (3x3x3x3)
    MATARIX <double>& t2, //either (6)   or (6x6)
    size_t iopt)//what operation to convert
{
    //
    //     voigt convention:
    //
    //cv      data ((ijv(n,m),m=1,2),n=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/
    //
    //  I really messed this one up :) R.M. 2021

    //size_t ijv_1d[12] = { 1, 1, 2, 2, 3, 3, 2, 3, 1, 3, 1, 2 };
    //MATARIX <size_t> ijv(ijv_1d, 6, 2);

    switch (iopt) {
    case 1:
        //for (size_t i = 1; i <= 6; i++) {
        //	const size_t i1 = ijv(i, 1);
        //	const size_t i2 = ijv(i, 2);
        //	t2(i1, i2) = t1(i);
        //	t2(i2, i1) = t1(i);
        //}

        t2(1) = t1(1);
        t2(2) = t1(6);
        t2(3) = t1(5);
        t2(4) = t1(6);
        t2(5) = t1(2);
        t2(6) = t1(4);
        t2(7) = t1(5);
        t2(8) = t1(4);
        t2(9) = t1(3);

        break;
    case 2:
        //for (size_t i = 1; i <= 6; i++) {
        //	const size_t i1 = ijv(i, 1);
        //	const size_t i2 = ijv(i, 2);
        //	t1(i) = t2(i1, i2);
        //}

        t1(1) = t2(1);
        t1(2) = t2(5);
        t1(3) = t2(9);
        t1(4) = t2(6);
        t1(5) = t2(3);
        t1(6) = t2(2);
        break;

    case 3:
        //for (size_t i = 1; i <= 6; i++) {
        //
        //	for (size_t j = 1; j <= 6; j++) {
        //		const size_t i1 = ijv(i, 1);
        //		const size_t i2 = ijv(i, 2);
        //		const size_t j1 = ijv(j, 1);
        //		const size_t j2 = ijv(j, 2);
        //		t2(i1, i2, j1, j2) = t1(i, j);
        //		t2(i2, i1, j1, j2) = t1(i, j);
        //		t2(i1, i2, j2, j1) = t1(i, j);
        //		t2(i2, i1, j2, j1) = t1(i, j);
        //	}
        //}


        t2(1) = t1(1);
        t2(2) = t1(6);
        t2(3) = t1(5);
        t2(4) = t1(6);
        t2(5) = t1(2);
        t2(6) = t1(4);
        t2(7) = t1(5);
        t2(8) = t1(4);
        t2(9) = t1(3);
        t2(10) = t1(31);
        t2(11) = t1(36);
        t2(12) = t1(35);
        t2(13) = t1(36);
        t2(14) = t1(32);
        t2(15) = t1(34);
        t2(16) = t1(35);
        t2(17) = t1(34);
        t2(18) = t1(33);
        t2(19) = t1(25);
        t2(20) = t1(30);
        t2(21) = t1(29);
        t2(22) = t1(30);
        t2(23) = t1(26);
        t2(24) = t1(28);
        t2(25) = t1(29);
        t2(26) = t1(28);
        t2(27) = t1(27);
        t2(28) = t1(31);
        t2(29) = t1(36);
        t2(30) = t1(35);
        t2(31) = t1(36);
        t2(32) = t1(32);
        t2(33) = t1(34);
        t2(34) = t1(35);
        t2(35) = t1(34);
        t2(36) = t1(33);
        t2(37) = t1(7);
        t2(38) = t1(12);
        t2(39) = t1(11);
        t2(40) = t1(12);
        t2(41) = t1(8);
        t2(42) = t1(10);
        t2(43) = t1(11);
        t2(44) = t1(10);
        t2(45) = t1(9);
        t2(46) = t1(19);
        t2(47) = t1(24);
        t2(48) = t1(23);
        t2(49) = t1(24);
        t2(50) = t1(20);
        t2(51) = t1(22);
        t2(52) = t1(23);
        t2(53) = t1(22);
        t2(54) = t1(21);
        t2(55) = t1(25);
        t2(56) = t1(30);
        t2(57) = t1(29);
        t2(58) = t1(30);
        t2(59) = t1(26);
        t2(60) = t1(28);
        t2(61) = t1(29);
        t2(62) = t1(28);
        t2(63) = t1(27);
        t2(64) = t1(19);
        t2(65) = t1(24);
        t2(66) = t1(23);
        t2(67) = t1(24);
        t2(68) = t1(20);
        t2(69) = t1(22);
        t2(70) = t1(23);
        t2(71) = t1(22);
        t2(72) = t1(21);
        t2(73) = t1(13);
        t2(74) = t1(18);
        t2(75) = t1(17);
        t2(76) = t1(18);
        t2(77) = t1(14);
        t2(78) = t1(16);
        t2(79) = t1(17);
        t2(80) = t1(16);
        t2(81) = t1(15);

        break;

    case 4:
        //for (size_t i = 1; i <= 6; i++) {
        //	for (size_t j = 1; j <= 6; j++) {
        //		t1(i, j) = t2(ijv(i, 1), ijv(i, 2), ijv(j, 1), ijv(j, 2));
        //	}
        //}


        t1(1) = t2(1);
        t1(2) = t2(5);
        t1(3) = t2(9);
        t1(4) = t2(6);
        t1(5) = t2(7);
        t1(6) = t2(4);
        t1(7) = t2(37);
        t1(8) = t2(41);
        t1(9) = t2(45);
        t1(10) = t2(42);
        t1(11) = t2(43);
        t1(12) = t2(40);
        t1(13) = t2(73);
        t1(14) = t2(77);
        t1(15) = t2(81);
        t1(16) = t2(78);
        t1(17) = t2(79);
        t1(18) = t2(74);
        t1(19) = t2(64);
        t1(20) = t2(68);
        t1(21) = t2(72);
        t1(22) = t2(71);
        t1(23) = t2(70);
        t1(24) = t2(47);
        t1(25) = t2(19);
        t1(26) = t2(23);
        t1(27) = t2(27);
        t1(28) = t2(26);
        t1(29) = t2(57);
        t1(30) = t2(58);
        t1(31) = t2(10);
        t1(32) = t2(14);
        t1(33) = t2(18);
        t1(34) = t2(17);
        t1(35) = t2(30);
        t1(36) = t2(31);

        break;
    }

    return;
}

//
// **********************************************************************
//
KOKKOS_FUNCTION void euler(size_t iopt, double& ph, double& th, double& tm, MATARIX <double>& a)
{
    double sph, cph, sth, cth, stm, ctm;
    if (iopt == 1) {
        th = acos(a(3, 3));
        if (fabs(a(3, 3)) >= 0.9999) {
            tm = 0.0;
            ph = atan2(a(1, 2), a(1, 1));
        }
        else {
            sth = sin(th);
            tm = atan2(a(1, 3) / sth, a(2, 3) / sth);
            ph = atan2(a(3, 1) / sth, -a(3, 2) / sth);
        }
    }
    else {

        sph = sin(ph);
        cph = cos(ph);
        sth = sin(th);
        cth = cos(th);
        stm = sin(tm);
        ctm = cos(tm);


        //a(1, 1) = ctm * cph - sph * stm * cth;
        //a(2, 1) = -stm * cph - sph * ctm * cth;
        //a(3, 1) = sph * sth;
        //a(1, 2) = ctm * sph + cph * stm * cth;
        //a(2, 2) = -sph * stm + cph * ctm * cth;
        //a(3, 2) = -sth * cph;
        //a(1, 3) = sth * stm;
        //a(2, 3) = ctm * sth;
        //a(3, 3) = cth;


        a(1) = ctm * cph - sph * stm * cth;
        a(2) = ctm * sph + cph * stm * cth;
        a(3) = sth * stm;
        a(4) = -stm * cph - sph * ctm * cth;
        a(5) = -sph * stm + cph * ctm * cth;
        a(6) = ctm * sth;
        a(7) = sph * sth;
        a(8) = -sth * cph;
        a(9) = cth;

    }
    return;
}
//
// **********************************************************************
//

// **********************************************************************
//     subroutine eshelby      --->      verson 07/oct/05
//
//     ioption=0: initialize arrays assoc. with gauss integration points.
//     ioption=1: calculate elastic eshelby tensor for elastic inclusion.
//     ioption=2: calculate incompressible eshelby tensors esim (strain-
//                rate) & escr (spin-rate) associated with the visco-
//                plastic inclusion.
//     ioption=3: calculate incompressible and hydrostatic eshelby tensor
//                pesh (deviatoric pressure), pdil (spherical pressure) &
//                desh (dilatation) for visco-plastic inclusion.
//     ioption=4: calculates d(s)/d(m) (term 1)
//     ioption=5: calculates d(s)/d(m) (term 2)
//
//     options 2-3-4-5 admit a non-zero visco-plastic bulk modulus keff.
//
//     algorithms are based in lebensohn et al, msmse 6 (1998) p.447.
//     uses explicit matrix inversion and explicit voigt notation (when
//     possible) to optimize computing time.
//
//     modified oct/2005 to adapt number of integration points to the sha
//     of the ellipsoid in order to keep eshelby tensor within a certain
//     tolerance (based on analysis done by gwenaelle proust).
//     aspect ratio criterion was adopted for the case when axis(2) is
//     largest and axis(3) is smallest.
//
//     inputs: axis, c4, keff, dldm, props, iopt
//     outputs: esim, escr, desh, pesh, ddsddm, pdil
// **********************************************************************
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
    size_t ioption)
{
    // **********************************************************************
    //     calculation of eshelby tensors
    //     for given stiffness 'c4' and ellipsoid axes 'axis'
    // **********************************************************************
    double abc, ratio1, ratio2, ro3, abcoro3;
    size_t cas;
	cas = 1;
    abc = axis(1) * axis(2) * axis(3);
    ratio1 = axis(2) / axis(3);
    ratio2 = axis(1) / axis(3);
    if (ratio1 < 25) cas = 1;
    if (ratio1 >= 25 && ratio2 <= 5) cas = 2;
    if (ratio1 >= 25 && ratio2 > 5) cas = 3;

    const size_t npoints = props.ngaussph(cas) * props.ngaussth(cas);

    fmat2(double, p, 3, 3);
    fmat2(double, c2, 6, 6); fmat2(double, gamma2, 6, 6); fmat2(double, dldm2, 6, 6);
    fmat4(double, gamma4, 3, 3, 3, 3);



    pdil = 0.0;

    p = 0.0;
    gamma2 = 0.0;

    // Calculate the start of index for the case
    size_t inds = 0;
    for (size_t i = 1; i < cas; i++) {
        inds += props.ngaussph(i) * props.ngaussth(i);
    }


    voigt(c2, c4, 4);

    //fullout(c4);

    if (ioption == 5) voigt(dldm2, dldm, 4);

    for (size_t ny = 1; ny <= npoints; ny++) {

        //   compute voigt components a1(1)-a(6) of tensor a(3,3) defined by eq.b
        //   --->  a(i,j)=l(i,j,k,l)*a(j)*a(l)


        //Used in loop
        fmat1(double, x1, 10); fmat1(double, a1, 10); fmat1(double, aa1, 6);


        //Copied from props
        fmat1(double, alpha, 3);
        double ww;

        slice(3, alpha, &props.alpha(inds + ny, 1));
        ww = props.ww(inds + ny);


        //create aa1 in loop rather than increase memory costs
        aa1(1) = alpha(1) * alpha(1);
        aa1(2) = alpha(2) * alpha(2);
        aa1(3) = alpha(3) * alpha(3);
        aa1(4) = alpha(2) * alpha(3);
        aa1(5) = alpha(1) * alpha(3);
        aa1(6) = alpha(1) * alpha(2);


        mult_voigt(c2, aa1, a1);
        //w
        if (ioption == 1) {

            //   if solving an elastic inclusion invert the system
            //   --> a(3,3) x x(3,3) = c(3,3)
            //   inverts a(3,3) using explicit voigt notation.
            //   uses explicit form of c(3,3) to calculate solution in voigt notation

            inv3_voigt(a1, x1);
            //for (size_t i = 1; i <= 6; i++) {
            //    x1(i) = a1inv(i);
            //}
        }
        else if (ioption >= 2) {

            //   if solving a visco-plastic inclusion defines components a1(7) to a1(
            //   solves the system given by eq.b4 --> a(4,4) x x(4,4) = c(4,4)
            //   inverts a(4,4) using explicit voigt notation.
            //   uses explicit form of c(4,4) to calculate solution in voigt notation
            //   the solution is symmetric. numerical deviation from symmetry is aver

            a1(7) = alpha(1);
            a1(8) = alpha(2);
            a1(9) = alpha(3);
            a1(10) = 0.0;
            if (keff > 0.0) a1(10) = -1.0 / keff;

            inv4_voigt(a1, x1);

            //for (size_t i = 1; i <= 10; i++) {
            //    x1(i) = a1inv(i);
            //}

        }
        //feb
        //fee

        double dum = alpha(1) * alpha(1) * axis(1) * axis(1) +
            alpha(2) * alpha(2) * axis(2) * axis(2) +
            alpha(3) * alpha(3) * axis(3) * axis(3);

        ro3 = dum * sqrt(dum);

        abcoro3 = abc / ro3;

        //   compute the eshelby integral eq.b11 defining:
        //         gamma(m,j,n,i)=t(m,n,i,j)=a(m)*a(j)*g(n,i)
        //   with the property:
        //         gamma(m,j,n,i)=gamma(j,m,n,i)=gamma(m,j,i,n)

        for (size_t i = 1; i <= 6; i++) {
            for (size_t j = 1; j <= 6; j++) {
                gamma2(i, j) += ww * aa1(i) * x1(j) * abcoro3;
            }
        }

        //   compute the pressure related eshelby integral eq.b14
        if (ioption == 3) {
            for (size_t i = 1; i <= 3; i++) {
                for (size_t j = 1; j <= 3; j++) {
                    p(i, j) += ww * alpha(j) * x1(i + 6) * abcoro3;
                }
            }
            pdil += ww * x1(10) * abcoro3;
        }

        // end of loop over double integration
    }

    // ********************************************************************
    //   go back to the 3*3*3*3 notation
    voigt(gamma2, gamma4, 3);

    //   compute symmetric (distortion) eshelby tensor from eq.b9.
    //       esim(n,m,k,l)=0.5*(gamma(m,j,n,i)+gamma(n,j,m,i))*c4(i,j,k,l)
    //   compute anti-symmetric (rotation) eshelby tensor from eq.b9.
    //       escr(n,m,k,l)=0.5*(gamma(m,j,n,i)-gamma(n,j,m,i))*c4(i,j,k,l)


    if (ioption < 4) {
        for (size_t n = 1; n <= 3; n++) {
            for (size_t m = 1; m <= 3; m++) {
                for (size_t k = 1; k <= 3; k++) {
                    for (size_t l = 1; l <= 3; l++) {
                        double dumsim = 0.0;
                        double dumcr = 0.0;

                        for (size_t i = 1; i <= 3; i++) {
                            for (size_t j = 1; j <= 3; j++) {
                                //
                                dumsim += (gamma4(m, j, n, i) + gamma4(n, j, m, i)) * c4(i, j, k, l);
                                dumcr += (gamma4(m, j, n, i) - gamma4(n, j, m, i)) * c4(i, j, k, l);
                                //
                            }
                        }

                        esim(n, m, k, l) = 0.5 * dumsim;
                        escr(n, m, k, l) = 0.5 * dumcr;
                    }
                }
            }
        }
    }
    else {
        for (size_t n = 1; n <= 3; n++) {
            for (size_t m = 1; m <= 3; m++) {
                for (size_t k = 1; k <= 3; k++) {
                    for (size_t l = 1; l <= 3; l++) {
                        double dum = 0.0;
                        for (size_t i = 1; i <= 3; i++) {
                            for (size_t j = 1; j <= 3; j++) {
                                dum += (gamma4(m, j, n, i) + gamma4(n, j, m, i)) * dldm(i, j, k, l);
                            }
                        }

                        dsddm(n, m, k, l) = 0.5 * dum;
                    }
                }
            }
        }
    }

    //   compute pressure & dilatation related eshelby tensors (eq.b13)

    if (ioption == 3) {

        pesh = 0.0;
        for (size_t k = 1; k <= 3; k++) {
            for (size_t l = 1; l <= 3; l++) {
                double dum = 0.0;
                for (size_t i = 1; i <= 3; i++) {
                    for (size_t j = 1; j <= 3; j++) {
                        dum += p(i, j) * c4(i, j, k, l);
                    }
                }
                pesh(k, l) = dum;
            }
        }
        for (size_t i = 1; i <= 3; i++) {
            for (size_t j = 1; j <= 3; j++) {

                desh(i, j) = (p(i, j) + p(j, i)) * 0.5;
            }
        }
    }
    return;
}

KOKKOS_FUNCTION void elsc(size_t ioption,
    StateVType& statv,
    PropType& props,
    IntVType& intv)
{
    fmat2(double, aux33, 3, 3);
    fmat4(double, aux3333, 3, 3, 3, 3); //auxiliary variable of dim (3x3x3x3)
    fmat2(double, cnew, 6, 6);
    fmat2(double, cub, 6, 6); // total stiffness
    fmat2(double, c2, 6, 6);
    fmat4(double, c4sa, 3, 3, 3, 3); // crystal stiffness rotated to sample orientation


    double rer, pdil;
    size_t it, ierrlu;

    cub = 0.0;

    // *** transforms single crystal stiffness to sample axes.
    // *** calculates upper bound (voigt) stiffness 'cub' for the polycrystal
    //kgx = 1;
    for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {

        MATARIX <double> c2ca(&props.c2ca(iph, 1, 1), 6, 6);

        voigt(c2ca, aux3333, 3);

        for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {

            MATARIX <double> ag(&statv.ag(kkk, 1, 1), 3, 3), cgr(&intv.cgr(kkk, 1, 1), 6, 6);

            double wgt = statv.wgt(kkk);

            //rotate_tens4(ag, aux3333, c4sa);
            rotate_tens4(ag, aux3333, c4sa);

            chg_basis46(cgr, c4sa);


            for (size_t i = 1; i <= 36; i++) {
                cub(i) += cgr(i) * wgt;
            }
            // end of do kkk
        }
        // end of do iph
    }

    chg_basis36(cub, c4sa);
    voigt(c2, c4sa, 4);

    if (ioption == 0) {

        intv.csc = cub;

        if (props.interaction == 0) {
            intv.ssc = intv.csc;
            mat_inverse(intv.ssc, 6, ierrlu);

            return;

        }
    }

    else if (ioption == 1) {
        //v
////$      print *, 'ioption=1 in elsc'
        //stop
        //v
        //v        do i=1,6
        //v        do j=1,6
        //v          csc(i,j)=cold(i,j)
        //v        enddo
        //v        enddo
    }

    // *** self-consistent calculation of the elastic stiffness
    //x      print *,
    it = 0;
    rer = 2 * 1.0e-05;
    while (rer > 1.0e-05 && it < props.itmaxext) {

        it++;

        cnew = 0.0;

        //printf("intv.csc\n");  matlabout(intv.csc);

        chg_basis36(intv.csc, c4sa);

        for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {

            fmat4(double, e4sa, 3, 3, 3, 3);
            fmat4(double, e4ga, 3, 3, 3, 3);
            fmat2(double, aux66, 6, 6);
            fmat2(double, ctilde, 6, 6);
            fmat2(double, ac2, 6, 6);
            fmat4(double, c4gab, 3, 3, 3, 3);
            fmat1(double, axb, 3);
            fmat2(double, eigb, 3, 3);

            for (size_t i = 1; i <= 3; i++) {
                axb(i) = statv.axisph(iph + 1, 1, i);
                for (size_t j = 1; j <= 3; j++) {
                    eigb(i, j) = statv.axisph(iph + 1, i + 1, j);

                    // printf("%.2f, ", statv.axisph(iph + 1, i + 1, j));
                }
                //printf("\n");
            }




            // *** rotation of stiffness 'c4' from sample to ellipsoid axes


            transpose(eigb, 3);
            rotate_tens4(eigb, c4sa, c4gab);
            transpose(eigb, 3);

            //v          call eshelby (axb,c4gab,0.,e4ga,aux3333,aux33,aux33,
            //v     #                 pdil,aux3333,aux3333,1)
            double zero = 0.0;



            eshelby(axb, c4gab, zero, e4ga, aux3333, aux33, aux33, pdil, aux3333, aux3333, props, 1);



            // *** rotates the eshelby tensor for the phase back into sample axes.

            rotate_tens4(eigb, e4ga, e4sa);

            chg_basis46(aux66, e4sa);
            mat_inverse(aux66, 6, ierrlu);

            //for (size_t i = 1; i <= 6; i++) {
            //    for (size_t j = 1; j <= 6; j++) {
            //        ctilde(i, j) = 0.0;
            //        for (size_t k = 1; k <= 6; k++) {
            //            ctilde(i, j) += intv.csc(i, k) * (aux66(k, j) - xid6(k, j));
            //        }
            //        ac2(i, j) = intv.csc(i, j) + ctilde(i, j);
            //    }
            //}

            peye(aux66, 6, -1.0);

            ctilde = 0.0;
            matMul<6>(intv.csc, aux66, ctilde);

            for (size_t i = 1; i <= 36; i++) {
                ac2(i) = intv.csc(i) + ctilde(i);
            }


            for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {
                fmat2(double, ac, 6, 6); fmat2(double, ac1, 6, 6); fmat2(double, cgr, 6, 6);

                slice(36, cgr, &intv.cgr(kkk, 1, 1));
                double wgt = statv.wgt(kkk);

                for (size_t i = 1; i <= 36; i++) {
                    ac1(i) = cgr(i) + ctilde(i);
                }

                mat_inverse(ac1, 6, ierrlu);

                ac = 0.0;
                matMul<6>(ac1, ac2, ac);

                aux66 = 0.0;
                matMul<6>(cgr, ac, aux66);


                for (size_t i = 1; i <= 36; i++) {
                    cnew(i) += aux66(i) * wgt;
                }
            }
            // end of loop over all grains in all phases
        }


        rer = tmismatch(36, intv.csc, cnew);

        // *** saves the symmetric tensor as the new guess
        for (size_t i = 1; i <= 6; i++) {
            for (size_t j = i; j <= 6; j++) {
                intv.csc(i, j) = 0.5 * (cnew(i, j) + cnew(j, i));
                intv.csc(j, i) = intv.csc(i, j);
            }
        }

        // end of (do..while)
    }

    //for (size_t i = 1; i <= 6; i++) {
    //    for (size_t j = 1; j <= 6; j++) {
    //        intv.ssc(i, j) = intv.csc(i, j);
    //    }
    //}
    intv.ssc = intv.csc;

    mat_inverse(intv.ssc, 6, ierrlu);

    return;
}

// **********************************************************************
//
//     subroutine newton_raphson (ex snlnr)   --->   version apr/2003
//
//     given an input grain stress 'x' solves the viscoplastic equation
//     using newton-raphson linearization and iterating until convergence
//     - fi(x)   : are the functions to minimize.
//     - fgradij : are the derivatives of fi with respect to xj.
//     corrections to 'x' and relaxation of tolerance eliminated (feb/200
//     passing all relevant crystal arrays through common (13/dec/02).
//     added relaxed constraints option (irc=1) (apr/03)
//
// **********************************************************************
KOKKOS_FUNCTION void newton_raphson(size_t irc,
    size_t kmax,
    double eps,
    double taulim,
    const size_t nsystx, // number of slip systems
    MATARIX <double>& x,
    MATARIX <double>& db,
    MATARIX <double>& xmastx,
    MATARIX <double>& scx,
    MATARIX <double>& itaux,
    MATARIX <double>& gamd0x,
    MATARIX <double>& nrsx,
    MATARIX <double>& isensex,
    size_t& ierror)
{
    //v
    xmat1(double, rss, NSYS_, nsystx);
    xmat1(double, sfactor, NSYS_, nsystx);
    xmat1(double, gd, NSYS_, nsystx);
    fmat1(double, f, 5);
    fmat2(double, fgrad, 5, 5);
    fmat1(double, xold, 5);

    double coef, rerror;
    size_t isingular;
    //v
    coef = 0.2;
    ierror = 0;
    bool breakflag;

    for (size_t k = 1; k <= kmax; k++) {

        // *** resolved shears calculation outside y.s. error managing block.
        // *** nrs may be even or odd:
        //     gd always > 0 --> gd is used to get derivatives to build the
        //     coefficient matrix for n-r method in f(i) calculation (independent
        //     term for n-r method)   rss*gd=gamdot -> sign(rss*gd)=sign(rss)

        breakflag = false;

        for (size_t is = 1; is <= nsystx; is++) {

            rss(is) = sumProd(5, nsystx, x, &scx(1, is)) * itaux(is);

            if (!(rss(is) > 0.0 || isensex(is) == 1.0) || (rss(is) < 1.0e-10 && rss(is) > -1.0e-10)) rss(is) = 0.0;

            if (rss(is) > taulim || rss(is) < -taulim) {
                breakflag = true;

                break;
            }


        }




        if (breakflag) {

            for (size_t i = 1; i <= 5; i++) { x(i) = xold(i) + coef * (x(i) - xold(i)); }
            breakflag = false;
            continue;
        }


        for (size_t is = 1; is <= nsystx; is++) {
            gd(is) = gamd0x(is) * powm(fabs(rss(is)), (size_t)(nrsx(is) - 1.0));
            sfactor(is) = gd(is) * itaux(is) * nrsx(is);
        }


        for (size_t i = 1; i <= 5; i++) {
            //// calc f
            f(i) = -db(i);
            f(i) += sumProd(5, x, &xmastx(i, 1));
            f(i) += sumProd(nsystx, &scx(i, 1), rss, gd);
        }

        for (size_t i = 1; i <= 5; i++) {
            //calc fgrad
            //enforce symmetry for stability
            double dum = 0.0;
            for (size_t is = 1; is <= nsystx; is++) {
                dum += sfactor(is) * scx(i, is) * scx(i, is);
            }
            fgrad(i, i) = -dum;// - sumProd(nsystx, sfactor, &scx(i, 1), &scx(i, 1));

            for (size_t j = i + 1; j <= 5; j++) {
                fgrad(i, j) = -sumProd(nsystx, sfactor, &scx(i, 1), &scx(j, 1));
                fgrad(j, i) = fgrad(i, j);
            }
        }

        for (size_t i = 1; i <= 25; i++) {
            fgrad(i) -= xmastx(i);
        }

        for (size_t i = 1; i <= 5; i++) xold(i) = x(i);


        //fullout(x);

        // *** solves linear system
        if (irc == 0) {
            fmat1(double, fo, 5);
            fo = f;

            mat_eqsystem(fgrad, f, 5, isingular);




            for (size_t i = 1; i <= 5; i++) {
                if (f(i) != f(i)) isingular = 1;
                x(i) = xold(i) + f(i);
            }

            if (isingular == 1) {

                printf("error in mat_eqsystem  (isingular=1) newton_raphson\n");
                printf("fgrad * f = x\n");
                for (size_t i = 1; i <= 5; i++) {
                    printf("[%12.4e, %12.4e, %12.4e, %12.4e, %12.4e] [%12.4e] = [%12.4e]\n", fgrad(i, 1), fgrad(i, 2), fgrad(i, 3), fgrad(i, 4), fgrad(i, 5), f(i), fo(i));
                }

                ierror = 1;
                return;
            }

            // *** bounds the stress correction to aKOKKOS_FUNCTION void large oscilations in converg

            //feb
            //fee




        }



        // relaxed constraints case
        else if (irc == 1) {
            fmat2(double, fgradx, 3, 3); fmat1(double, fx, 3);

            fx(1) = f(1);
            fx(2) = f(2);
            fx(3) = f(5);
            fgradx(1, 1) = fgrad(1, 1);
            fgradx(1, 2) = fgrad(1, 2);
            fgradx(1, 3) = fgrad(1, 5);
            fgradx(2, 1) = fgrad(2, 1);
            fgradx(2, 2) = fgrad(2, 2);
            fgradx(2, 3) = fgrad(2, 5);
            fgradx(3, 1) = fgrad(5, 1);
            fgradx(3, 2) = fgrad(5, 2);
            fgradx(3, 3) = fgrad(5, 5);
            mat_eqsystem(fgradx, fx, 3, isingular);
            if (isingular == 1) {
                ierror = 1;
                printf("error in mat_eqsystem  (isingular=1) newton_raphson\n");
                return;
            }
            x(1) = xold(1) + fx(1);
            x(2) = xold(2) + fx(2);
            x(3) = 0.;
            x(4) = 0.;
            x(5) = xold(5) + fx(3);
        }


        rerror = tmismatch(5, x, xold);
        if (rerror < eps) { return; } // keep this for now //printf("%i\n", k);  
        //if(rerror < eps && fnorm < eps) return

        // end of master do
    }
    // *******************************************************************
    ierror = 2;

    return;
}


//
// **********************************************************************
//     subroutine scale_3   ---->   version 26/may/2000
//
//     this subroutine is meant to be used in combination with mts model
//     where rate effects are accounted for in the functional form of crs
//
//     by redefining the factor 'gamd0' of the order of the strain-rate
//     the ratio (rss/crss) will be of order one and there are no rate
//     sensitivity effects associated with the n'th power.
//     the n'th power is only a way to aKOKKOS_FUNCTION void ambiguities and is unique.
// **********************************************************************
KOKKOS_FUNCTION void  scale_3(size_t istep, StateVType& statv, PropType& props, IntVType& intv)
{

    double refrate;

    refrate = tnorm(intv.dbar);

    for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {
        for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {
            statv.gamd0g(kkk) = refrate;
        }
    }

    return;
}

// **********************************************************************
//     subroutine twin_orientation    --->    version of 04/dec/02
//
//     given the burgers vector 'bur' of the twin system (in crystal axes
//     makes a rotation of 180 deg around it to define the orientation of
//     the twin related crystal.
//     the matrix 'atwin' transforms from twin to crystal axes.
// **********************************************************************
KOKKOS_FUNCTION void twin_orientation(MATARIX <double>& bur, MATARIX <double>& atwin)
{
    fmat2(double, hpi, 3, 3);
    fmat2(double, aux, 3, 3);
    double ang1, ang2, zero;
    //size_t i, j, k1, k2;

    //v      data hpi/-1.,0.,0.,0.,-1.,0.,0.,0.,1./
    //v
    hpi = 0.0;
    hpi(1, 1) = -1.0;
    hpi(2, 2) = -1.0;
    hpi(3, 3) = 1.0;

    ang1 = atan2(bur(2), bur(1)) + (PI / 2.0);
    ang2 = sqrt(bur(1) * bur(1) + bur(2) * bur(2));
    ang2 = atan2(ang2, bur(3));
    zero = 0.0;
    euler(2, ang1, ang2, zero, aux);

    for (size_t i = 1; i <= 3; i++) {
        for (size_t j = 1; j <= 3; j++) {
            double dum = 0.0;
            for (size_t k1 = 1; k1 <= 3; k1++) {
                for (size_t k2 = 1; k2 <= 3; k2++) {
                    dum += aux(k1, i) * hpi(k1, k2) * aux(k2, j);
                }
            }

            atwin(i, j) = dum;
        }
    }

    return;
}

// **********************************************************************
//     subroutine update_crss_voce     --->      version of 28/dec/00
//
//     a voce law function of the accumulated shear in each grain is adde
//     as a multiplicative factor that modulates the original linear hard
//     the units (and the strength) of the hardening are carried by the
//     multiplicative factor 'voce'.
//     the self latent coupling coefficients 'hard' are dimensionless
//     constants relative to the factor 'voce'.
//     the increment of crss is given by an analytic integral expression
//     the voce modulating function, rather than using forward extrapolat
//     as was done before 28/set/00.
//     the parameters in voce expression may adopt 'non-kosher' values (d
//***********************************************************************
KOKKOS_FUNCTION void update_crss_voce(size_t ioption, StateVType& statv, PropType& props, IntVType& intv)
{
    ////$  use vpsc_type_def

    //v      include 'vpsc7std.dim'

    double gamtotx, deltgam, dtau, tau0, tau1, thet0, thet1, tiny, voce, fact, expini, expdel;
    size_t iphel, kgx, cnt = 0;

    if (ioption == 1) {
        //z
        //v       crss=0.
        statv.crss = 0.;
        //z
        for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {
            for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {
                for (size_t is = 1; is <= props.nsyst(iph); is++) {
                    //v              crss(is,kkk)=tau(is,0,iphel)
                    statv.crss(kkk, is) = props.tau(iph, is, 1);
                }
            }
        }
    }
    else if (ioption == 2) {

        //for(is = 1; is <= props.nsyst(iphel); is++){
        //                  for(size_t js = 1; js <= props.nsyst(iphel); js++){
        //                  printf("%4.1F, ",props.hard(1, is, js));
        //                  }
        //                  printf("\n");
        //                  }

        kgx = 1;
        for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {
            //         nslsys=nsyst(iphel)-ntwsys(iphel)
            for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {

                gamtotx = statv.gtotgr(kkk);
                deltgam = 0.0;
                for (size_t is = 1; is <= props.nsyst(iph); is++) {
                    deltgam += fabs(statv.gamdot(kkk, is)) * intv.tincr;
                }



                if (deltgam > 0.0) {

                    for (size_t is = 1; is <= props.nsyst(iph); is++) {
                        dtau = 0.0;
                        for (size_t js = 1; js <= props.nsyst(iph); js++) {
                            cnt++;
                            dtau += props.hard(iph, is, js) * fabs(statv.gamdot(kkk, js)) * intv.tincr;
                            //printf("gamd : %u %16.8E\n",  cnt, statv.gamdot(kgx, js));
                        }
                        //printf("dtau : %16.8E\n", dtau);


                        tau0 = props.tau(iphel, is, 1);
                        tau1 = props.tau(iphel, is, 2);
                        thet0 = props.thet(iphel, is, 1);
                        thet1 = props.thet(iphel, is, 2);
                        tiny = tau0 * (0.0001);

                        voce = 0.0;
                        if (fabs(thet0) > tiny) {
                            voce = thet1 * deltgam;
                            if (fabs(tau1) > tiny) {
                                fact = fabs(thet0 / tau1);
                                expini = exp(-gamtotx * fact);
                                expdel = exp(-deltgam * fact);
                                voce = voce - (fact * tau1 - thet1) / fact * expini * (expdel - 1.0) - thet1 / fact * expini * (expdel * ((gamtotx + deltgam) * fact + 1.0) - (gamtotx * fact + 1.0));
                            }
                            //printf("%d  %2d  dcrss = %26.16E\n",kkk, is,  dtau*voce/deltgam);

                        }

                        statv.crss(kkk, is) += dtau * voce / deltgam;

                    }

                }
                //feb
                //fee
                statv.gtotgr(kkk) = gamtotx + deltgam;
                //printf("gtotgr %u : %12.6E\n", kkk, statv.gtotgr(kkk));
                kgx++;
            }
        }

    }



    return;
}

// **********************************************************************
//     subroutine update_crss_temp     --->      version of 10/1/2020
//
//     a temperature and strain-rate dependent law for updating crss
//     derived from temperature and strain-rate sensitive initial slip
//     resistance law from beyerlein 2008
//***********************************************************************
KOKKOS_FUNCTION void update_crss_temp(StateVType& statv, PropType& props, IntVType& intv)
{
    for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {

        for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {
            for (size_t is = 1; is <= props.nsyst(iph); is++) {

                //update crss directly without keeping a copy of crss_ini
                statv.crss(kkk, is) += props.tau0_mode(iph, is, 1) *
                    (exp(-(statv.temp - props.tau0_mode(iph, is, 3)) / props.tau0_mode(iph, is, 2)) -
                        exp(-(statv.tempold - props.tau0_mode(iph, is, 3)) / props.tau0_mode(iph, is, 2)));
            }
        }
    }



    return;
}

//
// **********************************************************************
//     subroutine update_fij      --->      version of nov/28/2005
//
//     uses the velocity gradient (average, phase or grain) in the step
//     to update incrementally the corresponding deformation tensor 'fij'
//
//     replaces previous subr. updfij_average subr. updfij_local  (sep/
// **********************************************************************
KOKKOS_FUNCTION void update_fij(size_t iph, StateVType& statv, IntVType& intv)
{
    fmat2(double, fnew, 3, 3);

    // *** updates the deform grad in the element 'fijph(i,j,0)' using the
    //     macroscopic velocity gradient 'udot'
    // *** updates the deform grad in each phase 'fijph(i,j,iph)' using the
    //     average velocity gradient for the phase 'xlijph' calculated inside
    //     subroutine update_orientation.
    // *** xlijph accounts for rotations but not for stretch when iflat(iph + 1)=
    // *** fijph coincides with fij of element if nph=1 and iflat(1 + 1)=0.

    for (size_t i = 1; i <= 3; i++) {
        for (size_t j = 1; j <= 3; j++) {
            fnew(i, j) = 0.0;
            if (iph == 0) {
                for (size_t k = 1; k <= 3; k++) {
                    fnew(i, j) += (intv.tincr * intv.udot(i, k) + xid(i, k)) * statv.fijph(1, k, j);
                }
            }
            else if (iph > 0) {
                for (size_t k = 1; k <= 3; k++) {
                    fnew(i, j) += (intv.tincr * intv.xlijph(iph, i, k) + xid(i, k)) * statv.fijph(iph + 1, k, j);
                }
            }
        }
    }
    for (size_t i = 1; i <= 3; i++) {
        for (size_t j = 1; j <= 3; j++) {
            statv.fijph(iph + 1, i, j) = fnew(i, j);
        }
    }
    //feb
    //fee
    //      print *,'in vpsc_update_fij, fijph0'
    //      print *, fijph(1,1,0),fijph(1,2,0),fijph(1,3,0)
    //      print *, fijph(2,1,0),fijph(2,2,0),fijph(2,3,0)
    //      print *, fijph(3,1,0),fijph(3,2,0),fijph(3,3,0)
    return;
}

//
// **********************************************************************
//     subroutine update_shape      --->      version 28/nov/2005
//
//     calculates the direction and length (eigenvectors and eigenvalues
//     deformation tensor fij) of the axes of the ellipsoid associated wi
//     average, phase and/or grain accumulated deformation.
//     calculates the euler angles of ellipsoid axes wrt sample axes.
//
//     replaces previous subr. graxes_average subr. graxes_local  (sep/
// *********************************************************************
KOKKOS_FUNCTION void update_shape(size_t iph, StateVType& statv, size_t interaction)
{
    // calls: eigen, eigensort, euler
    ////$  use vpsc_type_def

    //v      include 'vpsc7std.dim'
    fmat1(double, w, 3); fmat2(double, bx, 3, 3); fmat2(double, b, 3, 3); fmat2(double, bt, 3, 3);
    double sign, exchange, ang1, ang2, ang3;
    size_t nrot, ier;

    // *** if iph=0 ellipsoid represents average deformation in element
    // *** if iph>0 ellipsoid represents average deformation in phase 'iph'
    // *** calculates eigenvalues, eigenvectors euler angles of element gra
    //     phase grain, or individual grain
    // *** 'axisph' transforms from ellipsoid to sample axes.
	
	bt = 0.0;

    for (size_t i = 1; i <= 3; i++) {
        for (size_t j = 1; j <= 3; j++) {
            bx(i, j) = 0.0;
            for (size_t k = 1; k <= 3; k++) {
                //wx
                //wx          bx(i,j)=bx(i,j)+fijph(i,k,iph)*fijph(j,k,iph)
                //wx
                //wx     uses udot forall phases to prevent
                //wx     numerical instability in so procedure
                //wx
                if (interaction == 5) bx(i, j) += statv.fijph(1, i, k) * statv.fijph(1, j, k);


                else bx(i, j) += statv.fijph(iph + 1, i, k) * statv.fijph(iph + 1, j, k);

                //wx
            }
        }
    }

    eigen(bx, 3, 3, w, b, nrot, ier);

    eigsort(w, b, 3, 3);

    //printf("w : %16.8f, %16.8f, %16.8f\n", w(1), w(2), w(3));

    if (ier == 1) {
        printf("error in update_shape for phase ellipsoid %zu\n", iph);
    }

    // *** eigenvalues (and assoc eigenvectors) are ordered from larger to sm
    // *** redefine axis(2) to be the largest in order to improve accuracy in
    //     calculation of the eshelby tensor.
    // *** if det(b)<0 means that the system is left handed. it is made right
    //     handed by exchanging 1 and 2.

    sign = -1.0;
    if (det(b) <= 0.0) sign = 1.0;
    for (size_t i = 1; i <= 3; i++) {
        exchange = b(i, 1);
        b(i, 1) = b(i, 2);
        b(i, 2) = exchange * sign;
    }
    exchange = w(1);
    w(1) = w(2);
    w(2) = exchange;

    for (size_t i = 1; i <= 3; i++) {
        //v        axisph(0,i,iph)=sqrt(w(i))
        statv.axisph(iph + 1, 1, i) = sqrt(w(i));
        for (size_t j = 1; j <= 3; j++) {
            //v          axisph(i,j,iph)=b(i,j)
            statv.axisph(iph + 1, i + 1, j) = b(i, j);
            bt(i, j) = b(j, i);
        }
    }

    euler(1, ang1, ang2, ang3, bt);
    //v      eulerph(1,iph)=ang1
    //v      eulerph(2,iph)=ang2
    //v      eulerph(3,iph)=ang3
    statv.eulerph(iph + 1, 1) = ang1;
    statv.eulerph(iph + 1, 2) = ang2;
    statv.eulerph(iph + 1, 3) = ang3;

    // *** stops updating the ellipsoid for a limit aspect ratio

    if (sqrt(w(2) / w(3)) > 50 && statv.iflat(iph + 1) == 0.0 && iph != 0) {

        //v        iflat(iph + 1)=1
        statv.iflat(iph + 1) = 1;
    }
    //feb
    //fee
    return;
}

//
// **********************************************************************
//     subroutine update_orientation    --->   version 28/mar/2007
//
//     updates grain orientations due to crystallographic shear, but does
//     not perform twin reorientation.
//     updates grain and phase distortion tensors required for updating t
//     phase and grain ellipsoid in subroutines update_fij update_shape
//
//     cnt: modified to use only phase-associated tensors when ishape <=
//     ral: split in 3 do loops to deal with co-rotations (17/02/00)
//     cnt: moved last loop to update_twinning (05/06/02)
// **********************************************************************
KOKKOS_FUNCTION void update_orientation(StateVType& statv, PropType& props, IntVType& intv)
{
    if (props.interaction == 0) intv.sbarold = 0.0;

    double rslbar = 0.0, rlcbar = 0.0;

    //kgx = 1;
    for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {
        //iphel = iph - props.iphbot + 1;

        size_t nsystx = props.nsyst(iph);
        fmat4(double, as, 3, 3, 3, 3);

        for (size_t i = 1; i <= 3; i++) {
            for (size_t j = 1; j <= 3; j++) {
                intv.xlijph(iph, i, j) = 0.;
            }
        }

        if (props.ishape(iph + 1) <= 1.0) {
            for (size_t i = 1; i <= 3; i++) {
                for (size_t j = 1; j <= 3; j++) {
                    for (size_t k = 1; k <= 3; k++) {
                        for (size_t l = 1; l <= 3; l++) {
                            as(i, j, k, l) = intv.asph(iph, i, j, k, l);
                        }
                    }
                }
            }

            //For clarity, this is removed
            //slice(81, as, &intv.asph(iph, 1, 1, 1, 1));
        }

        // *** the following to aKOKKOS_FUNCTION void empty 'child' phase in comp grain model
        if (statv.wph(iph) < 1.e-6) continue;

        for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {

            // *** calculates local rotation rotloc=pi*s**(-1)*(dg-dav) for every gra
            // *** rotloc is zero for taylor calculation.
            fmat2(double, rot, 3, 3);
            fmat2(double, aa, 3, 3); fmat2(double, arot, 3, 3);
            fmat2(double, rotslip, 3, 3); fmat2(double, rotloc, 3, 3);
            fmat2(double, xlijgrx, 3, 3); fmat2(double, xlijgr0, 3, 3);

            double wgt = statv.wgt(kkk);


            if (props.interaction <= 0) {
                rotloc = 0.0;
            }
            else if (props.interaction > 0) {
                fmat2(double, east, 3, 3);
                fmat1(double, aux5, 5);

                for (size_t i = 1; i <= 5; i++) {
                    aux5(i) = intv.dg(kkk, i) - intv.dav(i);
                }

                chg_basis15(aux5, east);

                for (size_t i = 1; i <= 3; i++) {
                    for (size_t j = 1; j <= 3; j++) {
                        //rotloc(i, j) = 0.;
                        //for (size_t k = 1; k <= 3; k++) {
                        //	for (size_t l = 1; l <= 3; l++) {
                        //		rotloc(i, j) += as(i, j, k, l) * east(k, l);
                        //	}
                        //}

                        rotloc(i, j) = sumProd(9, &as(i, j, 1, 1), east);
                    }
                }

            }

            // *** calculates velocity gradient in each phase and each grain
            for (size_t i = 1; i <= 3; i++) {
                for (size_t j = 1; j <= 3; j++) {
                    aa(i, j) = statv.ag(kkk, i, j);
                    xlijgr0(i, j) = 0.;
                }
            }

            for (size_t is = 1; is <= nsystx; is++) {
                fmat1(double, dnsa, 3); fmat1(double, dbsa, 3);

                double gamdot = statv.gamdot(kkk, is);

                for (size_t i = 1; i <= 3; i++) {
                    dnsa(i) = 0.0;
                    dbsa(i) = 0.0;

                    for (size_t j = 1; j <= 3; j++) {
                        dnsa(i) += aa(i, j) * props.dnca(iph, is, j);
                        dbsa(i) += aa(i, j) * props.dbca(iph, is, j);
                    }
                }

                for (size_t i = 1; i <= 3; i++) {
                    for (size_t j = 1; j <= 3; j++) {
                        xlijgr0(i, j) += dbsa(i) * dnsa(j) * gamdot;
                    }
                }
            }

            for (size_t i = 1; i <= 3; i++) {
                for (size_t j = 1; j <= 3; j++) {
                    xlijgrx(i, j) = intv.rotbar(i, j) + rotloc(i, j);
                    if (statv.iflat(iph + 1) == 0.0) {
                        xlijgrx(i, j) += (xlijgr0(i, j) + xlijgr0(j, i)) / 2.0;
                    }
                    intv.xlijph(iph, i, j) += xlijgrx(i, j) * wgt / statv.wph(iph);
                }
            }
            //feb
            //fee
            // *** crystallographic grain rotation (rigid minus plastic)
            for (size_t i = 1; i <= 3; i++) {
                for (size_t j = 1; j <= 3; j++) {
                    rotslip(i, j) = (xlijgr0(i, j) - xlijgr0(j, i)) / 2.;
                }
            }
            for (size_t i = 1; i <= 3; i++) {
                for (size_t j = 1; j <= 3; j++) {
                    rot(i, j) = (intv.rotbar(i, j) + rotloc(i, j) - rotslip(i, j)) * intv.tincr;
                }
            }

            // *** average slip and local rotation (for statistical purposes only)
            rslbar = rslbar + sqrt(rotslip(3, 2) * rotslip(3, 2) + rotslip(1, 3) * rotslip(1, 3) + rotslip(2, 1) * rotslip(2, 1)) * wgt;
            rlcbar = rlcbar + sqrt(rotloc(3, 2) * rotloc(3, 2) + rotloc(1, 3) * rotloc(1, 3) + rotloc(2, 1) * rotloc(2, 1)) * wgt;


            MATARIX <double> ag(&statv.ag(kkk, 1, 1), 3, 3);

            aa = ag;
            // *** calculate the new trasformation matrix and update

            rodrigues(rot, arot);

            ag = 0.0;
            matMul<3>(arot, aa, ag);

            if (props.interaction == 0) {
                fmat2(double, aux33, 3, 3);
                fmat2(double, sg33, 3, 3);
                fmat1(double, aux5, 5);
                // rotate deviatoric part of stress
                //for (size_t i = 1; i <= 5; i++) {
                //    aux5(i) = statv.sg(kkk, i);
                //}

                slice(5, aux5, &statv.sg(kkk, 1));

                chg_basis15(aux5, sg33);

                for (size_t i = 1; i <= 3; i++) {
                    for (size_t j = 1; j <= 3; j++) {
                        double dum = 0.0;
                        for (size_t k = 1; k <= 3; k++) {
                            for (size_t l = 1; l <= 3; l++) {
                                dum += arot(i, k) * arot(j, l) * sg33(k, l);
                            }
                        }
                        aux33(i, j) = dum;
                    }
                }

                chg_basis25(aux5, aux33);

                //global sum

                for (size_t i = 1; i <= 5; i++) {
                    statv.sg(kkk, i) = aux5(i);
                    intv.sbarold(i) += aux5(i) * wgt;
                }

            }
            //kgx++;
            // end of do loop over grains
        }
        // end of do loop over phases
    }
    return;
}

//
// **********************************************************************
//     subroutine update_schmid
//
//     rotates schmid tensors of each grain from crystal to sample axes
// **********************************************************************
KOKKOS_FUNCTION void update_schmid(size_t ngr, MATARIX <double> nsyst, MATARIX <double> schca, MATARIX <double> ag, MATARIX <double> sch, MATARIX <double> giph)
{
    for (size_t kkk = 1; kkk <= ngr; kkk++) {

        fmat2(double, aa, 3, 3);
        slice(9, aa, &ag(kkk, 1, 1));

        size_t iph = giph(kkk);
        size_t nsystx = nsyst(iph);

        for (size_t is = 1; is <= nsystx; is++) {
            fmat1(double, aux5, 5); fmat2(double, aux33, 3, 3); fmat2(double, aux33r, 3, 3);

            slice(5, aux5, &schca(iph, is, 1));

            chg_basis15(aux5, aux33);

            for (size_t i = 1; i <= 3; i++) {
                for (size_t j = 1; j <= 3; j++) {
                    double dum = 0.0;
                    for (size_t i1 = 1; i1 <= 3; i1++) {
                        for (size_t j1 = 1; j1 <= 3; j1++) {
                            dum += aa(i, i1) * aa(j, j1) * aux33(i1, j1);
                        }
                    }

                    aux33r(i, j) = dum;
                }
            }

            chg_basis25(aux5, aux33r);

            for (size_t j = 1; j <= 5; j++) {
                sch(kkk, j, is) = aux5(j);
            }
        }
    }

    return;
}

//
// **********************************************************************
//    subroutine update_twinning      ---> version of 21/set/2005
//
// --> cnt modified on sept 2005: now the stats on twinning are done insi
//     1st do loop. the tagging for reorientation is done inside 2nd loop
// **********************************************************************
//
//    iptsgr(kgx):    index of predominant twin system in the grain.
//    iptmgr(kgx):    index of predominant twin mode in the grain.
//    twfrgr(twm):    accumulated twin fract in each mode in grain kkk
//                    in the step (relative to grain).
//
//    twfrsy(tws,kkk):accummulated twin volume fraction in system tws
//                    in grain kkk (relative to grain).
//    twfrph(twm,iph):accumulated twin fract in each mode of each phase
//                    (relative to phase volume).
//    eftwfr(twm,iph):accum twin fraction given by the reoriented grains
//                    in each mode of each phase (relative to phase volum
//    ktwsmx(kkk):index of the predominant twin system (max twin volume)
//                in the grain. used for reorienting grain.
//                   =0 if the grain is not tagged for reorientation.
//                   >0 if the grain is tagged for reorientation.
//                   <0 if the grain has been reoriented and is not to be
//                      retagged for secondary twinning
//
//    pritw(iph) :accummulated twin fract in each phase associated with
//                grains that were twin reoriented once.
//    sectw(iph) :accummulated twin fract in each phase associated with
//                grains that were twin reoriented more than once.
//    ntwevents(kkk) :number of times that grain 'kkk' has been reoriente
//                    0: none, 1: primary, 2:secondary, 3:tertiary
// **********************************************************************
KOKKOS_FUNCTION void update_twinning(size_t iph,
    StateVType& statv,
    PropType& props,
    IntVType& intv)
{
    fmat1(double, btwin, 3);
    fmat2(double, atwin, 3, 3);
    fmat2(double, acrys, 3, 3);

    //v
    double twfmax, gmode, gabs, rand, thres1, thres2, twflim, twshx;
    size_t kkk, its, kts, itm, kk, nsmx, ngphase, kleft, i, j, k, index, kpoint, kacum, igr, ipts, iptm, ktw, is;
    //v
    if (statv.wph(iph) == 0.0) return;

    for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {
        statv.iptsgr(kkk) = 0;
        statv.iptmgr(kkk) = 0;
        twfmax = 0.0;
        its = props.nsyst(iph) - props.ntwsys(iph);
        kts = 0;
        for (size_t itm = 1; itm <= props.ntwmod(iph); itm++) {
            kk = itm + props.nmodes(iph) - props.ntwmod(iph);
            gmode = 0.;
            statv.twfgr(itm, kkk) = 0.;
            twshx = props.twsh(iph, itm);

            for (size_t nsmx = 1; nsmx <= props.nsm(iph, kk); nsmx++) {
                // shifted counter for twin systems
                kts++;
                // absolute counter over all systems
                its++;
                if (statv.gamdot(kkk, its) > 0.) {
                    gabs = fabs(statv.gamdot(kkk, its)) * intv.tincr;
                    gmode += gabs;
                    statv.twfrsy(kkk, kts) += gabs / twshx;
                }
                statv.twfgr(itm, kkk) += statv.twfrsy(kkk, kts);
                rand = 1.0;
                //  rand=0.8+0.2*ran2(jran)   // adds stochasticity to selection
                if (statv.twfrsy(kkk, kts) * rand > twfmax) {
                    statv.iptsgr(kkk) = its;
                    statv.iptmgr(kkk) = itm;
                    twfmax = statv.twfrsy(kkk, kts);
                }
            }

            statv.twfrph(itm, iph) += (gmode / twshx) * (statv.wgt(kkk) / statv.wph(iph));
            // end of loop over twin modes in the phase
        }
        //kgx++;
        // end of loop over grains in the phase
    }

    // *********************************************************************
    //    this loop scans all grains in the phase/elem in a random way and
    //    tags them for checking against the twin-reorientation criterion.
    //    statv.imark(1:ngphase) tags with '1' the grains as they are PIcked.
    //    in each loop iterat a grain 'kkk' is randomly PIcked and reoriented
    //    in the pts if the 'real' twin fraction is not exceeded.
    // *********************************************************************

    ngphase = props.ngr(iph + 1) - props.ngr(iph);
    kleft = ngphase;
    for (size_t i = 1; i <= ngphase; i++) {
        statv.imark(i) = 0;
    }

    // tags one of the yet unchecked grains
    for (size_t index = 1; index <= ngphase; index++) {
        //v        rand=random2(jran)
        rand = 0.0;// random2(statv.jran);
        kpoint = size_t(rand * kleft) + 1;
        kleft = kleft - 1;
        kacum = 0;
        igr = 0;
        // identifies tagged
        while (igr < ngphase && kpoint != 0) {
            igr++;
            if (statv.imark(igr) == 0) {
                kacum++;
                if (kacum == kpoint) {
                    // absolute index of grain
                    kkk = igr + props.ngr(iph);
                    // index of grain relative to elemen
                    statv.imark(igr) = 1;
                    ipts = statv.iptsgr(kkk);
                    iptm = statv.iptmgr(kkk);

                    MATARIX <double> ag(&statv.ag(kkk, 1, 1), 3, 3);

                    //     labels the grain if the accumulated twin volume in the predominant
                    //     twinning system (pts) has exceeded a threshold value.
                    //     this grain will be completely reoriented by twinning.
                    //     reorientation is stopped when the effective reoriented fraction
                    //     reaches the 'real' twin fraction calculated from the shears.
                    //     thres1: min accum twin fraction in grain before twin reor is switc
                    //     thres1+thres2: eventual accum twin fr in gr required to twin reori

                    if (statv.ktwsmx(kkk) == 0.0 && ipts != 0) {
                        if (statv.eftwfr(iptm, iph) < statv.twfrph(iptm, iph)) {
                            thres1 = props.twthres(iph, iptm, 1);
                            thres2 = props.twthres(iph, iptm, 2);
                            twflim = thres1 + thres2 * statv.eftwfr(iptm, iph) / (statv.twfrph(iptm, iph) + 1.e-6);
                            if (statv.twfgr(iptm, kkk) > twflim) {
                                statv.eftwfr(iptm, iph) += statv.wgt(kkk) / statv.wph(iph);
                                statv.ktwsmx(kkk) = ipts;
                                //                   statv.ktwsmx(kkk)=0    // to suppress twin reorientn
                            }
                        }

                        // *** reorients by twinning grains tagged with ktwsmx>0 and updates the
                        // *** orientation matrix. this is done after having called update_orient
                        // *** which reorients the grain because of shear (including twin shears)

                        if (statv.ktwsmx(kkk) > 0) {
                            ktw = statv.ktwsmx(kkk);
                            for (size_t i = 1; i <= 3; i++) {
                                for (size_t j = 1; j <= 3; j++) {
                                    acrys(i, j) = ag(i, j);
                                }
                                btwin(i) = props.dbca(iph, ktw, i);
                            }

                            twin_orientation(btwin, atwin);


                            for (size_t i = 1; i <= 3; i++) {
                                for (size_t j = 1; j <= 3; j++) {
                                    ag(i, j) = 0.0;
                                    for (size_t k = 1; k <= 3; k++) {
                                        ag(i, j) += acrys(i, k) * atwin(k, j);
                                    }
                                }
                            }
                            // *** reset accumulated twinned fractions in the grain
                            for (size_t is = 1; is <= props.ntwsys(iph); is++) {
                                statv.twfrsy(kkk, is) = 0.0;
                            }
                            //feb
                            //fee
                            statv.ntwevents(kkk)++;
                            //v

                            //if(props.isectw(iphel, ktw) == 0  || statv.ntwevents(kkk) > 1) statv.ktwsmx(kkk)=-statv.ktwsmx(kkk); // supp; !changeme

                            if (props.isectw(iph, ktw) == 1.0 && statv.ntwevents(kkk) <= 1.0) statv.ktwsmx(kkk) = 0;// allow;


                            // *** reset crss in sl tw systems to initial (or other) values
                            //                 statv.gtotgr(kkk)=0
                            //             for(size_t is = 1; is <= nsyst(iphel)-ntwsys(iphel); is++){
                            //                   crss(is,kkk)=tau(is,0,iphel)
                            //                 }
                            //                 nslsys=nsyst(iphel)-ntwsys(iphel)
                            //             for(size_t is = nslsys+1; is <= nsyst(iphel); is++){
                            //                   crss(is,kkk)=5.0*crss(is,kkk)
                            //                 }

                            // end of if ktwsmx(kkk)>0
                        }
                        // end of if ktwsmx(kkk)=0
                    }
                    // end of if(kacum == kpoint)
                }
                // end of if(statv.imark(igr) == 0)
            }
            // end of do while(igr < ngphase)
        }
    }
    // end of do index=1,ngphase
}

//
// **********************************************************************
//     subroutine grain_stress (ex visc)   --->   version of 19/mar/2003
//
//     given a guess stress 'x(i)' and grain index 'kkk', solves interact
//     equation to find the stress 'x' compatible with macroscopic compli
//     if interx=-1 solves only the power law for a relaxed constraints c
//     if interx= 0 solves only the power law for a taylor case.
//     if interx> 0 solves interaction equation for a self-consistent cas
// **********************************************************************
// **********************************************************************
//     subroutine grain_rate_and_moduli (ex micro) -->  version 17/dec/02
//
//     given the stress 'x' in grain 'kkk', calculates strain-rate and
//     visco-plastic moduli using the rate sensitivity kinematic law.
// **********************************************************************
//
KOKKOS_FUNCTION void grain_stress_moduli(size_t interx, size_t jsc, size_t kcalmod, bool calcStress, StateVType& statv, PropType& props, IntVType& intv)
{
    //printf("Enter Grain Stress Moduli\n");
    for (size_t kkk = props.grbot; kkk <= props.grtop; kkk++) {

        size_t iph = statv.giph(kkk);

        //v      include 'vpsc7std.dim'
        const size_t nsystx = props.nsyst(iph);

        double taulim, taumax, rssx, eps;

        size_t irc, itmx, ierror, ierrlu;

        //inverse of taux
        xmat1(double, itaux, NSYS_, nsystx);
        xmat1(double, gamd0x, NSYS_, nsystx);
        xmat1(double, nrsx, NSYS_, nsystx);
        xmat1(double, isensex, NSYS_, nsystx);
        xmat2(double, scx, (5 * NSYS_), 5, nsystx);
        fmat2(double, xmastx, 5, 5);
        fmat1(double, x, 5);
        fmat1(double, xori, 5);
        fmat1(double, db, 5);


        //rate and moduli
        fmat1(double, sx, 5);
        fmat1(double, dx, 5);
        fmat1(double, dg, 5);
        xmat1(double, rss, NSYS_, nsystx);
        xmat1(double, gamdot, NSYS_, nsystx);


        gamdot = &statv.gamdot(kkk, 1);
        dg = &intv.dg(kkk, 1);

        ierrlu = 0;

        //     empiric algorithm to get taulim for nr subroutine (rl: 01/feb/00)

        taulim = 2. * pow((tnorm(intv.dbar) / statv.gamd0g(1)), (1.0 / (double)props.nrsmin));
        if (taulim < 2.0) taulim = 2.0;

        //     copy main arrays into auxiliar arrays for computational efficiency
        //     and to make 'newton_raphson' a 'stand-alone' subroutine

        isensex = &props.isense(iph, 1);

        if (interx == 1 && props.irsvar == 1) nrsx = intv.jjxrs;
        else                                  nrsx = &props.nrs(iph, 1);


        itaux = &statv.crss(kkk, 1);
        gamd0x = statv.gamd0g(kkk);
        scx = &statv.sch(kkk, 1, 1);

        for (size_t is = 1; is <= nsystx; is++) {
            itaux(is) = (1.0 / itaux(is));
        }

        if (calcStress) {

            x = &intv.stry(kkk, 1);

            irc = 0;

            // *** corrects stress 'x' if it exceeds the yield surface.
            taumax = 0.0;
            for (size_t is = 1; is <= nsystx; is++) {
                rssx = x(1) * scx(1, is) + x(2) * scx(2, is) + x(3) * scx(3, is) + x(4) * scx(4, is) + x(5) * scx(5, is);

                if (!(rssx > 0.0 || isensex(is) == 1.0)) rssx = 0.0;

                rssx *= itaux(is);

                if (fabs(rssx) > taumax) taumax = fabs(rssx);
            }



            if (taumax < 1.0e-10) {
                printf("taumax < 1.0e-10\n");
            }

            if (taumax > taulim || interx <= 0) {
                x /= taumax;
            }

            if (interx <= 0) {
                db = intv.dbar;
                xmastx = 0.0;

                if (props.grtop == 1) {
                    fmat1(double, aux5, 5); fmat2(double, aux55, 5, 5);
                    fmat1(double, deps6, 6); fmat2(double, deps33, 3, 3);

                    MATARIX <double> cgr(&intv.cgr(kkk, 1, 1), 6, 6);
                    MATARIX <double>  sg(&statv.sg(kkk, 1), 5);


                    // adjust xmastx and db for elasto-viscoplastic response
                    for (size_t i = 1; i <= 5; i++) {
                        for (size_t j = 1; j <= 5; j++) {
                            aux55(i, j) = cgr(i, j);
                        }
                    }

                    mat_inverse(aux55, 5, ierrlu); // deviatioric compliance

                    double itincr = 1.0 / intv.tincr;

                    for (size_t i = 1; i <= 25; i++) {
                        xmastx(i) = aux55(i) * itincr;
                    }

                    for (size_t i = 1; i <= 9; i++) {
                        deps33(i) = intv.dsimtot(i) * intv.tincr;
                    }

                    chg_basis26(deps6, deps33);

                    for (size_t i = 1; i <= 5; i++) {

                        //for (size_t j = 1; j <= 6; j++) {
                        //    aux5(i) += cgr(i, j) * deps6(j);
                        //}

                        aux5(i) = sumProd(6, &cgr(i, 1), deps6) + sg(i);
                    }

                    //matmul(cgr, deps6, aux5);



                    //for (size_t i = 1; i <= 5; i++) {
                    //    db(i) = 0.0;
                    //    for (size_t j = 1; j <= 5; j++) {
                    //        db(i) += aux55(i, j) / intv.tincr * (sg(j) + aux5(j));
                    //    }
                    //}

                    db = 0.0;
                    vecMul<5>(aux55, aux5, db);
                    db /= intv.tincr;

                }
            }
            else if (interx > 0) {
                //if (props.ishape(iph + 1) <= 1.0) {
                xmastx = &intv.xmastph(iph, 1, 1);

                //}
                //else {
                //	xmastx = 0.0;
                //}

                for (size_t i = 1; i <= 5; i++) {
                    db(i) = intv.dast(i) + statv.dzero(i);
                    db(i) += sumProd(5, &xmastx(i, 1), intv.sastav);
                }

            }

            // *** calls newton-raphson subroutine to calculate grain stress
            // *** internally it does a 5d (irc=0) or a 3d (irc=1) convergence.

            xori = x;

            itmx = 10000;
            eps = 5.0e-04;




            //
            //v      call newton_raphson (x,irc,itmx,eps,taulim,ierror ,iprx)
            newton_raphson(irc, itmx, eps, taulim, nsystx, x, db, xmastx, scx, itaux, gamd0x, nrsx, isensex, ierror);

            //v
            if (ierror > 0) {
                if (ierror == 1) printf(" singular system in newtraph -->\n");
                if (ierror == 1) printf(" cannot solve grain %zu in phase %zu\n", kkk, iph);
                if (ierror == 2) printf(" itmax was reach in newtraph -->\n");
                printf(" kinc %zu, k_cell: %zu, grain : %zu\n", intv.kinc, intv.k_cell, kkk);
                //std::cout << "x(1) : " << x(1) << " xmastx(1) : " << xmastx(1) << " " << std::endl;
                for (size_t i = 1; i <= 5; i++) {
                    x(i) = xori(i);
                }

                x = xori;
                //if (ierror == 1) exit(1);
                ierror = 0;

            }

            //for (size_t i = 1; i <= 5; i++) {
            //    statv.sg(kkk, i) = x(i);
            //}

            slice(5, &statv.sg(kkk, 1), x);
        }

        sx = &statv.sg(kkk, 1);
        //v
        //v
        for (size_t is = 1; is <= nsystx; is++) {
            double rssx = sx(1) * scx(1, is) + sx(2) * scx(2, is) + sx(3) * scx(3, is) + sx(4) * scx(4, is) + sx(5) * scx(5, is); //sumProd(5, nsystx, &sx(1), &scx(1, is)); //

            if (!(rssx > 0 || isensex(is) == 1.0)) rssx = 0.0;

            rssx *= itaux(is);

            double sign = sgn(rssx);

            rssx = fabs(rssx);

            rss(is) = gamd0x(is) * pow(rssx, (size_t)(nrsx(is) - 1.0)) * itaux(is);

            gamdot(is) = gamd0x(is) * pow(rssx, (size_t)nrsx(is)) * sign;
        }


        //     calculate strain-rate in grain 'dg'
        for (size_t i = 1; i <= 5; i++) {
            dg(i) = sumProd(nsystx, &scx(i, 1), gamdot);
        }

        dx = dg;

        if (kcalmod == 1) {
            MATARIX <double> dczero(&statv.dczero(kkk, 1), 5);
            fmat2(double, xmctg, 5, 5);

            //     calculate crystal compliance --> explain next lines//

            if (!(props.interaction == 2 || props.interaction == 3 || props.interaction == 4)) {
                for (size_t is = 1; is <= nsystx; is++) {
                    rss(is) *= nrsx(is);
                }
            }

            for (size_t i = 1; i <= 5; i++) {
                double dum = 0.0;
                for (size_t is = 1; is <= nsystx; is++) {
                    dum += rss(is) * scx(i, is) * scx(i, is);
                }
                xmctg(i, i) = dum; // sumProd(nsystx, rss, &scx(i, 1), &scx(i, 1));

                for (size_t j = i + 1; j <= 5; j++) {
                    xmctg(i, j) = sumProd(nsystx, rss, &scx(i, 1), &scx(j, 1));
                    xmctg(j, i) = xmctg(i, j);
                }
            }


            if (props.interaction == 1 || props.interaction == 5 || props.interaction == 0) {
                for (size_t i = 1; i <= 5; i++) {
                    dczero(i) = dx(i) - sumProd(5, &xmctg(i, 1), sx);
                }
            }
            else {
                dczero = 0.0;
            }

            slice(25, &statv.xmctg(kkk, 1, 1), xmctg);

            // kcalmod }
        }

        if (props.interaction == 0 && props.strain_control == 1) {
            fmat2(double, deps33, 3, 3); fmat1(double, deps6, 6); fmat1(double, aux6, 6);
            // pointer to the 6th row in this grain for cgr
            MATARIX <double> cgr6(&intv.cgr(kkk, 6, 1), 6);
            // calclulate hydrostatic stress in grains for interaction=0

            for (size_t i = 1; i <= 9; i++) deps33(i) = intv.dsimtot(i) * intv.tincr;


            chg_basis26(deps6, deps33);

            for (size_t i = 1; i <= 5; i++) {
                aux6(i) = dx(i);
            }
            aux6(6) = 0.0;

            double dum = 0.0;
            for (size_t i = 1; i <= 6; i++) {
                dum += cgr6(i) * (deps6(i) - aux6(i) * intv.tincr);
            }
            statv.sghyd(kkk) += dum / sqrt(3.0); // to cartesian from b-basis for hydrostatic stress;

        }

        slice(nsystx, &statv.gamdot(kkk, 1), gamdot);
        slice(5, &intv.dg(kkk, 1), dg);
    }
    return;
}


// ********************************************************************
//     subroutine elas_bc
// ********************************************************************
KOKKOS_FUNCTION void  elas_bc(double dt,
    MATARIX <double> xltg,
    MATARIX <double> dzero,
    MATARIX <double> dsimtot,
    MATARIX <double> csc,
    MATARIX <double> dbar,
    MATARIX <double> sbar,
    MATARIX <double> sbarold)
{
    fmat1(double, aux6, 6);
    fmat2(double, aux33, 3, 3);
    fmat2(double, c55, 5, 5);
    fmat1(double, deps6, 6);
    fmat2(double, deps33, 3, 3);
    fmat2(double, aux55, 5, 5);
    fmat1(double, aux5, 5);
    fmat2(double, aux33x, 3, 3);
    fmat1(double, aux5x, 5);

    size_t ierrlu;

    // strain inc
    //deps33 = intv.dsimtot; deps33 *= intv.tincr;
    for (size_t i = 1; i <= 9; i++) deps33(i) = dsimtot(i) * dt;

    chg_basis26(deps6, deps33);

    // strain rate calculation
    // --
    for (size_t i = 1; i <= 5; i++) {
        for (size_t j = 1; j <= 5; j++) {
            aux55(i, j) = xltg(i, j) + csc(i, j) * dt;
        }
    }

    mat_inverse(aux55, 5, ierrlu);

    // elastic predictor
    for (size_t i = 1; i <= 5; i++) {
        double dum = 0.0;
        for (size_t j = 1; j <= 6; j++) {
            dum += csc(i, j) * deps6(j);
        }
        aux5(i) = dum;
    }

    // auxiliary term
    for (size_t i = 1; i <= 5; i++) {
        aux5x(i) = 0.0;
        for (size_t j = 1; j <= 5; j++) {
            aux5x(i) += xltg(i, j) * dzero(j);
        }
        aux5x(i) += aux5(i) + sbarold(i);
    }

    // strain rate
    for (size_t i = 1; i <= 5; i++) {
        dbar(i) = 0.0;
        for (size_t j = 1; j <= 5; j++) {
            dbar(i) += aux55(i, j) * aux5x(j);
        }
    }
    // --

    // stress
    for (size_t i = 1; i <= 5; i++) {
        sbar(i) = 0.0;
        for (size_t j = 1; j <= 5; j++) {
            sbar(i) += xltg(i, j) * (dbar(j) - dzero(j));
        }
    }
    return;
}


KOKKOS_FUNCTION void vpsc_affine(size_t istep, StateVType& statv, PropType& props, IntVType& intv) {
    // record stress

        //auxiliary variables
    fmat1(double, aux5, 5);
    fmat1(double, aux6, 6);
    fmat2(double, aux33, 3, 3);


    double resratio, dummy, dbarn, dbarthresh, dbarnorm, resn, resoldn, sbaroldn, dum, xlambda;

    size_t it2, ierrlu;

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GPU TEST IN PROGRESS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    vmat2(double, sgold, props.grtop, 5); vmat1(double, sghydold, props.grtop); vmat2(double, sgtmp, props.grtop, 5);
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GPU TEST IN PROGRESS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fmat1(double, ddbarnr, 5); fmat1(double, ddbarsec, 5); fmat1(double, ddbar, 5);
    fmat1(double, res, 5); fmat2(double, dres_ddbar, 5, 5); fmat1(double, dsimtotb, 6);
    fmat1(double, szero, 5); fmat1(double, dbarold, 5);

    sgold = statv.sg;
    sgtmp = statv.sg;

    //matlabout(sgold);

    //for (size_t j = 1; j <= props.ngr(props.nph + 1); j++) {
    //    sghydold(j) = statv.sghyd(j);
    //}

    sghydold = statv.sghyd;

    // norm of plastic strain rate and total deviatoric strain rate
    dbarn = tnorm(intv.dbar);

    dummy = (intv.dsimtot(1, 1) + intv.dsimtot(2, 2) + intv.dsimtot(3, 3)) / 3.0;
    aux33 = intv.dsimtot;
    for (size_t i = 1; i <= 3; i++) {
        aux33(i, i) = aux33(i, i) - dummy;
    }

    dbarthresh = (tnorm(aux33)) * 1.0e-5;
    chg_basis26(dsimtotb, intv.dsimtot);

    // initialize for newton's method
    it2 = 0;
    resn = 2.0 * props.errd;
    resoldn = 0.0;
    ddbar = 0.0;
    dbarold = intv.dbar;
    sbaroldn = tnorm(intv.sbarold);

    if (sbaroldn == 0.0) sbaroldn = 1.0;

    if (props.ipr >= 1) printf("taylor solution:\n");

    // newton's loop
    while (resn > props.errd && it2 < props.itmaxext && dbarn > dbarthresh) {

        it2++;

        // find step length
        xlambda = 1.0;
        for (size_t it = 1; it <= 20; it++) {


            for (size_t i = 1; i <= 5; i++) {
                intv.dbar(i) = dbarold(i) + xlambda * ddbar(i);
            }

            // retrieve stress from previous increment
            // (grain_stress and grain_rate_and_moduli update stress to end of time increment for taylor)
            statv.sg = sgold;
            statv.sghyd = sghydold;

            // update stress and moduli

            dbarnorm = tnorm(intv.dbar);

            if (intv.kinc == 1) {

                for (size_t j = 1; j <= 5; j++) {
                    aux5(j) = intv.dbar(j) / dbarnorm;
                }

                for (size_t kkk = props.ngr(props.iphbot) + 1; kkk <= props.ngr(props.iphtop + 1); kkk++) {
                    //for (size_t j = 1; j <= 5; j++) {
                    //    intv.stry(kkk, j) = aux5(j);
                    //}

                    slice(5, &intv.stry(kkk, 1), aux5);
                }
            }
            else {
                intv.stry = sgtmp;
            }

            //grain_stress and rates
            //for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {
            //	for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {
            //
            //		grain_stress(0, kkk, iph, statv, props, intv);
            //
            //		grain_rate_and_moduli(0, 1, kkk, iph, statv, props, intv);
            //	}
            //}

            grain_stress_moduli(0, 0, 1, true, statv, props, intv);


            sgtmp = statv.sg;


            // *** calculate average stress and strain-rate

            for (size_t i = 1; i <= 5; i++) {
                intv.sav(i) = 0.0;
                intv.dav(i) = 0.0;
                for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {
                    for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {
                        intv.sav(i) += statv.sg(kkk, i) * statv.wgt(kkk);
                        intv.dav(i) += intv.dg(kkk, i) * statv.wgt(kkk);
                    }
                }
                // only required to calculate dvm
                intv.sbar(i) = intv.sav(i);
                if (props.iphtop == 1 && props.ngr(2) == 1.0) intv.dbar(i) = intv.dav(i);
            }

            chg_basis15(intv.sbar, intv.sdeviat);

            // residual
            for (size_t i = 1; i <= 5; i++) {
                aux6(i) = intv.dbar(i);
            }
            aux6(6) = 0.0;
            for (size_t i = 1; i <= 5; i++) {
                res(i) = 0.0;
                for (size_t j = 1; j <= 6; j++) {
                    res(i) -= intv.csc(i, j) * (dsimtotb(j) - aux6(j)) * intv.tincr;
                }
                res(i) += (intv.sbar(i) - intv.sbarold(i));
            }
            resn = tnorm(res) / sbaroldn;

            //printf(resn << "\n";

            resratio = resn / resoldn;

            // step cutting
            if (resoldn != 0.0 && resratio > 1.0) { xlambda *= 0.5; }

            else {
                resoldn = resn;
                dbarold = intv.dbar;
                //printf("recalc ddbar\n";
                break;
            }

        }

        //     calculate macroscopic moduli 'mtg' as the inverse
        //     of the average of the grain's stiffnesses

        //for (size_t i = 1; i <= 5; i++) {
        //    statv.dzero(i) = 0.0;
        //    szero(i) = 0.0;
        //
        //    for (size_t j = 1; j <= 5; j++) {
        //        statv.xmtg(i, j) = 0.0;
        //    }
        //}

        statv.dzero = 0.0;
        szero = 0.0;
        statv.xmtg = 0.0;



        if (props.iphtop == 1 && props.ngr(2) == 1.0) {

            for (size_t i = 1; i <= 5; i++) {
                statv.dzero(i) = statv.dczero(1, i);
                szero(i) = 0.0;
                for (size_t j = 1; j <= 5; j++) {
                    statv.xmtg(i, j) = statv.xmctg(1, i, j);
                }
            }

            statv.xltg = statv.xmtg;

            mat_inverse(statv.xltg, 5, ierrlu);
        }
        else {
            //more than one crystal

            //printf("ngr : " << props.ngr(1 + 1) << "\n";
            for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {
                for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {

                    fmat1(double, sczero, 5); fmat2(double, auxtan, 5, 5);
                    double wgt = statv.wgt(kkk);

                    MATARIX <double> dczero(&statv.dczero(kkk, 1), 5);

                    //for (size_t i = 1; i <= 5; i++) {
                    //    for (size_t j = 1; j <= 5; j++) {
                    //        auxtan(i, j) = statv.xmctg(kkk, i, j);
                    //    }
                    //}

                    slice(25, auxtan, &statv.xmctg(kkk, 1, 1));

                    mat_inverse(auxtan, 5, ierrlu);


                    //sczero = 0.0;
                    //vecMul<5>(auxtan, dczero, sczero);
                    //sczero *= -1;
                    for (size_t i = 1; i <= 5; i++) {
                        //sczero(i) = -sumProd(5, &auxtan(i, 1), dczero);
                        for (size_t j = 1; j <= 5; j++) {
                            sczero(i) -= auxtan(i, j) * dczero(j);
                        }
                    }



                    for (size_t i = 1; i <= 5; i++) {
                        statv.dzero(i) += dczero(i) * wgt;
                        szero(i) += sczero(i) * wgt;
                    }
                    for (size_t i = 1; i <= 25; i++) {
                        statv.xmtg(i) += auxtan(i) * wgt;
                    }

                }
            }
            //w
            statv.xltg = statv.xmtg;


            mat_inverse(statv.xmtg, 5, ierrlu);


            for (size_t i = 1; i <= 5; i++) {
                statv.dzero(i) = 0.0;
                for (size_t j = 1; j <= 5; j++) {
                    statv.dzero(i) -= statv.xmtg(i, j) * szero(j);
                }
            }

        }

        // true nr
        //dres_ddbar=intv.csc*intv.tincr+statv.xltg
        //call mat_inverse(dres_ddbar,5,ierrlu)
        fmat2(double, aux55, 5, 5);
        for (size_t i = 1; i <= 5; i++) {
            for (size_t j = 1; j <= 5; j++) {
                aux55(i, j) = intv.csc(i, j) * intv.tincr;
            }
        }
        //aux66 = intv.csc*intv.tincr;
        eye(dres_ddbar, 5, 1.0);
        matMul<5>(aux55, statv.xmtg, dres_ddbar);

        mat_inverse(dres_ddbar, 5, ierrlu);

        //matlabout(dres_ddbar);
        aux55 = 0.0;
        matMul<5>(statv.xmtg, dres_ddbar, aux55);

        dres_ddbar = aux55;

        for (size_t i = 1; i <= 5; i++) {
            ddbarnr(i) = 0.0;
            for (size_t j = 1; j <= 5; j++) {
                ddbarnr(i) -= dres_ddbar(i, j) * res(j);
            }
        }

        ddbar = ddbarnr;

        chg_basis15(intv.dbar, intv.dsim);
        chg_basis15(intv.sbar, intv.sdeviat);
        intv.scauchy = intv.sdeviat;

        // norm
        dbarn = tnorm(intv.dbar);

        // exit if we have only one grain
        if (props.iphtop == 1 && props.ngr(2) == 1.0) { break; }

    }

    // elastic behavior, recalculate sg
    if (dbarn > dbarthresh) {

        intv.sav = 0.0;
        for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {
            for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {
                MATARIX <double> sg(&statv.sg(kkk, 1), 5), cgr(&intv.cgr(kkk, 1, 1), 6, 6);
                //for (size_t i = 1; i <= 5; i++) {
                //    aux6(i) = intv.dbar(i);
                //}

                slice(5, &aux6(1), intv.dbar);
                aux6(6) = 0.0;

                for (size_t i = 1; i <= 5; i++) {
                    dum = 0.0;
                    for (size_t j = 1; j <= 6; j++) {
                        dum += cgr(i, j) * (dsimtotb(j) - aux6(j)) * intv.tincr; // elastic increment in stress
                    }
                    sg(i) = sgold(kkk, i) + dum;
                    intv.sav(i) += sg(i) * statv.wgt(kkk);//switch later
                }

                // kgx++;
            }
        }
        intv.sbar = intv.sav;
        chg_basis15(intv.sbar, intv.sdeviat);
        intv.scauchy = intv.sdeviat;

    }

    // von mises
    intv.svm = 0.0;
    intv.dvm = 0.;
    for (size_t i = 1; i <= 5; i++) {
        intv.svm += intv.sbar(i) * intv.sbar(i);
        intv.dvm += intv.sbar(i) * intv.dbar(i);
    }
    intv.svm = sqrt(intv.svm * 1.5);
    intv.dvm = intv.dvm / intv.svm;

    return;
}

//
//
//********************************************************************
//     subroutine vpsc     --->      version 29/nov/05
//********************************************************************
KOKKOS_FUNCTION void vpsc(size_t istep,
    StateVType& statv,
    PropType& props,
    IntVType& intv)
{

    fmat1(double, aux0, 1);// to trick eshelby


    double erreso, erraso, relsgr, rels, reld, rer1tg, dummy, pdil, sdev, sastmix, dbarn, dbarthresh;

    size_t jxrsini, jxrsfin, jxrstep, irs, itso, kso, it2, it1tg, iskip, ioption, ierrlu;

    // *** for a taylor case (interaction=0) solves vp equation for grain str
    // *** calculate strain-rate 'dg' and moduli 'xmctg' for every grain.
    pdil = 0;
    if (props.ipr > 1) printf("max sg : %f\n", max(statv.sg));
    jxrsini = props.nrs(1, 1);
    jxrsfin = props.nrs(1, 1);
    jxrstep = props.nrs(1, 1);
    //wx
    //wx   comment next line to perform
    //wx   rs loop at every step
    //wx
    if (istep > 1) jxrsini = jxrsfin;
    //
    irs = 0;
    for (size_t jxrs = jxrsini; jxrs <= jxrsfin; jxrs += jxrstep) {
        //v
        intv.jjxrs = jxrs;
        //v
        irs++;
        //feb
        //fee
        if (jxrsini != jxrsfin) {
            //        print *,
            //        print *, 'nrs iteration', irs
            //        print *,
        }
        //
        itso = 0;
        erreso = 2. * props.errso;
        erraso = 2. * props.errso;
        kso = 1;

        while ((kso == 1) && (erreso > props.errso || erraso > props.errso) && itso < props.itmaxso)
        {
            if (props.interaction != 5) kso = 0;
            itso++;
            if (props.interaction == 5) {
            }

            relsgr = 2 * props.errs;
            rels = 2 * props.errs;
            reld = 2 * props.errd;

            // *** outer loop: varies stress and compliance in the grains

            // norm of plastic strain rate and total deviatoric strain rate

            {
                fmat2(double, aux33, 3, 3);
                dbarn = tnorm(intv.dbar);
                dummy = (intv.dsimtot(1, 1) + intv.dsimtot(2, 2) + intv.dsimtot(3, 3)) / 3.0;
                aux33 = intv.dsimtot;
                for (size_t i = 1; i <= 3; i += 1) {
                    aux33(i, i) -= dummy;
                }
                dbarthresh = tnorm(aux33) * 1.0e-5;
            }


            it2 = 0;
            //
            while ((((relsgr > props.errs || reld > props.errd || rels > props.errs) && it2 < props.itmaxext)) && dbarn > dbarthresh)
            {
                it2++;

                fmat1(double, phave, 5); fmat2(double, bcinv, 5, 5);
                fmat2(double, xmtg, 5, 5); fmat2(double, xltg, 5, 5);

                xltg = statv.xltg;
                xmtg = statv.xmtg;

                // *** inner loop: varies overall tangent compliance

                it1tg = 0;
                rer1tg = 2 * props.errm;
                //
                while (rer1tg > props.errm && it1tg < props.itmaxint) {
                    it1tg++;


                    // **********************************************************************
                    //     calculate eshelby tensor s=s(mtg) and fs=[1/(i-s)]s.
                    //     s(mtg) will be the same for every grain in a given phase when
                    //     ishape.le.1.
                    //     it will be different for every grain when ishape.ge.2 because mtg
                    //     msec have to be rotated to ellipsoid axes and, in addition, the
                    //     size of the ellipsoid axes is different for every grain.
                    //     when mtg=n*msec then s(mtg)=s(msec), independent of ishape value.
                    // **********************************************************************

                    // *** skip every other calculation of mstar
                    iskip = 1;
                    if (iskip == 1) {
                        //
                        fmat4(double, c4sa, 3, 3, 3, 3);


                        chg_basis35(xltg, c4sa);

                        // **********************************************************************
                        // *** loop #1 over phases and grains
                        // **********************************************************************

                        //size_t kgx = 1;
                        for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {
                            fmat2(double, eigb, 3, 3);
                            fmat4(double, ega, 3, 3, 3, 3);
                            fmat4(double, rga, 3, 3, 3, 3);
                            fmat4(double, esa, 3, 3, 3, 3);
                            fmat4(double, rsa, 3, 3, 3, 3);
                            fmat2(double, e5, 5, 5);
                            fmat2(double, fs, 5, 5);
                            fmat2(double, xmast, 5, 5);

                            fmat1(double, axb, 3);

                            fmat2(double, pga, 3, 3);
                            fmat4(double, c4ga, 3, 3, 3, 3);
                            fmat2(double, ximsinv, 5, 5);
                            // if we're on the first iteration of phase loop: iph = 0 or 1.
                            // only done once in the phase loop.
                            if (props.ishape(iph + 1) <= 1.0) {
                                for (size_t i = 1; i <= 3; i++) {
                                    axb(i) = statv.axisph(iph + 1, 1, i);
                                    for (size_t j = 1; j <= 3; j++) {
                                        eigb(i, j) = statv.axisph(iph + 1, i + 1, j);//FCheck
                                    }
                                }
                            }


                            //for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {


                                // *** eshelby calculation for every phase or for every grain
                                // if we're in a single phase, on its first grain:
                                // (only needed once for each phase in the loop.)
                            //	if (props.ishape(iph + 1) >= 2.0 || kkk == props.ngr(iph) + 1) {


                                    //     rotation of stiffness 'c4sa' to ellipsoid principal axes

                                    //printf("xltg(1) = %e\n", xltg(1));
                                    //printf("c4sa(1) = %e\n", c4sa(1));

                            transpose(eigb, 3);
                            rotate_tens4(eigb, c4sa, c4ga);
                            transpose(eigb, 3);

                            if (props.icauchy == 0) ioption = 2;
                            if (props.icauchy == 1) ioption = 3;
                            //v      call eshelby(axb,c4ga,0.,ega,rga,aux33,pga,pdil,
                            //v     #             aux3333,aux3333,ioption



                            // compute eshelby tensor for each phase and grain.
                            // note: c4ga and axb are determined only by phase so
                            // they're the same for each grain of the same phase //////
                            eshelby(axb, c4ga, 0, ega, rga, aux0, pga, pdil, aux0, aux0, props, ioption);


                            // rotates the distortion, rotation and pressure eshelby tensors
                            // for the phase or for each grain back into sample axes.

                            rotate_tens4(eigb, ega, esa);
                            //

                            //rotate_tens4(eigb, ega, rga, esa, rsa);

                            chg_basis45(e5, esa);

                            for (size_t i = 1; i <= 5; i++) {
                                for (size_t j = 1; j <= 5; j++) {
                                    ximsinv(i, j) = xid(i, j) - e5(i, j);
                                }
                            }

                            mat_inverse(ximsinv, 5, ierrlu);

                            fs = 0.0;
                            matMul<5>(ximsinv, e5, fs);

                            //
                            //     calculate m* = inv(i-s) * s * mtg
                            //

                            xmast = 0.0;
                            matMul<5>(fs, xmtg, xmast); //sym

                            //for (size_t i = 1; i <= 5; i++) {
                            //	for (size_t j = i + 1; j <= 5; j++) {
                            //		xmast(i, j) = (xmast(i, j) + xmast(j, i)) * 0.5;
                            //		xmast(j, i) = xmast(i, j);
                            //	}
                            //}

                            //printf("xmtg(1) = %e\n", statv.xmtg(1));

                            if (props.interaction == 3) {
                                xmast *= 10;
                            }
                            else if (props.interaction == 4) {
                                if (props.irsvar == 1)      xmast *= intv.jjxrs;
                                else                        xmast *= props.nrs(1);
                            }


                            //
                            // *** copies tensors into phase arrays or into grain arrays.
                            // *** 'asph' or 'asgr' are used in subr. update_orientation to calculate
                            //     spin-rate deviations.
                            // *** 'd5ph' or 'd5gr' are used in subr. calc_cauchy to calculate
                            //     pressure deviations.

                            if (props.ishape(iph + 1) <= 1.0) {

                                fmat2(double, e5inv, 5, 5);  fmat4(double, einvsa, 3, 3, 3, 3);
                                rotate_tens4(eigb, rga, rsa);
                                e5inv = e5;
                                mat_inverse(e5inv, 5, ierrlu);
                                chg_basis35(e5inv, einvsa);

                                for (size_t i = 1; i <= 3; i++) {
                                    for (size_t j = 1; j <= 3; j++) {
                                        for (size_t k = 1; k <= 3; k++) {
                                            for (size_t l = 1; l <= 3; l++) {
                                                double dum = 0.0;
                                                for (size_t k1 = 1; k1 <= 3; k1++) {
                                                    for (size_t l1 = 1; l1 <= 3; l1++) {
                                                        dum += rsa(i, j, k1, l1) * einvsa(k1, l1, k, l);
                                                    }
                                                }
                                                intv.asph(iph, i, j, k, l) = dum;
                                            }
                                        }
                                    }
                                }
                            }
                            //feb
                            //fee
                            if (props.ishape(iph + 1) <= 1.0) {


                                slice(25, &intv.xmastph(iph, 1, 1), xmast);
                                //for (size_t i = 1; i <= 5; i++) {
                                //    for (size_t j = 1; j <= 5; j++) {
                                //        intv.xmastph(iph, i, j) = xmast(i, j);
                                //    }
                                //}
                            }
                            //feb
                            //fee
                            //  end of if(ishape(iph.ge.2  ||  kkk.eq.ngr(iph-1) - once for every phase combo.
                    //	}

                        //kgx++;

                    //}           //  end of do over grains   (loop #1)

                        }           //  end of do over phases   (loop #1)

                        //  end of eshelby calculation big if
                    }



                    // **********************************************************************
                    // *** loop #2 over phases and grains
                    // **********************************************************************

                    // *** solves iteratively sc equation for general case of different shape
                    // *** required for multi-phase when grains of each phase have different
                    //     or for ishape>1, when every grain has different shape, or both.

                    //for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {
                        //iphel = iph - props.iphbot + 1;


                    fmat2(double, xmtbave, 5, 5); fmat2(double, bcave, 5, 5); fmat1(double, xmtphave, 5);
                    fmat2(double, xmtnew, 5, 5); fmat1(double, dzero, 5);


                    dzero = statv.dzero;
                    xmtphave = 0.0;
                    bcave = 0.0;
                    xmtbave = 0.0;
                    phave = 0.0;

                    for (size_t kkk = props.grbot; kkk <= props.grtop; kkk++) {
                        // shouldn't this be done outs

                        fmat2(double, xmctg, 5, 5);
                        fmat1(double, dczero, 5);
                        xmctg = &statv.xmctg(kkk, 1, 1);
                        dczero = &statv.dczero(kkk, 1);


                        fmat2(double, bc, 5, 5);
                        fmat1(double, phic, 5);
                        fmat2(double, xmcbctg, 5, 5);
                        fmat1(double, xmcphc, 5);

                        fmat2(double, xmast, 5, 5);
                        fmat2(double, bc1, 5, 5);
                        fmat2(double, bc2, 5, 5);
                        fmat1(double, aux5, 5);
                        size_t iph = statv.giph(kkk);

                        double wgt = statv.wgt(kkk);

                        if (props.ishape(iph + 1) <= 1.0) {
                            xmast = &intv.xmastph(iph, 1, 1);


                            for (size_t i = 1; i <= 25; i++) {
                                bc2(i) = xmast(i) + xmtg(i);
                            }

                        }
                        //feb
                        //fee

                        for (size_t i = 1; i <= 25; i++) {
                            bc1(i) = xmast(i) + xmctg(i);
                        }


                        mat_inverse(bc1, 5, ierrlu);

                        //matrix multiplication A * B = C when B is symmetric 
                        bc = 0.0;
                        matMul<5>(bc1, bc2, bc);

                        for (size_t i = 1; i <= 5; i++) {
                            double dum = 0.0;
                            for (size_t j = 1; j <= 5; j++) {
                                dum += bc1(i, j) * (dzero(j) - dczero(j));
                            }

                            phic(i) = dum;
                        }

                        xmcbctg = 0.0;
                        matMul<5>(xmctg, bc, xmcbctg);

                        xmcphc = dczero;
                        vecMul<5>(xmctg, phic, xmcphc);

                        for (size_t i = 1; i <= 25; i++) {
                            bcave(i) += bc(i) * wgt;
                            xmtbave(i) += xmcbctg(i) * wgt;
                        }

                        for (size_t i = 1; i <= 5; i++) {
                            xmtphave(i) += xmcphc(i) * wgt;
                        }



                        if (props.ibcinv == 1) {
                            //phic *= statv.wgt(kkk);
                            //phave += phic;
                            for (size_t i = 1; i <= 5; i++) {
                                phave(i) += phic(i) * wgt;
                            }
                        }

                    }
                    //}

                    // ****************************************************************
                    // *** end of loop #2 over grains and phases
                    // ****************************************************************
                    if (props.ibcinv == 1)  bcinv = bcave;

                    else                    eye(bcinv, 5, 1.0);



                    mat_inverse(bcinv, 5, ierrlu);

                    xmtnew = 0.0;
                    matMul<5>(xmtbave, bcinv, xmtnew);


                    // **********************************************************************
                    //     mold-mnew (tangent) comparison

                    //z        rer1tg=tmismatch5(xmtg,xmtnew,5,5)
                    rer1tg = tmismatch(25, xmtg, xmtnew);
                    //feb
                    //fee

                    //printf("%f\n", rer1tg);

                    xmtg = xmtnew;
                    xltg = xmtnew;

                    mat_inverse(xltg, 5, ierrlu);

                    //
                    if (props.interaction == 1 || props.interaction == 5) {
                        fmat1(double, dzero1, 5);
                        dzero1 = 0.0;
                        vecMul<5>(xmtnew, phave, dzero1);

                        for (size_t i = 1; i <= 5; i++) {
                            statv.dzero(i) = xmtphave(i) - dzero1(i);
                        }
                    }
                    else {
                        statv.dzero = 0.0;
                    }
                }




                // **********************************************************************
                //     boundary conditions:
                //     * if diagonal imposed calls subroutine 'state_5x5' (rl: 26/1/00)
                //       and solves stress-strain components in deviatoric space
                //     * else, calls subrout 'state_6x6' with dsim(3,3) scauchy(3,3)
                //       plus indices of the known components. solves for unknown compone
                //       in 6-dim cauchy space.
                //


                //for (size_t i = 1; i <= 5; i++) {
                //    dbaux(i) = intv.dbar(i) - statv.dzero(i);
                //}

                //


                elas_bc(intv.tincr, xltg, statv.dzero, intv.dsimtot, intv.csc, intv.dbar, intv.sbar, intv.sbarold);

                chg_basis15(intv.dbar, intv.dsim);
                chg_basis15(intv.sbar, intv.sdeviat);
                intv.scauchy = intv.sdeviat;

                //
                intv.svm = sqrt(3. / 2.) * tnorm(intv.sbar);

                // **********************************************************************
                //     s* d* required for different grain shapes
                //     sastav(i) is used inside grain_stress
                //     sastbar(i) is used below to calculate dast(i)

                sastmix = 0.;
                if (props.interaction != 3 && props.interaction != 4) sastmix = 0.75;

                intv.sastav = 0.0;
                intv.sastbar = 0.0;

                {
                    fmat1(double, aux5, 5);
                    for (size_t i = 1; i <= 5; i++) aux5(i) = intv.sav(i) - phave(i);
                    vecMul<5>(bcinv, aux5, intv.sastav);

                    for (size_t i = 1; i <= 5; i++) aux5(i) = intv.sbar(i) - phave(i);
                    vecMul<5>(bcinv, aux5, intv.sastbar);
                }


                for (size_t i = 1; i <= 5; i++) {
                    double savex = intv.sastav(i);
                    double sbarx = intv.sastbar(i);
                    intv.sastav(i) = (1. - sastmix) * savex + sastmix * sbarx;
                    intv.sastbar(i) = sastmix * savex + (1. - sastmix) * sbarx;
                }

                //for (size_t i = 1; i <= 5; i++) {
                //    intv.dast(i) = 0.;
                //    for (size_t j = 1; j <= 5; j++) {
                //        intv.dast(i) += statv.xmtg(i, j) * intv.sastbar(j);
                //    }
                //}

                intv.dast = 0.0;
                vecMul<5>(xmtg, intv.sastbar, intv.dast);

                // *** save current stress in grain.
                // *** starting from converged 'xmtg' and 'stry' recalculates 'sg' using
                //     interaction equation.
                // *** uses new 'sg' to calculate strain-rate 'dg', shear rates 'gamdot'
                //     moduli 'xmctg' for every grain.

                intv.stry = statv.sg;
                statv.xmtg = xmtg;
                statv.xltg = xltg;

                grain_stress_moduli(1, 1, 1, true, statv, props, intv);


                // *** update average stress and average strain-rate.
                // *** calculates stress convergence in grains.

                sdev = 0.;


                for (size_t i = 1; i <= 5; i++) {
                    double dums = 0.;
                    double dumd = 0.;
                    //kgx = 1;
                    //for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {
                    for (size_t kkk = props.grbot; kkk <= props.grtop; kkk++) {
                        double wgt = statv.wgt(kkk);
                        dums += statv.sg(kkk, i) * wgt;
                        dumd += intv.dg(kkk, i) * wgt;
                        sdev += powm((statv.sg(kkk, i) - intv.stry(kkk, i)), 2) * wgt;
                        //kgx++;
                    }

                    intv.sav(i) = dums;
                    intv.dav(i) = dumd;
                    //}
                }

                //printf("sav dav\n");
                //fullout(intv.sav);
                //fullout(intv.dav);
                //printf("end sav dav\n");

                //printf(intv.sav;

                // *** check consistency of local, macroscopic and average magnitudes.
                //      a) <sg-stry> .lt.err
                //      b) /dbar-dav/.lt.err
                //      c) /sbar-sav/.lt.err
                relsgr = tnorm(intv.sbar);
                //if (relsgr != 0.0)
                relsgr = sqrt(sdev) / relsgr;
                reld = tmismatch(5, intv.dbar, intv.dav);
                rels = tmismatch(5, intv.sbar, intv.sav);

                // *** defines empirical mix coef for sastav sastbar based on convergen
                //     convnew=relsgr+reld+rels
                //     convdif=convnew-convold
                //     sastmix=0.
                //     if(it2.gt.10 .and. convdif.gt.0.) sastmix=0.25
                //     if(it2.gt.20 .and. convdif.gt.0.) sastmix=0.50
                //     if(it2.gt.30 .and. convdif.gt.0.) sastmix=0.75
                //     convold=convnew
                //feb
                //fee

                dbarn = tnorm(intv.dbar);
                // end of (do..while) for outer iteration (it2) on 'sg'
            }

            //printf("Inner\n");
            //
            intv.svm = 0.;
            intv.dvm = 0.;
            for (size_t i = 1; i <= 5; i++) {
                intv.svm += intv.sbar(i) * intv.sbar(i);
                intv.dvm += intv.sbar(i) * intv.dbar(i);
            }
            intv.svm = sqrt(intv.svm * 1.5);
            //if (intv.svm != 0.0)
            intv.dvm = intv.dvm / intv.svm;
            //feb
            //fee
            // end of (do..while) for so iteration (itso)
        }
        //feb
        //fee

        // end of do jxrs for nrs iteration
    }


    //feb
    //fee
    return;
}




KOKKOS_FUNCTION void update_temperature(double svm, double evmp, double& evmpold, double& temp, double& dtemp, double wrplastic, double& wplastic, double b_int)
{
    double rho, cpa0, cpa1, cpa2, eta, cp, devmp;


    eta = b_int; // 1 default for poly tantalum, adjusted based on number of active slips;
    rho = 16640; // kg/m^3;
    cpa0 = 145.5; // j/kgk;
    cpa1 = 0.009544; // j/kgk^2;
    cpa2 = -68900; // jk/kg;

    cp = cpa0 + cpa1 * temp + cpa2 / (temp * temp);

    devmp = evmp - evmpold;
    evmpold = evmp;

    wrplastic = svm * devmp;
    wplastic = wplastic + wrplastic;

    dtemp = eta * wrplastic * 1e6 / (rho * cp);

    temp = temp + dtemp;
    return;
}

//
// **********************************************************************
//     subroutine initial_state_guess     --->     version 07/dec/05
//
//     if strain is imposed, uses a first stress guess colinear with the
//     strain rate and calculates grain stress given by taylor.
//     if stress is imposed, sets the grain stress equal to the macroscop
//     calculates grain strain rates 'dg' for either case.
//     calculates average stress 'sav' and strain-rate 'dav'.
//     calculates visco-plastic moduli 'xmctg' for every grain.
//     calculates an initial guess for the macroscopic visco-plastic modu
// **********************************************************************
//
//v      subroutine initial_state_guess */
KOKKOS_FUNCTION void initial_state_guess(StateVType& statv, PropType& props, IntVType& intv)
{

    fmat1(double, szero, 5); fmat1(double, sczero, 5);
    //MATARIX <double> auxtan(5, 5), szero(5), sczero(5);

    double dbarnorm;
    size_t ierrlu;

    //fullout(intv.dbar);

    dbarnorm = tnorm(intv.dbar);

    //kgx = 1;
    if (props.strain_control == 1) {

        fmat1(double, aux5, 5);
        for (size_t j = 1; j <= 5; j++) aux5(j) = intv.dbar(j) / dbarnorm;

        for (size_t kkk = props.grbot; kkk <= props.grtop; kkk++) {
            for (size_t j = 1; j <= 5; j++) {
                intv.stry(kkk, j) = aux5(j);
            }
        }
    }
    else if (props.strain_control == 0) {

        fmat1(double, aux5, 5);
        for (size_t j = 1; j <= 5; j++) aux5(j) = intv.sbar(j);

        for (size_t kkk = props.grbot; kkk <= props.grtop; kkk++) {
            for (size_t j = 1; j <= 5; j++) {
                statv.sg(kkk, j) = aux5(j);
            }
        }
    }


    //Calculate stress and rates in each grain
    grain_stress_moduli(0, 0, 1, props.strain_control == 1, statv, props, intv);

    intv.dav = 0.0;
    intv.sav = 0.0;
    for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {
        for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {
            for (size_t i = 1; i <= 5; i++) {
                intv.dav(i) += intv.dg(kkk, i) * statv.wgt(kkk);
                intv.sav(i) += statv.sg(kkk, i) * statv.wgt(kkk);
            }
        }
    }

    //printf("sg1 : %e, %e, %e, %e, %e\n", statv.sg(1, 1), statv.sg(1, 2), statv.sg(1, 3), statv.sg(1, 4), statv.sg(1, 5));


    //     calculate initial guess for macroscopic moduli 'mtg' as the invers
    //     of the average of the grain's stiffnesses

    statv.dzero = 0.0;
    szero = 0.0;
    statv.xmtg = 0.0;

    for (size_t kkk = props.grbot; kkk <= props.grtop; kkk++) {

        fmat2(double, auxtan, 5, 5);
        auxtan = &statv.xmctg(kkk, 1, 1);

        double wgt = statv.wgt(kkk);

        MATARIX <double> dczero(&statv.dczero(kkk, 1), 5);

        mat_inverse(auxtan, 5, ierrlu);


        if (ierrlu == 1) {
            printf("auxtan singular\n");
            for (size_t i = 1; i <= 5; i++) {
                for (size_t j = 1; j <= 5; j++) {
                    printf("%f, ", auxtan(i, j));
                }
                printf("\n");
            }
            printf("\n");
        }

        for (size_t i = 1; i <= 5; i++) {
            sczero(i) = 0.0;
            for (size_t j = 1; j <= 5; j++) {
                sczero(i) -= auxtan(i, j) * dczero(j);
            }
        }

        for (size_t i = 1; i <= 5; i++) {
            statv.dzero(i) += dczero(i) * wgt;
            szero(i) += sczero(i) * wgt;
        }

        for (size_t i = 1; i <= 25; i++) {
            statv.xmtg(i) += auxtan(i) * wgt;
        }
    }


    //w
    statv.xltg = statv.xmtg;

    mat_inverse(statv.xmtg, 5, ierrlu);


    for (size_t i = 1; i <= 5; i++) {
        statv.dzero(i) = 0.0;
        for (size_t j = 1; j <= 5; j++) {
            statv.dzero(i) -= statv.xmtg(i, j) * szero(j);
        }
    }

    if (props.strain_control == 0) intv.dbar = intv.dav;
    if (props.strain_control == 1) intv.sbar = intv.sav;

    //fullout(intv.dbar);

    intv.svm = 0.0;
    intv.dvm = 0.0;
    for (size_t i = 1; i <= 5; i++) {
        intv.svm += intv.sbar(i) * intv.sbar(i);
        intv.dvm += intv.sbar(i) * intv.dbar(i);
    }
    intv.svm = sqrt(intv.svm * 1.5);
    //if (intv.svm != 0.0) 
    intv.dvm = intv.dvm / intv.svm;

    return;
}

// ----------------------------------------------------
// purpose: vpsc update of one cell over one increment.
// ----------------------------------------------------
KOKKOS_FUNCTION void evol_vpsc(StateVType& statv, PropType& props, IntVType& intv, const ViewCArrayKokkos <double> &vel_grad, double* stress_in, double dtime, size_t k_cell, ViewCMatrixKokkos <double>& ddsdde) {

    //auxilliary variables
    fmat1(double, aux5, 5);  fmat1(double, aux6, 6); fmat2(double, aux33, 3, 3); fmat2(double, aux55, 5, 5); fmat2(double, aux66, 6, 6);

    fmat2(double, dstrandev, 3, 3);
    //save previous increment for rates
    fmat2(double, xmtgold, 5, 5); fmat1(double, dzeroold, 5);
    //local stress array
    fmat2(double, dstress, 3, 3); fmat1(double, dstressb, 6); fmat1(double, stressb, 6);
    //local strain array
    fmat2(double, dstran, 3, 3); fmat1(double, dstranb, 6);
    //linear stress estimation
    fmat2(double, stress_lin, 3, 3);
    //local rotation matrices
    fmat2(double, drot, 3, 3); fmat2(double, drot_t, 3, 3);
	fmat2(double, D, 3, 3);  fmat2(double, W, 3, 3);  fmat2(double, L, 3, 3);  fmat2(double, dW, 3, 3);  


    //fmat2(double, ddsdde, 6, 6);
    double dum1, dum2, xerr, plratio, dbarthresh, dbartn;
    double dstran_norm, trace, deltaevmt;
    size_t ierrlu, iel, kinc;

    // assign passed-in structs to module structs
    // ---
    for (size_t i = 1; i <= 9; i++) L(i) = vel_grad(i-1);
	
    for (size_t i = 1; i <= 3; i++) {
        for (size_t j = 1; j <= 3; j++) {
			W(i, j) = (L(i, j) - L(j, i)) * 0.5;
			D(i, j) = (L(i, j) + L(j, i)) * 0.5;
			
			dW(i, j) = W(i, j)*dtime;
			dstran(i, j) = D(i, j) * dtime;
		}
	}
	
	rodrigues(dW, drot);
	
				
    //dtime = dtime_in*1.0e3;
	//for (size_t i = 1; i <= 9; i++)
	//{
	//	dstran(i) = dstran_in(i);
	//	drot(i) = drot_in(i);
	//}
	
	
	MATARIX stress(&stress_in[0], 3, 3);
	
	
	statv.kinc = statv.kinc + 1;
	kinc = statv.kinc;
	
	//printf("kinc : %zu; dstran(3, 3) : %16.4e; stress(3, 3) : %16.8f;\n", kinc, dstran_in(3, 3), stress(3, 3));

    if (props.ipr == 1) printf("\n*****************************\n");
    if (props.ipr == 1) printf("     enter evol_vpsc\n");
    if (props.ipr == 1) printf("*****************************\n");
    if (props.ipr == 1) printf("k_cell: %zu\n", k_cell);
    if (props.ipr == 1) printf("kinc: %zu\n", kinc);
    if (props.ipr == 1) printf("dstran : %f, %f, %f\n", dstran(1, 1), dstran(2, 2), dstran(3, 3));
    if (props.ipr == 1) printf("stress : %f, %f, %f\n", stress(1, 1), stress(2, 2), stress(3, 3));
	




    // perform linear stress update and determine
    // whether to solve the full problem
    // ---
    if (kinc > 1) {
        //for when we need more than one auxiliary matrix of size (6) or (6, 6)
        fmat1(double, auxx6, 6); fmat2(double, auxx66, 6, 6);
        fmat2(double, xmtg66, 6, 6); fmat1(double, stressoldb, 6);
        fmat1(double, dzero6, 6);

        // linear extrapolation
        // ---
        if (props.ipr == 1) printf("linear extrapolation\n");


        // rotate the accumulated strain increment for applied drot
        rotate_mat(statv.dstranacc, drot); // rotate;


        //fullout(statv.dstranacc);

        // rotate the stress for applied  drot
        rotate_mat(statv.stressold, drot); // rotate

        // update the accumulated strain increment, time increment and rotation increment
        statv.dstranacc += dstran;

        aux33 = dstran;
        trace = (aux33(1, 1) + aux33(2, 2) + aux33(3, 3));
        aux33(1, 1) = aux33(1, 1) - trace / 3.0;
        aux33(2, 2) = aux33(2, 2) - trace / 3.0;
        aux33(3, 3) = aux33(3, 3) - trace / 3.0;

        dum1 = tnorm(aux33);

        statv.dstranaccvm += sqrt(2.0 / 3.0) * dum1;
        statv.dtimeacc += dtime;

        aux33 = 0.0;
        matMul<3>(drot, statv.drotacc, aux33);
        statv.drotacc = aux33;

        // dstran to b-basis
        chg_basis(dstranb, statv.dstranacc, 2, 6);

        // stress to b-basis
        chg_basis(stressoldb, statv.stressold, 2, 6);


        // xmtg to 6x6 dim
        //aux55 = statv.xmtgdot;
        //aux55 *= statv.dtimeacc;
        //aux55 += statv.xmtg;
        for (size_t i = 1; i <= 25; i++) aux55(i) = statv.xmtg(i) + statv.dtimeacc * statv.xmtgdot(i);

        xmtg66 = 0.0;

        for (size_t i = 1; i <= 5; i++) {
            for (size_t j = 1; j <= 5; j++) {
                xmtg66(i, j) = aux55(i, j);
            }

            //aux5(i) = statv.dzero(i) + statv.dzerodot(i) * statv.dtimeacc;
        }


        // dzero to 6 dim
        //aux5 = statv.dzero + statv.dzerodot * statv.dtimeacc;
        for (size_t i = 1; i <= 5; i++) {
            aux5(i) = statv.dzero(i) + statv.dzerodot(i) * statv.dtimeacc;
        }


        dzero6(6) = 0.0;
        for (size_t i = 1; i <= 5; i++) {
            dzero6(i) = aux5(i);
        }

        // update stress
        aux66 = 0.0;
        matMul<6>(statv.cscold, xmtg66, aux66);

        aux66 *= statv.dtimeacc;

        //for (size_t i = 1; i <= 6; i++) {
        //    aux66(i, i) += 1.0;
        //}

        peye(aux66, 6, 1.0);


        mat_inverse(aux66, 6, ierrlu);

        for (size_t i = 1; i <= 6; i++) {
            aux6(i) = dstranb(i) - dzero6(i) * statv.dtimeacc;
        }

        auxx6 = 0.0;
        vecMul<6>(statv.cscold, aux6, auxx6);

        for (size_t i = 1; i <= 6; i++) aux6(i) = auxx6(i) + stressoldb(i);
        //aux6 = auxx6;
        //aux6 += stressoldb;

        stressb = 0.0;
        vecMul<6>(aux66, aux6, stressb);

        // stress to tensor
        chg_basis(stressb, stress_lin, 1, 6);

        // updated ddsdde
        auxx66 = 0.0;
        matMul<6>(aux66, statv.cscold, auxx66);
        aux66 = auxx66;
        for (size_t i = 1; i <= 36; i++) ddsdde(i) = aux66(i);
        if (props.ipr == 1) printf("dstranvmth: %f\n   ", statv.dstranvmth);



        // calcluate the plastic ratio
        aux6 = dzero6;
        vecMul<6>(xmtg66, stressb, aux6);
        dum2 = tnorm(aux6);

        for (size_t i = 1; i <= 5; i++) aux6(i) = dstranb(i) / statv.dtimeacc - aux6(i);

        aux6(6) = 0.0;
        dum1 = tnorm(aux6);
        plratio = dum2 / (dum1 + dum2);

        if (plratio < 0.9) {
            // write(props.ipru,'(5x,a,3e16.8)')'elastic or unloading, use el-pl thresh';
            if (props.ipr == 1) printf("elastic or unloading, use el-pl thresh\n");
            statv.dstranvmth = props.dstranvmth_elpl;
        }

        // determine whether to accept the linear update or solve the full problem
        if (statv.dstranaccvm < statv.dstranvmth) {//false) {//

            if (props.ipr == 1) printf("accept linear updated\n");
            // use the linear extrapolation and exit
            stress = stress_lin;
            if (props.ipr == 1) printf("stress : %f, %f, %f\n", stress(1, 1), stress(2, 2), stress(3, 3));


            return;
        }
        else {

            if (props.ipr == 1) printf("solve the full problem\n");

            // set dstran, drot and dtime, and execute full solution
            dstran = statv.dstranacc; // in current configuration, rotated;
            dtime = statv.dtimeacc;
            stress = statv.stressold; // in current configuration, rotated;
            drot = statv.drotacc;

            // reset accumulated values to zero
            statv.dstranacc = 0.0;
            statv.dstranaccvm = 0.0;
            statv.dtimeacc = 0.0;
            eye(statv.drotacc, 3, 1.0);

            if (props.ipr == 1) printf("stress : %f, %f, %f\n", stress(1, 1), stress(2, 2), stress(3, 3));


        }

    }

    // full solution starts here
    // ---

    // updates orientation, hardening and shape of every grain making
    // a linear extrapolation of the calculated rate from beginning of time increment
    // ---



    // calculate macroscopic spin from drot (to be used for rotation of quantities)
    drot2spin(drot, dtime, aux33);

    intv.rotbar = aux33;
    //call rodrigues(intv.rotbar*dtime, aux33) // verify the spin
    statv.udot += intv.rotbar; // statv.udot contains udot from t and here we add current rotation;


    // old stress to b-basis
    chg_basis(intv.sbarold, stress, 2, 5); // for taylor stress it will be rotated in update_orient;


    if (tnorm(statv.udot) > 1e-12) {
        if (props.ipr == 1) printf("update state vars\n");

        // reinitialize the intv with previous increment rates for updating
        intv.dvm = statv.dvm;
        intv.udot = statv.udot;
        intv.dav = statv.dav;
        intv.asph = statv.asph;
        intv.dg = statv.dg;
        intv.tincr = dtime;

        // update plastic von mises strain
        statv.evmp += intv.dvm * intv.tincr;
        statv.svm = intv.svm;

        update_fij(0, statv, intv); // updates deform tensor of element;
        update_shape(0, statv, props.interaction);     // updates shape of element;

        if (props.ihardlaw <= 1) {

            // update twinning (not verified)
            for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {
                //size_t iphel = iph - props.iphbot + 1;
                if (props.ntwmod(iph) != 0) update_twinning(iph, statv, props, intv);
            }

            // update orientation
            if (props.iupdori == 1) update_orientation(statv, props, intv);

            // update shape
            for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {
                update_fij(iph, statv, intv); // updates def tensor of phase grains;
                if (props.iupdshp == 1) update_shape(iph, statv, props.interaction);
            }

            if (props.iupdhar == 1) {
                // voce hardening plus predominant twin reorientation scheme
                if (props.ihardlaw == 0) update_crss_voce(2, statv, props, intv);

                if (props.itemphard == 1) update_crss_temp(statv, props, intv);

            }

        }
    }


    // reinitialize intv with zero values and store xmtg and dzero from beginning of time inc
    //if (.not. allocated(intv.sastav)) call init_intv()
    intv.udot = 0.;
    intv.dsim = 0.;
    intv.rotbar = 0.;
    intv.kinc = kinc;
    intv.k_cell = k_cell;
    xmtgold = statv.xmtg;
    dzeroold = statv.dzero;
    if (kinc == 1) xmtgold = 0.0;
    if (kinc == 1) dzeroold = 0.0;
    intv.tincr = dtime;

    if (props.ipr == 1) printf("updating of schmid tensors and elastic stiffness\n");

    // update schmid tensors
    update_schmid(props.grtop, props.nsyst, props.schca, statv.ag, statv.sch, statv.giph);


    // elastic self-consistent stiffness
    elsc(0, statv, props, intv);


    // save stiffness
    statv.cscold = intv.csc;


    // calculate acumulated rotation
    aux33 = 0.0;
    matMul<3>(drot, statv.rot_acum, aux33);
    statv.rot_acum = aux33;

    //// inverse rotation
    //transpose(statv.rot_acum, qg2l);
    //
    //// transformation matrices
    //transpose(qg2l, qg2l_t);


    // solve the full problem
    //---
    dstran_norm = tnorm(dstran);
    if (dstran_norm <= 1.0e-15 && kinc == 1) { //deal with small first kinc deformation;

        dstran = 1.0e-14; //small dummy strain inc for initialization of variables;
        dstran_norm = tnorm(dstran);
    }

    if (props.ipr == 1) printf("call vpsc\n");

    // dstran to rate
    //intv.dsimtot = dstran; intv.dsimtot /= dtime;
    for (size_t i = 1; i <= 9; i++) intv.dsimtot(i) = dstran(i) / dtime;

    //matlabout(intv.dbar);

    // initial guess for dbar
    dbartn = tnorm(statv.dbart);

    dum1 = (intv.dsimtot(1, 1) + intv.dsimtot(2, 2) + intv.dsimtot(3, 3)) / 3.0;
    aux33 = intv.dsimtot;

    aux33(1, 1) -= dum1;
    aux33(2, 2) -= dum1;
    aux33(3, 3) -= dum1;

    iel = 0;
    dum1 = tnorm(aux33);
    dbarthresh = 1.0e1 * dum1 * 1.0e-5; // threshold value scaled;


    if (props.ipr == 1) printf("dbartn : %f   dbarthresh : %f\n", dbartn, dbarthresh);
    //
    if ((kinc == 1 || dbartn < dbarthresh) && props.ngr(props.iphtop + 1) > 1) {//false) {//

        if (props.ipr == 1) printf("test elastic behavior\n");

        // test elastic behavior
        chg_basis(dstranb, dstran, 2, 6);

        for (size_t i = 1; i <= 5; i++) {
            for (size_t j = 1; j <= 5; j++) {
                aux55(i, j) = intv.csc(i, j);
            }
        }


        intv.sbar = intv.sbarold;
        vecMul<5>(aux55, dstranb, intv.sbar); // elastic predictor;

        props.strain_control = 0;
        initial_state_guess(statv, props, intv);

        if (props.ipr == 1) printf("dbar : %f, %f, %f, %f, %f\n", intv.dbar(1), intv.dbar(2), intv.dbar(3), intv.dbar(4), intv.dbar(5));

        if (false) {//tnorm(intv.dbar) < dbarthresh * 0.1) {//
            iel = 1; // adopt sachs response;

            if (props.ipr == 1) printf("adopt sachs response : dbarnorm = %f\n", tnorm(intv.dbar));
        }
        else {
            iel = 0;
            chg_basis(intv.dbar, intv.dsimtot, 2, 5);
            dum1 = tnorm(intv.dbar);
            intv.dbar *= (dbarthresh / dum1); // set the norm to scaled threshold;
        }
    }

    else if ((kinc == 1 || dbartn < dbarthresh) && props.ngr(props.iphtop) == 1.0) {

        iel = 0;
        chg_basis(intv.dbar, intv.dsimtot, 2, 5);
        dum1 = tnorm(intv.dbar);
        intv.dbar /= (dum1 * dbarthresh); // set the norm to scaled threshold;
    }
    else {

        iel = 0;
        intv.dbar = statv.dbart; // use previous solution;

    }

    if (props.ipr == 1) printf("iel = %zu\n", iel);

    if (iel == 0) {
        // vpsc
        props.strain_control = 1;
        if (props.interaction > 0) {
            //printf("dbar : %e, %e, %e, %e, %e\n", intv.dbar(1), intv.dbar(2), intv.dbar(3), intv.dbar(4), intv.dbar(5));
            initial_state_guess(statv, props, intv);
            vpsc(0, statv, props, intv); // vpsc;
        }
        else {
            vpsc_affine(0, statv, props, intv); // taylor;
        }
        props.strain_control = 0;
    }

    // total stress
    chg_basis(dstranb, dstran, 2, 6);

    for (size_t i = 1; i <= 5; i++) aux6(i) = intv.dbar(i);


    aux6(6) = 0.0;

    //dstressb = dstranb - aux6 * dtime;
    for (size_t i = 1; i <= 6; i++) dstressb(i) = dstranb(i) - dtime * aux6(i);

    aux6 = 0.0;
    vecMul<6>(intv.csc, dstressb, aux6);
    dstressb = aux6;



    chg_basis(dstressb, dstress, 1, 6);
    if (props.interaction == 0) {
        dum1 = (stress(1, 1) + stress(2, 2) + stress(3, 3)) / 3.0;// hydrostatic stress;
        chg_basis(intv.sbarold, stress, 1, 5); // for taylor stress it will be rotated in update_orient;

        //printf("hydrostatic stress : " << dum1 << "\n");

        for (size_t i = 1; i <= 3; i++) {
            stress(i, i) += dum1;
        }
    }


    if (props.interaction == 0) intv.dsim = intv.dsimtot;


    if (props.ipr == 1) printf("stress = (%f, %f, %f)\n", stress(1, 1), stress(2, 2), stress(3, 3));
    if (props.ipr == 1) printf("dstress = (%f, %f, %f)\n", dstress(1, 1), dstress(2, 2), dstress(3, 3));

    // update stress
    stress += dstress;

    //fullout(intv.dbar);

    chg_basis(intv.dbar, aux33, 1, 5);

    //fullout(intv.dbar);

    //if(props.ipr==1)write(props.ipru,'(x,a,9e16.8)')'plastic strain inc:',dtime*aux33;
    //if(props.ipr==1)write(props.ipru,'(x,a,9e16.8)')'stress:',stress;

    // jacobian
    aux66 = intv.ssc;
    for (size_t i = 1; i <= 5; i++) {
        for (size_t j = 1; j <= 5; j++) {
            aux66(i, j) += statv.xmtg(i, j) * dtime;
        }
    }

    mat_inverse(aux66, 6, ierrlu);
    //chg_basis(aux66, ddsdde, 3, 6);
    //ddsdde = aux3333;
    for (size_t i = 1; i <= 36; i++) ddsdde(i) = aux66(i);
    //ddsdde = aux66;

    // copy intv variables to corresponding statev and update rate of moduli change
    intv.udot = intv.dsim;
    statv.dbart = intv.dbar;
    statv.dvm = intv.dvm;
    statv.udot = intv.udot;
    statv.dav = intv.dav;
    statv.asph = intv.asph;
    statv.dg = intv.dg;

    for (size_t i = 1; i <= 25; i++) {
        statv.xmtgdot(i) = (statv.xmtg(i) - xmtgold(i)) / dtime;
    }


    for (size_t i = 1; i <= 5; i++) {
        statv.dzerodot(i) = (statv.dzero(i) - dzeroold(i)) / dtime;
    }


    // calculate deviatoric increment
    dstrandev = dstran;
    trace = dstran(1, 1) + dstran(2, 2) + dstran(3, 3);
    dstrandev(1, 1) -= trace / 3.0;
    dstrandev(2, 2) -= trace / 3.0;
    dstrandev(3, 3) -= trace / 3.0;

    // von mises quantities
    deltaevmt = sqrt(2.0 / 3.0) * tnorm(dstrandev);
    //dvmt = deltaevmt / intv.tincr;
    statv.evmt += deltaevmt;
    intv.epsacu = statv.evmt;

    // define the threshold vm strain increment
    //if(props.ipr==1)write(props.ipru,'(x,a)')'determine the dstranvmth';

    statv.stressold = stress;
    dum2 = tnorm(intv.dbar);
    chg_basis(dstranb, dstran, 2, 6);
    //aux6 = 0.0;

    for (size_t i = 1; i <= 5; i++) {
        aux5(i) = dstranb(i) / dtime - intv.dbar(i);
    }



    dum1 = tnorm(aux5);
    plratio = dum2 / (dum1 + dum2);

    for (size_t i = 1; i <= 9; i++) {
        aux33 = stress_lin(i) - stress(i);
    }




    if ((kinc == 1) || (statv.evmp <= props.stranvm_pl) || (plratio < 0.9)) {//

        // in the elasto-plastic transition use small strain inc
        statv.dstranvmth = props.dstranvmth_elpl;

        //if(props.ipr==1)write(props.ipru,'(3x,a)')'elasto-plastic region, use dstranvmth_elpl';
    }
    else {

        // define the difference between linear update and true stress update


        for (size_t i = 1; i <= 9; i++) aux33(i) = stress_lin(i) - stress(i);

        xerr = tnorm(aux33) / tnorm(stress);



        // increase the threshold if the difference is small
        if (xerr < 0.01) {
            statv.dstranvmth *= 1.25;
        }
        else {
            statv.dstranvmth *= 0.75;
        }

        // limit the threshold values
        if (statv.dstranvmth > props.dstranvmth_max) statv.dstranvmth = props.dstranvmth_max;
        if (statv.dstranvmth < props.dstranvmth_min) statv.dstranvmth = props.dstranvmth_min;

        //if(props.ipr==1)write(props.ipru,'(3x,a)')'plastic region';
        //if(props.ipr==1)write(props.ipru,'(3x,a,e16.8)')'xerr:',xerr;
        //
    }

    return;
}

// ---- BEGINNING OF HOST CALCULATIONS ----//

template <typename T> KOKKOS_FUNCTION bool ismember(T a, MATARIX <T>& b) {

    size_t sz = b.size();
    bool out = false;
    for (size_t i = 1; i <= sz; i++) {
        if (b(i) == a) {
            out = true;
            return out;
        }
    }

    return out;
}


KOKKOS_FUNCTION void nextline(std::fstream& in, size_t nline) {
    char next;
    size_t line = 0;
    while (in.get(next))
    {
        if (next == '\n')  line++;

        if (line >= nline) break;
    }
}

KOKKOS_FUNCTION void gauss_legendre(double x1, double x2, MATARIX <double>& x, MATARIX <double>& w, size_t n)
{
    //      parameter(eps=3.d-14)
    //
    //     changed by r.l. -  8/2/97
    //
    const double eps = 1.0e-07;

    double xm, xl, xn, xi, z = 100.0, z1 = 0.0, p1, p2, xj, p3, pp;
    size_t m, iter;

    m = (n + 1) / 2;
    xm = 0.5 * (x1 + x2);
    xl = 0.5 * (x2 - x1);
    xn = n;
    pp = 0.0;
    for (size_t i = 1; i <= m; i++) {
        xi = i;
        z = cos(PI * (xi - 0.25) / (xn + 0.5));
        //
        iter = 0;
        while (abs(z - z1) > eps) {
            iter = iter + 1;
            //
            //     r.l. 8/2/97
            //
            if (iter > 10000) {
                ////$        write(*,*)'gauleg warning: tol 1.e-07 never reached - err = ', &
                ////$             abs(z-z1)
                return;
            }
            //
            p1 = 1.0;
            p2 = 0.0;
            for (size_t j = 1; j <= n; j++) {
                xj = j;
                p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (xj - 1.0) * p3) / xj;
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;

        }
        x(i) = xm - xl * z;
        x(n + 1 - i) = xm + xl * z;
        w(i) = 2.0 * xl / ((1.0 - z * z) * pp * pp);
        w(n + 1 - i) = w(i);
    }
    return;
}


KOKKOS_FUNCTION void eshelby_init(PropType& props) {

    // **********************************************************************
    //     initialization run
    //     calculates gauss-legendre integration points and weights in the
    //     interval [0,pi].
    //     initializes arrays associated with each point to aKOKKOS_FUNCTION void repeating
    //     its calculation at every call.
    // **********************************************************************
    //MATARIX <double> xph(NGAUMX_), xth(NGAUMX_), wph(NGAUMX_), wth(NGAUMX_);
    //MATARIX <double> aa2x(3, 3), aaww2x(3, 3);

    fmat1(double, xph, NGAUMX_); fmat1(double, xth, NGAUMX_); fmat1(double, wph, NGAUMX_); fmat1(double, wth, NGAUMX_); fmat2(double, aa2x, 3, 3); fmat2(double, aaww2x, 3, 3);

    props.ngaussph(1) = 16;
    props.ngaussth(1) = 16;
    props.ngaussph(2) = 32;
    props.ngaussth(2) = 32;
    props.ngaussph(3) = 16;
    props.ngaussth(3) = 32;

    for (size_t cas = 1; cas <= 3; cas++) {
        if (props.ngaussph(cas) > NGAUMX_ || props.ngaussth(cas) > NGAUMX_) {
            ////$            print *, ' dimension ngaumx exceeded in subr eshelby ////'
            //stop;
        }

        // Calculate the start of index for the case
        size_t inds = 0;
        for (size_t i = 1; i < cas; i++) {
            inds += props.ngaussph(i) * props.ngaussth(i);
        }

        //v
        double zero = 0.0;
        //v
        gauss_legendre(zero, PI, xph, wph, props.ngaussph(cas));
        gauss_legendre(zero, PI, xth, wth, props.ngaussth(cas));

        // *** integration [0,pi][0,pi] adds a factor 2 in eqs. b11 & b14.

        for (size_t ith = 1; ith <= props.ngaussth(cas); ith++) {
            double sinth, costh, simbtet;
            sinth = sin(xth(ith));
            costh = cos(xth(ith));
            simbtet = wth(ith) * sinth / (2.0 * PI);

            for (size_t iph = 1; iph <= props.ngaussph(cas); iph++) {
                size_t ny = iph + (ith - 1) * props.ngaussph(cas);

                MATARIX <double> alpha(&props.alpha(inds + ny, 1), 3);

                props.ww(inds + ny) = simbtet * wph(iph);
                alpha(1) = sinth * cos(xph(iph));
                alpha(2) = sinth * sin(xph(iph));
                alpha(3) = costh;




                //for (short i = 1; i <= 3; i++) {
                //	for (short j = 1; j <= 3; j++) {
                //		//aa2x(i, j) = props.alpha(cas, ny, i) * props.alpha(cas, ny, j);
                //		//aaww2x(i, j) = aa2x(i, j) * props.ww(cas, ny);
                //
                //		aa2x(i, j) = alpha(i) * alpha(j);
                //	}
                //}
                //
                //voigt(aa1, aa2x, 2);
                //for (size_t i = 1; i <= 6; i++) {
                //    props.aa1(cas, ny, i) = aa1x(i);
                //    props.aaww1(cas, ny, i) = aaww1x(i);
                //}
                // *** array aww is used only if icauchy=1.
               //for (short i = 1; i <= 3; i++) {
               //    props.aww(cas, ny, i) = props.alpha(cas, ny, i) * props.ww(cas, ny);
               //}
            }
        }
        // end of do cas=1,3
    }

    return;
    // endif for ioption=0
}


KOKKOS_FUNCTION void init_crss_voce(size_t ioption, StateVType& statv, PropType& props)
{

    //v      include 'vpsc7std.dim'

    //size_t iph, iphel, kkk, is;

    if (ioption == 1) {
        //z
        for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {
            //iphel = iph - props.iphbot + 1;

            for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {
                for (size_t is = 1; is <= props.nsyst(iph); is++) {
                    statv.crss(kkk, is) = props.tau(iph, is, 1);
                }
            }
        }

    }

    return;
}

// zf add temperature-derived initial crss
KOKKOS_FUNCTION void init_crss_temp(StateVType& statv, PropType& props)
{
    //size_t iph, iphel, kkk, is;

    statv.crss = 0.0;
    //z
    for (size_t iph = props.iphbot; iph <= props.iphtop; iph++) {
        //iphel = iph - props.iphbot + 1;
        for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {
            for (size_t is = 1; is <= props.nsyst(iph); is++) {
                //statv.crss(kkk, is) = props.tau0_mode(iph, is, 1) + 
                //	props.tau0_mode(iph, is, 2) * exp(-statv.temp_ini / props.tau0_mode(iph, is, 3));

                statv.crss(kkk, is) = props.tau0_mode(iph, is, 1) * exp(-(statv.temp - props.tau0_mode(iph, is, 3)) / props.tau0_mode(iph, is, 2));
            }
        }
    }




    return;
}


KOKKOS_FUNCTION void write_tex(size_t iph, std::string filename, PropType& props, StateVType& statv) {

    //size_t i, j, kkk, n;
    fmat2(double, aa, 3, 3);
    double eul1, eul2, eul3, rad2deg = 180.0 / PI;

    char line[800];

    std::ofstream texfile(filename.c_str());

    if (texfile.is_open()) {
        texfile << "\n\n\nB   " << props.ngr(iph + 1) - props.ngr(iph) << "\n";

        for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {
            MATARIX <double> ag(&statv.ag(kkk, 1, 1), 3, 3);
            //for (size_t i = 1; i <= 3; i++) {
            //    for (size_t j = 1; j <= 3; j++) {
            //        aa(j, i) = statv.ag(kkk, i, j);
            //    }
            //}

            euler(1, eul1, eul2, eul3, ag);

            eul1 *= rad2deg;
            eul2 *= rad2deg;
            eul3 *= rad2deg;

            sprintf(line, "%16.12F  %16.12F  %16.12F  %20.12E", eul1, eul2, eul3, statv.wgt(kkk));

            texfile << line << "\n";
        }

        texfile.close();
    }

    else {
        printf("error opening texture file\n");
    }

    return;

}

//
//  Reads Texture File and writes it to props and statv
//
KOKKOS_FUNCTION void data_texture(size_t iph, std::string filename, PropType& props, StateVType& statv)
{
    double ncrys, totwgt, eul1, eul2, eul3, deg2rad = PI / 180.0;
    double da, db, dc;
    fmat2(double, fijx, 3, 3); fmat2(double, fnew, 3, 3); fmat2(double, aa, 3, 3);
    //MATARIX <double> fijx(3, 3), fnew(3, 3), aa(3, 3);

    std::string line, nomen;
    std::ifstream texfile(filename.c_str());

    // dummy lines for strain data
    std::getline(texfile, line);
    std::getline(texfile, line);
    std::getline(texfile, line);


    texfile >> nomen >> ncrys;

    props.ngr(iph + 1) = props.ngr(iph) + ncrys;

    //printf("Number of Crystals in Phase %zu : %f\n", iph, ncrys);

    // read in euler angles
    totwgt = 0.0;
    for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {

        MATARIX <double> ag(&statv.ag(kkk, 1, 1), 3, 3);

        texfile >> eul1 >> eul2 >> eul3 >> statv.wgt(kkk);

        eul1 *= deg2rad;
        eul2 *= deg2rad;
        eul3 *= deg2rad;

        euler(2, eul1, eul2, eul3, ag);

        transpose(ag, 3);

        totwgt += statv.wgt(kkk);
    }

    texfile.close();


    for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {
        statv.wgt(kkk) = statv.wgt(kkk) / totwgt * statv.wph(iph);
    }


    // **********************************************************************
    //     initial f tensor and eigenvectors:
    //
    //     * if ishape=0 assumes same initial shape (axisph) and orientation
    //       (eulerph) for all the grains in the phase.
    //       calculates eshelby tensor with the average grain shape and does
    //       not keep track of local grain shape.
    //
    //     * if ishape=1 assumes same initial shape (axisph) and orientation
    //       (eulerph) for all the grains in the phase.
    //       calculates eshelby tensor with the average grain shape and keeps
    //       track of local grain shape.
    //
    //     * if ishape=2 assumes same initial shape (axisph) and orientation
    //       (eulerph) for all the grains in the phase.
    //       calculates eshelby tensor with the individual grain shape and
    //       keeps track of local grain shape.
    //
    //     * if ishape=3 reads individual initial grain axes orientation
    //       from fileaxes (da,db,dc), assumes same shape (axisph) for all
    //       the grains in the phase.
    //       calculates eshelby tensor with the individual grain shape and
    //       keeps track of local grain shape.
    //
    //     * if ishape=4 reads individual initial grain axes orientation
    //       and shape from fileaxes (da,db,dc,ax(1),ax(2),ax(3)).
    //       calculates eshelby tensor with individual grain shapes and
    //       keeps track of local grain shape
    //
    //       'eulerph' angles of (g) wrt (s).
    //       'aa'    transforms from (s) to (g)
    //       'fijx'  columns are grain axes expressed in grain system
    //       'fnew'  columns are grain axes expressed in sample system

    da = statv.eulerph(iph + 1, 1);
    db = statv.eulerph(iph + 1, 2);
    dc = statv.eulerph(iph + 1, 3);

    euler(2, da, db, dc, aa);

    fijx = 0.0;
    fnew = 0.0;

    for (size_t i = 1; i <= 3; i++) {
        fijx(i, i) = statv.axisph(iph + 1, 1, i);
    }

    for (size_t i = 1; i <= 3; i++) {
        for (size_t j = 1; j <= 3; j++) {
            for (size_t m = 1; m <= 3; m++) {
                fnew(i, j) += aa(m, i) * fijx(m, j);
            }
        }
    }

    for (size_t i = 1; i <= 3; i++) {
        for (size_t j = 1; j <= 3; j++) {
            statv.fijph(iph + 1, i, j) = fnew(i, j);
        }
    }


    return;
}

//
// **********************************************************************
//     subroutine crystal_symmetry   --->   version 04/sep/06
//
// *** if ioption=1:
//     reads crystal symmetry 'icrysym' and unit cell parameters.
//     generates vectors 'cvec(i,n)' of the unit cell.
//     generates symmetry operators 'h(i,j,nsymop)' for crystal symmetry.
// *** if ioption=2:
//     reads miller indices of systems in 3 or 4-index notation 'isn(i)'
//     'isb(i)'. calculates normal burgers vectors 'sn(i)' 'sb(i)'
// *** if ioption=3:
//     reads miller indices of diffraction planes 'isn(i)' and spherical
//     angles 'chi , eta'of diffraction direction.
//     generates crystallographically equivalent orientations sneq(i,n) o
//     a sn(i) by applying all the symmetry operations to it.
//     discards multiplicity and defines 'npol'
// *** simmetry parameter icrysym:
//        1: cubic
//        2: hexagonal
//        3: trigonal
//        4: tetragonal
//        5: orthorhombic
//        6: monoclinic
//        7: triclinic
// **********************************************************************
KOKKOS_FUNCTION void crystal_symmetry(size_t ioption, std::fstream& file, PropType& props, size_t& icrysym, MATARIX <double>& sn, MATARIX <double>& sneq, MATARIX <double>& sb, size_t& npol) {
    fmat3(double, hx, 6, 3, 3); fmat1(double, itag, 24);
    fmat1(double, isn, 4); fmat1(double, isb, 4);
    fmat1(double, cdim, 3); fmat1(double, cang, 3);
    //MATARIX <double>  hx(6, 3, 3);
    //MATARIX <double>  itag(24), isn(4), isb(4);
    //MATARIX <double>  cdim(3), cang(3);

    double ang, chi, eta, snnor, sbnor, snpro, prod;
    size_t mn, nrot, nind;

    int isign;

    isn = 0.0;



    // -----------------------------------------------------------------------

    if (ioption == 1) {
        std::string crysym;
        std::string crycmp[7] = { "CUBIC", "HEXAG", "TRIGO", "TETRA", "ORTHO", "MONOC", "TRICL" };
        nextline(file, 1);
        file >> crysym;
        nextline(file, 1);
        icrysym = 0;

        //printf(crysym;

        for (size_t i = 1; i <= 7; i++) {
            if (crysym == crycmp[i - 1])
                icrysym = i;
        }
        if (icrysym == 0) {
            printf("cannot recognize the crystal symmetry\n");
            exit(1);
        }

        //char test[50]; for(size_t ii = 0; ii < 50; ii++) file.get(test[ii]);


        file >> cdim(1) >> cdim(2) >> cdim(3) >> cang(1) >> cang(2) >> cang(3);
        nextline(file, 1);


        for (size_t i = 1; i <= 3; i++) {
            cang(i) *= (PI / 180.0);
        }
        props.cvec(1, 1) = 1.;
        props.cvec(2, 1) = 0.;
        props.cvec(3, 1) = 0.;
        props.cvec(1, 2) = cos(cang(3));
        props.cvec(2, 2) = sin(cang(3));
        props.cvec(3, 2) = 0.;
        props.cvec(1, 3) = cos(cang(2));
        props.cvec(2, 3) = (cos(cang(1)) - cos(cang(2)) * cos(cang(3))) / sin(cang(3));
        props.cvec(3, 3) = sqrt(1. - pow(props.cvec(1, 3), 2) - pow(props.cvec(2, 3), 2));

        for (size_t i = 1; i <= 3; i++) {
            for (size_t j = 1; j <= 3; j++) {
                props.cvec(i, j) = cdim(j) * props.cvec(i, j);
            }
        }


        props.h = 0.0;

        hx = 0.0;

        // --- identity operation ---> triclinic all symmetries
        for (size_t i = 1; i <= 3; i++) {
            props.h(i, 1, i) = 1.0;
        }
        props.nsymop = 1;

        // --- 180 deg rotation around (001) ---> orthorhombic, monoclinic
        if (icrysym == 5 || icrysym == 6) {
            props.h(2, 1, 1) = cos(PI);
            props.h(2, 2, 2) = cos(PI);
            props.h(2, 3, 3) = 1.0;
            props.h(2, 1, 2) = -sin(PI);
            props.h(2, 2, 1) = sin(PI);
            props.nsymop = 2;
        }

        // --- x-mirror y-mirror ---> orthorhombic
        if (icrysym == 5) {
            props.h(3, 1, 1) = -1.0;
            props.h(3, 2, 2) = 1.0;
            props.h(3, 3, 3) = 1.0;

            props.h(4, 1, 1) = 1.0;
            props.h(4, 2, 2) = -1.0;
            props.h(4, 3, 3) = 1.0;
            props.nsymop = 4;
        }

        // --- cubic symmetry
        if (icrysym == 1) {

            // --- rotations of (pi/3) (2*pi/3) around <111>
            hx(1, 1, 3) = 1.0;
            hx(1, 2, 1) = 1.0;
            hx(1, 3, 2) = 1.0;

            hx(2, 1, 2) = 1.0;
            hx(2, 2, 3) = 1.0;
            hx(2, 3, 1) = 1.0;

            for (size_t m = 1; m <= 2; m++) {
                for (size_t n = 1; n <= props.nsymop; n++) {
                    mn = m * props.nsymop + n;
                    for (size_t i = 1; i <= 3; i++) {
                        for (size_t j = 1; j <= 3; j++) {
                            for (size_t k = 1; k <= 3; k++) {
                                props.h(mn, i, j) = props.h(mn, i, j) + hx(m, i, k) * props.h(n, k, j);
                            }
                        }
                    }
                }
            }
            props.nsymop = mn;

            // --- mirror across the plane (110)
            hx(3, 1, 2) = 1.0;
            hx(3, 2, 1) = 1.0;
            hx(3, 3, 3) = 1.0;

            for (size_t n = 1; n <= props.nsymop; n++) {
                mn = props.nsymop + n;
                for (size_t i = 1; i <= 3; i++) {
                    for (size_t j = 1; j <= 3; j++) {
                        for (size_t k = 1; k <= 3; k++) {
                            props.h(mn, i, j) = props.h(mn, i, j) + hx(3, i, k) * props.h(n, k, j);
                        }
                    }
                }
            }
            props.nsymop = mn;

            // --- rotations of 90, 180, 270 around x3

            for (size_t m = 1; m <= 3; m++) {
                ang = PI / 2.0 * float(m);
                hx(m, 1, 1) = cos(ang);
                hx(m, 2, 2) = cos(ang);
                hx(m, 3, 3) = 1.0;
                hx(m, 1, 2) = -sin(ang);
                hx(m, 2, 1) = sin(ang);
                hx(m, 1, 3) = 0.0;
                hx(m, 3, 1) = 0.0;
                hx(m, 2, 3) = 0.0;
                hx(m, 3, 2) = 0.0;
            }

            for (size_t m = 1; m <= 3; m++) {
                for (size_t n = 1; n <= props.nsymop; n++) {
                    mn = m * props.nsymop + n;
                    for (size_t i = 1; i <= 3; i++) {
                        for (size_t j = 1; j <= 3; j++) {
                            for (size_t k = 1; k <= 3; k++) {
                                props.h(mn, i, j) = props.h(mn, i, j) + hx(m, i, k) * props.h(n, k, j);
                            }
                        }
                    }
                }
            }
            props.nsymop = mn;

            //end of condition for icrysym=1
        }

        // --- hexagonal, trigonal and tetragonal symmetry

        if (icrysym >= 2 && icrysym <= 4) {
            if (icrysym == 2) nrot = 6;
            if (icrysym == 3) nrot = 3;
            if (icrysym == 4) nrot = 4;

            // *** mirror plane at 30 deg or 60 deg or 45 deg with respect to x1
            ang = PI / double(nrot);
            props.h(2, 1, 1) = pow(cos(ang), 2) - pow(sin(ang), 2);
            props.h(2, 2, 2) = -props.h(2, 1, 1);
            props.h(2, 3, 3) = 1.0;
            props.h(2, 1, 2) = 2. * cos(ang) * sin(ang);
            props.h(2, 2, 1) = props.h(2, 1, 2);
            props.nsymop = 2;

            // --- rotations of 2*pi/6 around axis <001> for hexagonals.
            // --- rotations of 2*pi/3 around axis <001> for trigonals.
            // --- rotations of 2*pi/8 around axis <001> for trigonals.
            for (size_t nr = 1; nr <= nrot - 1; nr++) {
                ang = nr * 2.0 * PI / nrot;
                hx(nr, 1, 1) = cos(ang);
                hx(nr, 2, 2) = cos(ang);
                hx(nr, 3, 3) = 1.0;
                hx(nr, 1, 2) = -sin(ang);
                hx(nr, 2, 1) = sin(ang);
            }

            for (size_t m = 1; m <= nrot - 1; m++) {
                for (size_t n = 1; n <= props.nsymop; n++) {
                    mn = m * props.nsymop + n;
                    for (size_t i = 1; i <= 3; i++) {
                        for (size_t j = 1; j <= 3; j++) {
                            for (size_t k = 1; k <= 3; k++) {
                                props.h(mn, i, j) = props.h(mn, i, j) + hx(m, i, k) * props.h(n, k, j);
                            }
                        }
                    }
                }
            }
            props.nsymop = mn;

            //end of condition for icrysym= 2,3,4
        }

        //     write(4,*)
        //     write(4,'(''  # of symmetry operations='',i4)') nsymop
        //     write(4,'(''  symmetry matrices'')')
        //     write(4,'(i3,9f7.3)') (n,((h(i,j,n),j=1,3),i=1,3),n=1,nsymop)

        //end of condition for ioption=1

        return;
    }

    // ---------------------------------------------------------------------
    //   reads miller-bravais indices for cubic (1), tetragonal (4), ortho-
    //   rhombic (5), monoclinic (6) triclinic (7) systems in 3-index notat
    //   for hexagonal (2) trigonal (3) systems reads 4-index notation.
    //   converts indices of plane normal and slip direction into normalized
    //   vectors sn(i) and sb(i), respectively.
    // ---------------------------------------------------------------------

    if (ioption == 2 || ioption == 3) {
        //char test[50]; for(size_t ii = 0; ii < 50; ii++) file.get(test[ii]);
        nind = 3;
        if (icrysym == 2 || icrysym == 3) nind = 4;
        if (ioption == 2) {

            for (size_t i = 1; i <= nind; i++) {
                isn(i) = 10000;
                isb(i) = 10000;
            }

            for (size_t i = 1; i <= nind; i++) {
                file >> isn(i);
            }
            for (size_t i = 1; i <= nind; i++) {
                file >> isb(i);
            }

            nextline(file, 1);
        }
        else if (ioption == 3) {
            for (size_t i = 1; i <= nind; i++) {
                file >> isn(i);
            }
            file >> chi >> eta;
            //feb
            //fee
            eta = eta * PI / 180.0;
            chi = chi * PI / 180.0;
            sb(1) = cos(eta) * sin(chi);
            sb(2) = sin(eta) * sin(chi);
            sb(3) = cos(chi);
        }
        if (nind == 4) isn(3) = isn(4);

        sn(1) = isn(1) / props.cvec(1, 1);
        sn(2) = (isn(2) - props.cvec(1, 2) * sn(1)) / props.cvec(2, 2);
        sn(3) = (isn(3) - props.cvec(1, 3) * sn(1) - props.cvec(2, 3) * sn(2)) / props.cvec(3, 3);

        snnor = (tnorm(sn));

        for (size_t j = 1; j <= 3; j++) {
            sn(j) = sn(j) / snnor;
            if (abs(sn(j)) < 1.e-03) sn(j) = 0.;
        }

        if (ioption == 2) {

            if (icrysym == 2 || icrysym == 3) {
                isb(1) -= isb(3);
                isb(2) -= isb(3);
                isb(3) = isb(4);
            }
            for (size_t i = 1; i <= 3; i++) {
                sb(i) = isb(1) * props.cvec(i, 1) + isb(2) * props.cvec(i, 2) + isb(3) * props.cvec(i, 3);
            }
            sbnor = (tnorm(sb));
            for (size_t j = 1; j <= 3; j++) {
                sb(j) = sb(j) / sbnor;
                if (abs(sb(j)) < 1.e-03) sb(j) = 0.;
            }

            prod = sn(1) * sb(1) + sn(2) * sb(2) + sn(3) * sb(3);
            if (prod >= 1.e-3) {
                printf("system is not orthogonal\n");
                exit(1);
            }

            // end of if(ioption == 2)
        }

    }

    // ----------------------------------------------------------------------
    // --- generates all symmetry related vectors sneq(i,n) with z>0.
    // --- eliminates redundant poles: coincidents and opposites
    // ----------------------------------------------------------------------

    if (ioption == 3) {

        for (size_t n = 1; n <= props.nsymop; n++) {
            itag(n) = 0;
            for (size_t i = 1; i <= 3; i++) {
                sneq(i, n) = 0.0;
                for (size_t j = 1; j <= 3; j++) {
                    sneq(i, n) += props.h(n, i, j) * sn(j);
                }
            }
        }

        for (size_t m = 1; m <= props.nsymop - 1; m++) {
            if (itag(m) == 0) {
                for (size_t n = m + 1; n <= props.nsymop; n++) {
                    snpro = sneq(1, m) * sneq(1, n) + sneq(2, m) * sneq(2, n) + sneq(3, m) * sneq(3, n);
                    if (snpro <= 1.001 && snpro >= 0.999) itag(n) = 1;
                    if (snpro >= -1.001 && snpro <= -0.999) itag(n) = 1;
                }
            }
        }

        npol = 0;
        for (size_t n = 1; n <= props.nsymop; n++) {
            if (itag(n) == 0) {
                npol++;
                isign = 1;
                if (sneq(3, n) < 0.) isign = -1;
                sneq(1, npol) = isign * sneq(1, n);
                sneq(2, npol) = isign * sneq(2, n);
                sneq(3, npol) = isign * sneq(3, n);
            }
        }

        //end of ioption=3
    }



    // **********************************************************************
    return;
}

// **********************************************************************
//     subroutine data_crystal        --->      version 04/sep/2006
// **********************************************************************
KOKKOS_FUNCTION void data_crystal(size_t iph, std::string filename, PropType& props, StateVType& statv)
{
    //x
    std::string prosa;
    //x
    fmat1(double, sn, 3); fmat1(double, sb, 3); fmat1(double, aux5, 5); fmat2(double, aux33, 3, 3); fmat2(double, aux55, 5, 5); fmat4(double, aux3333, 3, 3, 3, 3); fmat2(double, aux324, 3, 24); fmat1(size_t, mode, 20); fmat2(double, hard_mode, 20, 20);
    // MATARIX <double>  sn(3), sb(3);
    //
    // MATARIX <double>  aux5(5), aux33(3, 3), aux55(5, 5), aux3333(3, 3, 3, 3);
    // MATARIX <double>  aux324(3, 24);
    //
    // MATARIX <size_t>  mode(20);
    // MATARIX <double>  hard_mode(20, 20);

    std::fstream  file;

    double twshx, thres1x, thres2x, tau0x, tau1x, thet0x, thet1x, dummy;
    double tau0_mode_a, tau0_mode_b, tau0_mode_c;
    size_t icrysym, npoles, i, j, nmodesx, kount, nsmx, nrsx, isensex, isectwx, nslsys, iclosepos, icloseneg, modex, nsysx, jount;
    //z
    for (size_t i = 1; i <= 20; i++) {
        mode(i) = 0;
    }

    // *** reads crystal symmetry and unit cell parameters.
    // *** generates all symmetry operations associated with crysym.



    file.std::fstream::open(filename, std::fstream::in);

    //q     crystal_symmetry (1,ur1,icrysym,sn,aux3,sb,npoles)
    crystal_symmetry(1, file, props, icrysym, sn, aux324, sb, npoles);


    // *** reads single crystal elastic stiffness
    nextline(file, 1);

    for (size_t i = 1; i <= 6; i++) {
        for (size_t j = 1; j <= 6; j++) {
            file >> props.c2ca(iph, i, j);
        }
        nextline(file, 1);
    }

    nextline(file, 1);

    // *** reads single crystal thermal expansion coefficients
    nextline(file, 1);


    // *** reads information about slip and twinning systems
    nextline(file, 1);

    file >> nmodesx; nextline(file, 1);

    file >> props.nmodes(iph); nextline(file, 1);
    ////std::cout  << props.nmodes(iph) << std::endl;

    //v      props.nmodes(iph)=nmodes(iph)
    for (size_t i = 1; i <= props.nmodes(iph); i++) {
        file >> mode(i);
    }

    nextline(file, 1);

    props.ntwmod(iph) = 0;
    props.nsyst(iph) = 0;
    props.ntwsys(iph) = 0;
    kount = 1;

    // *** reads deformation modes and associated parameters from filecrys

    for (size_t nm = 1; nm <= nmodesx; nm++) {

        // skip header
        nextline(file, 1);

        file >> modex >> nsmx >> nrsx >> isensex; nextline(file, 1);



        if (modex != nm) {
            printf("mode numbers must be sequential in crystal file\n");

            exit(1);
        }

        // *** skips non-active mode if it is not in the list.
        //printf("modex, mode(kount), %zu, %zu\n", modex, mode(kount));
        if (modex == mode(kount)) {

            if (props.ipr == 1) printf("MODE NUMBER : %zu; SLIP SYSTEM NUMBER : %zu\n", modex, nsmx);

            //x        if(isensex == 0) ics=0      // ics=0 --> non-centro-symm scys

            file >> twshx >> isectwx >> thres1x >> thres2x; nextline(file, 1);

            file >> tau0x >> tau1x >> thet0x >> thet1x; nextline(file, 1);

            file >> tau0_mode_a >> tau0_mode_b >> tau0_mode_c >> statv.temp_ini >> statv.temp_fact; nextline(file, 1);

            jount = 0;
            for (size_t jm = 1; jm <= nmodesx; jm++) {

                file >> dummy;//props.hard(iph, kount, jm);
                //props.hard(iph, jm, jm) = 1.0;

                if (ismember(jm, mode)) {
                    jount++;
                    hard_mode(kount, jount) = dummy;
                    //printf("latent hard %f\n", dummy);
                }
            }

            nextline(file, 1);

            //WRITEOUT(file)
            //hlatez -> hard

            // *** checks whether voce parameters are kosher:
            //       tau0>0 , tau1 >= 0 , thet0 >= thet1 >= 0
            //       tau1=0   corresponds to linear hardening.
            //       theta0=0 forces no-hardening.
            // *** if voce parameters are non-kosher checks for ill-posed hardening.

            //if(props.ihardlaw != 1) data_crystal_voce(kount,iph,tau0x,tau1x,thet0x,thet1x);

            props.nsm(iph, kount) = nsmx;

            // *** verification of twinning data to be sure program will run properly

            if (twshx == 0. && props.ntwmod(iph) != 0) {
                printf("twinning modes must follow slip modes\n");

                exit(1);
            }

            if (twshx != 0.) {
                //v        ntwmod(iph)=ntwmod(iph)+1
                props.ntwmod(iph)++;

                //x        if(ntwmod(iph) > ntwmmx) {
                //x         write(*,'('' ntwmod in phase'',i3,'' is'',i3)') iph,ntwmod(ip
                //x         write(*,'('' change parameter ntwmmx in vpsc6.dim'')')
                //x         stop
                //x        }

                size_t ntwmod = props.ntwmod(iph);

                props.twsh(iph, ntwmod) = twshx;
                props.twthres(iph, ntwmod, 1) = thres1x;
                props.twthres(iph, ntwmod, 2) = thres2x;

            }

            for (size_t js = 1; js <= props.nsm(iph, kount); js++) {

                //v        nsyst(iph)=nsyst(iph)+1
                props.nsyst(iph)++;
                nsysx = props.nsyst(iph);

                //printf("NSYST : %zu\n", nsysx);



                if (twshx != 0.) {
                    //v         ntwsys(iph)=ntwsys(iph)+1
                    props.ntwsys(iph)++;
                }

                //   initializes parameters associated with each system in the mode.

                props.nrs(iph, nsysx) = nrsx;
                props.isense(iph, nsysx) = isensex;
                props.isectw(iph, nsysx) = isectwx;
                props.tau(iph, nsysx, 1) = tau0x;
                props.tau(iph, nsysx, 2) = tau1x;
                props.tau0_mode(iph, nsysx, 1) = tau0_mode_a;
                props.tau0_mode(iph, nsysx, 2) = tau0_mode_b;
                props.tau0_mode(iph, nsysx, 3) = tau0_mode_c;
                props.thet(iph, nsysx, 1) = thet0x;
                props.thet(iph, nsysx, 2) = thet1x;

                // *** calculates cartesian components of slip and normal vectors

                //q       crystal_symmetry (2,ur1,icrysym,sn,aux3,sb,npoles)

                crystal_symmetry(2, file, props, icrysym, sn, aux324, sb, npoles);
                //char test[50]; for(size_t ii = 0; ii < 50; ii++) file.get(test[ii]);
                for (size_t j = 1; j <= 3; j++) {
                    props.dnca(iph, nsysx, j) = sn(j);
                    props.dbca(iph, nsysx, j) = sb(j);
                }



                // ***  defines schmid vector in crystal axes for each system

                for (size_t i = 1; i <= 3; i++) {
                    for (size_t j = 1; j <= 3; j++) {
                        aux33(i, j) = (props.dnca(iph, nsysx, i) * props.dbca(iph, nsysx, j) + props.dnca(iph, nsysx, j) * props.dbca(iph, nsysx, i)) / 2.;

                        //printf("%26.16E\n", aux33(i, j));
                    }
                }

                chg_basis(&aux5(1), &aux33(1), 2, 5);

                //fullout(aux5);

                for (size_t i = 1; i <= 5; i++) {
                    props.schca(iph, nsysx, i) = aux5(i);
                }

                // end of loop over deformation modes
            }


            kount++;
            // end of loop over all modes in phase 'iph'
        }
        else {

            nextline(file, nsmx + 4);

        }



    }

    // assign hardening for each slip mode in phase
    i = 0;
    for (size_t im = 1; im <= props.nmodes(iph); im++) {
        for (size_t is = 1; is <= props.nsm(iph, im); is++) {
            i++;
            j = 0;
            for (size_t jm = 1; jm <= props.nmodes(iph); jm++) {
                for (size_t js = 1; js <= props.nsm(iph, jm); js++) {
                    j++;
                    props.hard(iph, i, j) = hard_mode(im, jm);
                }
            }
            props.hard(iph, i, i) = 1.0;
        }
    }
    // *** checks whether the single crystal yield surface is open

    //printf("IPH : %zu\nNSYST : %f\nNTWS : %f\n", iph, props.nsyst(iph), props.ntwsys(iph));

    nslsys = props.nsyst(iph) - props.ntwsys(iph);
    for (size_t icomp = 1; icomp <= 5; icomp++) {
        iclosepos = 0;
        icloseneg = 0;
        for (size_t ns = 1; ns <= nslsys; ns++) {
            //z          if(abs(schca(icomp,ns,iph)) > 1.e-3) {

            //printf("%f, ", props.schca(iph, ns, icomp));

            if (abs(props.schca(iph, ns, icomp)) > 1.e-3) {
                iclosepos = 1;
                icloseneg = 1;
            }
        }
        //std::cout << std::endl;
        if (props.ntwsys(iph) != 0) {
            for (size_t ns = nslsys + 1; ns <= props.nsyst(iph); ns++) {
                //z            if(schca(icomp,ns,iph) >  1.e-3) iclosepos=1
                //z            if(schca(icomp,ns,iph) < -1.e-3) icloseneg=1
                if (props.schca(iph, ns, icomp) > 1.e-3) iclosepos = 1;
                if (props.schca(iph, ns, icomp) < -1.e-3) icloseneg = 1;
            }
        }
        if (iclosepos != 1 || icloseneg != 1) {
            ////$          write(*,'('' warning // the scys is open for phase'',i5, ////$                    '' along direction'',i5)') iph,icomp
            //z          pause
            printf("the scys is open\n");
            //exit(1);
        }
    }

    // *** initialize self latent hardening coefs for each system of the ph
    //     absolute units are accounted for by modulating factor in hardening

    file.close();

    return;
}


KOKKOS_FUNCTION void user_mat_init(double* statev1D)
{
    fmat1(double, aux5, 5); fmat2(double, aux55, 5, 5); fmat1(double, aux3, 3); fmat2(double, aux33, 3, 3); fmat4(double, aux3333, 3, 3, 3, 3);
    //MATARIX <double> aux5(5), aux55(5, 5), aux3(3), aux33(3, 3), aux3333(3, 3, 3, 3);
    size_t ncompa;

    double wphtot;

    std::fstream vpsc7in;

    std::string line;

    std::string filetex, filecrys;

    PropType props;
    StateVType statv;
    IntVType intv;

    size_t cnt = 0;
    props.init(&statev1D[0], cnt);
    statv.init(&statev1D[0], cnt);
    intv.init(&statev1D[0], cnt);

    //eye(props.xid3, 3, 1.0);
    //eye(props.xid5, 5, 1.0);


    //init_vpsc(statv, props);

    //CHANGEME
    props.cv = 140.0;


    vpsc7in.std::fstream::open("vpsc7.txt", std::fstream::in);

    vpsc7in >> props.nph;

    nextline(vpsc7in, 1);

    for (size_t i = 1; i <= props.nph; i++) {
        vpsc7in >> statv.wph(i);
    }

    nextline(vpsc7in, 1);



    props.iphbot = 1;
    props.iphtop = props.nph;
    props.ngr(1) = 0.0;
    props.ishape(1) = 0.0;
    statv.iflat(1) = 0.0;

    wphtot = 0.0;

    for (size_t iph = 1; iph <= props.nph; iph++) wphtot += statv.wph(iph);

    for (size_t iph = 1; iph <= props.nph; iph++) statv.wph(iph) /= wphtot;

    for (size_t i = 1; i <= 3; i++) {
        for (size_t j = 1; j <= 3; j++) {
            statv.fijph(1, i, j) = xid(i, j);
        }
    }

    update_shape(0, statv, props.interaction);

    for (size_t iph = 1; iph <= props.nph; iph++) {

        nextline(vpsc7in, 1);


        vpsc7in >> statv.axisph(iph + 1, 1, 1) >> statv.axisph(iph + 1, 1, 2) >> statv.axisph(iph + 1, 1, 3);

        nextline(vpsc7in, 1);

        vpsc7in >> statv.eulerph(iph + 1, 1) >> statv.eulerph(iph + 1, 2) >> statv.eulerph(iph + 1, 3);

        nextline(vpsc7in, 2);

        vpsc7in >> filetex;

        nextline(vpsc7in, 2);

        vpsc7in >> filecrys;

        nextline(vpsc7in, 1);

        if (props.ipr == 1) printf("Texture File : %s\nCrystal Properties : %s\n", filetex.c_str(), filecrys.c_str());

        data_crystal(iph, filecrys, props, statv);

        data_texture(iph, filetex, props, statv);

        update_shape(iph, statv, props.interaction);

        //WARNING DANGER AHEAD
        props.ishape(iph + 1) = 0.0;
        //LEAVING DANGER
    }


    // **********************************************************************
    // *** reads settings for convergence procedures.

    nextline(vpsc7in, 1);
    vpsc7in >> props.errs >> props.errd >> props.errm >> props.errso; nextline(vpsc7in, 1);


    vpsc7in >> props.itmaxext >> props.itmaxint >> props.itmaxso; nextline(vpsc7in, 1);

    //printf("props.itmaxext = %zu\n", props.itmaxext);

    props.irsvar = 0;

    vpsc7in >> props.ibcinv;
    //v      props.ibcinv=ibcinv


    // **********************************************************************
    // *** reads input/output settings for the run.

    //v      irecover=0
    //v      isave=0
    props.irecover = 0;
    props.isave = 0;

    props.icauchy = 0;

    // **********************************************************************
    // *** reads modeling conditions.

    nextline(vpsc7in, 2);
    //char test[50]; for(size_t ii = 0; ii < 50; ii++) vpsc7in.get(test[ii]);
    vpsc7in >> props.ihardlaw; nextline(vpsc7in, 1);

    vpsc7in >> props.iratesens; nextline(vpsc7in, 1);

    vpsc7in >> props.interaction; nextline(vpsc7in, 1);
    vpsc7in >> props.iupdori >> props.iupdshp >> props.iupdhar >> props.itemphard; nextline(vpsc7in, 1);
    vpsc7in >> props.ipr; nextline(vpsc7in, 1);
    vpsc7in >> props.stranvm_pl >> props.dstranvmth_elpl >> props.dstranvmth_max >> props.dstranvmth_min; nextline(vpsc7in, 1);
    vpsc7in >> props.kinctex; nextline(vpsc7in, 1);
    vpsc7in >> props.neltex; nextline(vpsc7in, 1);

    //props.eltex = MATARIX <size_t> (props.neltex);

    for (size_t i = 1; i <= props.neltex; i++)
    {
        vpsc7in >> props.eltex(i);
    }

    // closes sx file after reading hardening

    vpsc7in.close();

    props.ipru = 789;

    //v      iflu=0
    //v      props.iflu=0

    // *** checks if rate sensitivity is the same for all systems.
    // *** search for nrsmin (needed to get taumax inside nr subroutine)

    props.nunique = 1;
    ncompa = props.nrs(1, 1);
    props.nrsmin = props.nrs(1, 1);
    for (size_t iph = 1; iph <= props.nph; iph++) {
        for (size_t is = 1; is <= props.nsyst(iph); is++) {
            if (props.nrs(iph, is) != ncompa) props.nunique = 0;
            if (props.nrs(iph, is) < props.nrsmin) props.nrsmin = props.nrs(iph, is);
        }
    }

    if (props.iupdori == 0 && props.iupdshp == 1) {
        printf("error in vpsc_input\n"); exit(1);
    }

    // *******************************************************************
    // *** initialize arrays associated with grains:
    //     hardening, accumulated shear, power, twinning parameters, etc


    statv.jran = 0;
    statv.gamd0g = 1.0;


    // *** initialize crss in each system of each grain.
    //     for 'voce' the parameters have been read from sx file.
    //     for 'mts' read hardening parameters here (option=0). will also
    //     initialize the structure/temp/rate independent array 'taue'.

    //v      if(ihardlaw == 0) call update_crss_voce (1,statev,props)
    if (props.ihardlaw == 0)  init_crss_voce(1, statv, props);
    if (props.itemphard == 1) init_crss_temp(statv, props);


    // *** initialize gauss-legendre coordinates and arrays used in the doubl
    // *** integrals giving the eshelby tensors

    props.grbot = 1;
    props.grtop = props.ngr(props.nph + 1);

    for (size_t iph = 1; iph <= props.nph; iph++) {
        for (size_t kkk = props.ngr(iph) + 1; kkk <= props.ngr(iph + 1); kkk++) {
            statv.giph(kkk) = iph;
        }
    }

    eshelby_init(props);

    for (size_t i = 1; i <= 9; i++) {
       statv.stressold(i) = 0.0;
    }
    statv.temp = statv.temp_ini;
    statv.tempold = statv.temp;
	statv.kinc = 0;


    props.write(&statev1D[0]);
    statv.write(&statev1D[0]);

    return;
}


KOKKOS_FUNCTION size_t user_mat_number_vars() {

    size_t SZ_STATV = NSTATV_;

    return SZ_STATV;
}

// ---- END OF HOST CALCULATIONS ----//
//
//KOKKOS_FUNCTION void user_mat_model(double* statev1D,
//    double* vel_grad, double dtime, size_t kinc, size_t k_cell, double den, double ie,
//    double* stress1D, double& sspd, double& temp, double& pressure) {
//
//
//
//    fmat2(double, ddsdde, 6, 6);
//
//
//
//    //Allocate change in strain and rotation
//    fmat2(double, deps, 3, 3);
//    fmat2(double, drot, 3, 3);
//
//    //Allocate asymmetric rotation tensor
//    fmat2(double, dW, 3, 3);
//
//    //pressure
//    double p;
//
//    // create the property classes 
//    PropType props;
//    StateVType statv;
//    IntVType intv;
//
//    //Interpret velocity gradient and stress 
//    MATARIX <double> L(&vel_grad[0], 3, 3),
//        stress(&stress1D[0], 3, 3);
//
//
//    //check for nans in vel_grad
//    //if (L(1, 1) != L(1, 1)) {
//    //	printf("L : %16.8f, %16.8f, %16.8f, %16.8f, %16.8f, %16.8f;\n", L(1, 1), L(2, 2), L(3, 3), L(1, 2), L(1, 3), L(2, 3));
//    //	exit(1);
//    //}
//
//
//    if (tnorm(L) > 1e-16) {
//
//
//        for (size_t i = 1; i <= 3; i++) {
//            for (size_t j = 1; j <= 3; j++) {
//                deps(i, j) = 0.5 * (L(i, j) + L(j, i)) * dtime;
//                dW(i, j) = 0.5 * (L(i, j) - L(j, i)) * dtime;
//            }
//        }
//
//        //for (size_t i = 1; i <= 3; i++) {
//        //	for (size_t j = i+1; j <= 3; j++) {
//        //		deps(i, j) *= 2.0;
//        //		deps(j, i) *= 2.0;
//        //	}
//        //}
//
//        //Calculate rotation matrix
//        rodrigues(dW, drot);
//
//        //interpret the state variable array
//        size_t cnt = 0;
//        size_t NSTATV = user_mat_number_vars();
//
//
//        props.init(&statev1D[0], cnt);
//        statv.init(&statev1D[0], cnt);
//        intv.init(&statev1D[0], cnt);
//
//        if (cnt != NSTATV) {
//            printf("NSTATV wrong size; NSTATV = %zu; cnt = %zu;\n", NSTATV, cnt);
//        }
//
//
//        //supply statv with defined vars
//        // update the temperature
//        update_temperature(statv.svm, statv.evmp, statv.evmpold, statv.temp, statv.dtemp, statv.wrplastic, statv.wplastic, props.b_int);
//
//        // temp = statv.temp;
//
//
//        // solve the problem
//        double dtimevpsc = dtime * 1.0e3;
//        evol_vpsc(statv, props, intv, deps, drot, stress, kinc, dtimevpsc, k_cell, ddsdde);
//
//
//
//
//        statv.tempold = statv.temp;
//
//        // write the scalars to the state variable array
//        props.write(&statev1D[0]);
//        statv.write(&statev1D[0]);
//
//        //printf("stress : %16.8f, %16.8f, %16.8f\n", stress(1, 1), stress(2, 2), stress(3, 3));
//
//        if (stress(1, 1) != stress(1, 1)) {
//            printf("deps(1, 1) : %16.4e; stress(1, 1) : %16.8f;\n", deps(1, 1), stress(1, 1));
//            exit(1);
//        }
//
//
//    }
//
//    // gruneisen
//
//    // should be in inputs.cpp
//    //double cv       = 1.0;         // specific heat
//    //double csmin    = 0.01;        // minimum sound speed
//    //double C0       = 0.394;       //
//    //double s1       = 1.49;        //
//    //double g0       = 2.0;         //
//    //double b        = 0.0;         //
//    //double den_ref  = 8.94;
//
//
//    double cv       = 0.14e3;      // specific heat
//    double csmin    = 1000;        // minimum sound speed
//    double C0       = 3400.0;      // bulk speed of sound
//    double s1       = 1.49;        //
//    double g0       = 2.0;         //
//    double b        = 0.0;         //
//    double den_ref  = 16.640e-3;
//
//    
//    
//    double b1       = 1.3333;      // linear slope of UsUp for Riemann solver
//    
//    
//
//    // end inputs
//
//    double dPde = 0.0;
//    double dPdr = 0.0;
//    double P = 0.0;
//
//    double mu = (den/den_ref) - 1.0;
//
//    if (mu < 0.0) {
//        dPde = g0 * den;
//        dPdr = C0 * C0 + g0 * ie;
//        P = den_ref * C0 * C0 * mu + g0 * den * ie;
//    }
//    else {
//        dPde = (g0 + b * mu) * den;
//
//        // two helper variables
//        double det = (1.0 - (s1 - 1.0) * mu) * (1.0 - (s1 - 1.0) * mu);
//        double fa = (1.0 + (1.0 - 0.5 * g0) * mu - 0.5 * b * mu * mu) * mu;
//
//        dPdr = C0 * C0 * ((1.0 + (1.0 - 0.5 * g0) * 2.0 * mu - 0.5 * b * 3.0 * mu * mu) * det -
//            fa * 2.0 * (s1 - 1.0) * (mu * (s1 - 1.0) - 1.0)) / (det * det) +
//            ie * ((2.0 * mu + 1.0) * b + g0);
//
//        P = den_ref * C0 * C0 * fa / det + (g0 + b * mu) * den * ie;
//    } // end if
//
//
//    // calculate sound speed
//    sspd = sqrt(fmax(P * dPde / (den * den) + dPdr, 0.0));
//
//    //printf("sspd = %.0f\n", sspd);
//    if (sspd < csmin) sspd = csmin;
//
//    //sspd = 3400.0e0;
//    pressure = 0.0;
//    return;
//}

