#include "matar.h"
using namespace mtr;
using real_t = double;

// class ChgBasis
// {
//   private:
//     const real_t SQR2_ = 1.41421356237309;
//     const real_t RSQ2_ = 0.70710678118654744;
//     const real_t RSQ3_ = 0.57735026918962584;
//     const real_t RSQ6_ = 0.40824829046386304;

//     DCMatrixKokkos <double> B_basis_;
        
//   public:
//     ChgBasis();

//     void set_b_basis();

// using namespace to avoid collision of functions that are the same in other models like EVPFFT or LS-EVPFFT
namespace EVP
{

    KOKKOS_FUNCTION
    void chg_basis_1(real_t *CE2_, real_t *C2_, int IOPT, int KDIM, real_t *B_basis_);//ptr) const;

    KOKKOS_FUNCTION
    void chg_basis_2(real_t *CE2_, real_t *C2_, int IOPT, int KDIM, real_t *B_basis_);//ptr) const;

    KOKKOS_FUNCTION
    void chg_basis_3(real_t *CE4_, real_t *C4_, int IOPT, int KDIM, real_t *B_basis_);//ptr) const;

    KOKKOS_FUNCTION
    void chg_basis_4(real_t *CE4_, real_t *C4_, int IOPT, int KDIM, real_t *B_basis_);//ptr) const;

//     KOKKOS_FUNCTION
//     real_t* B_basis_device_pointer() const;

//     real_t* B_basis_host_pointer() const;

// };


KOKKOS_FUNCTION
void inverse_gj(real_t *a_, int n);

KOKKOS_FUNCTION
void lu_inverse(real_t *a_, int n);
KOKKOS_FUNCTION
void ludcmp(real_t *a_, int n, int np, int *indx_, real_t d, int &isingular);
KOKKOS_FUNCTION
void lubksb(real_t *a_, int n, int np, int *indx_, real_t *b_);

KOKKOS_FUNCTION
void euler(int iopt, real_t &ph, real_t &th, real_t &tm, real_t *a_);

KOKKOS_FUNCTION
void update_strain(real_t *strain_, real_t *strain_n_, real_t *rate_, const real_t dt);

KOKKOS_FUNCTION
void update_stiff(real_t *cg66_, real_t *cc66_, real_t *aa_, real_t *B_basis_);

KOKKOS_FUNCTION
void update_schmid(real_t *sch_, real_t *schca_, real_t *aa_, const int nsm, real_t *B_basis_);

KOKKOS_FUNCTION
void evpal(real_t *stress_, real_t *edotp_, real_t *gamdot_, real_t *stress_n_, real_t *strain_,
	real_t* ept_n_, real_t *cg66_, real_t *sc_, real_t *crss_, const real_t gamd0, const real_t nrs, 
	const int nsmx, const real_t dt, const size_t cycle, real_t *B_basis_, const size_t elnum); //const real_t gacumgr, real_t *voceParam, real_t *hard, 


KOKKOS_FUNCTION
void update_orient(real_t *ag_, real_t *vel_grad_, real_t *gamdot_, real_t *dnca_, real_t *dbca_, const int nsm, const real_t dt);

KOKKOS_FUNCTION
void orient(real_t *a_, real_t *c_);

KOKKOS_FUNCTION
void harden(real_t *crss_, real_t &gacumgr, real_t *gamdot_, real_t *voceParam_, real_t *hard_, const int nsm, const real_t dt);

} // end namespace EVP