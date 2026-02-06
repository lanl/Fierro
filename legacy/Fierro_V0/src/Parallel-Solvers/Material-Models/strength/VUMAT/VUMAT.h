#pragma once

#include "Simulation_Parameters/Material.h"
#include "matar.h"
using namespace mtr;


/* VUMAT */
class VUMAT {
public:
  
  VUMAT(
    const DCArrayKokkos <material_t> &material,
    const DCArrayKokkos <double> &eos_state_vars,
    const DCArrayKokkos <double> &strength_state_vars,
    const DCArrayKokkos <double> &eos_global_vars,
    const DCArrayKokkos <double> &strength_global_vars,
    const DCArrayKokkos <double> &elem_user_output_vars,
    const size_t mat_id,
    const size_t elem_gid);

  
  // ~VUMAT();

  int calc_stress(
    const DViewCArrayKokkos <double> &elem_pres,
    const DViewCArrayKokkos <double> &elem_stress,
    const size_t elem_gid,
    const size_t mat_id,
    const DCArrayKokkos <double> &eos_state_vars,
    const DCArrayKokkos <double> &strength_state_vars,
    const DCArrayKokkos <double> &eos_global_vars,
    const DCArrayKokkos <double> &strength_global_vars,
    const DCArrayKokkos <double> &elem_user_output_vars,
    const DViewCArrayKokkos <double> &elem_sspd,
    const double den,
    const double sie,
    const ViewCArrayKokkos <double> &vel_grad,
    const ViewCArrayKokkos <size_t> &elem_node_gids,
    const DViewCArrayKokkos <double> &node_coords,
    const DViewCArrayKokkos <double> &node_vel,
    const double vol,
    const double dt, 
    const double rk_alpha,
    const size_t cycle,
    const size_t rk_level,
    const double time);


private:

  // Abaqus VUMAT parameters
  int nblock;
  int ndir;
  int nshr;
  int nstatev;
  int nfieldv;
  int nprops;
  int lanneal;
  double stepTime = 0.0; // reduction is done on this
  double totalTime = 0.0; // reduction is done on this
  double dt;
  char cmname[80]; // fortran is expecting 80 char length

  // Abaqus VUMAT arrays parameters
  FMatrix <double> coordMp;
  FMatrix <double> charLength;
  ViewFMatrix <double> props;
  FMatrix <double> density;
  FMatrix <double> strainInc;
  FMatrix <double> relSpinInc;
  FMatrix <double> tempOld;
  FMatrix <double> stretchOld;
  FMatrix <double> defgradOld;
  FMatrix <double> fieldOld;
  FMatrix <double> stressOld;
  FMatrix <double> stateOld;
  FMatrix <double> enerInternOld;
  FMatrix <double> enerInelasOld;
  FMatrix <double> tempNew;
  FMatrix <double> stretchNew;
  FMatrix <double> defgradNew;
  FMatrix <double> fieldNew;
  FMatrix <double> stressNew;
  ViewFMatrix <double> stateNew;
  FMatrix <double> enerInternNew;
  FMatrix <double> enerInelasNew;

  // Other parameters for intermediate calculations
  int ind_[6*2] = {1,2,3,1,2,3,1,2,3,2,3,1};
  ViewFMatrixKokkos<int> ind {ind_,6,2};

  FMatrix <double> Fvel_grad {3,3}; // vel_grad in Fortran array layout
  FMatrix <double> D {3,3}; // symmetric part (D) of the velocity gradient (L)
  FMatrix <double> W {3,3}; // skew symmetric part (W) of the velocity gradient (L)

  FMatrix <double> F_old {3,3}; // old deformation gradient
  FMatrix <double> R_old {3,3}; // old rotation matrix
  FMatrix <double> U_old {3,3}; // old right stretch tensor

  FMatrix <double> F_new {3,3}; // new deformation gradient
  FMatrix <double> R_new {3,3}; // new rotation matrix
  FMatrix <double> U_new {3,3}; // new right stretch tensor


};


