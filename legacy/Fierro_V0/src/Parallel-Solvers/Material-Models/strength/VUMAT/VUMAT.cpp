#include "VUMAT.h"
#include "vumat_function.h"


// Below are the list of variables calculated or that have values.
// Variables which are "dummy" are not currently calculated yet. If your vumat
// needs these dummy variables please calculated them in calc_stress before the 
// call to vumat_ and update this list;

//    nblock = 1
//    ndir = 3
//    nshr = 3
//    nstatev = strength_state_vars.dims(1) // from fierro yaml
//    nfieldv = 0; // Fierro does not currently support user-defined external field variables
//    nprops = strength_global_vars.dims(1); // from fierro yaml
//    lanneal - 0; // Fierro does not currently support annealing process
//    stepTime - calculated in calc_stress before calling vumat_
//    totalTime - calculated in calc_stress before calling vumat_
//    dt - calculated in calc_stress before calling vumat_
//    cmname - string(mat_id) // mat_id from fierro
//    coordMp - dummy
//    charLength - dummy
//    props - strength_global_vars.host(mat_id,0); // from fierro yaml
//    density - dummy
//    strainInc - calculated in calc_stress before calling vumat_
//    relSpinInc - calculated in calc_stress before calling vumat_
//    tempOld - dummy
//    stretchOld - dummy
//    defgradOld - dummy
//    fieldOld - dummy
//    stressOld - assigned in calc_stress before calling vumat_
//    stateOld - assigned in calc_stress before calling vumat_
//    enerInternOld - dummy
//    enerInelasOld - dummy
//    tempNew - dummy
//    stretchNew - dummy
//    defgradNew - dummy
//    fieldNew - dummy
//    stressNew - updated stress after call to vumat_
//    stateNew - updated state variables after call to vumat_
//    enerInternNew - dummy
//    enerInelasNew - dummy

VUMAT::VUMAT(
  const DCArrayKokkos <material_t> &material,
  const DCArrayKokkos <double> &eos_state_vars,
  const DCArrayKokkos <double> &strength_state_vars,
  const DCArrayKokkos <double> &eos_global_vars,
  const DCArrayKokkos <double> &strength_global_vars,
  const DCArrayKokkos <double> &elem_user_output_vars,
  const size_t mat_id,
  const size_t elem_gid)
{

  auto run_loc = material.host(mat_id).strength_run_location;
  if (run_loc == RUN_LOCATION::device) {
    throw std::runtime_error("VUMAT only runs on Host");
  }

  // check that the user remembered to enter number of nstatev and nprops in the yaml
  if (strength_state_vars.dims(1) == 0) {
    throw std::runtime_error("VUMAT: num_strength_state_vars must be specified = abaqus nstatev");
  }

  nblock = 1;  // num of gauss integration points in the element
  ndir = 3;    // for 3D (ndir: direct stress components)
  nshr = 3;    // for 3D (nshr: shear stress components)
  nstatev = strength_state_vars.dims(1);
  nfieldv = 0; // Fierro does not currently support user-defined external field variables
  nprops = strength_global_vars.dims(1);
  lanneal = 0; // Fierro does not currently support annealing process

  // This long process in needed to pass string to fortran
  std::string mat_id_str = std::to_string(mat_id);
  std::strncpy(cmname, mat_id_str.c_str(), mat_id_str.size());
  cmname[79] = '\0'; // Ensure null termination

  // Allocate Arrays
  coordMp = FMatrix <double> (nblock, ndir);
  charLength = FMatrix <double> (nblock);
  props = ViewFMatrix <double> (&strength_global_vars.host(mat_id,0), nprops);
  density = FMatrix <double> (nblock);
  strainInc = FMatrix <double> (nblock, ndir+nshr);
  relSpinInc = FMatrix <double> (nblock, nshr);
  tempOld = FMatrix <double> (nblock);
  stretchOld = FMatrix <double> (nblock, ndir+nshr);
  defgradOld = FMatrix <double> (nblock,ndir+2*nshr);
  fieldOld = FMatrix <double> (nblock, nfieldv);
  stressOld = FMatrix <double> (nblock, ndir+nshr);
  stateOld = FMatrix <double> (nblock, nstatev);
  enerInternOld = FMatrix <double> (nblock);
  enerInelasOld = FMatrix <double> (nblock);
  tempNew = FMatrix <double> (nblock);
  stretchNew = FMatrix <double> (nblock, ndir+nshr);
  defgradNew = FMatrix <double> (nblock,ndir+2*nshr);
  fieldNew = FMatrix <double> (nblock, nfieldv);
  stressNew = FMatrix <double> (nblock, ndir+nshr);
  stateNew = ViewFMatrix <double> (&strength_state_vars.host(elem_gid,0), nblock, nstatev);
  enerInternNew = FMatrix <double> (nblock);
  enerInelasNew = FMatrix <double> (nblock);

}

// VUMAT::~VUMAT()
// {}

int VUMAT::calc_stress(
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
  const double dt_in,
  const double rk_alpha,
  const size_t cycle,
  const size_t rk_level,
  const double time)
{

  //set dt
  dt = dt_in;

  // slice stress view from elem_stress
  ViewFMatrix<double> stress_view (&elem_stress.host(rk_level,elem_gid,0,0), 3, 3);

  // stateOld:
  for (int i = 1; i <= nblock; i++) {
    for (int j = 1; j <= nstatev; j++) {
      stateOld(i,j) = stateNew(i,j);
    }
  }

  // write vel_grad to Fortran array layout.
  for (size_t i = 1; i <= 3; i++) {
    for (size_t j = 1; j <= 3; j++) {
      Fvel_grad(i,j) = vel_grad(i-1, j-1);
    }
  }

  // decompose vel_grad
  for(size_t i = 1; i <= 3; i++) {
    for(size_t j = 1; j <= 3; j++) {
      D(i,j) = 0.5*(Fvel_grad(i,j) + Fvel_grad(j,i));
      W(i,j) = 0.5*(Fvel_grad(i,j) - Fvel_grad(j,i));
    }   
  }

  // calculate strainInc
  for (int i = 1; i <= (ndir+nshr); i++) {
    strainInc(1,i) = D(ind(i,1), ind(i,2)) * dt;
  }

  // calculate relSpinInc
  relSpinInc(1,1) = W(3,2) * dt;
  relSpinInc(1,2) = W(1,3) * dt;
  relSpinInc(1,3) = W(2,1) * dt;

  // stressOld:
  
  for (int i = 1; i <= (ndir+nshr); i++) {
    stressOld(1,i) = stress_view(ind(i,1), ind(i,2));
  }

  //---------------------------------------------------------------------------------------------------
  //--------------------------------------------call vumat---------------------------------------------
  //---------------------------------------------------------------------------------------------------
  vumat_(
        &nblock, &ndir, &nshr, &nstatev, &nfieldv, &nprops, &lanneal,
        &stepTime, &totalTime, &dt, cmname, coordMp.pointer(), charLength.pointer(),
        props.pointer(), density.pointer(), strainInc.pointer(), relSpinInc.pointer(),
        tempOld.pointer(), stretchOld.pointer(), defgradOld.pointer(), fieldOld.pointer(),
        stressOld.pointer(), stateOld.pointer(), enerInternOld.pointer(), enerInelasOld.pointer(),
        tempNew.pointer(), stretchNew.pointer(), defgradNew.pointer(), fieldNew.pointer(),
        stressNew.pointer(), stateNew.pointer(), enerInternNew.pointer(), enerInelasNew.pointer() );
  //---------------------------------------------------------------------------------------------------


  // stressNew:
  for (int i = 1; i <= (ndir+nshr); i++) {
    stress_view(ind(i,1), ind(i,2)) = stressNew(1,i);
    stress_view(ind(i,2), ind(i,1)) = stressNew(1,i);
  }

  // update stepTime, totalTime
  stepTime += dt;
  totalTime += dt;

  return 0;
}

