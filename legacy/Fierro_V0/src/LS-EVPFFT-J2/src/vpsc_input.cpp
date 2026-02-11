// ********************************************************************
//     SUBROUTINE VPSC_INPUT      --->      VERSION 31/jan/99
//
//     READS CHARACTERISTICS OF THE RUN: # OF PHASES, NAMES OF INPUT FILES,
//     DEFORMATION TO BE IMPOSED, CONVERGENCE PARAMETERS, ETC.
//     READS SINGLE CRYSTAL PROPERTIES: DEFORMATION MODES, CRSS, HARDENING
//     READS CRYSTAL AND MORPHOLOGIC TEXTURES.
//     INITIALIZES ARRAYS REQUIRED TO RUN VPSC.
//     OPENS AND CLOSES INPUT FILES.   OPENS OUTPUT FILES.
//
//     MODIFIED 21/07/98 by CNT:
//     INITIALIZATION RELATED TO 'ELEMENTS' IS DONE INSIDE A SINGLE BLOCK.
// *****************************************************************************



#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <fstream>
#include <cmath>
#include <string>
#include <fstream>

#include "evpfft.h"
#include "definitions.h"
#include "utilities.h"

#include "math_functions.h"

using namespace utils;


void EVPFFT::vpsc_input()
{
  real_t dum;
  real_t dummy;
  real_t xmi;
  real_t sigma;
  int iph;

  std::ifstream ur0;
  std::string filetext;
  std::string filecryspl;
  std::string filecrysel;
  std::string prosa;
  std::string filethermo;

// *********   INITIALIZATION BLOCK   ***************************
  
  ur0.open(cmd.input_filename);
  check_that_file_is_open(ur0, cmd.input_filename.c_str());

  ur0 >> NPHMX >> NMODMX >> NTWMMX >> NSYSMX; CLEAR_LINE(ur0);
  ur0 >> npts1_g >> npts2_g >> npts3_g; CLEAR_LINE(ur0);

  allocate_memory();

  ur0 >> nph; CLEAR_LINE(ur0);

  // THE FOLLOWING REQUIRED FOR SEVERAL ROUTINES WITH 'do iph=iphbot,iphtop'
  iphbot = 1;
  iphtop = nph;

  // RVE DIMENSIONS
  ur0 >> delt(1) >> delt(2) >> delt(3); CLEAR_LINE(ur0);
  deltvol3 = pow(delt(1) * delt(2) * delt(3), 1.0/3.0);

  ur0 >> prosa; CLEAR_LINE(ur0);
  ur0 >> filetext; CLEAR_LINE(ur0);

// ***************************************************************************
//     LOOP OVER PHASES
// ***************************************************************************
  for (int iph = 1; iph <= nph; iph++) {

    ur0 >> prosa; CLEAR_LINE(ur0);
    ur0 >> igas.host(iph); CLEAR_LINE(ur0);
    ur0 >> prosa; CLEAR_LINE(ur0);
    ur0 >> filecryspl; CLEAR_LINE(ur0);
    ur0 >> filecrysel; CLEAR_LINE(ur0);
    
    // READS SLIP AND TWINNING MODES FOR THE PHASE

    if (igas.host(iph) == 0) {
      data_crystal(iph, filecryspl);
      data_crystal_elast(iph, filecrysel);
    } // end if (igas.host(iph) == 0)

  } // end for (int iph = 1; iph <= nph; iph++)
 
  // update device (these variables were modified in data_crystal)
  igas.update_device();
  nsyst.update_device();
  gamd0.update_device();
  nrs.update_device();
  tau.update_device();
  thet.update_device();
  dnca.update_device();
  dbca.update_device();
  schca.update_device();
  sigma0.update_device();
  sigma1.update_device();
  thet0_j2.update_device();
  thet1_j2.update_device();
  nrs_j2.update_device();
  edotp0_j2.update_device();
  iJ2.update_device();
#ifdef NON_SCHMID_EFFECTS
  schcnon.update_device();
#endif
  hard.update_device();

  // READS INITIAL TEXTURE FROM FILETEXT
  // and calculates the local elastic stiffness
  data_grain(filetext);
 

// ****************************************************************************
//     READ BOUNDARY CONDITIONS ON OVERALL STRESS AND STRAIN-RATE
// ****************************************************************************

  ur0 >> prosa; CLEAR_LINE(ur0);
  ur0 >> prosa; CLEAR_LINE(ur0);

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      ur0 >> iudot(i,j);
    }
    CLEAR_LINE(ur0);
  }
 
  // check that the components of iudot are correct
  check_iudot();

  CLEAR_LINE(ur0);
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      ur0 >> udot.host(i,j);
    }
    CLEAR_LINE(ur0);
  }
  //update device
  udot.update_device();
  Kokkos::fence();

  // calculate strain-rate, rotation-rate, and dvm from udot
  decompose_vel_grad(udot.host_pointer());

  // initialize dvm used to calculate evm which is the denominatior for TOL_strain_field
  init_dvm();

  CLEAR_LINE(ur0);
  ur0 >> iscau(1) >> iscau(6) >> iscau(5); CLEAR_LINE(ur0);
  ur0 >> iscau(2) >> iscau(4); CLEAR_LINE(ur0);
  ur0 >> iscau(3); CLEAR_LINE(ur0);

  // check that mixed BC are correct
  check_mixed_bc();

  CLEAR_LINE(ur0);
  ur0 >> scauchy(1,1) >> scauchy(1,2) >> scauchy(1,3); CLEAR_LINE(ur0);
  ur0 >> scauchy(2,2) >> scauchy(2,3); CLEAR_LINE(ur0);
  ur0 >> scauchy(3,3); CLEAR_LINE(ur0);

  scauchy(3,2) = scauchy(2,3);
  scauchy(3,1) = scauchy(1,3);
  scauchy(2,1) = scauchy(1,2);

  CLEAR_LINE(ur0);
  ur0 >> dummy; CLEAR_LINE(ur0);
  ur0 >> ictrl; CLEAR_LINE(ur0);

  if (ictrl == -1) tdot=dummy;
  if (ictrl == 0)  tdot=dummy/dvm;
  if (ictrl > 0) {
    printf("ICTRL>0 NOT IMPLEMENTED\n");
    exit(1);
  }

  CLEAR_LINE(ur0);
  ur0 >> nsteps; CLEAR_LINE(ur0);
  ur0 >> error;  CLEAR_LINE(ur0);
  ur0 >> itmax;  CLEAR_LINE(ur0);
  ur0 >> irecover; CLEAR_LINE(ur0);

  if (irecover == 1) {
    std::ifstream f50;
    f50.open("stress.in");
    check_that_file_is_open(f50, "stress.in");
  }

  ur0 >> isave;   CLEAR_LINE(ur0);
  ur0 >> iupdate; CLEAR_LINE(ur0);
  ur0 >> iuphard; CLEAR_LINE(ur0);
  ur0 >> iwtex;   CLEAR_LINE(ur0);
  ur0 >> iwfields >> iwstep; CLEAR_LINE(ur0);

  init_crss_voce();

  ur0 >> ithermo; CLEAR_LINE(ur0);
  if (ithermo == 1) {
    int num_crystals = 200000; // maximum number of crystals/grains
    eth = MatrixTypeRealHost (3,3,num_crystals);

    ur0 >> filethermo; CLEAR_LINE(ur0);
    std::ifstream f49;
    f49.open(filethermo);
    check_that_file_is_open(f49, filethermo.c_str());

    for (int ig = 1; ig <= num_crystals; ig++) {
      f49 >> eth(1,1,ig) >> eth(2,2,ig)
          >> eth(3,3,ig) >> eth(2,3,ig)
          >> eth(3,1,ig) >> eth(1,2,ig); CLEAR_LINE(ur0);

      eth(3,2,ig) = eth(2,3,ig);
      eth(1,3,ig) = eth(3,1,ig);
      eth(2,1,ig) = eth(1,2,ig);
    }

    f49.close();

  //} 
  //else {
  //  for (int ig = 1; ig <= num_crystals; ig++) {
  //    for (int j = 1; j <= 3; j++) {
  //      for (int i = 1; i <= 3; i++) {
  //        eth(i,j,ig) = 0.0;
  //      }
  //    }
  //  }
  } // end if (ithermo == 1)

  ur0.close();
}

void EVPFFT::check_iudot()
{
  if (iudot(1,1) + iudot(2,2) + iudot(3,3) == 2) {
    printf("\nCHECK DIAGONAL BOUNDARY CONDITIONS IUDOT\n");
    printf("CANNOT ENFORCE ONLY TWO DEVIATORIC COMPONENTS\n");
    exit(1);
  }

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      if ((i != j) && (iudot(i,j) + iudot(j,i) == 0)) {
        printf("CHECK OFF-DIAGONAL BOUNDARY CONDITIONS IUDOT\n");
        exit(1);
      }
    }
  }
}

void EVPFFT::decompose_vel_grad(real_t* vel_grad_ptr)
{
//     SYMMETRIC STRAIN-RATE, ANTISYMMETRIC ROTATION-RATE TENSORS
//     AND INDICES OF IMPOSED COMPONENTS

  ViewMatrixTypeReal vel_grad(vel_grad_ptr,3,3);

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      dsim(i,j) = (vel_grad(i,j) + vel_grad(j,i)) / 2.0;
      tomtot.host(i,j) = (vel_grad(i,j) - vel_grad(j,i)) / 2.0;
    }
  }
  // update device
  tomtot.update_device();
  Kokkos::fence();
}

void EVPFFT::init_dvm()
{
//     WRITES STRAIN RATE DSIM(I,J) IN b-BASIS AS A 5-DIM VECTOR DBAR(K)
  MatrixTypeRealHost dbar5 (5);
  MatrixTypeRealHost dbar6 (6);
  MatrixTypeRealHost dsimdev (3,3);
  cb.chg_basis_2(dbar6.pointer(), dsim.pointer(), 2, 6, cb.B_basis_host_pointer());

  for (int i = 1; i <= 5; i++) {
    dbar5(i) = dbar6(i);
  }

  cb.chg_basis_1(dbar5.pointer(), dsimdev.pointer(), 1, 5, cb.B_basis_host_pointer());

  dvm = 0.0;
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      dvm += POW2(dsimdev(i,j));
    }
  }
  dvm = sqrt(2.0/3.0 * dvm);
}

void EVPFFT::check_mixed_bc()
{
  for (int i = 1; i <= 3; i++) {
    idsim(i) = iudot(i,i);
  }

  idsim(4) = 0;
  if ((iudot(2,3) == 1) && (iudot(3,2) == 1)) idsim(4) = 1;
  idsim(5) = 0;
  if ((iudot(1,3) == 1) && (iudot(3,1) == 1)) idsim(5) = 1;
  idsim(6) = 0;
  if ((iudot(1,2) == 1) && (iudot(2,1) == 1)) idsim(6) = 1;

  for (int i = 1; i <= 6; i++) {
    if ((iscau(i) * idsim(i) != 0) || (iscau(i) + idsim(i) != 1)) {
      printf("CHECK BOUNDARY CONDITS ON STRAIN-RATE AND STRESS\n");
      printf("IDSIM = ");
      for (int j = 1; j <= 6; j++) printf("%d ", idsim(j));
      printf("\n");

      printf("ISCAU = ");
      for (int j = 1; j <= 6; j++) printf("%d ", iscau(j));
      printf("\n");

      exit(1);
    }
  }
}

void EVPFFT::init_crss_voce()
{
  for (int kk = 1; kk <= npts3; kk++) {
    for (int jj = 1; jj <= npts2; jj++) {
      for (int ii = 1; ii <= npts1; ii++) {

        int iph = jphase.host(ii,jj,kk);

        if (igas.host(iph) == 0) {
          if (iJ2.host(iph) == 0) {
            for (int i = 1; i <= nsyst.host(iph); i++) {
              crss.host(i,1,ii,jj,kk) = tau.host(i,1,iph);
              crss.host(i,2,ii,jj,kk) = tau.host(i,2,iph);
              //trialtau.host(i,1,ii,jj,kk)  = tau.host(i,1,iph);
              //trialtau.host(i,2,ii,jj,kk)  = tau.host(i,2,iph);
            }
          } else {
            sigma0gr.host(ii,jj,kk) = sigma0.host(iph);
          }
        } // end if (igas.host(iph) == 0)

      }  // end for ii
    }  // end for jj
  } // end for kk

  crss.update_device();
  sigma0gr.update_device();
}

