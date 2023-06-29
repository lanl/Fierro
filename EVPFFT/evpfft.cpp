#include "evpfft.h"
#include "utilities.h"
#include <math.h>
#include <algorithm>

//#ifndef NDEBUG
#include <fenv.h> // for debugging floating point exceptions
//#endif

EVPFFT::EVPFFT(const CommandLineArgs cmd_, const real_t stress_scale_, const real_t time_scale_)
//-------------------------------------------------
// Data Members needed for EVPFFT Calculations
//-------------------------------------------------
  : cmd (cmd_)
  , npts1 (cmd.nn[0])
  , npts2 (cmd.nn[1])
  , npts3 (cmd.nn[2])
  , wgt (real_t(1.0) / real_t(npts1*npts2*npts3))
  , cb ()
  , fft (npts1,npts2,npts3)

  , npts1_cmplx (npts1/2 + 1)
  , npts2_cmplx (npts2)
  , npts3_cmplx (npts3)

  , udot (3,3)
  , dsim (3,3)
  , scauchy (3,3)
  , sdeviat (3,3)
  , tomtot (3,3)
  //, dbar5 (5)
  //, dbar6 (6)
  , delt (3)

  , disgradmacro (3,3)
  , ddisgradmacro (3,3)
  , scauav (3,3)
  , ddisgradmacroacum (3,3)
  , velmax (3)
  , iudot (3,3)
  , idsim (6)
  , iscau (6)
  , nsteps (0)

  , dnca (3,NSYSMX,NPHMX)
  , dbca (3,NSYSMX,NPHMX)
  , schca (5,NSYSMX,NPHMX)
  , tau (NSYSMX,3,NPHMX)
  , tau0_mode (NSYSMX,3,NPHMX)
  , hard (NSYSMX,NSYSMX,NPHMX)
  , thet (NSYSMX,2,NPHMX)
  , nrs (NSYSMX,NPHMX)
  , gamd0 (NSYSMX,NPHMX)
#ifdef NON_SCHMID_EFFECTS
  , dtca (3,NSYSMX,NPHMX)
  , cns (5,NMODMX,NPHMX)
  , schcnon (5,NSYSMX,NPHMX)
#endif

  , twsh (NSYSMX,NPHMX)
  , nsm (NMODMX,NPHMX)
  , nmodes (NPHMX)
  , nsyst (NPHMX)
  , ntwmod (NPHMX)
  , ntwsys (NPHMX)
  , isectw (NSYSMX,NPHMX)
  , icryst (NPHMX)

  , pi (4.0*atan(1.0))

  , igas (NPHMX)

  , cc (3,3,3,3,NPHMX)
  , c0 (3,3,3,3)
  , s0 (3,3,3,3)
  , c066 (6,6)

  , scauav1 (3,3)

  , num_crystals (200000)
  //, eth (3,3,num_crystals) // eth will be allocated if ithermo==1

  , xk_gb (npts1)
  , yk_gb (npts2)
  , zk_gb (npts3)

  , sg (3, 3, npts1, npts2, npts3)
  , disgrad (3, 3, npts1, npts2, npts3)
  , velgrad (3, 3, npts1, npts2, npts3)
  , edotp (3, 3, npts1, npts2, npts3)
  , cg66 (6, 6, npts1, npts2, npts3)
  , ept (3, 3, npts1, npts2, npts3)
  , ag (3, 3, npts1, npts2, npts3)
  , crss (NSYSMX, 2, npts1, npts2, npts3)
  , sch (5, NSYSMX, npts1, npts2, npts3)
#ifdef NON_SCHMID_EFFECTS
  , schnon (5, NSYSMX, npts1, npts2, npts3)
#endif

  , gamdot (NSYSMX, npts1, npts2, npts3)
  , gacumgr (npts1, npts2, npts3)
  //, trialtau (NSYSMX, 2, npts1, npts2, npts3)
  , xkin (NSYSMX, npts1, npts2, npts3)

  , ph_array (npts1, npts2, npts3)
  , th_array (npts1, npts2, npts3)
  , om_array (npts1, npts2, npts3)
  , jphase (npts1, npts2, npts3)
  , jgrain (npts1, npts2, npts3)

  , work (3, 3, npts1, npts2, npts3)
  , workim (3, 3, npts1_cmplx, npts2_cmplx, npts3_cmplx)
  , data (npts1, npts2, npts3)
  , data_cmplx (2, npts1_cmplx, npts2_cmplx, npts3_cmplx)

  , imicro(0)
  , epav (3,3)
  , edotpav (3,3)
  , disgradmacroactual (3,3)
  , disgradmacrot (3,3)
  , velgradmacro (3,3)

  , active(false)
  , stress_scale(stress_scale_)
  , time_scale(time_scale_)
  , M66 (6,6)
  , udotAcc(3,3)
  , dtAcc(0.0)

  , ofile_mgr ()
//-----------------------------------------------
// End Of EVPFFT Data Members
//-----------------------------------------------

{
  fft.make_plan(data.host_pointer(), data_cmplx.host_pointer());
  set_some_voxels_arrays_to_zero();
  // Read input parameters and do some setup
  vpsc_input();
  init_after_reading_input_data();
}

//EVPFFT::~EVPFFT(){}

void EVPFFT::set_some_voxels_arrays_to_zero()
{
  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    for (int jj = 1; jj <= 6; jj++) {
      for (int ii = 1; ii <= 6; ii++) {
        cg66(ii,jj,i,j,k) = 0.0;
      }
    }

    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        ag(ii,jj,i,j,k) = 0.0;
        velgrad(ii,jj,i,j,k) = 0.0;
      }
    }

    for (int ii = 1; ii <= NSYSMX; ii++) {     
      xkin(ii,i,j,k) = 0.0;
      gamdot(ii,i,j,k) = 0.0;
    }

    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        disgrad(ii,jj,i,j,k) = 0.0;
        edotp(ii,jj,i,j,k) = 0.0;
        sg(ii,jj,i,j,k) = 0.0;
      }
    }

    gacumgr(i,j,k) = 0.0;

  }); // end FOR_ALL_CLASS
  Kokkos::fence();



  // update host
  cg66.update_host();
  ag.update_host();
  velgrad.update_host();
  gamdot.update_host();
  disgrad.update_host();
  edotp.update_host();
  sg.update_host();
  gacumgr.update_host();

  for (int jj = 1; jj <= 3; jj++) {
    for (int ii = 1; ii <= 3; ii++) {
      scauav(ii,jj) = 0.0;  // macroscopic stress
      udotAcc(ii,jj) = 0.0;
    }
  }
}

void EVPFFT::init_after_reading_input_data()
{
    init_xk_gb();
    init_disgradmacro();
    init_ept();
    init_sg();
    init_evm();
}

void EVPFFT::init_xk_gb()
{
  real_t delt1 = delt(1);
  real_t delt2 = delt(2);
  real_t delt3 = delt(3);

  FOR_ALL_CLASS(kxx, 0, npts1, {
        int kx;
        kx = kxx;
        if (kx >  npts1/2) kx = kx-npts1;
        xk_gb(kxx+1) = kx/(delt1*npts1);
  });

  FOR_ALL_CLASS(kyy, 0, npts2, {
        int ky;
        ky = kyy;
        if (ky >  npts2/2) ky = ky-npts2;
        yk_gb(kyy+1) = ky/(delt2*npts2);
  });

  FOR_ALL_CLASS(kzz, 0, npts3, {
        int kz;
        kz = kzz;
        if (kz >  npts3/2) kz = kz-npts3;
        zk_gb(kzz+1) = kz/(delt3*npts3);
  });

  // update host
  xk_gb.update_host();
  yk_gb.update_host();
  zk_gb.update_host();
  Kokkos::fence();
}

void EVPFFT::init_disgradmacro()
{
  for (int ii = 1; ii <= 3; ii++) {
    for (int jj = 1; jj <= 3; jj++) {
      disgradmacrot(ii,jj) = 0.0;
      disgradmacro.host(ii,jj)  = udot.host(ii,jj) * tdot;
    }
  }

  // update_device
  disgradmacro.update_device();
  Kokkos::fence();
}

void EVPFFT::init_ept()
{

  for (int k = 1; k <= npts3; k++) {
  for (int j = 1; j <= npts2; j++) {
  for (int i = 1; i <= npts1; i++) {

    int jgr;
    jgr = jgrain(i,j,k);
    if (ithermo == 1 and jgr > 0) {
      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          ept.host(ii,jj,i,j,k) = eth(ii,jj,jgr);
        }
      }
    } else {
      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          ept.host(ii,jj,i,j,k) = 0.0;
        }
      }
    }

  } // end for i
  } // end for j
  } // end for k

  // update device
  ept.update_device();
  Kokkos::fence();

}

void EVPFFT::init_sg()
{
  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        sg(ii,jj,i,j,k) = 0.0;
      }
    }
  }); // end FOR_ALL_CLASS


  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    int jph;
    real_t cg66aux_[6*6];
    real_t cg_[3*3*3*3];
    ViewMatrixTypeReal cg66aux(cg66aux_,6,6);
    ViewMatrixTypeReal cg(cg_,3,3,3,3);

    jph = jphase(i,j,k);

    if (igas(jph) == 0) {
      for (int ii = 1; ii <= 6; ii++) {
        for (int jj = 1; jj <= 6; jj++) {
          cg66aux(ii,jj) = cg66(ii,jj,i,j,k);
        }
      }

      cb.chg_basis_3(cg66aux.pointer(), cg.pointer(), 3, 6, cb.B_basis_device_pointer());

      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {

          sg(ii,jj,i,j,k) = 0.0;

          for (int kk = 1; kk <= 3; kk++) {
            for (int ll = 1; ll <= 3; ll++) {
              sg(ii,jj,i,j,k) += cg(ii,jj,kk,ll) * disgradmacro(kk,ll);
            }
          }
        }
      }

    } else {
      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          sg(ii,jj,i,j,k) = 0.0;
        }
      }
    }
  }); // end FOR_ALL_CLASS

  // update host
  sg.update_host();
  Kokkos::fence(); 
}

void EVPFFT::init_evm()
{
  evm = dvm*tdot;
}

void EVPFFT::evolve()
{

    //for (imicro = 1; imicro <= nsteps; imicro++) {
#ifndef AbsoluteNoOutput
      printf(" ***************************************************\n");
      printf("       Current  Time  STEP = %d\n", imicro);
      if (nsteps != 1) {
        fprintf(ofile_mgr.err_file.get(), "STEP = %d\n", imicro);
      }
      fprintf(ofile_mgr.conv_file.get(), "STEP = %d\n", imicro);
#endif

      step_update_disgrad();

      if (imicro == 1 || iupdate == 1) {
        update_schmid();
      }

      step_set_ddisgradmacro_and_ddisgradmacroacum_to_zero();

      iter = 0;
      erre = 2.0*error;
      errs = 2.0*error;

      while (iter < itmax && (errs > error || erre > error)) {

        iter += 1;

#ifndef AbsoluteNoOutput
        printf(" --------------------\n");
        printf(" STEP = %d\n", imicro);
        printf(" ITER = %d\n", iter);
        printf(" DIRECT FFT OF STRESS FIELD\n");
        fprintf(ofile_mgr.conv_file.get(), "ITER = %d\n", iter);
#endif

        // forward fft of sg
        forward_fft();

#ifndef AbsoluteNoOutput
        printf(" CALCULATING G^pq,ij : SG^ij ...\n");
#endif

        // nr_lu_minval takes local part of work,
        // updates local work, then allgathers
        inverse_the_greens();

#ifndef AbsoluteNoOutput
        printf(" INVERSE FFT TO GET STRAIN FIELD\n");
#endif

        // backward fft of work/workim
        backward_fft();

        initialize_disgrad();

#ifndef AbsoluteNoOutput
        printf(" UPDATE STRESS FIELD\n\n");
#endif

        // Newton-Raphson iteration
        // evpal takes local disgrad, trialtau, xkin, crss, gacumgr
        // updates local sg, edotp, trialtau, gamdot
        evpal(imicro);

        get_smacro();

        if (evm > 0.0) erre = erre / evm;
        if (svm > 0.0) errs = errs / svm;

#ifndef AbsoluteNoOutput
        printf(" STRAIN FIELD ERROR = %24.14E\n", erre);
        printf(" STRESS FIELD ERROR = %24.14E\n", errs);
        fprintf(ofile_mgr.err_file.get(), "%d,%10.4E,%10.4E,%10.4E\n", 
                iter, erre, errs, svm);
#endif
      } // end while loop

      // kinematic hardening (ithermo == 1)
      if (ithermo == 1 && imicro == 1) {
        kinhard_param();
      }

      step_update_velgrad_etc();

      step_vm_calc();

      if (iupdate == 1 && (ithermo != 1 || imicro > 1)) {
        step_texture_rve_update();
      }

      if (iuphard == 1 && (ithermo != 1 || imicro > 2)) {
        harden(imicro);
      }

#ifndef AbsoluteNoOutput
      write_macro_state();
      write_micro_state(imicro);
#endif
    //} // end for imicro

#ifndef AbsoluteNoOutput
  if (imicro == nsteps) write_texture();
#endif
}


void EVPFFT::check_macrostress()
{
  bool nan_stress = false;

  for (int j = 1; j <= 3; j++) {
    for (int i = 1; i <= 3; i++) {
      if (isnan(scauav(i,j))) {
        nan_stress = true;
      }
    } // end for i
  } // end for j

  if (nan_stress == true) {
    printf("NaN stress at elem %d, cycle %d, dt %24.14E \n", elem_id, fierro_cycle, tdot);
    print_vel_grad();
    exit(1);
  }

}

void EVPFFT::print_vel_grad()
{
    printf("vel_grad at elem_gid %d, cycle %d, dt %24.14E is :\n", elem_id, fierro_cycle, tdot);

    for (int j = 1; j <= 3; j++) {
      for (int i = 1; i <= 3; i++) {
        printf("  velgrad_%d%d : %24.14E", i, j, udot.host(i,j));
      } // end for i
      printf("\n");
    } // end for j

}
