/**********************************************************************************************
 © 2020. Triad National Security, LLC. All rights reserved.
 This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
 National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
 Department of Energy/National Nuclear Security Administration. All rights in the program are
 reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
 Security Administration. The Government is granted for itself and others acting on its behalf a
 nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
 derivative works, distribute copies to the public, perform publicly and display publicly, and
 to permit others to do so.
 This program is open source under the BSD-3 License.
 Redistribution and use in source and binary forms, with or without modification, are permitted
 provided that the following conditions are met:
 
 1.  Redistributions of source code must retain the above copyright notice, this list of
 conditions and the following disclaimer.
 
 2.  Redistributions in binary form must reproduce the above copyright notice, this list of
 conditions and the following disclaimer in the documentation and/or other materials
 provided with the distribution.
 
 3.  Neither the name of the copyright holder nor the names of its contributors may be used
 to endorse or promote products derived from this software without specific prior
 written permission.
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **********************************************************************************************/


#include "evpfft.h"
#include "utilities.h"
#include <math.h>
#include <algorithm>
#include "Profiler.h"

EVPFFT::EVPFFT(const MPI_Comm mpi_comm_, const CommandLineArgs cmd_, const real_t stress_scale_, const real_t time_scale_)
//-------------------------------------------------
// Data Members needed for EVPFFT Calculations
//-------------------------------------------------

  : cmd (cmd_)
  , mpi_comm (mpi_comm_)
  , root (0)
  , my_rank (get_mpi_comm_rank(mpi_comm))
  , num_ranks (get_mpi_comm_size(mpi_comm))
  , cb()

  , nsteps (0)
  , pi (4.0*atan(1.0))
  , imicro(0)

  , active(false)
  , stress_scale(stress_scale_)
  , time_scale(time_scale_)
  , dtAcc(0.0)

  , ofile_mgr ()
  , hdf5_filename ("micro_state_evpfft.h5")
//-----------------------------------------------
// End Of EVPFFT Data Members
//-----------------------------------------------
{
  Profiler profiler(__FUNCTION__);

  // Reading and initialization of arrays and parameters
  vpsc_input();

  //set_some_voxels_arrays_to_zero();
  init_after_reading_input_data();

  //... For file management
  if (0 == my_rank) {
    ofile_mgr.open_files();
  }

  //...
}

EVPFFT::~EVPFFT(){}

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
  defgradp.update_host();
  gacumgr.update_host();

  for (int jj = 1; jj <= 3; jj++) {
    for (int ii = 1; ii <= 3; ii++) {
      scauav(ii,jj) = 0.0;  // macroscopic stress
      udotAcc(ii,jj) = 0.0;
      // if (ii == jj) {
      //   defgradinvavgc_inv(ii,jj) = 1.0;
      // } else {
      //   defgradinvavgc_inv(ii,jj) = 0.0;
      // }
    }
  }
}

void EVPFFT::init_after_reading_input_data()
{
    init_xk_gb();
    calc_wfhat();
    init_disgradmacro_velgradmacro();
    init_ept();
    init_disgrad();
    init_evm();
    init_defgrad();
    if (ibc == 1) {
      init_Xvert();
    }

#ifndef ABSOLUTE_NO_OUTPUT
    if (iwfields == 1) {
      int imicro = 0;
      //write_micro_state_xdmf();
      write_micro_state_pvtu();
    }
#endif
// the variables initialized in the funcitons below are reduced into
// and should be done once, hence the need for #if guard since the variables
// needs to be initialized after udot and dt are know from fierro
#ifndef BUILD_LSNP_EVPFFT_FIERRO
    if (ibc == 0) init_sg();
    // init_c0_s0();
#endif

}

void EVPFFT::init_xk_gb()
{

  FOR_ALL_CLASS(kxx, 0, npts1, {
        int kx;
        kx = kxx + local_start1_cmplx;
        if (kx >  npts1_g/2) kx = kx-npts1_g;
        xk_gb(kxx+1) = 1.0*kx/npts1_g;
  });

  FOR_ALL_CLASS(kyy, 0, npts2, {
        int ky;
        ky = kyy + local_start2_cmplx;
        if (ky >  npts2_g/2) ky = ky-npts2_g;
        yk_gb(kyy+1) = 1.0*ky/npts2_g;
  });

  FOR_ALL_CLASS(kzz, 0, npts3, {
        int kz;
        kz = kzz + local_start3_cmplx;
        if (kz >  npts3_g/2) kz = kz-npts3_g;
        zk_gb(kzz+1) = 1.0*kz/npts3_g;
  });
  Kokkos::fence();

  // update host
  xk_gb.update_host();
  yk_gb.update_host();
  zk_gb.update_host();
}

void EVPFFT::init_disgradmacro_velgradmacro()
{
  for (int ii = 1; ii <= 3; ii++) {
    for (int jj = 1; jj <= 3; jj++) {
      disgradmacrot(ii,jj) = 0.0;
      //disgradmacro.host(ii,jj)  = udot.host(ii,jj) * tdot;
      disgradmacro.host(ii,jj)  = 0.0;
      velgradmacro.host(ii,jj)  = udot.host(ii,jj);
    }
  }

  // update_device
  disgradmacro.update_device();
  velgradmacro.update_device();
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

void EVPFFT::init_disgrad()
{

  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        disgrad(ii,jj,i,j,k) = 0.0;
      }
    }
  }); // end FOR_ALL_CLASS
  Kokkos::fence();

  // update device
  disgrad.update_host();
}

void EVPFFT::init_sg()
{
  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        sgPK1(ii,jj,i,j,k) = 0.0;
        sg(ii,jj,i,j,k) = 0.0;
        sgt(ii,jj,i,j,k) = 0.0;
      }
    }
  }); // end FOR_ALL_CLASS


  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {


    if (iframe(i,j,k) == 0) {

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

      if (ibc == 0) {
        for (int ii = 1; ii <= 3; ii++) {
          for (int jj = 1; jj <= 3; jj++) {
  
            sg(ii,jj,i,j,k) = 0.0;
  
            for (int kk = 1; kk <= 3; kk++) {
              for (int ll = 1; ll <= 3; ll++) {
                // sg(ii,jj,i,j,k) += cg(ii,jj,kk,ll) * disgradmacro(kk,ll);
                sg(ii,jj,i,j,k) += cg(ii,jj,kk,ll) * udot(kk,ll) * tdot;
              }
            }
          }
        }
      } else if (ibc == 1){
        for (int ii = 1; ii <= 3; ii++) {
          for (int jj = 1; jj <= 3; jj++) {
  
            sg(ii,jj,i,j,k) = 0.0;
  
            for (int kk = 1; kk <= 3; kk++) {
              for (int ll = 1; ll <= 3; ll++) {
                // sg(ii,jj,i,j,k) += cg(ii,jj,kk,ll) * disgradmacro(kk,ll);
                sg(ii,jj,i,j,k) += cg(ii,jj,kk,ll) * eigenvelgradref(kk,ll,i,j,k) * tdot;
              }
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
    }
  }); // end FOR_ALL_CLASS
  Kokkos::fence();

  // update host
  sg.update_host();
  sgt.update_host();
}

void EVPFFT::init_evm()
{
  evm = dvm*tdot;
}

void EVPFFT::init_c0_s0()
{
  for (int ii = 1; ii <= 6; ii++) {
    for (int jj = 1; jj <= 6; jj++) {
      c066.host(ii,jj) *= tdot;
      dF6_dP6.host(ii,jj) /= tdot;
    }
  }

  for (int ii = 1; ii <= 3; ii++) {
    for (int jj = 1; jj <= 3; jj++) {
      for (int kk = 1; kk <= 3; kk ++) {
        for (int ll = 1; ll <= 3; ll ++) {
          c0.host(ii,jj,kk,ll)  *= tdot;
          s0.host(ii,jj,kk,ll)  /= tdot;
        }
      }
    }
  }

  // update_device
  c066.update_device();
  c0.update_device();
  s0.update_device();
  Kokkos::fence();
}

void EVPFFT::deinit_c0_s0()
{
  for (int ii = 1; ii <= 6; ii++) {
    for (int jj = 1; jj <= 6; jj++) {
      c066.host(ii,jj) /= tdot;
      dF6_dP6.host(ii,jj) *= tdot;
    }
  }

  for (int ii = 1; ii <= 3; ii++) {
    for (int jj = 1; jj <= 3; jj++) {
      for (int kk = 1; kk <= 3; kk ++) {
        for (int ll = 1; ll <= 3; ll ++) {
          c0.host(ii,jj,kk,ll)  /= tdot;
          s0.host(ii,jj,kk,ll)  *= tdot;
        }
      }
    }
  }

  // update_device
  c066.update_device();
  c0.update_device();
  s0.update_device();
  Kokkos::fence();
}

void EVPFFT::init_defgrad() {

  for (int jj = 1; jj <= 3; jj++) {
    for (int ii = 1; ii <= 3; ii++) {
      defgradavg(jj,ii) = 0.0;
      defgradinvavgc_inv(jj,ii) = 0.0;
    }
    defgradavg(jj,jj) = delt(jj);
    defgradinvavgc_inv(jj,jj) = delt(jj);
  } 

  DViewFMatrixKokkos <real_t> defgradavg_kokkos(&defgradavg(1,1),3,3);
  DViewFMatrixKokkos <real_t> delt_kokkos(&delt(1),3);
  defgradavg_kokkos.update_device();
  delt_kokkos.update_device();

  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    for (int jj = 1; jj <= 3; jj++) {
      for (int ii = 1; ii <= 3; ii++) {
        defgradp (ii,jj,i,j,k) = 0.0;
        defgrade (ii,jj,i,j,k) = 0.0;
        defgrad (ii,jj,i,j,k) = 0.0;
        defgradinv (ii,jj,i,j,k) = 0.0;
        defgradini (ii,jj,i,j,k) = 0.0;
      }
      defgradp (jj,jj,i,j,k) = 1.0;
      defgrade (jj,jj,i,j,k) = 1.0;
      defgrad (jj,jj,i,j,k) = defgradavg_kokkos(jj,jj);
      defgradinv (jj,jj,i,j,k) = 1.0/defgradavg_kokkos(jj,jj);
      defgradini (jj,jj,i,j,k) = defgradavg_kokkos(jj,jj);
    }
    detF(i,j,k) = delt_kokkos(1)*delt_kokkos(2)*delt_kokkos(3);

    wgtc(i,j,k) = wgt;
  }); // end FOR_ALL_CLASS

  FOR_ALL_CLASS(k, 1, npts3+1,
                j, 1, npts2+1,
                i, 1, npts1+1, {

    // thread private arrays
    real_t x_[3];

    // create views of thread private arrays
    ViewMatrixTypeReal x(x_,3);

    if (igamma == 0) {
      x(1) = 1.0*(i + local_start1 - dnpts1_g);
      x(2) = 1.0*(j + local_start2 - dnpts2_g);
      x(3) = 1.0*(k + local_start3 - dnpts3_g);
    } else if (igamma == 1) {
      x(1) = 1.0*(i + local_start1 - dnpts1_g) - 0.5; // nodal
      x(2) = 1.0*(j + local_start2 - dnpts2_g) - 0.5;
      x(3) = 1.0*(k + local_start3 - dnpts3_g) - 0.5;
    }
    for (int ii = 1; ii <= 3; ii++) {
      x_grid(ii,i,j,k) = 0.0;
      for (int jj = 1; jj <= 3; jj++) {
        x_grid(ii,i,j,k) = x_grid(ii,i,j,k) + defgradavg_kokkos(ii,jj)*x(jj);
      }
    }

  }); // end FOR_ALL_CLASS

  FOR_ALL_CLASS(k, 1, npts3+2,
                j, 1, npts2+2,
                i, 1, npts1+2, {

    // thread private arrays
    real_t x_[3];

    // create views of thread private arrays
    ViewMatrixTypeReal x(x_,3);

    x(1) = 1.0*(i + local_start1 - dnpts1_g) - 0.5;
    x(2) = 1.0*(j + local_start2 - dnpts2_g) - 0.5;
    x(3) = 1.0*(k + local_start3 - dnpts3_g) - 0.5;

    for (int ii = 1; ii <= 3; ii++) {
      xnode(ii,i,j,k) = 0.0;
      for (int jj = 1; jj <= 3; jj++) {
        xnode(ii,i,j,k) = xnode(ii,i,j,k) + defgradavg_kokkos(ii,jj)*x(jj);
      }
    }

  }); // end FOR_ALL_CLASS

}

void EVPFFT::init_Xvert() {

  Xvert(1,1) = 0.5; // 1
  Xvert(2,1) = 0.5; // 1
  Xvert(3,1) = 0.5; // 1

  Xvert(1,2) = 1.0*(npts1_g - dnpts1_g*2) + 0.5; // 2
  Xvert(2,2) = 0.5; // 2
  Xvert(3,2) = 0.5; // 2

  Xvert(1,3) = 1.0*(npts1_g - dnpts1_g*2) + 0.5; // 3
  Xvert(2,3) = 1.0*(npts2_g - dnpts2_g*2) + 0.5; // 3
  Xvert(3,3) = 0.5; // 3

  Xvert(1,4) = 0.5; // 4
  Xvert(2,4) = 1.0*(npts2_g - dnpts2_g*2) + 0.5; // 4
  Xvert(3,4) = 0.5; // 4

  Xvert(1,5) = 0.5; // 5
  Xvert(2,5) = 0.5; // 5
  Xvert(3,5) = 1.0*(npts3_g - dnpts3_g*2) + 0.5; // 5

  Xvert(1,6) = 1.0*(npts1_g - dnpts1_g*2) + 0.5; // 6
  Xvert(2,6) = 0.5; // 6
  Xvert(3,6) = 1.0*(npts3_g - dnpts3_g*2) + 0.5; // 6

  Xvert(1,7) = 1.0*(npts1_g - dnpts1_g*2) + 0.5; // 7
  Xvert(2,7) = 1.0*(npts2_g - dnpts2_g*2) + 0.5; // 7
  Xvert(3,7) = 1.0*(npts3_g - dnpts3_g*2) + 0.5; // 7

  Xvert(1,8) = 0.5; // 8
  Xvert(2,8) = 1.0*(npts2_g - dnpts2_g*2) + 0.5; // 8
  Xvert(3,8) = 1.0*(npts3_g - dnpts3_g*2) + 0.5; // 8

  xvert_current = Xvert;

}

void EVPFFT::evolve()
{

    //for (imicro = 1; imicro <= nsteps; imicro++) {
#ifndef ABSOLUTE_NO_OUTPUT
      if (0 == my_rank) {
        printf(" ***************************************************\n");
        printf("       Current  Time  STEP = %d\n", imicro);
        fprintf(ofile_mgr.conv_file.get(), "STEP = %d\n", imicro);
      }
#endif
      init_c0_s0();
      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          velgradmacro.host(ii,jj)  = udot.host(ii,jj);
        }
      }
      if (imicro == 1) {
        step_update_velgrad();
      }
      // step_update_velgrad();

      if (imicro == 1 || iupdate == 1) {
        update_schmid();
      }

      calc_Goperr0();
      calc_c066mod();

      step_set_dvelgradmacro_and_dvelgradmacroacum_to_zero();

      int itertot;

      err_disp_bc_vel = 2.0*error_bc;
      iter_out = 0;
      itertot = 0;
      while ((iter_out <= itmax_out && err_disp_bc_vel > error_bc) || iter_out < itmin_out) {
      
      iter_out += 1;

      if (ibc == 1) {
        if (iter_out > 1 && mult_incr_accum < mult_incr_mx) frame_c = frame_c*mult_incr;
        if (iter_out > 1 && mult_incr_accum < mult_incr_mx) mult_incr_accum = mult_incr_accum*mult_incr;
        calc_IC0a_inv();
      }

      iter = 0;
      erre = 2.0*error;
      errs = 2.0*error;

      while ((iter < itmax && (errs > error || erre > error)) || iter < 2) {

        iter += 1;
        itertot += 1;

#ifndef ABSOLUTE_NO_OUTPUT
        if (0 == my_rank) {
          printf(" --------------------\n");
          printf(" STEP = %d\n", imicro);
          printf(" ITER = %d %d\n", iter, itertot);
          printf(" DIRECT FFT OF STRESS FIELD\n");
          fprintf(ofile_mgr.conv_file.get(), "ITER = %d\n", iter);
        }
#endif

        // transform to PK1
        Cauchy_to_PK1();

        // forward fft of sgPK1
        forward_fft();

#ifndef ABSOLUTE_NO_OUTPUT
        if (0 == my_rank) {
          printf(" CALCULATING G^pq,ij : SG^ij ...\n");
        }
#endif

        // nr_lu_minval takes local part of work,
        // updates local work, then allgathers
        inverse_the_greens();

#ifndef ABSOLUTE_NO_OUTPUT
        if (0 == my_rank) {
          printf(" INVERSE FFT TO GET STRAIN FIELD\n");
        }
#endif

        // backward fft of work/workim
        backward_fft();

        initialize_velgrad();

#ifndef ABSOLUTE_NO_OUTPUT
        if (0 == my_rank) {
          printf(" UPDATE STRESS FIELD\n\n");
        }
#endif

        // Newton-Raphson iteration
        // evpal takes local disgrad, trialtau, xkin, crss, gacumgr
        // updates local sg, edotp, trialtau, gamdot
        evpal(imicro);

        get_smacro();

        if (dvm > 0.0) erre = sqrt(erre) / dvm;
        if (svm > 0.0) errs = sqrt(errs) / svm;

#ifndef ABSOLUTE_NO_OUTPUT
        if (0 == my_rank) {
          printf(" STRAIN FIELD ERROR = %24.14E\n", erre);
          printf(" STRESS FIELD ERROR = %24.14E\n", errs);
          printf(" AVG NUM OF NR ITER = %24.14E\n", avg_nr_iter);
          fprintf(ofile_mgr.err_file.get(), "%d,%d,%10.4E,%10.4E,%10.4E,%10.4E,%10.4E\n", 
                  imicro, itertot, erre, errs, err_disp_bc_vel, svm, avg_nr_iter);
        }
#endif
      } // end while loop

      calc_velocity();

      if (ibc == 0) {
        err_disp_bc_vel = 0.0;
      } else if (ibc == 1) {
        calc_vel_bc_error();
        if (0 == my_rank) {
#ifndef ABSOLUTE_NO_OUTPUT
          fprintf(ofile_mgr.err_file.get(), "%d,%d,%10.4E,%10.4E,%10.4E,%10.4E,%10.4E\n", 
                  imicro, itertot, erre, errs, err_disp_bc_vel, svm, avg_nr_iter);
#endif
        }
      } 
      

      }

      // kinematic hardening (ithermo == 1)
      if (ithermo == 1 && imicro == 1) {
        kinhard_param();
      }

      step_update();

      step_vm_calc();

      // if (iupdate == 1 && (ithermo != 1 || imicro > 1)) {
      if (iupdate == 1 ) { // do not skip first increment
        // step_texture_rve_update(); // ag now stores inital crystallographic orientation
        update_defgradp();
        // update_grid_velgrad();      
        update_grid();
        update_defgrad();
        update_defgrade();
        update_el_stiff();
        update_xvert();
      }

      if (iuphard == 1 && (ithermo != 1 || imicro > 2)) {
        harden(imicro);
      }

      calc_c0();
      deinit_c0_s0();

#ifndef ABSOLUTE_NO_OUTPUT
      write_macro_state();
      if (iwfields == 1 and imicro % iwstep == 0) {
        //write_micro_state_xdmf();
        write_micro_state_pvtu();
      }
      if (iwtex == 1 and imicro == nsteps) {
        write_texture();
      }
#endif
    //} // end for imicro

  return;
}

void EVPFFT::solve() {
  for (imicro = 1; imicro <= nsteps; imicro++) {
    if (ibc == 1) {
      // velvert
      calc_velvert();
      calc_vel_boundary_lin_el();
      calc_eigenvelgradref();
      if (imicro == 1) {
        init_sg();
      }
    }
    evolve();
  }
  return;
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
