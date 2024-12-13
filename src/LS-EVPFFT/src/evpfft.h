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


#pragma once

#include "mpi.h"
#include "command_line_args.h"
#include "chg_basis.h"
#include "heffte_fft.h"
#include "output_file_manager.h"
#include "micro_output_writer.h"
#include "manager_mpi_io.h"
#include <memory>

class EVPFFT
{
public:
//-----------------------------------------------
// EVPFFT Data Members
//-----------------------------------------------
  
  const MPI_Comm mpi_comm;
  const int root;
  const int my_rank;
  const int num_ranks;
  CommandLineArgs cmd;
  ChgBasis cb;
  std::shared_ptr<FFT3D_R2C<heffte_backend,real_t>> fft;

  int npts1_g;   // global array size in x dim.
  int npts2_g;   // global array size in y dim.
  int npts3_g;   // global array size in z dim.
  int npts1;     // local per rank array size in x dim.
  int npts2;     // local per rank array size in y dim.
  int npts3;     // local per rank array size in z dim.
  int local_start1; // start index of local array in x dim.
  int local_start2; // start index of local array in y dim.
  int local_start3; // start index of local array in z dim.
  int local_end1;   // end index of local array in x dim.
  int local_end2;   // end index of local array in y dim.
  int local_end3;   // end index of local array in z dim.
  real_t wgt;

  int npts1_g_cmplx;   // complex global array size in x dim.
  int npts2_g_cmplx;   // complex global array size in y dim.
  int npts3_g_cmplx;   // complex global array size in z dim.
  int npts1_cmplx;     // complex local per rank array size in x dim.
  int npts2_cmplx;     // complex local per rank array size in y dim.
  int npts3_cmplx;     // complex local per rank array size in z dim.
  int local_start1_cmplx; // complex start index of local array in x dim.
  int local_start2_cmplx; // complex start index of local array in y dim.
  int local_start3_cmplx; // complex start index of local array in z dim.
  int local_end1_cmplx;   // complex end index of local array in x dim.
  int local_end2_cmplx;   // complex end index of local array in y dim.
  int local_end3_cmplx;   // complex end index of local array in z dim.

  int NPHMX;    //Maximum number of phases
  int NMODMX;   //Maximum number of active SL+TW modes in any phase
  int NTWMMX;   //Maximum number of active twin modes in any phase
  int NSYSMX;   //Maximum number of active SL+TW systems in any phase

  MatrixTypeRealDual udot;
  MatrixTypeRealHost dsim;
  MatrixTypeRealHost scauchy;
  MatrixTypeRealHost sdeviat;
  MatrixTypeRealDual tomtot;
  //MatrixTypeRealHost dbar5;
  //MatrixTypeRealHost dbar6;
  MatrixTypeRealHost delt;
  real_t deltvol3;
  real_t tdot;

  MatrixTypeRealDual disgradmacro;
  MatrixTypeRealHost dvelgradmacro;
  MatrixTypeRealHost scauav;
  MatrixTypeRealHost dvelgradmacroacum;
  MatrixTypeRealHost velmax;
  MatrixTypeIntHost iudot;
  MatrixTypeIntHost idsim;
  MatrixTypeIntHost iscau;
  MatrixTypeRealHost defgradavg;
  MatrixTypeRealHost defgradinvavgc_inv;
  int ictrl;
  int nsteps;

  MatrixTypeRealDual dnca;
  MatrixTypeRealDual dbca;
  MatrixTypeRealDual schca;
  MatrixTypeRealDual tau;
  MatrixTypeRealDual hard;
  MatrixTypeRealDual thet;
  MatrixTypeIntDual  nrs;
  MatrixTypeRealDual gamd0;
#ifdef NON_SCHMID_EFFECTS
  MatrixTypeRealHost dtca;
  MatrixTypeRealHost cns;
  MatrixTypeRealDual schcnon;
#endif

  MatrixTypeRealHost twsh;
  MatrixTypeIntHost nsm;
  MatrixTypeIntHost nmodes;
  MatrixTypeIntDual nsyst;
  MatrixTypeIntHost ntwmod;
  MatrixTypeIntHost ntwsys;
  MatrixTypeIntHost isectw;
  MatrixTypeStringHost icryst;

  int nph;
  int nelem;
  int iphbot;
  int iphtop;
  
  const real_t pi;

  real_t error;
  real_t fact2;
  int itmax;
  int irecover;
  int isave;
  int iwtex;
  int iwdeq;

  real_t svm;
  real_t dvm;
  real_t evm;
  real_t errs;
  real_t erre;
  real_t erre2;
  real_t avg_nr_iter;
  int iupdate;
  int iuphard;
  MatrixTypeIntDual igas;
  int iwfields;
  int iwstep;

  MatrixTypeRealHost cc;
  MatrixTypeRealDual c0;
  MatrixTypeRealDual s0;
  MatrixTypeRealDual c066;

  real_t wph1;
  real_t svm1;
  MatrixTypeRealHost scauav1;

  MatrixTypeRealHost eth;
  int ithermo;

  MatrixTypeRealDual xk_gb;
  MatrixTypeRealDual yk_gb;
  MatrixTypeRealDual zk_gb;

  MatrixTypeRealDual sg;
  MatrixTypeRealDual disgrad;
  MatrixTypeRealDual velgrad;
  MatrixTypeRealDual edotp;
  MatrixTypeRealDual cg66;
  MatrixTypeRealDual ept;
  MatrixTypeRealDual ag;
  MatrixTypeRealDual crss;
  MatrixTypeRealDual sch;
  MatrixTypeRealDual sgt;
  MatrixTypeRealDual defgradp;
  MatrixTypeRealDual defgrade;
  MatrixTypeRealDual defgrad;
  MatrixTypeRealDual defgradinv;
  MatrixTypeRealDual defgradini;
  MatrixTypeRealDual wgtc;
  MatrixTypeRealDual detF;
  MatrixTypeRealDual sgPK1;
  MatrixTypeRealDual c066mod;
  MatrixTypeRealDual velgradref;
  MatrixTypeRealDual xnode;
#ifdef NON_SCHMID_EFFECTS
  MatrixTypeRealDual schnon;
#endif

  MatrixTypeRealDual gamdot;
  MatrixTypeRealDual gacumgr;
  //MatrixTypeRealDual trialtau;
  MatrixTypeRealDual xkin;

  MatrixTypeRealHost ph_array;
  MatrixTypeRealHost th_array;
  MatrixTypeRealHost om_array;
  MatrixTypeIntDual  jphase;
  MatrixTypeIntHost  jgrain;

  MatrixTypeRealDual work;
  MatrixTypeRealDual workim;
  MatrixTypeRealDual data;
  MatrixTypeRealDual data_cmplx;

  int elem_id;
  int imicro;
  int iter;
  real_t evmp;
  real_t dvmp;
  MatrixTypeRealHost epav;
  MatrixTypeRealHost edotpav;
  MatrixTypeRealHost velgradmacroactual;
  MatrixTypeRealHost disgradmacrot;
  MatrixTypeRealDual velgradmacro;

  bool active;
  size_t fierro_cycle;
  const real_t stress_scale;
  const real_t time_scale;
  MatrixTypeRealHost M66;
  MatrixTypeRealHost edotp_avg;
  MatrixTypeRealHost dedotp66_avg;
  MatrixTypeRealHost cg66_avg;
  MatrixTypeRealHost sg66_avg;
  MatrixTypeRealHost udotAcc; 
  double dtAcc;

  // For file management
  std::string hdf5_filename;
  OutputFileManager ofile_mgr;
  std::shared_ptr<Manager_MPI_IO<real_t,MPI_ORDER_FORTRAN>> mpi_io_real_t;
  std::shared_ptr<Manager_MPI_IO<int,MPI_ORDER_FORTRAN>> mpi_io_int;
  std::shared_ptr<MicroOutputWriter<MPI_ORDER_FORTRAN>> micro_writer;

  // Device memory space
  MatrixTypeRealDevice  rss;
  MatrixTypeRealDevice rss1;
  MatrixTypeRealDevice rss2;

//-----------------------------------------------
// End Of EVPFFT Data Members
//-----------------------------------------------


//-----------------------------------------------
// EVPFFT Functions
//-----------------------------------------------
  EVPFFT(const MPI_Comm mpi_comm_, const CommandLineArgs cmd_, const real_t stress_scale_=1.0, const real_t time_scale_=1.0);
  ~EVPFFT();
  void vpsc_input();
  void check_iudot();
  void decompose_vel_grad(real_t* vel_grad);
  void init_dvm();
  void check_mixed_bc();
  void init_after_reading_input_data();
  void solve();
  void solve(real_t* vel_grad, real_t* stress, real_t dt, size_t cycle, size_t elem_gid, real_t udotAccThIn);
  void evolve();
  void check_macrostress();
  void print_vel_grad();
  void allocate_memory();

  void set_some_voxels_arrays_to_zero();
  void init_ept();
  void init_disgrad();
  void init_sg();
  void init_xk_gb();
  void init_disgradmacro_velgradmacro();
  void init_evm();
  void init_c0_s0();
  void deinit_c0_s0();
  void init_defgrad();
  void data_crystal(int iph, const std::string & filecryspl);
  void data_crystal_elast(int iph, const std::string & filecrysel);
  void data_grain(const std::string & filetext);
  void step_update_velgrad();
  void update_schmid();
  void step_set_dvelgradmacro_and_dvelgradmacroacum_to_zero();
  void forward_fft();
  void backward_fft();
  void inverse_the_greens();
  void initialize_velgrad();
  void evpal(int imicro);
  void get_smacro();
  void kinhard_param();
  void step_update();
  void step_vm_calc();
  void step_texture_rve_update();
  void update_orient();
  void update_defgradp();
  void update_defgrad();
  void update_defgrade();
  void update_el_stiff();
  void calc_c0();
  void Cauchy_to_PK1();
  void calc_c066mod();
  void update_grid_velgrad();
  void harden(int imicro);
  void read_classic_los_alamos_texture_file(const std::string & filetext);
  void read_hdf5_texture_file(const std::string & filetext);
  void read_vtk_lattice_structure(const std::string & filetext);
  void calculate_eel(MatrixTypeRealDual &eel);
  void write_macro_state();
  void write_micro_state_xdmf();
  void write_micro_state_pvtu();
  void write_texture();

  void init_crss_voce();
  void update_crss_voce();
};
