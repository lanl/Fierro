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

  const int npts1_g;   // global array size in x dim.
  const int npts2_g;   // global array size in y dim.
  const int npts3_g;   // global array size in z dim.
  const int npts1;     // local per rank array size in x dim.
  const int npts2;     // local per rank array size in y dim.
  const int npts3;     // local per rank array size in z dim.
  const int local_start1; // start index of local array in x dim.
  const int local_start2; // start index of local array in y dim.
  const int local_start3; // start index of local array in z dim.
  const int local_end1;   // end index of local array in x dim.
  const int local_end2;   // end index of local array in y dim.
  const int local_end3;   // end index of local array in z dim.
  const real_t wgt;

  const int npts1_g_cmplx;   // complex global array size in x dim.
  const int npts2_g_cmplx;   // complex global array size in y dim.
  const int npts3_g_cmplx;   // complex global array size in z dim.
  const int npts1_cmplx;     // complex local per rank array size in x dim.
  const int npts2_cmplx;     // complex local per rank array size in y dim.
  const int npts3_cmplx;     // complex local per rank array size in z dim.
  const int local_start1_cmplx; // complex start index of local array in x dim.
  const int local_start2_cmplx; // complex start index of local array in y dim.
  const int local_start3_cmplx; // complex start index of local array in z dim.
  const int local_end1_cmplx;   // complex end index of local array in x dim.
  const int local_end2_cmplx;   // complex end index of local array in y dim.
  const int local_end3_cmplx;   // complex end index of local array in z dim.

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
  MatrixTypeRealDual ddisgradmacro;
  MatrixTypeRealHost scauav;
  MatrixTypeRealHost ddisgradmacroacum;
  MatrixTypeRealHost velmax;
  MatrixTypeIntHost iudot;
  MatrixTypeIntHost idsim;
  MatrixTypeIntHost iscau;
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
  int ngr;
  
  const real_t pi;
  int jran;

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

  const int num_crystals; //# OF Crystals (grains) in the input texture
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
  MatrixTypeRealHost disgradmacroactual;
  MatrixTypeRealHost disgradmacrot;
  MatrixTypeRealHost velgradmacro;

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
  OutputFileManager ofile_mgr;
  std::shared_ptr<Manager_MPI_IO<real_t,MPI_ORDER_FORTRAN>> mpi_io_real_t;
  std::shared_ptr<Manager_MPI_IO<int,MPI_ORDER_FORTRAN>> mpi_io_int;
  std::shared_ptr<MicroOutputWriter<MPI_ORDER_FORTRAN>> micro_writer;
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

  void set_some_voxels_arrays_to_zero();
  void init_ept();
  void init_sg();
  void init_xk_gb();
  void init_disgradmacro();
  void init_evm();
  void data_crystal(int iph, const std::string & filecryspl);
  void data_crystal_elast(int iph, const std::string & filecrysel);
  void data_grain(const std::string & filetext);
  void step_update_disgrad();
  void update_schmid();
  void step_set_ddisgradmacro_and_ddisgradmacroacum_to_zero();
  void forward_fft();
  void backward_fft();
  void inverse_the_greens();
  void initialize_disgrad();
  void evpal(int imicro);
  void get_smacro();
  void kinhard_param();
  void step_update_velgrad_etc();
  void step_vm_calc();
  void step_texture_rve_update();
  void update_orient();
  void harden(int imicro);
  void read_classic_los_alamos_texture_file(const std::string & filetext, int & nph1);
  void read_hdf5_texture_file(const std::string & filetext, int & nph1);
  void calculate_eel(MatrixTypeRealDual &eel);
  void write_macro_state();
  void write_micro_state();
  void write_texture();

  void init_crss_voce();
  void update_crss_voce();
};
