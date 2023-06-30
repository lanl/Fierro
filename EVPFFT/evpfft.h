#pragma once

#include "command_line_args.h"
#include "chg_basis.h"
#include "FFT3D_R2C.h"
#include "output_file_manager.h"

class EVPFFT
{
public:
//-----------------------------------------------
// EVPFFT Data Members
//-----------------------------------------------
  const CommandLineArgs cmd;
  const int npts1;     // RVE size in x dim.
  const int npts2;     // RVE size in y dim.
  const int npts3;     // RVE size in z dim.
  const real_t wgt;
  ChgBasis cb;
  FFT3D_R2C <real_t> fft;

  const int npts1_cmplx;     // RVE complex size in x dim.
  const int npts2_cmplx;     // RVE complex size in y dim.
  const int npts3_cmplx;     // RVE complex size in z dim.

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
  MatrixTypeRealDual tau0_mode;
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
  int iupdate;
  int iuphard;
  int itemphard;
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

  real_t temp_ini;
  real_t temp;
  real_t tempold;
  real_t temp_fact;
  real_t evmpold;
  real_t wplastic;

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
//-----------------------------------------------
// End Of EVPFFT Data Members
//-----------------------------------------------


//-----------------------------------------------
// EVPFFT Functions
//-----------------------------------------------
  EVPFFT(const CommandLineArgs cmd_, const real_t stress_scale_=1.0, const real_t time_scale_=1.0);
  //~EVPFFT();
  void vpsc_input();
  void check_iudot();
  void decompose_vel_grad(real_t* vel_grad);
  void init_dvm();
  void check_mixed_bc();
  void init_after_reading_input_data();
  void solve(real_t* vel_grad, real_t* stress, real_t dt, size_t cycle, size_t elem_gid);
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
  void write_micro_state(int imicro);
  void write_texture();

  void init_crss_voce();
  void init_crss_temp();
  void update_crss_voce();
  void update_crss_temp();
  void update_temperature();
};
