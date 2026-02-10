#include "evpfft.h"

void EVPFFT::allocate_memory()
{

  fft = std::make_shared<FFT3D_R2C<heffte_backend,real_t>>(mpi_comm, std::array<int,3>{npts1_g,npts2_g,npts3_g});

  npts1 = fft->localRealBoxSizes[my_rank][0];
  npts2 = fft->localRealBoxSizes[my_rank][1];
  npts3 = fft->localRealBoxSizes[my_rank][2];
  local_start1 = fft->localRealBoxes[my_rank].low[0];
  local_start2 = fft->localRealBoxes[my_rank].low[1];
  local_start3 = fft->localRealBoxes[my_rank].low[2];
  local_end1 = fft->localRealBoxes[my_rank].high[0];
  local_end2 = fft->localRealBoxes[my_rank].high[1];
  local_end3 = fft->localRealBoxes[my_rank].high[2];
  wgt = real_t(1.0) / real_t(npts1_g*npts2_g*npts3_g);

  npts1_g_cmplx = fft->globalComplexBoxSize[0];
  npts2_g_cmplx = fft->globalComplexBoxSize[0];
  npts3_g_cmplx = fft->globalComplexBoxSize[0];
  npts1_cmplx = fft->localComplexBoxSizes[my_rank][0];
  npts2_cmplx = fft->localComplexBoxSizes[my_rank][1];
  npts3_cmplx = fft->localComplexBoxSizes[my_rank][2];
  local_start1_cmplx = fft->localComplexBoxes[my_rank].low[0];
  local_start2_cmplx = fft->localComplexBoxes[my_rank].low[1];
  local_start3_cmplx = fft->localComplexBoxes[my_rank].low[2];
  local_end1_cmplx = fft->localComplexBoxes[my_rank].high[0];
  local_end2_cmplx = fft->localComplexBoxes[my_rank].high[1];
  local_end3_cmplx = fft->localComplexBoxes[my_rank].high[2];

  udot  = MatrixTypeRealDual (3,3);
  dsim = MatrixTypeRealHost (3,3);
  scauchy = MatrixTypeRealHost (3,3);
  sdeviat = MatrixTypeRealHost (3,3);
  tomtot = MatrixTypeRealDual (3,3);
  delt = MatrixTypeRealHost (3);

  disgradmacro = MatrixTypeRealDual (3,3);
  dvelgradmacro = MatrixTypeRealHost (3,3);
  scauav = MatrixTypeRealHost (3,3);
  dvelgradmacroacum = MatrixTypeRealHost (3,3);
  velmax = MatrixTypeRealHost (3);
  iudot = MatrixTypeIntHost (3,3);
  idsim = MatrixTypeIntHost (6);
  iscau = MatrixTypeIntHost (6);
  defgradinvavgc_inv = MatrixTypeRealHost (3,3);

  dnca = MatrixTypeRealDual (3,NSYSMX,NPHMX);
  dbca = MatrixTypeRealDual (3,NSYSMX,NPHMX);
  schca = MatrixTypeRealDual (5,NSYSMX,NPHMX);
  tau = MatrixTypeRealDual (NSYSMX,3,NPHMX);
  hard = MatrixTypeRealDual (NSYSMX,NSYSMX,NPHMX);
  thet = MatrixTypeRealDual (NSYSMX,2,NPHMX);
  nrs = MatrixTypeIntDual (NSYSMX,NPHMX);
  gamd0 = MatrixTypeRealDual (NSYSMX,NPHMX);
#ifdef NON_SCHMID_EFFECTS
  dtca = MatrixTypeRealHost (3,NSYSMX,NPHMX);
  cns = MatrixTypeRealHost (5,NMODMX,NPHMX);
  schcnon = MatrixTypeRealDual (5,NSYSMX,NPHMX);
#endif
  sigma0 = MatrixTypeRealDual (NPHMX);
  sigma1 = MatrixTypeRealDual (NPHMX);
  thet0_j2 = MatrixTypeRealDual (NPHMX);
  thet1_j2 = MatrixTypeRealDual (NPHMX);
  nrs_j2 = MatrixTypeIntDual (NPHMX);
  edotp0_j2 = MatrixTypeRealDual (NPHMX);

  twsh = MatrixTypeRealHost (NSYSMX,NPHMX);
  nsm = MatrixTypeIntHost (NMODMX,NPHMX);
  nmodes = MatrixTypeIntHost (NPHMX);
  nsyst = MatrixTypeIntDual (NPHMX);
  ntwmod = MatrixTypeIntHost (NPHMX);
  ntwsys = MatrixTypeIntHost (NPHMX);
  isectw = MatrixTypeIntHost (NSYSMX,NPHMX);
  icryst = MatrixTypeStringHost (NPHMX);

  igas = MatrixTypeIntDual (NPHMX);
  iJ2 = MatrixTypeIntDual (NPHMX);

  cc = MatrixTypeRealHost (3,3,3,3,NPHMX);
  c0 = MatrixTypeRealDual (3,3,3,3);
  s0 = MatrixTypeRealDual (3,3,3,3);
  c066 = MatrixTypeRealDual (6,6);

  scauav1 = MatrixTypeRealHost (3,3);

  xk_gb = MatrixTypeRealDual (npts1);
  yk_gb = MatrixTypeRealDual (npts2);
  zk_gb = MatrixTypeRealDual (npts3);

  sg = MatrixTypeRealDual (3, 3, npts1, npts2, npts3);
  disgrad = MatrixTypeRealDual (3, 3, npts1, npts2, npts3);
  velgrad = MatrixTypeRealDual (3, 3, npts1, npts2, npts3);
  edotp = MatrixTypeRealDual (3, 3, npts1, npts2, npts3);
  cg66 = MatrixTypeRealDual (6, 6, npts1, npts2, npts3);
  ept = MatrixTypeRealDual (3, 3, npts1, npts2, npts3);
  ag = MatrixTypeRealDual (3, 3, npts1, npts2, npts3);
  crss = MatrixTypeRealDual (NSYSMX, 2, npts1, npts2, npts3);
  sch = MatrixTypeRealDual (5, NSYSMX, npts1, npts2, npts3);
  sgt = MatrixTypeRealDual (3, 3, npts1, npts2, npts3);
  defgradp = MatrixTypeRealDual (3, 3, npts1, npts2, npts3);
  defgrade = MatrixTypeRealDual (3, 3, npts1, npts2, npts3);
  defgrad = MatrixTypeRealDual (3, 3, npts1, npts2, npts3);
  defgradinv = MatrixTypeRealDual (3, 3, npts1, npts2, npts3);
  defgradini = MatrixTypeRealDual (3, 3, npts1, npts2, npts3);
  wgtc = MatrixTypeRealDual (npts1, npts2, npts3);
  detF = MatrixTypeRealDual (npts1, npts2, npts3);
  sgPK1 = MatrixTypeRealDual (3, 3, npts1, npts2, npts3);
  c066mod = MatrixTypeRealDual (6, 6, npts1, npts2, npts3);
  velgradref = MatrixTypeRealDual (3, 3, npts1, npts2, npts3);
  xnode = MatrixTypeRealDual (3, npts1 + 1, npts2 + 1, npts3 + 1);
  sigma0gr = MatrixTypeRealDual (npts1, npts2, npts3);
#ifdef NON_SCHMID_EFFECTS
  schnon = MatrixTypeRealDual (5, NSYSMX, npts1, npts2, npts3);
#endif

  gamdot = MatrixTypeRealDual (NSYSMX, npts1, npts2, npts3);
  gacumgr = MatrixTypeRealDual (npts1, npts2, npts3);
  //trialtau = MatrixTypeRealDual (NSYSMX, 2, npts1, npts2, npts3);
  xkin = MatrixTypeRealDual (NSYSMX, npts1, npts2, npts3);

  ph_array = MatrixTypeRealHost (npts1, npts2, npts3);
  th_array = MatrixTypeRealHost (npts1, npts2, npts3);
  om_array = MatrixTypeRealHost (npts1, npts2, npts3);
  jphase = MatrixTypeIntDual (npts1, npts2, npts3);
  jgrain = MatrixTypeIntHost (npts1, npts2, npts3);

  work = MatrixTypeRealDual (3, 3, npts1, npts2, npts3);
  workim = MatrixTypeRealDual (3, 3, npts1_cmplx, npts2_cmplx, npts3_cmplx);
  data = MatrixTypeRealDual (npts1, npts2, npts3);
  data_cmplx = MatrixTypeRealDual (2, npts1_cmplx, npts2_cmplx, npts3_cmplx);

  epav = MatrixTypeRealHost (3,3);
  edotpav = MatrixTypeRealHost (3,3);
  velgradmacroactual = MatrixTypeRealHost (3,3);
  disgradmacrot = MatrixTypeRealHost (3,3);
  velgradmacro = MatrixTypeRealDual (3,3);

  // for fierro linking
  M66 = MatrixTypeRealHost (6,6);
  edotp_avg = MatrixTypeRealHost (3,3);
  dedotp66_avg = MatrixTypeRealHost (6,6);
  cg66_avg = MatrixTypeRealHost (6,6);
  sg66_avg = MatrixTypeRealHost (6,6);
  udotAcc = MatrixTypeRealHost (3,3);

  //... For file management
  mpi_io_real_t = std::make_shared<Manager_MPI_IO<real_t,MPI_ORDER_FORTRAN>>
    (mpi_comm, fft->globalRealBoxSize, fft->localRealBoxSizes[my_rank], fft->localRealBoxes[my_rank].low, 3, 0);

  mpi_io_int = std::make_shared<Manager_MPI_IO<int,MPI_ORDER_FORTRAN>>
    (mpi_comm, fft->globalRealBoxSize, fft->localRealBoxSizes[my_rank], fft->localRealBoxes[my_rank].low, 3, 0);

  micro_writer = std::make_shared<MicroOutputWriter<MPI_ORDER_FORTRAN>>
    (mpi_comm, hdf5_filename.c_str(), 3, fft->globalRealBoxSize.data(), fft->localRealBoxSizes[my_rank].data(),
     fft->localRealBoxes[my_rank].low.data());
  //...


  // Device memory space
  rss = MatrixTypeRealDevice (NSYSMX, npts1, npts2, npts3);
  rss1 = MatrixTypeRealDevice (NSYSMX, npts1, npts2, npts3);
  rss2 = MatrixTypeRealDevice (NSYSMX, npts1, npts2, npts3);

  set_some_voxels_arrays_to_zero();
  return;
}
