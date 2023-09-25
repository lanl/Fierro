#include "evpfft.h"
#include "utilities.h"
#include "vm.h"
#include "math_functions.h"
#include "Profiler.h"
#ifndef NDEBUG
  #include <cfenv>
#endif


#ifndef BUILD_EVPFFT_FIERRO
int main(int argc, char *argv[])
{
#ifndef NDEBUG
    // These are used for Floating point exception debuging.
    // However, comment out hdf5 output (`micro_writer`) or compile 
    // EVPFFT with ABSOLUTE_NO_OUTPUT before using these
    // because `H5Screate_simple` function has a Floating point exception bug
    // when run with these on. Although, it does not affect anything in the EVPFFT code.

    //feenableexcept (FE_DIVBYZERO); 
    //feenableexcept (FE_INVALID);
    //feenableexcept (FE_OVERFLOW);
#endif 

  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  { // kokkos scope

    int my_rank = get_mpi_comm_rank(MPI_COMM_WORLD);
    int num_ranks = get_mpi_comm_size(MPI_COMM_WORLD);

    { Profiler profiler("Total");

    // Command line arguments
    CommandLineArgs cmd;
    cmd.parse_command_line(argc, argv);

    // EVPFFT
    EVPFFT evpfft(MPI_COMM_WORLD, cmd);
    evpfft.solve();

    } // end profiler("Total") scope

    Profiler::print(MPI_COMM_WORLD);

  } // end of kokkos scope
  Kokkos::finalize();
  MPI_Finalize();

  return 0;
}
#endif




#if BUILD_EVPFFT_FIERRO
void EVPFFT::solve(real_t* vel_grad, real_t* stress, real_t dt, size_t cycle, size_t elem_gid, real_t udotAccThIn)
{
#ifndef NDEBUG
    feenableexcept (FE_DIVBYZERO); 
    feenableexcept (FE_INVALID);
    feenableexcept (FE_OVERFLOW);
#endif 

  /* All tensors must come in in a F-layout */

  ViewFMatrix vel_grad_view (vel_grad,3,3);
  ViewFMatrix stress_view (stress,3,3);
  double udotAccTh = udotAccThIn;

  // calculate L2 norm of vel_grad
  real_t L2norm = 0.0;
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      L2norm += vel_grad_view(i,j) * vel_grad_view(i,j);
    }
  }
  L2norm = sqrt(L2norm);

  // return if L2norm is too small
  if (L2norm < 1.0e-16) {
    return;
  }

  // Accumulate udot and dt
  dtAcc += dt;
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      udotAcc(i,j) += vel_grad_view(i,j) * dt;
    }
  }
  // Symmetrize udotAcc
  MatrixTypeRealHost udotAccSymm(3,3);
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      udotAccSymm(i,j) = 0.5*(udotAcc(i,j) + udotAcc(j,i));
    }
  }
  // Calculate VonMises of udotAccSymm
  double udotAccVm = vm(udotAccSymm.pointer());

  // calculate strain increament
  MatrixTypeRealHost dstran(3,3);
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      dstran(i,j) = 0.5*(vel_grad_view(i,j) + vel_grad_view(j,i)) * dt;
    }
  }

  // Linear extrapolation
  if (active == true and udotAccVm < udotAccTh) {

    // calculate M66
    for (int ii = 1; ii <= 6; ii++) {
      for (int jj = 1; jj <= 6; jj++) {
        M66(ii,jj) = sg66_avg(ii,jj) + dedotp66_avg(ii,jj) * dt;
      }
    }
    invert_matrix(M66.pointer(), 6);

    MatrixTypeRealHost M3333(3,3,3,3);
    cb.chg_basis_3(M66.pointer(), M3333.pointer(), 3, 6, cb.B_basis_host_pointer());
    MatrixTypeRealHost dstress(3,3);
    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        dstress(ii,jj) = 0.0;
        for (int kk = 1; kk <= 3; kk++) {
          for (int ll = 1; ll <= 3; ll++) {
            dstress(ii,jj) += M3333(ii,jj,kk,ll) * (dstran(kk,ll) - edotp_avg(kk,ll)*dt);
          }
        }
        stress_view(ii,jj) += dstress(ii,jj);
      }
    }
    return;
  }

  // set dt
  dt = dtAcc;

  real_t dtAcc_inv = 1.0 / dtAcc;
  // copy udotAcc into udot
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      udot.host(i,j) = udotAcc(i,j) * dtAcc_inv;
    }
  }
  // update device
  udot.update_device();

  // Reset Acculated variables
  dtAcc = 0.0;
  for (int i = 0; i < 9; i++) udotAcc.pointer()[i] = 0.0;

  // nsteps should be set to 1 either in the input file or here
  nsteps = 1;

  // calculate strain-rate and rotation-rate
  decompose_vel_grad(udot.host_pointer());

  // assign recieved variables to evpfft global variables
  tdot = dt;
  imicro += 1; // imicro starts at 1 while cycle starts at 0 in fierro.
  elem_id = elem_gid;
  fierro_cycle = cycle;

  // do some initialization if first load step
  if (active == false) {
    init_dvm();
    init_evm();
    init_disgradmacro();
    init_sg();
  }
  active = true;

  //--------------------------------
  // EVPFFT evolve
  evolve();

  // check macrostress for NaN
  check_macrostress();

  // copy scauav into stress
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      stress_view(i,j) = scauav(i,j);
    }
  }

  return;
}
#endif
