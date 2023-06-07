#include "evpfft.h"

#ifndef BUILD_EVPFFT_FIERRO
int main(int argc, char *argv[])
{

  Kokkos::initialize(argc, argv);
  { // kokkos scope

    auto begin = std::chrono::high_resolution_clock::now();


    CommandLineArgs cmd;
    cmd.parse_command_line(argc, argv);

    // EVPFFT
    EVPFFT evpfft(cmd);

    // Do calculation
    for (evpfft.imicro = 1; evpfft.imicro <= evpfft.nsteps; evpfft.imicro++) {
      evpfft.evolve();
    }


    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf("Elapsed time: %12.4E seconds\n", elapsed_time.count() * 1e-9);

  } // end of kokkos scope
  Kokkos::finalize();

  return 0;
}
#endif




#if BUILD_EVPFFT_FIERRO
void EVPFFT::solve(real_t* vel_grad, real_t* stress, real_t dt, size_t cycle, size_t elem_gid)
{
#if 0
//#ifndef NDEBUG
    feenableexcept (FE_DIVBYZERO); 
    feenableexcept (FE_INVALID);
    feenableexcept (FE_OVERFLOW);
//#endif 
#endif

  // nsteps should be set to 1 either in the input file or here
  nsteps = 1;

  // copy vel_grad into udot
  for (int i = 0; i < 9; i++) {
    udot.host_pointer()[i] = vel_grad[i];
  } 
  // update device
  udot.update_device();

  // calculate L2 norm of udot
  real_t ddnorm = 0.0;
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      ddnorm += udot.host(i,j)*udot.host(i,j);
    }
  }
  ddnorm = sqrt(ddnorm);


  if (ddnorm > 1.0e-16) {

    //if (fierro_cycle%1000 == 0) print_vel_grad();
    
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

  } // end if ddnorm

  // copy scauav into stress
  for (int i = 0; i < 9; i++) {
    stress[i] = scauav.pointer()[i];
  }

}
#endif
