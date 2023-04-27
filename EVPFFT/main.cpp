#include "evpfft.h"

#if 0
int main(int argc, char *argv[])
{

  Kokkos::initialize(argc, argv);
  { // kokkos scope

    auto begin = std::chrono::high_resolution_clock::now();


    CommandLineArgs cmd;
    cmd.nn = {16,16,16};
    cmd.input_filename = "fft_library.in";
    cmd.micro_filetype = 0;
    cmd.check_cmd_args();
 
    // EVPFFT
    EVPFFT evpfft(cmd);

    // vel_grad
    MatrixTypeRealHost vel_grad (3,3);
    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
        vel_grad(i,j) = 0.0;
      }
    }
    vel_grad(3,3) = 1.0;

    // dt
    real_t dt = 0.00005;

    // stress
    MatrixTypeRealHost stress (3,3);

    size_t max_cycle = 30;
    for (size_t cycle = 0; cycle < max_cycle; cycle++) {
      // Do calculation
      evpfft.solve(vel_grad.pointer(), stress.pointer(), dt, cycle, 0);
    }


    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf("Elapsed time: %12.4E seconds\n", elapsed_time.count() * 1e-9);

  } // end of kokkos scope
  Kokkos::finalize();

  return 0;
}
#endif
