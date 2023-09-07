#include <stdio.h>
#include "matar.h"
#include <Kokkos_Random.hpp> // for Kokkos random number generator

#define SEED  5374857

// Kokkos provides two random number generator pools one for 64bit states and one for 1024 bit states.
// Choose one.
//using gen_t = Kokkos::Random_XorShift64_Pool<DefaultExecSpace>;
using gen_t = Kokkos::Random_XorShift1024_Pool<DefaultExecSpace>;

int main ()
{
  Kokkos::initialize();
  { // kokkos scope

  // Seed random number generator
  gen_t rand_pool(SEED);

  // DCArrayKokkos type to store the random numbers generated on the device
  // and print out on the host
  const int N = 100;
  DCArrayKokkos<int> arr(N);

  // Generate random numbers
  FOR_ALL(i, 0, N, {

    // Get a random number state from the pool for the active thread
    gen_t::generator_type rand_gen = rand_pool.get_state();
 
    // rand_gen.rand() generates integers from (0,MAX_RAND]
    // rand_gen.rand(END) generates integers from (0,END]
    // rand_gen.rand(START, END) generates integers from (START,END]
    // Note, frand() or drand() can be used in place of rand() to generate floats and
    // doubles, respectively. Please check out Kokkos_Random.hpp for all the other type of
    // scalars that are supported.
    
    // generate random numbers in the range (0,10]
    arr(i) = rand_gen.rand(10);

    // Give the state back, which will allow another thread to acquire it
    rand_pool.free_state(rand_gen);
  }); // end FOR_ALL

  // update host
  arr.update_host();

  for (int i = 0; i < N; i++) {
    printf(" %d", arr.host(i));
  }
  printf("\n");

  } // end kokkos scope
  Kokkos::finalize();

  return 0;
}
