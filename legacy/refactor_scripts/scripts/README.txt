all scripts can be built fully by entering the appropriate directory and doing
source build-it.sh <arg1 arg2 ...>
build-it.sh can take up to 3 arguments (minimum of 2)
source build-it.sh <environment type> <parallelism> <build directory name (optional)>
    environment has two options: 'hpc' or 'macos'
        hpc: builds by loading modules and can perform parallel builds (make -j)
        macos: does not load anything externally and expects the environment to be set up on your mac. Additionally, the builds will all be serial (make)
    parallelism has four options: 'cuda', 'hip', 'openmp', 'none'
    Note - all builds use Kokkos. The 'none' option will still use Kokkos, however only utilizing the Kokkos Serial build
        cuda: loads cuda module and a working gcc module pairing (these can be changed, more info later)
        hip: loads hip module and a working clang module pairing
        openmp: loads gcc module and sets openmp environment variables
        none: loads gcc module
    build directory is an optional argument which will create the name of the build directory

All other scripts will be called with the appropriate arguments as a result of running build-it.

If you need to simply rebuild the app and not get a new kokkos installation, simply
source cmake_build.sh <args> 
with the same arguments you would with build-it.sh

If you log onto a machine for the first time (or get a new allocation) you will need to run
source setup-env.sh <args>
with the same arguments you would with build-it.sh
 
