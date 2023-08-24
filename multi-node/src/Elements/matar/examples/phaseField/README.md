## Phase Field Simulation of Phase Seperation
This folder uses the MATAR library to handle the array data structures needed for the implementation of phase field method for phase separation in a binary alloy. 
The implementation allows for both 2D and 3D simulations.
The Cahn-Hilliard equation in phase field modeling is solved using the semi-implicit Fourier spectral algorithm [[1]](#1). The [FFTW](http://www.fftw.org/download.html) or [cuFFT](https://docs.nvidia.com/cuda/cufft/index.html) libraries are used to perform the Fourier transforms on the CPU or GPU, respectively. Performance on the CPU and different GPU architectures are detailed in `report.pdf`.

## Usage
### Before build
In the `.../test/CMakeLists.txt` file uncomment:
```
include_directories(phaseField/src)
add_subdirectory(phaseField/src)
```
Here, `src` is either `srcMacros` or `srcKokkosVerbose`. Both produce the same results. The difference between `srcMacros` or `srcKokkosVerbose` is to illustrate the use of kokkos parallel region versus MATAR Macros. <br />
If compiling for the CPU with OpenMP backend, ensure that [FFTW](http://www.fftw.org/download.html) is installed. Specify the link directory for `fftw3_threads` and `fftw3` libraries in the `.../phaseField/src/CMakeLists.txt` file: 
```
link_directories("/usr/local/lib")
```
Nothing needs to be done if compiling for the GPU with CUDA backend. The `find_package(CUDA REQUIRED)` command in the `.../phaseField/src/CMakeLists.txt` takes care of that. <br />
Use `-DOUT_OF_PLACE_FFT=1` or `-DIN_PLACE_FFT=1` to experiment with out-of-place or in-place FFT, respectively. <br />
Note that `-DOUT_OF_PLACE_FFT=1` is set by default.

### Simulation parameters
Default simulation parameters are provided in the `.../phaseField/src/sim_parameters.cpp` file and can be changed accordingly. 

### Simulation outputs
The simulation outputs vtk files which can be visualized using [paraview](https://www.paraview.org/). Also, the evolution of the total free energy of the system is outputted in the `total_free_energy.csv` file. 

## References
<a id="1">[1]</a> 
Chen, L.Q. and Shen, J., 1998. Applications of semi-implicit Fourier-spectral method to phase field equations. 
*Computer Physics Communications*, 108(2-3), pp.147-158.
