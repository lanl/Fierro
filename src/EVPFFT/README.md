# Elasto-Viscoplastic Fast Fourier Transform-based (EVPFFT) Micromechanical Solver

EVPFFT is a parallel implementation of an elasto-viscoplastic fast Fourier transform-based (EVPFFT) micromechanical solver to facilitate computationally efficient crystal plasticit modeling of polycrystalline materials. Emphasis is placed on performance portability to allow for a single implementation of EVPFFT to run optimally on variety of homogeneous architectures that are comprised of multi-core central processing units(CPUs), and on heterogeneous computer architectures comprised of multi-core CPUs and graphics processing units (GPUs) from different vendors. To achieve this, we utilize MATAR, a C++ software library that allows the straightforward creation and usage of multidimensional dense or sparse matrix and array data structures that are also portable across disparate architectures using Kokkos, a performance portable library. Also, we utilize the message passing interface (MPI) to distribute the computational workload amongst the processors. The heFFTe (Highly Efficient FFT for Exascale) library is used to facilitate the performance portability of the fast Fourier transforms (FFTs) computation.

# Building EVPFFT as a standalone program

EVPFFT depends on the following to build:

1. MPI
2. [MATAR](https://github.com/lanl/MATAR)
3. [Kokkos](https://github.com/kokkos/kokkos)
4. [HeFFTe](https://github.com/icl-utk-edu/heffte) (FFTW,CUFFT,ROCFFT)
5. [HDF5](https://www.hdfgroup.org/solutions/hdf5/)

To make it easy to build EVPFFT, we have included `scripts/build-scripts/build_evpfft.sh` which when executed will download and install all required dependencies and EVPFFT. The script assumes that the user already has MPI installed.

Run the build script as:

```
source build_evpfft.sh --help
```

Which outputs:

```
Usage: source build_evpfft.sh [OPTION]
Valid options:
  --heffte_build_type=<fftw|cufft|rocfft>
  --kokkos_build_type=<serial|cuda|hip>
  --help: Display this help message
```

To build EVPFFT you would need to provide both the `--heffte_build_type` and `--kokkos_build_type` options. The command below build EVPFFT using FFTW and Serial version of Kokkos:

```
source build_evpfft.sh --heffte_build_type=fftw --kokkos_build_type=serial
```

This will build EVPFFT in the folder `evpfft_heffte_{fftw}_kokkos_{serial}`. The binary, `evpfft` is found in that folder

# Using EVPFFT as a standalone program

To get help on how to run `evpfft` use the `-h` or `--help` opton:

```
./evpfft --help
```

Which gives the following output:

```
Required arguments:
  -f, --infile   <path>     input file (full path)
...
```

We have provided example input files in the `example_input_files` directory. In summary, EVPFFT needs these files:

1. input file (see `example_input_files/example_evpfft_standalone_inputfile.txt`)
2. plastic parameters file (see `example_input_files/example_plastic_parameters.txt`)
3. elastic parameters file (see `example_input_files/example_elastic_parameters.txt`)
4. microstructure RVE file (see `example_input_files/random_microstructure_8x8x8.txt`)

With all these files available you can now run EVPFFT as:

```
mpirun -np 4 evpfft --infile=example_evpfft_standalone_inputfile.txt
```

# Using EVPFFT with Fierro

To use EVPFFT as a user material model for Fierro, build Fierro with the following cmake flags `-DBUILD_EVPFFT_FIERRO=ON`, choose the FFT backend with `-D<USE_FFTW,USE_CUFFT,USE_ROCFFT>=ON`, and use `-DABSOLUTE_NO_OUTPUT=ON` to specify that EVPFFT does not provide output per element which can be cumbersome and memory intensive. It is best to modify EVPFFT to output for a few elements of interst.

The yaml file for Fierro should include:

```
material_options:
  - eos_model: user_eos_model
    eos_run_location: device
    strength_model: user_strength_model
    strength_type: hypo
    strength_run_location: host
    global_vars:
      - 8 #N1 (EVPFFT RVE X dimension)
      - 8 #N2 (EVPFFT RVE Y dimension)
      - 8 #N3 (EVPFFT RVE Z dimension)
      - 0.0001 #udotAccTh (EVPFFT accumulated strain threshold before full solve)
      - 3000000 #[mm/s] (EVPFFT material sound speed)
```

Also, input files for each material should be defined as `evpfft1.in`, `evpfft2.in`, etc. and should be provided in the directory Fierro is being run from.

# Additional Information

These information are for users interested in more customization of the build process. You can add to the `cmake` options in `build_evpfft.sh` as needed.

Scripts for building each dependency are available in the `install-scripts` folder and can be run individually.

Note that MATAR is used as an include header file.

If you already have an installation of HDF5 (parallel version) on you system you can use your installation, just comment out the part of the `build_evpfft.sh` that calls `install_hdf5.sh`.

During the installation of HeFFTe libray we assume the required FFT backend (FFTW, CUFFT, ROCM) is installed on the system. If this is not the case then look at the [HeFFTe installation page](https://github.com/icl-utk-edu/heffte) for more information.

The EVPFFT `CMakeLists.txt` uses the following default values to build EVPFFT:

```
  TWO_SIGN_SLIP_SYSTEMS: OFF
  NON_SCHMID_EFFECTS: OFF
  ABSOLUTE_NO_OUTPUT: OFF
  ENABLE_PROFILE: OFF
```

To change these default options include the `-D OPTION=<value>` in the `cmake` option, E.g. `-D ABSOLUTE_NO_OUTPUT=ON` in the `build_evpfft.sh`.

# Using EVPFFT for Lattice Structure Homogenization

Example for input files needed to run EVPFFT for lattice structure homogenization is shown in `example_input_files/lattice_input_files`. In that file you will see how to set up evpft input file, elastic and plastic parameter files.

Provide a `STRUCTURED_POINTS` vtk file type that contains information about which grid point is solid (1) or void (0), example is shown in `example_input_files/lattice_input_files/void_in_solid.vtk`.

Run EVPFFT as:
```
  mpirun -n 1 -f evpfft_lattice_input.txt -m 2
```
the `-m 2` option tells EVPFFT to use the vtk lattice file microstructure file type.

