# Large Strain Elasto-Viscoplastic Fast Fourier Transform-based (LS-EVPFFT) Micromechanical Solver

LS-EVPFFT is a parallel implementation of an elasto-viscoplastic fast Fourier transform-based (EVPFFT) micromechanical solver to facilitate computationally efficient crystal plasticit modeling of polycrystalline materials. Emphasis is placed on performance portability to allow for a single implementation of LS-EVPFFT to run optimally on variety of homogeneous architectures that are comprised of multi-core central processing units(CPUs), and on heterogeneous computer architectures comprised of multi-core CPUs and graphics processing units (GPUs) from different vendors. To achieve this, we utilize MATAR, a C++ software library that allows the straightforward creation and usage of multidimensional dense or sparse matrix and array data structures that are also portable across disparate architectures using Kokkos, a performance portable library. Also, we utilize the message passing interface (MPI) to distribute the computational workload amongst the processors. The heFFTe (Highly Efficient FFT for Exascale) library is used to facilitate the performance portability of the fast Fourier transforms (FFTs) computation.

# Building LS-EVPFFT as a standalone program

LS-EVPFFT depends on the following to build:

1. MPI
2. [MATAR](https://github.com/lanl/MATAR)
3. [Kokkos](https://github.com/kokkos/kokkos)
4. [HeFFTe](https://github.com/icl-utk-edu/heffte) (FFTW,CUFFT,ROCFFT)
5. [HDF5](https://www.hdfgroup.org/solutions/hdf5/)

To make it easy to build LS-EVPFFT, we have included `scripts/build-scripts/build_ls-evpfft.sh` which when executed will download and install all required dependencies and LS-EVPFFT. It is important to note that LS-EVPFFT depends on three main libraries which can be time consuming to download from scratch namely, MPI, FFTW (MPI version), and HDF5 (MPI version). Therefore the user is advised to download these packages themselves either using `sudo apt install` of Anaconda package manager.

## Building LS-EVPFFT with Anaconda
It is advised to use Anaconda package manager to build LS-EVPFFT as follows:

1. Create an environment and activate:
```
conda create --name ls-evpfftEnv
conda activate ls-evpfftEnv
```

2. Install needed packages:
```
conda install cxx-compiler -c conda-forge
conda install cmake
conda install "fftw=*=mpi_openmpi*" -c conda-forge
conda install "hdf5=*=mpi_openmpi*" -c conda-forge
conda install openmpi-mpicxx -c conda-forge 
```

3. Run the build script as:
```
source build_ls-evpfft.sh --help
```

Which outputs:

```
Usage: source build_ls-evpfft.sh [OPTION]
Required arguments:
  --heffte_build_type=<fftw|cufft|rocfft>
  --kokkos_build_type=<serial|openmp|pthreads|cuda|hip>

Optional arguments:
  --build_fftw: builds fftw from scratch
  --build_hdf5: builds hdf5 from scratch
  --machine=<darwin|chicoma|linux|mac> (default: none)
  --num_jobs=<number>: Number of jobs for 'make' (default: 1, on Mac use 1)
  --help: Display this help message
```

To build LS-EVPFFT you would need to provide both the `--heffte_build_type` and `--kokkos_build_type` options. The command below build LS-EVPFFT using FFTW and Serial version of Kokkos:

```
source build_ls-evpfft.sh --heffte_build_type=fftw --kokkos_build_type=serial
```

This will build LS-EVPFFT in the folder `ls-evpfft_{fftw}_{serial}`. The binary, `ls-evpfft` is found in that folder.

## Building LS-EVPFFT without Anaconda
Install HDF5 (MPI version) and FFTW (MPI version) with `sudo apt install` which will install the libraries in the default location. Then run the build script as above.

If you would like to build HDF5 and FFTW from scratch then use `--build_hdf5` and `--build_fftw` options of the build script.

If you do not want to build HDF5 and FFTW from scratch because you already have a version on your system remove the `--build_hdf5` and `--build_fftw` option from the command. If your installed HDF5 and FFTW are located in a directory other than the default linux directories for libraries, then specify the following before running the build script.

```
export HDF5_ROOT=/path/to/where/HDF5/lib/and/include/folders/are/located
export FFTW_ROOT=/path/to/where/FFTW/lib/and/include/folders/are/located
```

# Using LS-EVPFFT as a standalone program

To get help on how to run `ls-evpfft` use the `-h` or `--help` opton:

```
./ls-evpfft --help
```

Which gives the following output:

```
Required arguments:
  -f, --infile   <path>     input file (full path)
...
```

We have provided example input files in the `example_input_files` directory. In summary, LS-EVPFFT needs these files:

1. input file (see `example_input_files/example_ls-evpfft_standalone_inputfile.txt`)
2. plastic parameters file (see `example_input_files/example_plastic_parameters.txt`)
3. elastic parameters file (see `example_input_files/example_elastic_parameters.txt`)
4. microstructure RVE file (see `example_input_files/random_microstructure_8x8x8.txt`)

With all these files available you can now run LS-EVPFFT as:

```
mpirun -np 4 ls-evpfft --infile=example_ls-evpfft_standalone_inputfile.txt
```

# Using LS-EVPFFT with Fierro

To use LS-EVPFFT as a user material model for Fierro, build Fierro with the following cmake flags `-DBUILD_LS_EVPFFT_FIERRO=ON`, choose the FFT backend with `-D<USE_FFTW,USE_CUFFT,USE_ROCFFT>=ON`, and use `-DABSOLUTE_NO_OUTPUT=ON` to specify that LS-EVPFFT does not provide output per element which can be cumbersome and memory intensive. It is best to modify LS-EVPFFT to output for a few elements of interst.

See `example_input_files/taylor_anvil.yaml` for the yaml file setup when using EVPFFT with Fierro.

Also, input files for each material should be defined as `evpfft1.in`, `evpfft2.in`, etc. and should be provided in the directory Fierro is being run from.

# Additional Information

These information are for users interested in more customization of the build process. You can add to the `cmake` options in `build_ls-evpfft.sh` as needed.

Scripts for building each dependency are available in the `install-scripts` folder and can be run individually.

Note that MATAR is used as an include header file.

If you already have an installation of HDF5 (parallel version) on you system you can use your installation, just comment out the part of the `build_ls-evpfft.sh` that calls `install_hdf5.sh`.

During the installation of HeFFTe libray we assume the required FFT backend (FFTW, CUFFT, ROCM) is installed on the system. If this is not the case then look at the [HeFFTe installation page](https://github.com/icl-utk-edu/heffte) for more information.

The LS-EVPFFT `CMakeLists.txt` uses the following default values to build LS-EVPFFT:

```
  TWO_SIGN_SLIP_SYSTEMS: OFF
  NON_SCHMID_EFFECTS: OFF
  ABSOLUTE_NO_OUTPUT: OFF
  ENABLE_PROFILE: OFF
```

To change these default options include the `-D OPTION=<value>` in the `cmake` option, E.g. `-D ABSOLUTE_NO_OUTPUT=ON` in the `build_ls-evpfft.sh`.

# Using LS-EVPFFT for Lattice Structure Homogenization

Example for input files needed to run LS-EVPFFT for lattice structure homogenization is shown in `example_input_files/lattice_input_files`. In that file you will see how to set up evpft input file, elastic and plastic parameter files.

Provide a vtk file type that contains information about which grid point is solid (1) or void (0), example is shown in `example_input_files/lattice_input_files/void_in_solid.vtk`.

Run LS-EVPFFT as:
```
  mpirun -n 1 ls-evpfft -f tension_11.txt -m 2
```
the `-m 2` option tells LS-EVPFFT to use the vtk lattice file microstructure file type.

