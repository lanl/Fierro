# Fierro
**Fierro** (LANL code number C21030) is a modern C++ code designed to simulate quasi-static solid mechanics problems and transient, compressible material dynamic problems with Lagrangian methods, which have meshes with constant mass elements that move with the material, or with Eulerian methods, which have stationary meshes.  **Fierro** is designed to aid a) material model research that has historically been done using commercial implicit and explicit finite element codes, b) numerical methods research, and c) computer science research.  The linear Lagrangian finite element methods in **Fierro** supports user developed material models.  **Fierro** is built on the **ELEMENTS** library that supports a diverse suite of element types, including high-order elements, and quadrature rules. The mesh class within the **ELEMENTS** library is designed for efficient calculations on unstructured meshes and to minimize memory usage.  **Fierro** is designed to readily accommodate a range of numerical methods including continuous finite element, finite volume, and discontinuous Galerkin methods.  **Fierro** is designed to support explicit and implicit time integration methods as well as implicit optimization methods.  


## Computer implementation
**Fierro** is implemented in C++ following modern programming practices.  **Fierro** leverages the unique features of the **ELEMENTS** library, as such, this code serves as an example for solving a system of partial differential equations using the mesh class and geometric functions within the **ELEMENTS** library.  **Fierro** registers state at material points within the element, registers polynomial fields in the element, and registers kinematic variables at element vertices.  The routines for the state are designed for highly efficient computations and to minimize memory usage.  The loops are written to aid fine-grained parallelization and to allow vectorization. **Fierro** is a light-weight software application, and cleanly written following modern programming practices, so it useful for researching computer science based technologies for software performances portability.  

## Spatial discretization methods 
**Fierro** has an established conservative low-order Lagrangian finite element method, a low-order Lagrangian cell-centered finite volume method, and an arbitrary-order Lagrangian Discontinuous Galerkin method for solving the governing equations (e.g., mass, momentum, and energy evolution equations) for compressible material dynamics using unstructured hexahedral meshes.  These methods are combined with a multidirectional approximate Riemann solver (MARS) for improved accuracy on smooth flows and stable solutions near velocity discontinuities and large gradients in a flow. **Fierro** is designed for both low and high-order Lagrangian methods research and development but other types of numerical methods can be readily added to the code.  Likewise, other high-order methods can be studied within the code because it is built upon the **ELEMENTS** library that supports high-order elements and high-order quadrature rules.  Numerical methods are being added to **Fierro** to simulate quasi-static solid mechanics problems.  Likewise, direct Eulerian hydrodynamic methods can be investigated using **Fierro**.

## Temporal discretization methods 
**Fierro** supports a range of published multi-step time integration methods. The code has an explicit multi-step Runge Kutta time integration method. Implicit time integration methods can be implemented in **Fierro**.

# Usage
## Anaconda
The recommended way to use **Fierro** is through the provided Anaconda package and command line utility. To use the anaconda package, follow the steps for your platform to install [anaconda](https://docs.anaconda.com/free/anaconda/install/index.html)/[miniconda](https://docs.conda.io/en/latest/miniconda.html)/[mamba](https://mamba.readthedocs.io/en/latest/installation.html). Then run the following command in your desired Anaconda environment:

```
conda install fierro-cpu -c kwelsh-lanl -c conda-forge
```

This will give you access to the **Fierro** command line interface. You can run the following to check that the package was installed correctly:
```
fierro -h
```

## Material models  
The classical ideal gas model is the only material model implemented in the code, and it is useful for verification tests of the software and simple gas dynamic simulations. The linear Lagrangian finite element methods for explicit material dynamics have an interface to user developed material models. The interface is to enable **Fierro** to be used for model research and development that has historically been done with commercial explicit finite element codes. 

To include your own custom material models, you need to implement them under `Fierro/Parallel-Solvers/User-Material-Interface` and re-build the project.
Steps:
1. Create an anaconda environment for your build
2. Install **Fierro** anaconda dependencies with `conda install elements fierro-trilinos-cpu elements -c kwelsh-lanl`
3. [Clone the code](#cloning-the-code)
4. [Build the code](#building-the-code)

Now, if all went correctly, you should be able to run your custom **Fierro** build by calling the executable located at `Fierro/build/bin/fierro`. The executables can be installed into your system directories with the `make install` command as well.

# Cloning the code
If the user has set up ssh keys with GitHub, type
```
git clone --recursive ssh://git@github.com/lanl/Fierro.git
```
The code can also be cloned using
```
git clone --recursive https://github.com/lanl/Fierro.git
```

# Building the code
Building the code from source allows you to compile with more targeted hardware optimizations that could offer a potentially faster executable. 
To build it yourself, run the following from the root directory. The native CPU architecture will automatically be taken into account. 
```
mkdir build
cd build
cmake .. -DBUILD_PARALLEL_EXPLICIT_SOLVER=ON -DBUILD_IMPLICIT_SOLVER=ON
make -j
```

GPU hardware will be leveraged according to the distribution of Trilinos that **Fierro** is built against.
You are welcome to only compile one solver or the other, and the one(s) that you did compile will be available through the CLI.

## Building dependencies
**Fierro** depends on both ELEMENTS and Trilinos. If you are building **Fierro** for hardware optimizations, you should also also build ELEMENTS from source. ELEMENTS is included in the **Fierro** source distribution and building ELEMENTS is enabled by default when building **Fierro**.

As for Trilinos, we recommend installing the Anaconda package for the desired build into a new Anaconda environment to satisfy **Fierro**'s dependency rather than building it from source. If you do wish to build it from source, however, sample build scripts for Trilinos can be found in `Fierro/Trilinos-Build-Scripts`. Build scripts for all Anaconda packages can be found in `Fierro/Anaconda-Packages/`.

## Alternative Build Workflows
In addition to the primary build workflow described above, there are build scripts for a variety of alternative workflows. These scripts can be found under `Fierro/scripts`.
### Building the explicit Lagrangian methods with Kokkos
Explicit Lagrangian codes are being added to the repository that are written using MATAR+Kokkos and run with fine-grained parallellism on multi-core CPUs and GPUs.  Build scripts are provided for each Lagrangian code, and those scripts follow those used in the [MATAR](https://github.com/lanl/MATAR/) GitHub repository. The scripts to build the Lagrangian codes (that use MATAR+Kokkos) are in the scripts folder.  The user must update the modules loaded by the build scripts (for the compiler etc.), and then type
```
source build-it.sh
```
The build-it.sh script sources the other scripts in the folder.  The compiled code will be in a folder (named after the explicit Lagrangian method) in the Fierro directory.  A range of scripts are provided for many architectures; however, they might not be correctly configured for the user's hardware.  The CPU architecture information needs to be listed if running with the Kokkos serial, OpenMP, and pthreads backends; GPU architecture information must be listed if using a Kokkos GPU backend. We refer the user to Kokkos compiling page to see the large list of compilation options,
``` 
https://github.com/kokkos/kokkos/wiki/Compiling
```
If the scripts fail to build a Lagrangian code, then carefully review the modules used and the computer architecture settings.  A more lenghtly discussion of the build scripts is provided in the MATAR GitHub repository. 


### Updating submodules
The ELEMENTS library and MATAR library can be updated to the newest release using
```
git submodule update --remote --merge
```







