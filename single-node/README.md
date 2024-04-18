# Fierro
**Fierro** (LANL code number C21030) is a modern C++ code designed to simulate quasi-static solid mechanics problems and transient, compressible material dynamic problems with Lagrangian methods, which have meshes with constant mass elements that move with the material, or with Eulerian methods, which have stationary meshes.  **Fierro** is designed to aid a) material model research that has historically been done using commercial implicit and explicit finite element codes, b) numerical methods research, and c) computer science research.  The linear Lagrangian finite element methods in **Fierro** supports user developed material models.  **Fierro** is built on the **ELEMENTS** library that supports a diverse suite of element types, including high-order elements, and quadrature rules. The mesh class within the **ELEMENTS** library is designed for efficient calculations on unstructured meshes and to minimize memory usage.  **Fierro** is designed to readily accommodate a range of numerical methods including continuous finite element, finite volume, and discontinuous Galerkin methods.  **Fierro** is designed to support explicit and implicit time integration methods as well as implicit optimization methods.  


## Computer implementation
**Fierro** is implemented in C++ following modern programming practices.  **Fierro** leverages the unique features of the **ELEMENTS** library, as such, this code serves as an example for solving a system of partial differential equations using the mesh class and geometric functions within the **ELEMENTS** library.  **Fierro** registers state at material points within the element, registers polynomial fields in the element, and registers kinematic variables at element vertices.  The routines for the state are designed for highly efficient computations and to minimize memory usage.  The loops are written to aid fine-grained parallelization and to allow vectorization. **Fierro** is a light-weight software application, and cleanly written following modern programming practices, so it useful for researching computer science based technologies for software performances portability.  

## Spatial discretization methods 
**Fierro** has an established conservative low-order Lagrangian finite element method, a low-order Lagrangian cell-centered finite volume method, and an arbitrary-order Lagrangian Discontinuous Galerkin method for solving the governing equations (e.g., mass, momentum, and energy evolution equations) for compressible material dynamics using unstructured hexahedral meshes.  These methods are combined with a multidirectional approximate Riemann solver (MARS) for improved accuracy on smooth flows and stable solutions near velocity discontinuities and large gradients in a flow. **Fierro** is designed for both low and high-order Lagrangian methods research and development but other types of numerical methods can be readily added to the code.  Likewise, other high-order methods can be studied within the code because it is built upon the **ELEMENTS** library that supports high-order elements and high-order quadrature rules.  Numerical methods are being added to **Fierro** to simulate quasi-static solid mechanics problems.  Likewise, direct Eulerian hydrodynamic methods can be investigated using **Fierro**.

## Temporal discretization methods 
**Fierro** supports a range of published multi-step time integration methods. The code has an explicit multi-step Runge Kutta time integration method. Implicit time integration methods can be implemented in **Fierro**.

## Material models  
The classical ideal gas model is the only material model implemented in the code, and it is useful for verification tests of the software and simple gas dynamic simulations. The Lagrangian finite element methods for explicit material dynamics have an interface to user developed material models. The interface is to enable **Fierro** to be used for model research and development that has historically been done with commercial explicit finite element codes. 


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
Building the code from source allows you to compile with more targeted hardware optimizations that could offer a potentially faster executable. Explicit Lagrangian codes are being added to the repository that are written using MATAR+Kokkos and run with fine-grained parallellism on multi-core CPUs and GPUs.  Build scripts are provided in,  
```
Fierro/single-node/scripts/
```

## Building Fierro with Anaconda
It is advised to use Anaconda package manager to build Fierro as follows:

1. Create an environment and activate:
```
conda create --name fierroEnv
conda activate fierroEnv
```

2. Install needed packages:
```
conda install cxx-compiler -c conda-forge
conda install cmake
```

3. Run the build script as:
```
source build-it.sh --help
```

Which outputs:

```
Usage: source build_evpfft.sh [OPTION]
Valid options:
  --solver=<1DSGH|SGH>. Default is 'SGH'
  --kokkos_build_type=<serial|openmp|pthreads|cuda|hip>. Default is 'serial'
  --build_action=<full-app|set-env|install-kokkos|fierro>. Default is 'full-app'
  --machine=<darwin|chicoma|linux|mac>. Default is 'linux'
  --build_cores=<Integers greater than 0>. Default is set 1
  --help: Display this help message
```

The build-it.sh script sources the other scripts in the folder.  A range of options are provided for computer architectures and machines; however, they might not be correctly configured for the user's hardware or computing environment. If the scripts fail to build the Lagrangian code, then carefully review the modules, environment variables used, and the computer architecture settings.  Taking the defaults as the example, the executable will be located in 
```
Fierro/single-node/build-SGH-serial/bin
```


### Updating submodules
The ELEMENTS library and MATAR library can be updated to the newest release using
```
git submodule update --remote --merge
```







