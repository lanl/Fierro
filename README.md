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
conda install fierro-cpu -c fierromechanics -c conda-forge
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
2. Install **Fierro** anaconda dependencies with `conda install elements fierro-trilinos-cpu elements -c fierromechanics`
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

GPU hardware will be leveraged according to the distribution of Trilinos that **Fierro** is built against.
You are welcome to only compile one solver or the other, and the one(s) that you did compile will be available through the CLI.

## Building dependencies
**Fierro** depends on both ELEMENTS and Trilinos. If you are building **Fierro** for hardware optimizations, you should also build ELEMENTS from source. ELEMENTS is included in the **Fierro** source distribution and building ELEMENTS is enabled by default when building **Fierro**.

As for Trilinos, we recommend installing the Anaconda package for the desired build into a new Anaconda environment to satisfy **Fierro**'s dependency rather than building it from source. If you do wish to build it from source, however, sample build scripts for Trilinos can be found in `Fierro/Trilinos-Build-Scripts`. Build scripts for all Anaconda packages can be found in `Fierro/Anaconda-Packages/`.

## Alternative Build Workflows
In addition to the primary build workflow described above, there are build scripts for a variety of alternative workflows. These scripts can be found under `Fierro/scripts`.
### Building the explicit and implicit Lagrangian methods with Trilinos+Kokkos
Explicit Lagrangian codes are being added to the repository that are written using MATAR+Kokkos and run with fine-grained parallellism on multi-core CPUs and GPUs.  Build scripts are provided for each Lagrangian code, and those scripts follow those used in the [MATAR](https://github.com/lanl/MATAR/) GitHub repository. The scripts to build the Lagrangian codes (that use MATAR+Kokkos) are in the scripts folder.  The user must update the modules loaded by the build scripts (for the compiler etc.), and then type
Immediate help with all scripts can be had running
```
<script> --help
```
The build-it script can take up to 3 arguments (with a minimum of 2)
```
source build-it.sh --machine=<arg> --kokkos_build_type=<arg> (optional)--build_cores=<arg>
```
machine currently has three options: darwin,linux,mc
```
    darwin: builds by loading modules on the darwin cluster. and can perform parallel builds (make -j). Changes can be made to scripts/[solver]/machines/darwin-setup.sh to reflect your machine.
    linux: does not load any modules and instead sources some paths related to your loaded software. Additionally, the builds can be parallel (make -j <build_cores>)
    macos: does not load any modules and instead sources some paths related to your loaded software. Additionally, the builds will be forced to beserial (make)
For Mac builds, please see the ***BUILD.md*** file for more info
```
kokkos_build_type  has four options: 'cuda', 'hip', 'openmp', 'serial'
***Note*** - all builds use Trilinos with Kokkos. The 'serial' option will utilize the Kokkos serial build
```
    cuda: loads cuda module and a working gcc module pairing
    hip: loads hip module and a working clang module pairing
    openmp: loads gcc module and sets openmp environment variables
    serial: loads gcc module
```
***Note*** - compiler can be changed with the appropriate variables in *setup-env.sh*, the ones provided are simply known to work together
If running on Mac, only the openmp and serial options are enabled

All other scripts will be called with the appropriate arguments as a result of running build-it.

If you need to simply rebuild Fierro without making a new Trilinos installation, simply
```
source cmake_build.sh
```
If you are getting back on to a machine or allocation to continue development, you will need to run
```
source setup-env.sh <same args you passed to build-it>
```
If the scripts fail to build a Lagrangian code, then carefully review the modules used and the computer architecture settings.  
#A more lenghtly discussion of the build scripts is provided in the MATAR GitHub repository. 

For help with [Trilinos](https://github.com/trilinos/Trilinos/wiki/New-Trilinos-Developers)
For help with [Kokkos](https://github.com/kokkos/kokkos/wiki/Compiling) compilation

### Updating submodules
The ELEMENTS library and MATAR library can be updated to the newest release using
```
git submodule update --remote --merge
```

# Contributing to Fierro
As an open source software team, we greatly appreciate any and all community involvement. There are many ways to contribute to Fierro, from tidying up code sections to implementing novel solvers. 
To streamline the integration of community contributions, we follow the following guidelines.

## Writing commit messages

Write your commit messages using these standard prefixes:
- `BUG:` Fix for runtime crash or incorrect result
- `COMP:` Compiler error or warning fix
- `DOC:` Documentation change
- `ENH:` New functionality
- `PERF:` Performance improvement
- `STYLE:` No logic impact (indentation, comments)
- `WIP:` Work In Progress not ready for merge

The commit message should assist in the review process and allow future developers to understand the change.  To that end, please follow these few rules:

1. Try to keep the subject line below 72 characters, ideally 50.
2. Be concise, but honor the change.
3. If the change is uncharacteristically large (or diverse), consider breaking the change into multiple distinct commits.
4. If the commit references a specific issue, link it in the message

[This](https://cbea.ms/git-commit/) is a great post on how to write a good commit message. 

## The PR Process

If you are new to Fierro development and don't have push access to the repository, here are the steps:
1. [Fork and clone](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo) the repository.
2. Create a branch for your work.
3. [Push](https://docs.github.com/en/get-started/using-git/pushing-commits-to-a-remote-repository) the branch to your GitHub fork.
4. Build and test your changes.
5. Update any necessary documentation.
6. Create a [Pull Request](https://github.com/lanl/Fierro/pulls).

This corresponds to the **Fork & Pull Model** described in the [GitHub collaborative development](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/getting-started/about-collaborative-development-models) documentation.

## Integrating a PR

Integrating your contributions to Fierro is relatively straightforward; here is the checklist:
- All test pass
- The changes build with no new compiler warnings/errors
- All feedback has been addressed
- Consensus is reached. This usually means that at least two reviewers approved the changes (or added a `LGTM` comment) and at least one business day passed without anyone objecting. `LGTM` is an acronym for Looks Good to Me.
- If you do NOT have push access, a Fierro core developer will integrate your PR. 

## Benevolent dictators for life

The [benevolent dictators](https://en.wikipedia.org/wiki/Benevolent_dictator_for_life) can integrate changes to keep the platform healthy and help interpret or address conflict related to the contribution guidelines and the platform as a whole.

These currently include:
- [Nathaniel Morgan](https://github.com/nathanielmorgan)
