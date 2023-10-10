# Fierro
These instructions will focus primarily on building **Fierro** on a Mac. Builds on Linux and Linux-based HPC systems are more straightforward, so the example instructions listed in the ***README*** are relevant for those builds.

## Environment Setup
The general steps to create a working Mac build will be
    - Install software packages using Homebrew
    - Setting environment variables to point to the appropriate software
    - Run the approriate build script(s)

### Homebrew Installations
Install [Homebrew](https://brew.sh) on your your computer. Once installed you will need to install several packages using
```
brew install <package name>
```
The minimum packages needed to be installed (using separate install commands) are
```
git
cmake
gcc # If you want to use gcc. This will install gcc-13
llvm # If you want to use llvm and clang. These are the default softwares so they will be needed without change to build scripts
libomp # This is only necessary if you are using llvm/clang. gcc installs it's own OpenMP
openmpi 
```

### Environment Variables
We need to have several environment variables point to these homebrew installations so that our build system knows where to look. These variables can either be set individually, or put into your ***~/.bashrc*** file and sourced. 
The current setup for the build scripts is to add the below lines to your environment automatically. For user that are using a different software stack such as GCC, these path's would need to point to those paths in the script file they are pulled from. This file is
```
scripts/<solver>/machines/mac-env.sh
```
If you are changing the environment paths in your ***~/.bashrc*** then the file could look like this
```
export CC=/opt/homebrew/opt/llvm/bin/clang
export CXX=/opt/homebrew/opt/llvm/bin/clang++
export PATH="/opt/homebrew/opt/llvm/bin:$PATH"
export LDFLAGS="-L/opt/homebrew/opt/llvm/lib"
export CPPFLAGS="-I/opt/homebrew/opt/llvm/include"
export OMPI_CC=${CC}
export OMPI_CXX=${CXX}
```
and then you would source those variables with 
```
source ~/.bashrc
```
For those that noticed, this looks exactly like the ***mac-env.sh*** file mentioned above.
There are two things of note with this example. First, it is using the llvm/clang compiler. This is one of the possible options. If you prefer to use GCC as your base compiler, all of the variables (CC,CXX,PATH,LDFLAGS,CPPFLAGS) would need to be modified to point to the corresonding spot in the gcc directory. The second thing to note with this example is the path that homebrew has installed your software. By default on Apple silicon hardware (M1/M2), it should be the one listed as in the above example
```
/opt/homebrew/opt/<software>
```
Intel-based Macs have homebrew install software in
```
/usr/local/opt/<software>
```
If you explicitally tell homebrew to install them elsewhere, or it does it in a different default locatation on your machine, you will have to have it point to the appropriate path.

### Running the build scripts
Now the assumption is that you have the software installed and the system knows where to look for that software. The code can be built using the scripts in the ***scripts*** directory. There are several options here depending on the problems you'd like to run, but the build process is identical, so we will use the *Explicit* folder for our example. The script can be run from anywhere in the Fierro repository. If we were at the top of the directory would run
```
./scripts/Explicit/build-it.sh --machine=mac --kokkos_build_type=openmp
```
***Note*** you can see more information on the script arguments in the *README.md* or running the script with the argument ***--help***. However this is the standard script you would run on a Mac. You could subistitue 'serial' in for 'openmp', but since the build supports OpenMP that is the suggestion option for better runtime performance.

## Running the example
Assuming a smooth build, you will be placed at the top of the repository again. You can move yourself into the build directory
```
cd build-explicit
```
From there you can run this example with
```
./bin/fierro-parallel-explicit ../src/Parallel-Solvers/Parallel-Explicit/example_simple.yaml
```
***Note*** For this example the mesh file that is passed into *example_simple.yaml* needs to be present in your build directory. 

## Important Notes
Hopefully the above steps result in a smooth build process. When making changes to code that require another cmake configuration, we suggest running the full ***build-it.sh <args>*** script again. This is because it sources the necessaryfiles to setup environment paths. Trilinos will *not* get re-built as a result of this.
The reason for this suggestion is ease of environment. It is possible to simply source the cmake script
```
source cmake_build.sh
```
However you may end up in an unexpected directory since all internal paths were not set. It is possible to
```
source setup-env.sh --machine=mac --kokkos_build_type=openmp
```
before the cmake step. This would set all path's globally, instead of self-contained inside the script run. However any software changes made after this needs to have environment-correcting commands.
Also, it is not advised to source the ***build-it.sh*** script from anywhere as the initial step, because of how Macs handle that source. In bash shells the behavior would be more or less expected. However on others (for instance Mac's default zsh shell) the paths would not be expected.
Therefore, we suggest one of the two options
- run the script with ***./<path to script>/build-it.sh*** from anywhere or
- ***source build-it.sh*** from *inside* the scripts/[solver] directory
