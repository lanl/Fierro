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
The paths are set in the **scripts/machine/mac.sh** script.
There are two things of note. First, it is using the llvm/clang compiler. This is one of the possible options. If you prefer to use GCC as your base compiler, all of the variables (CC,CXX,PATH,LDFLAGS,CPPFLAGS) would need to be modified to point to the corresonding spot in the gcc directory. The second thing to note with this example is the path that homebrew has installed your software. By default on Apple silicon hardware (M1/M2), it should be the one listed as in the above example
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
Now the assumption is that you have the software installed and the system knows where to look for that software. The code can be built as normal now, with the exception that we always need to pass in the machine argument and set it to mac (also of note, unlike linux systems, we do not ***source*** script, we instead have to just run the script with ./)
```
./{path-to-repo}/scripts/build-fierro.sh --machine=mac [any other parameters we need to pass]
```
***Note*** you can still see more information on the script arguments by running the script with the argument ***--help***. However this is the standard script you would run on a Mac.

## Running the example
Assuming a smooth build, you will be placed at the top of the repository again. You can move yourself into the build directory
```
cd build-explicit
```
From there you can run this example with
```
./bin/fierro-parallel-explicit ../src/Parallel-Solvers/Parallel-Explicit/example_simple.yaml
```
***Note*** For this example the mesh file that is passed into *example_simple.yaml* needs to be present in your build directory. For the examples, you can copy the one in {path-to-repo}/single-node/src/Explicit-Lagrange/meshes/mesh_Sedov_32.geo 
