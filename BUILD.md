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
The minimum packages need to be installed (using separate install commands) are
```
git
cmake
gcc # If you want to use gcc. This will install gcc-13
llvm # If you want to use llvm and clang.
libomp # This is only necessary if you are using llvm/clang. gcc installs it's own OpenMP
openmpi 
```
