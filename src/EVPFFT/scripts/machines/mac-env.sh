# Assign path's to the current software you want to use.
# These path's will need to be changed based on your computer
# The default for homebrew should place packages as shown below

export CC=/opt/homebrew/opt/llvm/bin/clang
export CXX=/opt/homebrew/opt/llvm/bin/clang++
export PATH="/opt/homebrew/opt/llvm/bin:$PATH"
export LDFLAGS="-L/opt/homebrew/opt/llvm/lib"
export CPPFLAGS="-I/opt/homebrew/opt/llvm/include"
export OMPI_CC=${CC}
export OMPI_CXX=${CXX}
#export HDF5_ROOT=/path/to/where/HDF5/lib/and/include/folders/are/located
#export FFTW_ROOT=/path/to/where/FFTW/lib/and/include/folders/are/located
