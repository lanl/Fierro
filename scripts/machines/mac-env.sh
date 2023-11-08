# Assign path's to the current software you want to use.
# These path's will need to be changed based on your computer
# The default for homebrew should place packages as shown below

export CC=/usr/local/opt/llvm/bin/clang
export CXX=/usr/local/opt/llvm/bin/clang++
export PATH="/usr/local/opt/llvm/bin:$PATH"
export LDFLAGS="-L/usr/local/opt/llvm/lib"
export CPPFLAGS="-I/usr/local/opt/llvm/include"
export OMPI_CC=${CC}
export OMPI_CXX=${CXX}
