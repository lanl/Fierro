### Load environment modules here
### Listed modules are for Darwin
module load cmake
module load openmpi

export scriptdir=`pwd`

cd ../../..
export basedir=`pwd`

export matardir=${basedir}/Elements/matar

export srcdir=${basedir}/Parallel-Solvers
export builddir=${basedir}/build-parallel-explicit-macosx
#export installdir=${basedir}/build-trilinos

export FIERRO_BASE_DIR=${basedir}
export FIERRO_EXPLICIT_SOURCE_DIR=${srcdir}/Parallel-Explicit
export FIERRO_BUILD_DIR=${builddir}

# Do this differently (in src tree) than other libs because
# of compile time
export TRILINOS_SOURCE_DIR=${basedir}/Trilinos
export TRILINOS_BUILD_DIR=${TRILINOS_SOURCE_DIR}/build
export TRILINOS_INSTALL_DIR=${TRILINOS_BUILD_DIR}

cd $scriptdir


