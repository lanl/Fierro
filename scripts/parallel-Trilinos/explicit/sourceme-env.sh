### Load environment modules here
### Listed modules are for Darwin
module load cmake
module load openmpi

export scriptdir=`pwd`

cd ../../..
export basedir=`pwd`
export srcdir=${basedir}/src
export matardir=${srcdir}/Elements/matar
export trilinosdir=${srcdir}

export builddir=${basedir}/build-parallel-explicit
# Only used for non-Trilinos installs
export installdir=${basedir}/install

export FIERRO_BASE_DIR=${basedir}
export FIERRO_EXPLICIT_SOURCE_DIR=${srcdir}/Parallel-Solvers/Parallel-Explicit
export FIERRO_BUILD_DIR=${builddir}

# Do this differently (in src tree) than other libs because
# of compile time
export TRILINOS_SOURCE_DIR=${trilinosdir}/Trilinos
export TRILINOS_BUILD_DIR=${TRILINOS_SOURCE_DIR}/build
export TRILINOS_INSTALL_DIR=${TRILINOS_BUILD_DIR}

cd $scriptdir


