module purge
### Load environment modules here
module load cmake
module load gcc/9.4.0
module list

export scriptdir=`pwd`
cd ../..

export basedir=`pwd`
export srcdir=${basedir}/src
export builddir=${basedir}/build-serial
export installdir=${srcdir}/install-serial

export MATAR_BASE_DIR=${basedir}
export MATAR_SOURCE_DIR=${srcdir}
export MATAR_BUILD_DIR=${builddir}

cd $scriptdir