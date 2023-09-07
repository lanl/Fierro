#!/bin/bash -e
SYSTEM=$1
PARALLEL=$2
CUSTOM_BUILD=$3

source setup-env.sh ${SYSTEM} ${PARALLEL} ${CUSTOM_BUILD}
source trilinos-install.sh ${SYSTEM} ${PARALLEL}
source cmake_build.sh ${SYSTEM} ${PARALLEL}
