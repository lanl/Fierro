#!/bin/bash
# Rudimentary GPU affinity setter for Summit supercomputer
# >$ jsrun -rs_per_host 1 -gpu_per_rs 6 <task/cpu option> ./gpu_setter.sh <your app>

# This script assumes your code does not attempt to set its own
# GPU affinity (e.g., with cudaSetDevice). Using this affinity script
# with a code that does its own internal GPU selection probably won't work

# Compute device number from OpenMPI local rank environment variable
# Keeping in mind Summit has 6 GPUs per node

mydevice=$((${OMPI_COMM_WORLD_LOCAL_RANK} % 6))

# CUDA_VISIBLE_DEVICES controls both what GPUs are visible to your process
# and the other they appear in. By putting "mydevice" first the in list, we
# make sure it shows up as device "0" to the process so it's automatically selected.
# The order of the other devices doesn't matter, only that all devices (0-5) are present

CUDA_VISIBLE_DEVICES="${mydevice},0,1,2,3,4,5"

# Process with sed to remove the duplicate and reform the list, keeping the order we set
CUDA_VISIBLE_DEVICES=$(sed -r ':a; s/\b([[:alnum:]]+)\b(.*)\b\1\b/\1\2/g; ta; s/(,,)+/,/g; s/, *$//' <<< $CUDA_VISIBLE_DEVICES)

export CUDA_VISIBLE_DEVICES

# Launch the application we were given
exec "$@"
