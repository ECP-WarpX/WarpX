#!/bin/bash

#SBATCH -A <project id>
#SBATCH -J warpx
#SBATCH -o %x-%j.out
#SBATCH -t 00:10:00
# Early access to the GPU partition
#SBATCH -p eap
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=8
#SBATCH --gpu-bind=closest

export MPICH_GPU_SUPPORT_ENABLED=1

# note (12-12-22)
# this environment setting is currently needed on LUMI to work-around a
# known issue with Libfabric
#export FI_MR_CACHE_MAX_COUNT=0  # libfabric disable caching
# or, less invasive:
export FI_MR_CACHE_MONITOR=memhooks  # alternative cache monitor

# note (9-2-22, OLCFDEV-1079)
# this environment setting is needed to avoid that rocFFT writes a cache in
# the home directory, which does not scale.
export ROCFFT_RTC_CACHE_PATH=/dev/null

export OMP_NUM_THREADS=1
srun ../warpx inputs > outputs
