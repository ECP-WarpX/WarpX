#!/bin/bash -l

# Copyright 2022 The WarpX Community
#
# This file is part of WarpX.
#
# Author: Axel Huebl
# License: BSD-3-Clause-LBNL

#SBATCH -t 01:00:00
#SBATCH --ntasks=1
#SBATCH -J WarpX
#SBATCH --gres=gpu:a100:1
#SBATCH --constraint="gpu"
#SBATCH --cpus-per-task=18
#SBATCH --mem=125000
#SBATCH --mail-type=none
#SBATCH --mail-user=userid@example.mpg.de
#SBATCH -o WarpX.o%j
#SBATCH -e WarpX.e%j

# exit immediately at any error
set -e

# put TMPDIR onto a RAMDISK for fast builds and cleanup, moreover the
# default /tmp would be too small to hold all temporary files
export TMPDIR=$JOB_SHMTMPDIR

# GPU-aware MPI
#export MPICH_GPU_SUPPORT_ENABLED=1

EXE=./warpx
#EXE=../WarpX/build/bin/warpx.3d.MPI.CUDA.DP.OPMD.QED
INPUTS=inputs_small

srun ${EXE} ${INPUTS} \
  > output.txt
