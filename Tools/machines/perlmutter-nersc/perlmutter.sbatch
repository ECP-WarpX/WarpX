#!/bin/bash -l

# Copyright 2021 Axel Huebl, Kevin Gott
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

#SBATCH -t 01:00:00
#SBATCH -N 4
#SBATCH -J WarpX
#    note: <proj> must end on _g
#SBATCH -A <proj>
#SBATCH -q regular
#SBATCH -C gpu
#SBATCH -c 32
#SBATCH --ntasks-per-gpu=1
#SBATCH --gpus-per-node=4
#SBATCH -o WarpX.o%j
#SBATCH -e WarpX.e%j

# GPU-aware MPI
export MPICH_GPU_SUPPORT_ENABLED=1

# expose one GPU per MPI rank
export CUDA_VISIBLE_DEVICES=$SLURM_LOCALID

EXE=./warpx
#EXE=../WarpX/build/bin/warpx.3d.MPI.CUDA.DP.OPMD.QED
#EXE=./main3d.gnu.TPROF.MPI.CUDA.ex
INPUTS=inputs_small

srun ${EXE} ${INPUTS} \
  > output.txt
