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
# for m3906_g LBNL/AMP users: for large runs, comment in
##SBATCH -q early_science
#SBATCH -C gpu
#SBATCH -c 32
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=single:1
#SBATCH -o WarpX.o%j
#SBATCH -e WarpX.e%j

# ============
# -N =                 nodes
# -n =                 tasks (MPI ranks, usually = G)
# -G =                 GPUs (full Perlmutter node, 4)
# -c =                 CPU per task (128 total threads on CPU, 32 per GPU)
#
# --ntasks-per-node=   number of tasks (MPI ranks) per node (full node, 4)
# --gpus-per-task=     number of GPUs per task (MPI rank) (full node, 4)
# --gpus-per-node=     number of GPUs per node (full node, 4)
#
# --gpu-bind=single:1  sets only one GPU to be visible to each MPI rank
#                         (quiets AMReX init warnings)
#
# Recommend using --ntasks-per-node=4, --gpus-per-task=1 and --gpu-bind=single:1,
# as they are fixed values and allow for easy scaling with less adjustments.
#
# ============

# GPU-aware MPI
export MPICH_GPU_SUPPORT_ENABLED=1

EXE=./warpx
#EXE=../WarpX/build/bin/warpx.3d.MPI.CUDA.DP.OPMD.QED
#EXE=./main3d.gnu.TPROF.MPI.CUDA.ex
INPUTS=inputs_small

srun ${EXE} ${INPUTS} \
  > output.txt
