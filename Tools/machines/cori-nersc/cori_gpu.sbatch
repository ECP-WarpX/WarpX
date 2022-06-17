#!/bin/bash -l

# Copyright 2021 Axel Huebl
# This file is part of WarpX.
# License: BSD-3-Clause-LBNL
#
# Ref:
# - https://docs-dev.nersc.gov/cgpu/hardware/
# - https://docs-dev.nersc.gov/cgpu/access/
# - https://docs-dev.nersc.gov/cgpu/usage/#controlling-task-and-gpu-binding

# Just increase this number of you need more nodes.
#SBATCH -N 2
#SBATCH -t 03:00:00
#SBATCH -J <job name>
#SBATCH -A m1759
#SBATCH -q regular
#SBATCH -C gpu
# 8 V100 GPUs (16 GB) per node
#SBATCH --gres=gpu:8
#SBATCH --exclusive
# one MPI rank per GPU (a quarter-socket)
#SBATCH --tasks-per-node=8
# request all logical (virtual) cores per quarter-socket
#SBATCH --cpus-per-task=10
#SBATCH -e WarpX.e%j
#SBATCH -o WarpX.o%j


# each Cori GPU node has 2 sockets of Intel Xeon Gold 6148 ('Skylake') @ 2.40 GHz
export WARPX_NMPI_PER_NODE=8

# each MPI rank per half-socket has 10 physical cores
#   or 20 logical (virtual) cores
# we split half-sockets again by 2 to have one MPI rank per GPU
# over-subscribing each physical core with 2x
#   hyperthreading leads to often to slight speedup on Intel
# the settings below make sure threads are close to the
#   controlling MPI rank (process) per half socket and
#   distribute equally over close-by physical cores and,
#   for N>20, also equally over close-by logical cores
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=10

# for async_io support: (optional)
export MPICH_MAX_THREAD_SAFETY=multiple

EXE="<path/to/executable>"

srun --cpu_bind=cores --gpus-per-task=1 --gpu-bind=map_gpu:0,1,2,3,4,5,6,7 \
  -n $(( ${SLURM_JOB_NUM_NODES} * ${WARPX_NMPI_PER_NODE} )) \
  ${EXE} <input file> \
  > output.txt
