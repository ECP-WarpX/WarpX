#!/bin/bash -l

# Copyright 2021 Axel Huebl
# This file is part of WarpX.
# License: BSD-3-Clause-LBNL
#
# Ref:
# - https://docs-dev.nersc.gov/cgpu/hardware/
# - https://docs-dev.nersc.gov/cgpu/access/

# Just increase this number of you need more nodes.
#SBATCH -N 2
#SBATCH -t 03:00:00
#SBATCH -J <job name>
#SBATCH -A <allocation ID>
#SBATCH -q regular
#SBATCH -C gpu
# 8 V100 GPUs (16 GB) per node
#SBATCH --gres=gpu:8
#SBATCH --exclusive
# one MPI rank per half-socket (see below)
#SBATCH --tasks-per-node=4
# request all logical (virtual) cores per half-socket
#SBATCH --cpus-per-task=20
#SBATCH -e WarpX.e%j
#SBATCH -o WarpX.o%j


# each Cori GPU node has 2 sockets of Intel Xeon Gold 6148 ('Skylake') @ 2.40 GHz
# each Xeon CPU is divided into 2 bus rings that each have direct L3 access (TODO: double-check)
export WARPX_NMPI_PER_NODE=4

# each MPI rank per half-socket has 10 physical cores
#   or 20 logical (virtual) cores
# over-subscribing each physical core with 2x
#   hyperthreading leads to often to slight speedup on Intel
# the settings below make sure threads are close to the
#   controlling MPI rank (process) per half socket and
#   distribute equally over close-by physical cores and,
#   for N>20, also equally over close-by logical cores
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=20

# for async_io support: (optional)
export MPICH_MAX_THREAD_SAFETY=multiple

EXE="<path/to/executable>"

srun --cpu_bind=cores -n $(( ${SLURM_JOB_NUM_NODES} * ${WARPX_NMPI_PER_NODE} )) \
  ${EXE} <input file> \
  > output.txt
