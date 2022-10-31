#!/bin/bash -l

# Copyright 2021-2022 Axel Huebl, Kevin Gott
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
#SBATCH --exclusive
#SBATCH --ntasks-per-gpu=1
#SBATCH --gpus-per-node=4
#SBATCH -o WarpX.o%j
#SBATCH -e WarpX.e%j

# GPU-aware MPI
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=GPU

# threads for OpenMP and threaded compressors per MPI rank
export SRUN_CPUS_PER_TASK=32

EXE=./warpx
#EXE=../WarpX/build/bin/warpx.3d.MPI.CUDA.DP.OPMD.QED
#EXE=./main3d.gnu.TPROF.MPI.CUDA.ex
INPUTS=inputs_small

# CUDA visible devices are ordered inverse to local task IDs
srun /bin/bash -l -c "  \
    export CUDA_VISIBLE_DEVICES=$((3-SLURM_LOCALID));
    ${EXE} ${INPUTS} \
      amrex.the_arena_is_managed=0 \
      amrex.use_gpu_aware_mpi=1"   \
  > output.txt
