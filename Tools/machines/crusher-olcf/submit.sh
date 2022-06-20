#!/usr/bin/env bash

#SBATCH -A <project id>
#SBATCH -J warpx
#SBATCH -o %x-%j.out
#SBATCH -t 00:10:00
#SBATCH -p batch
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=closest
#SBATCH -N 1

# From the documentation:
# Each Crusher compute node consists of [1x] 64-core AMD EPYC 7A53
# "Optimized 3rd Gen EPYC" CPU (with 2 hardware threads per physical core) with
# access to 512 GB of DDR4 memory.
# Each node also contains [4x] AMD MI250X, each with 2 Graphics Compute Dies
# (GCDs) for a total of 8 GCDs per node. The programmer can think of the 8 GCDs
# as 8 separate GPUs, each having 64 GB of high-bandwidth memory (HBM2E).

# note (5-16-22)
# this environment setting is currently needed on Crusher to work-around a
# known issue with Libfabric
export FI_MR_CACHE_MAX_COUNT=0

export OMP_NUM_THREADS=8
srun ./warpx inputs > output.txt
