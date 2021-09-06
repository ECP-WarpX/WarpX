#!/bin/bash -l

# Just increase this number of you need more nodes.
#SBATCH -N 1
#SBATCH -t 03:00:00
#SBATCH -q regular
#SBATCH -C haswell
#SBATCH -J <job name>
#SBATCH -A <allocation ID>
#SBATCH -e error.txt
#SBATCH -o output.txt
# one MPI rank per half-socket (see below)
#SBATCH --tasks-per-node=4
# request all logical (virtual) cores per half-socket
#SBATCH --cpus-per-task=16


# each Cori Haswell node has 2 sockets of Intel Xeon E5-2698 v3
# each Xeon CPU is divided into 2 bus rings that each have direct L3 access
export WARPX_NMPI_PER_NODE=4

# each MPI rank per half-socket has 8 physical cores
#   or 16 logical (virtual) cores
# over-subscribing each physical core with 2x
#   hyperthreading leads to a slight (3.5%) speedup
# the settings below make sure threads are close to the
#   controlling MPI rank (process) per half socket and
#   distribute equally over close-by physical cores and,
#   for N>8, also equally over close-by logical cores
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=16

# for async_io support: (optional)
export MPICH_MAX_THREAD_SAFETY=multiple

EXE="<path/to/executable>"

srun --cpu_bind=cores -n $(( ${SLURM_JOB_NUM_NODES} * ${WARPX_NMPI_PER_NODE} )) \
  ${EXE} <input file> \
  > output.txt
