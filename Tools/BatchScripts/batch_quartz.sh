#!/bin/bash -l

# Just increase this number of you need more nodes.
#SBATCH -N 2
#SBATCH -t 24:00:00
#SBATCH -A <allocation ID>

#SBATCH -J WarpX
#SBATCH -q pbatch
#SBATCH --qos=normal
#SBATCH --license=lustre1,lustre2
#SBATCH --export=ALL
#SBATCH -e error.txt
#SBATCH -o output.txt
# one MPI rank per half-socket (see below)
#SBATCH --tasks-per-node=2
# request all logical (virtual) cores per half-socket
#SBATCH --cpus-per-task=18


# each Quartz node has 1 socket of Intel Xeon E5-2695 v4
# each Xeon CPU is divided into 2 bus rings that each have direct L3 access
export WARPX_NMPI_PER_NODE=2

# each MPI rank per half-socket has 9 physical cores
#   or 18 logical (virtual) cores
# over-subscribing each physical core with 2x
#   hyperthreading led to a slight (3.5%) speedup on Cori's Intel Xeon E5-2698 v3,
#   so we do the same here
# the settings below make sure threads are close to the
#   controlling MPI rank (process) per half socket and
#   distribute equally over close-by physical cores and,
#   for N>9, also equally over close-by logical cores
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=18

EXE="<path/to/executable>"  # e.g. ./warpx

srun --cpu_bind=cores -n $(( ${SLURM_JOB_NUM_NODES} * ${WARPX_NMPI_PER_NODE} )) ${EXE} <input file>
