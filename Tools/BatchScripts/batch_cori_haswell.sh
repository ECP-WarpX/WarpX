#!/bin/bash -l

#SBATCH -N 1
#SBATCH -t 03:00:00
#SBATCH -q regular
#SBATCH -C haswell
#SBATCH -J <job name>
#SBATCH -A <allocation ID>
#SBATCH -e error.txt
#SBATCH -o output.txt
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=16

export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=16

export WARPX_NMPI_PER_NODE=4

EXE="<path/to/executable>"

srun --cpu_bind=cores -n $(( ${SLURM_JOB_NUM_NODES} * ${WARPX_NMPI_PER_NODE} )) ${EXE} <input file>
