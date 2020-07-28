#!/bin/bash -l

#SBATCH -N 1
#SBATCH -t 03:00:00
#SBATCH -q regular
#SBATCH -C haswell
#SBATCH -J n1m4o8
#SBATCH -A m2852
#SBATCH -e error_spread.txt
#SBATCH -o output_spread.txt
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=16

export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=16

export WARPX_NMPI_PER_NODE=4

EXE="$HOME/WarpX/Bin/main3d.gnu.haswell.TPROF.MPI.OMP.ex"

srun --cpu_bind=cores -n $(( ${SLURM_JOB_NUM_NODES} * ${WARPX_NMPI_PER_NODE} )) ${EXE} inputs
