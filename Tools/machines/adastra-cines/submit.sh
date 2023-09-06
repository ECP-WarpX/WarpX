#!/bin/bash
#SBATCH --job-name=warpx
#SBATCH --account=<account_to_charge>
#SBATCH --constraint=MI250
#SBATCH --ntasks-per-node=8 --cpus-per-task=8 --gpus-per-node=8
#SBATCH --threads-per-core=1 # --hint=nomultithread
#SBATCH --exclusive
#SBATCH --output=%x-%j.out
#SBATCH --time=00:10:00
#SBATCH --nodes=2

module purge

# Architecture
module load craype-accel-amd-gfx90a craype-x86-trento
# A compiler to target the architecture
module load PrgEnv-cray
# Some architecture related libraries and tools
module load amd-mixed

export MPICH_GPU_SUPPORT_ENABLED=1

# note
# this environment setting is currently needed to work-around a
# known issue with Libfabric
#export FI_MR_CACHE_MAX_COUNT=0  # libfabric disable caching
# or, less invasive:
export FI_MR_CACHE_MONITOR=memhooks  # alternative cache monitor

# note
# this environment setting is needed to avoid that rocFFT writes a cache in
# the home directory, which does not scale.
export ROCFFT_RTC_CACHE_PATH=/dev/null

export OMP_NUM_THREADS=1
export WARPX_NMPI_PER_NODE=8
export TOTAL_NMPI=$(( ${SLURM_JOB_NUM_NODES} * ${WARPX_NMPI_PER_NODE} ))
srun -N${SLURM_JOB_NUM_NODES} -n${TOTAL_NMPI} --ntasks-per-node=${WARPX_NMPI_PER_NODE} \
    ./warpx inputs > output.txt
