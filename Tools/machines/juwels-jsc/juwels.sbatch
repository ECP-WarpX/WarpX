#!/bin/bash -l

#SBATCH -A $proj
#SBATCH --partition=booster
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --time=00:05:00
#SBATCH --job-name=warpx
#SBATCH --output=warpx-%j-%N.txt
#SBATCH --error=warpx-%j-%N.err

export OMP_NUM_THREADS=1
export OMPI_MCA_io=romio321  # for HDF5 support in openPMD

# you can comment this out if you sourced the warpx.profile
# files before running sbatch:
module load GCC
module load OpenMPI
module load CUDA/11.3
module load HDF5
module load Python

srun -n 8 --cpu_bind=sockets $HOME/src/warpx/build/bin/warpx.3d.MPI.CUDA.DP.OPMD.QED inputs
