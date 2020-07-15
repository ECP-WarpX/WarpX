#!/bin/bash -l

#SBATCH -A $proj
#SBATCH --partition=gpus
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --time=00:05:00
#SBATCH --job-name=warpx
#SBATCH --output=warpx-%j-%N.txt
#SBATCH --error=warpx-%j-%N.err

export OMP_NUM_THREADS=1

module load GCC
module load OpenMPI
module load CUDA

srun -n 8 --cpu_bind=sockets $HOME/src/warpx/Bin/main3d.gnu.TPROF.MPI.CUDA.ex inputs
