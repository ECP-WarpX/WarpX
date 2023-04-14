#!/bin/bash

#SBATCH -A <project id>
#SBATCH -J warpx
#SBATCH -o %x-%j.out
#SBATCH -t 00:10:00
#SBATCH -p ecp
#SBATCH -N 1

export OMP_NUM_THREADS=1
srun -n 4 -c 2 --ntasks-per-node=4 ./warpx inputs > output.txt
