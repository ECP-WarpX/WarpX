#!/bin/bash

# Copyright 2019-2020 Maxence Thevenet, Axel Huebl
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL
#
# Refs.:
#   https://jsrunvisualizer.olcf.ornl.gov/?s4f0o11n6c7g1r11d1b1l0=
#   https://docs.olcf.ornl.gov/systems/summit_user_guide.html#cuda-aware-mpi

#BSUB -P <allocation ID>
#BSUB -W 00:10
#BSUB -nnodes 2
#BSUB -alloc_flags smt4
#BSUB -J WarpX
#BSUB -o WarpXo.%J
#BSUB -e WarpXe.%J

module load gcc
module load cuda

export OMP_NUM_THREADS=1
jsrun -r 6 -a 1 -g 1 -c 7 -l GPU-CPU -d packed -b rs --smpiargs="-gpu" <path/to/executable> <input file> > output.txt
