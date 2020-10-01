#!/bin/bash

# Copyright 2020 Axel Huebl
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL
#
# Refs.:
#   https://jsrunvisualizer.olcf.ornl.gov/?s4f0o11n6c7g1r11d1b1l0=
#   https://hpc.llnl.gov/training/tutorials/using-lcs-sierra-system#quick16

#BSUB -G <allocation ID>
#BSUB -W 00:10
#BSUB -nnodes 2
#BSUB -alloc_flags smt4
#BSUB -J WarpX
#BSUB -o WarpXo.%J
#BSUB -e WarpXe.%J

export OMP_NUM_THREADS=1
jsrun -r 4 -a 1 -g 1 -c 7 -l GPU-CPU -d packed -b rs -M "-gpu" <path/to/executable> <input file> > output.txt
