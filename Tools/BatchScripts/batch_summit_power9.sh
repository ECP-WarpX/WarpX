#!/bin/bash

# Copyright 2019-2020 Maxence Thevenet, Axel Huebl, Michael Rowan
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL
#
# Refs.:
#   https://jsrunvisualizer.olcf.ornl.gov/?s1f0o121n2c21g0r11d1b1l0=

#BSUB -P <allocation ID>
#BSUB -W 00:10
#BSUB -nnodes 1
#BSUB -alloc_flags "smt1"
#BSUB -J WarpX
#BSUB -o WarpXo.%J
#BSUB -e WarpXe.%J

# make output group-readable by default
umask 0027

# fix problems with collectives since RHEL8 update: OLCFHELP-3545
# disable all the IBM optimized barriers and drop back to HCOLL or OMPI's barrier implementations
export OMPI_MCA_coll_ibm_skip_barrier=true

# libfabric 1.6+: limit the visible devices
# Needed for ADIOS2 SST staging/streaming workflows since RHEL8 update
#   https://github.com/ornladios/ADIOS2/issues/2887
#export FABRIC_IFACE=mlx5_0   # ADIOS SST: select interface (1 NIC on Summit)
#export FI_OFI_RXM_USE_SRX=1  # libfabric: use shared receive context from MSG provider

# OpenMP: 21 threads per MPI rank
export OMP_NUM_THREADS=21

# run WarpX
jsrun -n 2 -a 1 -c 21 -r 2 -l CPU-CPU -d packed -b rs <path/to/executable> <input file> > output.txt
