#!/bin/bash

set -e
set -x

apt-get install -y \
    libfftw3-mpi-dev    \
    libhdf5-openmpi-dev \
    libopenmpi-dev \
    openmpi-bin \
    ;
