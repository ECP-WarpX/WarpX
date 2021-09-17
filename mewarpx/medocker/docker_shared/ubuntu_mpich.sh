#!/bin/bash

set -e
set -x

echo "Missing fftw package for mpich"
exit 1

apt-get install -y \
    libhdf5-mpich-dev \
    libmpich-dev \
    mpich \
    ;
