#!/usr/bin/env bash
#
# Copyright 2023 The WarpX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Luca Fedeli

set -eu -o pipefail

# `man apt.conf`:
#   Number of retries to perform. If this is non-zero APT will retry
#   failed files the given number of times.
echo 'Acquire::Retries "3";' | sudo tee /etc/apt/apt.conf.d/80-retries

sudo apt-get -qqq update
sudo apt-get install -y \
    cmake               \
    clang-14            \
    clang-tidy-14       \
    libc++-14-dev       \
    libboost-math-dev   \
    libfftw3-dev        \
    libfftw3-mpi-dev    \
    libhdf5-openmpi-dev \
    libopenmpi-dev      \
    ninja-build
