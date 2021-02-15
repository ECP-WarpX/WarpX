#!/usr/bin/env bash
#
# Copyright 2021 The WarpX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

sudo apt-get -qqq update
sudo apt-get install -y \
    build-essential     \
    ca-certificates     \
    ccache              \
    clang               \
    cmake               \
    gnupg               \
    libblas-dev         \
    libboost-math-dev   \
    libfftw3-dev        \
    libfftw3-mpi-dev    \
    libhdf5-openmpi-dev \
    liblapack-dev       \
    libopenmpi-dev      \
    make                \
    pkg-config          \
    python3             \
    python3-pip         \
    python3-setuptools  \
    wget

# cmake-easyinstall
#
sudo curl -L -o /usr/local/bin/cmake-easyinstall https://git.io/JvLxY
sudo chmod a+x /usr/local/bin/cmake-easyinstall
export CEI_SUDO="sudo"

# BLAS++ & LAPACK++
cmake-easyinstall \
  --prefix=/usr/local                      \
  git+https://bitbucket.org/icl/blaspp.git \
  -Duse_openmp=OFF                         \
  -Dbuild_tests=OFF                        \
  -DCMAKE_VERBOSE_MAKEFILE=ON

cmake-easyinstall \
  --prefix=/usr/local                        \
  git+https://bitbucket.org/icl/lapackpp.git \
  -Duse_cmake_find_lapack=ON                 \
  -Dbuild_tests=OFF                          \
  -DCMAKE_VERBOSE_MAKEFILE=ON
