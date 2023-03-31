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
    ninja-build         \
    pkg-config          \
    python3             \
    python3-pip         \
    python3-setuptools  \
    wget

# cmake-easyinstall
#
sudo curl -L -o /usr/local/bin/cmake-easyinstall https://raw.githubusercontent.com/ax3l/cmake-easyinstall/main/cmake-easyinstall
sudo chmod a+x /usr/local/bin/cmake-easyinstall
export CEI_SUDO="sudo"
export CEI_TMP="/tmp/cei"

# BLAS++ & LAPACK++
cmake-easyinstall \
  --prefix=/usr/local                      \
  git+https://github.com/icl-utk-edu/blaspp.git \
  -Duse_openmp=OFF                         \
  -Dbuild_tests=OFF                        \
  -DCMAKE_CXX_COMPILER_LAUNCHER=$(which ccache) \
  -DCMAKE_VERBOSE_MAKEFILE=ON

cmake-easyinstall \
  --prefix=/usr/local                        \
  git+https://github.com/icl-utk-edu/lapackpp.git \
  -Duse_cmake_find_lapack=ON                 \
  -Dbuild_tests=OFF                          \
  -DCMAKE_CXX_COMPILER_LAUNCHER=$(which ccache) \
  -DCMAKE_VERBOSE_MAKEFILE=ON
