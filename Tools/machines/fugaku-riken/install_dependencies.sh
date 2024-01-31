#!/bin/bash
#
# Copyright 2023 The WarpX Community
#
# This file is part of WarpX.
#
# Author: Axel Huebl, Luca Fedeli
# License: BSD-3-Clause-LBNL

# Exit on first error encountered #############################################
#
set -eu -o pipefail

# Remove old dependencies #####################################################
#
SW_DIR="${HOME}/sw/fugaku/a64fx"
rm -rf ${SW_DIR}
mkdir -p ${SW_DIR}

# General extra dependencies ##################################################
#

# BLAS++ (for PSATD+RZ)
if [ -d $HOME/src/blaspp ]
then
  cd $HOME/src/blaspp
  git fetch --prune
  git checkout master
  git pull
  cd -
else
  git clone https://github.com/icl-utk-edu/blaspp.git $HOME/src/blaspp
fi
rm -rf $HOME/src/blaspp-frontier-gpu-build
CXX=$(which CC) cmake -S $HOME/src/blaspp -B $HOME/src/blaspp-frontier-gpu-build -Duse_openmp=OFF -Dgpu_backend=hip -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=${SW_DIR}/blaspp-master
cmake --build $HOME/src/blaspp-frontier-gpu-build --target install --parallel 16
rm -rf $HOME/src/blaspp-frontier-gpu-build

# LAPACK++ (for PSATD+RZ)
if [ -d $HOME/src/lapackpp ]
then
  cd $HOME/src/lapackpp
  git fetch --prune
  git checkout master
  git pull
  cd -
else
  git clone https://github.com/icl-utk-edu/lapackpp.git $HOME/src/lapackpp
fi
rm -rf $HOME/src/lapackpp-frontier-gpu-build
CXX=$(which CC) CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S $HOME/src/lapackpp -B $HOME/src/lapackpp-frontier-gpu-build -DCMAKE_CXX_STANDARD=17 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=${SW_DIR}/lapackpp-master
cmake --build $HOME/src/lapackpp-frontier-gpu-build --target install --parallel 16
rm -rf $HOME/src/lapackpp-frontier-gpu-build

# c-blosc (I/O compression)
if [ -d $HOME/src/c-blosc ]
then
  cd $HOME/src/c-blosc
  git fetch --prune
  git checkout master
  git pull
  cd -
else
  git clone -b v1.21.1 https://github.com/Blosc/c-blosc.git $HOME/src/c-blosc
fi
  rm -rf $HOME/src/c-blosc-fugaku-build
  cmake -S $HOME/src/c-blosc -B $HOME/src/c-blosc-fugaku-build -DBUILD_SHARED_LIBS=OFF -DBUILD_SHARED=OFF -DBUILD_STATIC=ON -DBUILD_TESTS=OFF -DBUILD_FUZZERS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=${HOME}/sw/a64fx/c-blosc-1.21.1-install
  cmake --build $HOME/src/c-blosc-fugaku-build --target install --parallel 48

# ADIOS2 (I/O)
if [ -d $HOME/src/c-blosc ]
then
  cd $HOME/src/adios2
  git fetch --prune
  git checkout master
  git pull
  cd -
else
  git clone -b v2.8.3 https://github.com/ornladios/ADIOS2.git $HOME/src/adios2
fi
rm -rf $HOME/src/adios2-fugaku-build
cmake -S src/adios2 -B $HOME/src/adios2-fugaku-build -DBUILD_SHARED_LIBS=OFF -DADIOS2_USE_Blosc=ON -DBUILD_TESTING=OFF -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=${HOME}/sw/a64fx/adios2-2.8.3-install
cmake --build $HOME/src/adios2-fugaku-build --target install -j 48
