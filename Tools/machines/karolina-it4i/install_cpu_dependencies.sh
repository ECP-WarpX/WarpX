#!/bin/bash
#
# Copyright 2023 The WarpX Community
#
# This file is part of WarpX.
#
# Author: Axel Huebl
# License: BSD-3-Clause-LBNL

# Exit on first error encountered #############################################
#
set -eu -o pipefail


# Check: ######################################################################
#
#   Was karolina_cpu_warpx.profile sourced and configured correctly?
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your karolina_cpu_warpx.profile file! Please edit its line 2 to continue!"; exit 1; fi


# Remove old dependencies #####################################################
#
SW_DIR="${HOME}/sw/karolina/cpu"
rm -rf ${SW_DIR}
mkdir -p ${SW_DIR}

# remove common user mistakes in python, located in .local instead of a venv
python3 -m pip uninstall -qq -y pywarpx
python3 -m pip uninstall -qq -y warpx
python3 -m pip uninstall -qqq -y mpi4py 2>/dev/null || true


# General extra dependencies ##################################################
#

# c-blosc (I/O compression)
if [ -d $HOME/src/c-blosc ]
then
  cd $HOME/src/c-blosc
  git fetch
  git checkout v1.21.1
  cd -
else
  git clone -b v1.21.1 https://github.com/Blosc/c-blosc.git $HOME/src/c-blosc
fi
rm -rf $HOME/src/c-blosc-cpu-build
cmake -S $HOME/src/c-blosc -B $HOME/src/c-blosc-cpu-build -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/c-blosc-1.21.1
cmake --build $HOME/src/c-blosc-cpu-build --target install --parallel 16
rm -rf $HOME/src/c-blosc-cpu-build

# HDF5
if [ -d $HOME/src/hdf5 ]
then
  cd $HOME/src/hdf5
  git fetch
  git checkout hdf5-1_14_1-2
  cd -
else
  git clone -b hdf5-1_14_1-2 https://github.com/HDFGroup/hdf5.git src/hdf5
fi
rm -rf $HOME/src/hdf5-build
cmake -S $HOME/src/hdf5 -B $HOME/src/hdf5-build -DBUILD_TESTING=OFF -DHDF5_ENABLE_PARALLEL=ON -DCMAKE_INSTALL_PREFIX=${SW_DIR}/hdf5-1.14.1.2
cmake --build $HOME/src/hdf5-build --target install --parallel 16
rm -rf $HOME/src/hdf5-build

# ADIOS2
if [ -d $HOME/src/adios2 ]
then
  cd $HOME/src/adios2
  git fetch
  git checkout v2.8.3
  cd -
else
  git clone -b v2.8.3 https://github.com/ornladios/ADIOS2.git $HOME/src/adios2
fi
rm -rf $HOME/src/adios2-cpu-build
cmake -S $HOME/src/adios2 -B $HOME/src/adios2-cpu-build -DADIOS2_USE_Blosc=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/adios2-2.8.3
cmake --build $HOME/src/adios2-cpu-build --target install -j 16
rm -rf $HOME/src/adios2-cpu-build

# BLAS++ (for PSATD+RZ)
if [ -d $HOME/src/blaspp ]
then
  cd $HOME/src/blaspp
  git fetch
  git checkout master
  git pull
  cd -
else
  git clone https://github.com/icl-utk-edu/blaspp.git $HOME/src/blaspp
fi
rm -rf $HOME/src/blaspp-cpu-build
cmake -S $HOME/src/blaspp -B $HOME/src/blaspp-cpu-build -Duse_openmp=ON -Dcpu_backend=OFF -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=${SW_DIR}/blaspp-master
cmake --build $HOME/src/blaspp-cpu-build --target install --parallel 16
rm -rf $HOME/src/blaspp-cpu-build

# LAPACK++ (for PSATD+RZ)
if [ -d $HOME/src/lapackpp ]
then
  cd $HOME/src/lapackpp
  git fetch
  git checkout master
  git pull
  cd -
else
  git clone https://github.com/icl-utk-edu/lapackpp.git $HOME/src/lapackpp
fi
rm -rf $HOME/src/lapackpp-cpu-build
CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S $HOME/src/lapackpp -B $HOME/src/lapackpp-cpu-build -DCMAKE_CXX_STANDARD=17 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=${SW_DIR}/lapackpp-master
cmake --build $HOME/src/lapackpp-cpu-build --target install --parallel 16
rm -rf $HOME/src/lapackpp-cpu-build


# Python ######################################################################
#
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade virtualenv
python3 -m pip cache purge
rm -rf ${SW_DIR}/venvs/warpx-cpu
python3 -m venv ${SW_DIR}/venvs/warpx-cpu
source ${SW_DIR}/venvs/warpx-cpu/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade wheel
python3 -m pip install --upgrade cython
python3 -m pip install --upgrade numpy
python3 -m pip install --upgrade pandas
python3 -m pip install --upgrade scipy
python3 -m pip install --upgrade mpi4py --no-cache-dir --no-build-isolation --no-binary mpi4py
python3 -m pip install --upgrade openpmd-api
python3 -m pip install --upgrade matplotlib
python3 -m pip install --upgrade yt
# install or update WarpX dependencies such as picmistandard
python3 -m pip install --upgrade -r $HOME/src/warpx/requirements.txt
# optional: for libEnsemble
python3 -m pip install -r $HOME/src/warpx/Tools/LibEnsemble/requirements.txt
# optional: for optimas (based on libEnsemble & ax->botorch->gpytorch->pytorch)
python3 -m pip install --upgrade torch --index-url https://download.pytorch.org/whl/cpu
python3 -m pip install -r $HOME/src/warpx/Tools/optimas/requirements.txt
