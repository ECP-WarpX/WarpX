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
#   Was frontier_warpx.profile sourced and configured correctly?
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your frontier_warpx.profile file! Please edit its line 2 to continue!"; exit 1; fi


# Check $proj variable is correct and has a corresponding CFS directory #######
#
if [ ! -d "${PROJWORK}/${proj}/" ]
then
    echo "WARNING: The directory $PROJWORK/$proj/ does not exist!"
    echo "Is the \$proj environment variable of value \"$proj\" correctly set? "
    echo "Please edit line 2 of your frontier_warpx.profile file to continue!"
    exit
fi


# Remove old dependencies #####################################################
#
SW_DIR="${HOME}/sw/frontier/gpu"
rm -rf ${SW_DIR}
mkdir -p ${SW_DIR}

# remove common user mistakes in python, located in .local instead of a venv
python3 -m pip uninstall -qq -y pywarpx
python3 -m pip uninstall -qq -y warpx
python3 -m pip uninstall -qqq -y mpi4py 2>/dev/null || true


# General extra dependencies ##################################################
#

# BLAS++ (for PSATD+RZ)
if [ -d $HOME/src/blaspp ]
then
  cd $HOME/src/blaspp
  git fetch --prune
  git checkout v2024.05.31
  cd -
else
  git clone -b v2024.05.31 https://github.com/icl-utk-edu/blaspp.git $HOME/src/blaspp
fi
rm -rf $HOME/src/blaspp-frontier-gpu-build
CXX=$(which CC) cmake -S $HOME/src/blaspp -B $HOME/src/blaspp-frontier-gpu-build -Duse_openmp=OFF -Dgpu_backend=hip -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=${SW_DIR}/blaspp-2024.05.31
cmake --build $HOME/src/blaspp-frontier-gpu-build --target install --parallel 16
rm -rf $HOME/src/blaspp-frontier-gpu-build

# LAPACK++ (for PSATD+RZ)
if [ -d $HOME/src/lapackpp ]
then
  cd $HOME/src/lapackpp
  git fetch --prune
  git checkout v2024.05.31
  cd -
else
  git clone -b v2024.05.31 https://github.com/icl-utk-edu/lapackpp.git $HOME/src/lapackpp
fi
rm -rf $HOME/src/lapackpp-frontier-gpu-build
CXX=$(which CC) CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S $HOME/src/lapackpp -B $HOME/src/lapackpp-frontier-gpu-build -DCMAKE_CXX_STANDARD=17 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=${SW_DIR}/lapackpp-2024.05.31
cmake --build $HOME/src/lapackpp-frontier-gpu-build --target install --parallel 16
rm -rf $HOME/src/lapackpp-frontier-gpu-build

# c-blosc (I/O compression, for OpenPMD)
if [ -d $HOME/src/c-blosc ]
then
  # git repository is already there
  :
else
  git clone -b v1.21.1 https://github.com/Blosc/c-blosc.git $HOME/src/c-blosc
fi
rm -rf $HOME/src/c-blosc-frontier-build
cmake -S $HOME/src/c-blosc -B $HOME/src/c-blosc-frontier-build -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/c-blosc-1.21.1
cmake --build $HOME/src/c-blosc-frontier-build --target install --parallel 16
rm -rf $HOME/src/c-blosc-frontier-build

# HDF5 (for openPMD)
if [ -d $HOME/src/hdf5 ]
then
  cd $HOME/src/hdf5
  git fetch --prune
  git checkout hdf5-1_14_1-2
  cd -
else
  git clone -b hdf5-1_14_1-2 https://github.com/HDFGroup/hdf5.git $HOME/src/hdf5
fi
rm -rf $HOME/src/hdf5-frontier-build
cmake -S $HOME/src/hdf5          \
      -B $HOME/src/hdf5-frontier-build  \
      -DBUILD_TESTING=OFF         \
      -DHDF5_ENABLE_PARALLEL=ON   \
      -DCMAKE_INSTALL_PREFIX=${SW_DIR}/hdf5-1.14.1.2
cmake --build $HOME/src/hdf5-frontier-build --target install --parallel 10
rm -rf $HOME/src/hdf5-frontier-build

# ADIOS2 v. 2.9.2 (for OpenPMD)
if [ -d $HOME/src/adios2 ]
then
  # git repository is already there
  :
else
  git clone -b v2.9.2 https://github.com/ornladios/ADIOS2.git $HOME/src/adios2
fi
rm -rf $HOME/src/adios2-frontier-build
cmake -S $HOME/src/adios2 -B $HOME/src/adios2-frontier-build -DADIOS2_USE_Blosc=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/adios2-2.9.2
cmake --build $HOME/src/adios2-frontier-build --target install -j 16
rm -rf $HOME/src/adios2-frontier-build


# Python ######################################################################
#
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade virtualenv
python3 -m pip cache purge
rm -rf ${SW_DIR}/venvs/warpx-frontier
python3 -m venv ${SW_DIR}/venvs/warpx-frontier
source ${SW_DIR}/venvs/warpx-frontier/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build
python3 -m pip install --upgrade packaging
python3 -m pip install --upgrade wheel
python3 -m pip install --upgrade setuptools
# cupy and h5py need an older Cython
# https://github.com/cupy/cupy/issues/4610
# https://github.com/h5py/h5py/issues/2268
python3 -m pip install --upgrade "cython<3.0"
python3 -m pip install --upgrade numpy
python3 -m pip install --upgrade pandas
python3 -m pip install --upgrade scipy
MPICC="cc -shared" python3 -m pip install --upgrade mpi4py --no-cache-dir --no-build-isolation --no-binary mpi4py
python3 -m pip install --upgrade openpmd-api
python3 -m pip install --upgrade matplotlib
python3 -m pip install --upgrade yt
# install or update WarpX dependencies such as picmistandard
python3 -m pip install --upgrade -r $HOME/src/warpx/requirements.txt
# cupy for ROCm
#   https://docs.cupy.dev/en/stable/install.html#building-cupy-for-rocm-from-source
#   https://github.com/cupy/cupy/issues/7830
CC=cc CXX=CC \
CUPY_INSTALL_USE_HIP=1  \
ROCM_HOME=${ROCM_PATH}  \
HCC_AMDGPU_TARGET=${AMREX_AMD_ARCH}  \
  python3 -m pip install -v cupy
# optional: for libEnsemble
python3 -m pip install -r $HOME/src/warpx/Tools/LibEnsemble/requirements.txt
# optional: for optimas (based on libEnsemble & ax->botorch->gpytorch->pytorch)
#python3 -m pip install --upgrade torch --index-url https://download.pytorch.org/whl/rocm5.4.2
#python3 -m pip install -r $HOME/src/warpx/Tools/optimas/requirements.txt
