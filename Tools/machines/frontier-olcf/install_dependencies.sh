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
#   Was perlmutter_gpu_warpx.profile sourced and configured correctly?
if [ -z "${proj}" ]; then echo "WARNING: The 'proj' variable is not yet set in your frontier_warpx.profile file! Please edit its line 2 to continue!"; return; fi


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
rm -rf ${HOME}/sw/frontier/gpu

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
  git fetch
  git pull
  cd -
else
  git clone https://github.com/icl-utk-edu/blaspp.git $HOME/src/blaspp
fi
rm -rf $HOME/src/blaspp-frontier-gpu-build
CXX=$(which CC) cmake -S $HOME/src/blaspp -B $HOME/src/blaspp-frontier-gpu-build -Duse_openmp=OFF -Dgpu_backend=hip -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=${HOME}/sw/frontier/gpu/blaspp-master
cmake --build $HOME/src/blaspp-frontier-gpu-build --target install --parallel 16

# LAPACK++ (for PSATD+RZ)
if [ -d $HOME/src/lapackpp ]
then
  cd $HOME/src/lapackpp
  git fetch
  git pull
  cd -
else
  git clone https://github.com/icl-utk-edu/lapackpp.git $HOME/src/lapackpp
fi
rm -rf $HOME/src/lapackpp-frontier-gpu-build
CXX=$(which CC) CXXFLAGS="-DLAPACK_FORTRAN_ADD_" cmake -S $HOME/src/lapackpp -B $HOME/src/lapackpp-frontier-gpu-build -DCMAKE_CXX_STANDARD=17 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=${HOME}/sw/frontier/gpu/lapackpp-master
cmake --build $HOME/src/lapackpp-frontier-gpu-build --target install --parallel 16


# Python ######################################################################
#
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade virtualenv
python3 -m pip cache purge
rm -rf ${HOME}/sw/frontier/gpu/venvs/warpx
python3 -m venv ${HOME}/sw/frontier/gpu/venvs/warpx
source ${HOME}/sw/frontier/gpu/venvs/warpx/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade wheel
python3 -m pip install --upgrade cython
python3 -m pip install --upgrade numpy
python3 -m pip install --upgrade pandas
python3 -m pip install --upgrade scipy
MPICC="cc -shared" python3 -m pip install --upgrade mpi4py --no-cache-dir --no-build-isolation --no-binary mpi4py
python3 -m pip install --upgrade openpmd-api
python3 -m pip install --upgrade matplotlib
python3 -m pip install --upgrade yt
# install or update WarpX dependencies such as picmistandard
python3 -m pip install --upgrade -r $HOME/src/warpx/requirements.txt
# optional: for libEnsemble
python3 -m pip install -r $HOME/src/warpx/Tools/LibEnsemble/requirements.txt
# optional: for optimas (based on libEnsemble & ax->botorch->gpytorch->pytorch)
#python3 -m pip install --upgrade torch --index-url https://download.pytorch.org/whl/rocm5.4.2
#python3 -m pip install -r $HOME/src/warpx/Tools/optimas/requirements.txt
