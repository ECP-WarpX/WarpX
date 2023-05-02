#!/usr/bin/env bash
#
# Copyright 2023 The WarpX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

conda update -y -n base conda
conda install -y -n base conda-libmamba-solver
conda config --set solver libmamba

conda install -y -n testing -c conda-forge blaspp ccache cmake compilers git lapackpp "openpmd-api=*=mpi_mpich*" python mpi4py numpy pandas scipy yt "fftw=*=mpi_mpich*" pkg-config setuptools matplotlib mpich pip virtualenv wheel
