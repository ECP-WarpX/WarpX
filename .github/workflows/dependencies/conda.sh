#!/usr/bin/env bash
#
# Copyright 2023 The WarpX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

conda update -n base conda
conda install -n base conda-libmamba-solver
conda config --set solver libmamba

conda install -c conda-forge blaspp ccache cmake compilers git lapackpp "openpmd-api=*=mpi_mpich*" python "mpi4py=*=mpi_mpich*" numpy pandas scipy yt "fftw=*=mpi_mpich*" pkg-config setuptools matplotlib mpich pip virtualenv wheel
