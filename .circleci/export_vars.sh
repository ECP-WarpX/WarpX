#!/bin/bash

export WARPX_DIMS='1;2;RZ'
export WARPX_MPI=ON
export WARPX_OPENPMD=OFF
export WARPX_PSATD=OFF
export WARPX_QED=OFF
export WARPX_EB=ON
export CC=$(which clang)
export CXX=$(which clang++)

export SETUPTOOLS_USE_DISTUTILS=stdlib
