#!/bin/bash

set -e
set -x

# Install mewarpx
pushd WarpX/mewarpx

make devel-all

popd

# Build and install warpx
pushd WarpX

export WarpX_MPI=ON
export WarpX_DIMS=2
export WarpX_LIBS=ON
export WarpX_EB=ON
export BUILD_PARALLEL=8
export CC=$(which clang)
export WarpX_COMPUTE=${WARPXCOMPUTE}

python3 -m pip install -v .

popd
