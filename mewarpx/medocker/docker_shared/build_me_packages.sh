#!/bin/bash

set -e
set -x

# Install mewarpx
pushd WarpX/mewarpx

make devel-all

popd

# Build and install warpx
pushd WarpX

export WARPX_MPI=ON
export WARPX_DIMS='1;2;RZ'
export BUILD_SHARED_LIBS=ON
export WARPX_QED=OFF
export WARPX_EB=ON
export BUILD_PARALLEL=8
export CC=$(which clang)
export WARPX_COMPUTE=${WARPXCOMPUTE}

# compile binaries for V100 and A100 GPUs
# also compile with single precision particles if on GPU
if [ "${WARPX_COMPUTE}" == "CUDA" ]; then
    export AMREX_CUDA_ARCH="7.0;8.0"
    # export WARPX_PARTICLE_PRECISION=SINGLE (turn on after restart issue is fixed)
fi

python3 -m pip install -v .

popd

# install minerva
pushd minerva

python3 -m pip install .

popd
