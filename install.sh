#!/bin/bash

# uninstall old versions
#rm -rf build
#rm -rf *.whl

# Need to use LLVM's clang instead of apple's clang to use openmp
export CC=/opt/homebrew/Cellar/llvm/18.1.8/bin/clang
export CXX=/opt/homebrew/Cellar/llvm/18.1.8/bin/clang++

# Activate python virtual environemnt
source /Users/archermarks/venvs/warpx/bin/activate

# Build warpx
cmake -S . -B build -Wno-dev \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DWarpX_LIB=ON \
    -DWarpX_APP=ON \
    -DWarpX_MPI=ON \
    -DWarpX_DIMS="1;2;3" \
    -DWarpX_PRECISION=DOUBLE \
    -DWarpX_PARTICLE_PRECISION=SINGLE \
    -DWarpX_COMPUTE=OMP

    #-DWarpX_PYTHON=ON \
cmake --build build -j 20
#cmake --build build --target pip_install -j 20
