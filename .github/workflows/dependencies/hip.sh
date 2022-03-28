#!/usr/bin/env bash
#
# Copyright 2020 The WarpX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

# Ref.: https://rocmdocs.amd.com/en/latest/Installation_Guide/Installation-Guide.html#ubuntu
wget -q -O - http://repo.radeon.com/rocm/rocm.gpg.key \
  | sudo apt-key add -
echo 'deb [arch=amd64] http://repo.radeon.com/rocm/apt/debian/ ubuntu main' \
  | sudo tee /etc/apt/sources.list.d/rocm.list

echo 'export PATH=/opt/rocm/llvm/bin:/opt/rocm/bin:/opt/rocm/profiler/bin:/opt/rocm/opencl/bin:$PATH' \
  | sudo tee -a /etc/profile.d/rocm.sh
# we should not need to export HIP_PATH=/opt/rocm/hip with those installs

sudo apt-get update

# Ref.: https://rocmdocs.amd.com/en/latest/Installation_Guide/Installation-Guide.html#installing-development-packages-for-cross-compilation
# meta-package: rocm-dkms
# OpenCL: rocm-opencl
# other: rocm-dev rocm-utils
sudo apt-get install -y --no-install-recommends \
    build-essential \
    gfortran        \
    libhiredis-dev  \
    libnuma-dev     \
    libopenmpi-dev  \
    libzstd-dev     \
    ninja-build     \
    openmpi-bin     \
    rocm-dev        \
    rocfft-dev      \
    rocprim-dev     \
    rocrand-dev

# activate
#
source /etc/profile.d/rocm.sh
hipcc --version
which clang
which clang++

# cmake-easyinstall
#
sudo curl -L -o /usr/local/bin/cmake-easyinstall https://git.io/JvLxY
sudo chmod a+x /usr/local/bin/cmake-easyinstall
export CEI_SUDO="sudo"
export CEI_TMP="/tmp/cei"

# ccache 4.2+
#
CXXFLAGS="" cmake-easyinstall --prefix=/usr/local \
    git+https://github.com/ccache/ccache.git@v4.6 \
    -DCMAKE_BUILD_TYPE=Release        \
    -DENABLE_DOCUMENTATION=OFF        \
    -DENABLE_TESTING=OFF              \
    -DWARNINGS_AS_ERRORS=OFF
