#!/usr/bin/env bash
#
# Copyright 2020-2022 The WarpX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

sudo apt-get -qqq update
sudo apt-get install -y \
    build-essential     \
    ca-certificates     \
    cmake               \
    gnupg               \
    libhiredis-dev      \
    libopenmpi-dev      \
    libzstd-dev         \
    ninja-build         \
    openmpi-bin         \
    pkg-config          \
    wget

wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-keyring_1.0-1_all.deb
sudo dpkg -i cuda-keyring_1.0-1_all.deb

sudo apt-get update
sudo apt-get install -y          \
    cuda-command-line-tools-11-8 \
    cuda-compiler-11-8           \
    cuda-cupti-dev-11-8          \
    cuda-minimal-build-11-8      \
    cuda-nvml-dev-11-8           \
    cuda-nvtx-11-8               \
    libcufft-dev-11-8            \
    libcurand-dev-11-8
sudo ln -s cuda-11.8 /usr/local/cuda

# if we run out of temporary storage in CI:
#du -sh /usr/local/cuda-11.8
#echo "+++ REDUCING CUDA Toolkit install size +++"
#sudo rm -rf /usr/local/cuda-11.8/targets/x86_64-linux/lib/libcu{fft,pti,rand}_static.a
#sudo rm -rf /usr/local/cuda-11.8/targets/x86_64-linux/lib/libnvperf_host_static.a
#du -sh /usr/local/cuda-11.8/
#df -h

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
