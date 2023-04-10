#!/usr/bin/env bash
#
# Copyright 2020-2022 The WarpX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

# `man apt.conf`:
#   Number of retries to perform. If this is non-zero APT will retry
#   failed files the given number of times.
echo 'Acquire::Retries "3";' | sudo tee /etc/apt/apt.conf.d/80-retries

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

wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/3bf863cc.pub
echo "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/ /" \
    | sudo tee /etc/apt/sources.list.d/cuda.list

sudo apt-get update
sudo apt-get install -y          \
    cuda-command-line-tools-11-0 \
    cuda-compiler-11-0           \
    cuda-cupti-dev-11-0          \
    cuda-minimal-build-11-0      \
    cuda-nvml-dev-11-0           \
    cuda-nvtx-11-0               \
    libcufft-dev-11-0            \
    libcurand-dev-11-0
sudo ln -s cuda-11.0 /usr/local/cuda

# if we run out of temporary storage in CI:
#du -sh /usr/local/cuda-11.0
#echo "+++ REDUCING CUDA Toolkit install size +++"
#sudo rm -rf /usr/local/cuda-11.0/targets/x86_64-linux/lib/libcu{fft,pti,rand}_static.a
#sudo rm -rf /usr/local/cuda-11.0/targets/x86_64-linux/lib/libnvperf_host_static.a
#du -sh /usr/local/cuda-11.0/
#df -h

# cmake-easyinstall
#
sudo curl -L -o /usr/local/bin/cmake-easyinstall https://raw.githubusercontent.com/ax3l/cmake-easyinstall/main/cmake-easyinstall
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
