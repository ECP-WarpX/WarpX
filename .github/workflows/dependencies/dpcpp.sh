#!/usr/bin/env bash
#
# Copyright 2020 The WarpX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

# Ref.: https://github.com/rscohn2/oneapi-ci
# intel-basekit intel-hpckit are too large in size
wget -q -O - https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB \
  | sudo apt-key add -
echo "deb https://apt.repos.intel.com/oneapi all main" \
  | sudo tee /etc/apt/sources.list.d/oneAPI.list

sudo apt-get update

df -h
# Install and reduce disk space
# https://github.com/ECP-WarpX/WarpX/pull/1566#issuecomment-790934878
sudo apt-get install -y --no-install-recommends \
    build-essential \
    cmake           \
    intel-oneapi-dpcpp-cpp-compiler intel-oneapi-mkl-devel \
    g++ gfortran    \
    libopenmpi-dev  \
    openmpi-bin     && \
du -sh /opt/intel/oneapi/
du -sh /opt/intel/oneapi/*/*
sudo rm -rf /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_sycl.a \
            /opt/intel/oneapi/compiler/latest/linux/lib/oclfpga    \
            /opt/intel/oneapi/compiler/latest/linux/lib/emu
du -sh /opt/intel/oneapi/
du -sh /opt/intel/oneapi/*/*
sudo rm -rf /opt/intel/oneapi/mkl/latest/lib/intel64/*.a
du -sh /opt/intel/oneapi/
du -sh /opt/intel/oneapi/*/*
df -h
