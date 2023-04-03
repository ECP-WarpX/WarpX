#!/usr/bin/env bash
#
# Copyright 2020-2021 The WarpX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

# `man apt.conf`:
#   Number of retries to perform. If this is non-zero APT will retry
#   failed files the given number of times.
echo 'Acquire::Retries "3";' | sudo tee /etc/apt/apt.conf.d/80-retries

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
tries=0
while [[ ${tries} -lt 5 ]]
do
    sudo apt-get install -y --no-install-recommends \
        build-essential \
        ccache          \
        cmake           \
        intel-oneapi-dpcpp-cpp-compiler intel-oneapi-mkl-devel \
        g++ gfortran    \
        libopenmpi-dev  \
        openmpi-bin     \
        && { sudo apt-get clean; tries=6; }     \
        || { sleep 10; tries=$(( tries + 1 )); }
done
if [[ ${tries} -eq 5 ]]; then exit 1; fi

du -sh /opt/intel/oneapi/
du -sh /opt/intel/oneapi/*/*
echo "+++ REDUCING oneAPI install size +++"
sudo rm -rf /opt/intel/oneapi/mkl/latest/lib/intel64/*.a           \
            /opt/intel/oneapi/compiler/latest/linux/lib/oclfpga    \
            /opt/intel/oneapi/compiler/latest/linux/lib/emu        \
            /opt/intel/oneapi/compiler/latest/linux/bin/intel64    \
            /opt/intel/oneapi/compiler/latest/linux/bin/lld        \
            /opt/intel/oneapi/compiler/latest/linux/bin/lld-link   \
            /opt/intel/oneapi/compiler/latest/linux/bin/wasm-ld
du -sh /opt/intel/oneapi/
du -sh /opt/intel/oneapi/*/*
df -h
