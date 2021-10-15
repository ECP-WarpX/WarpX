#!/usr/bin/env bash
#
# Copyright 2020-2021 The WarpX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

export DEBIAN_FRONTEND=noninteractive
sudo apt-get -qqq update
sudo apt-get install -y \
  build-essential \
  cmake           \
  ca-certificates \
  gnupg           \
  pkg-config      \
  wget

# Ref.: https://github.com/rscohn2/oneapi-ci
sudo wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB
echo "deb https://apt.repos.intel.com/oneapi all main" | \
    sudo tee /etc/apt/sources.list.d/oneAPI.list
sudo apt-get update
sudo apt-get install -y intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic

# activate now via
set +e
source /opt/intel/oneapi/setvars.sh
set -e

# cmake-easyinstall
sudo curl -L -o /usr/local/bin/cmake-easyinstall https://git.io/JvLxY
sudo chmod a+x /usr/local/bin/cmake-easyinstall
export CEI_SUDO="sudo"

# openPMD-api
CXX=$(which icpc) CC=$(which icc) \
  cmake-easyinstall               \
  --prefix=/usr/local             \
  git+https://github.com/openPMD/openPMD-api.git@0.14.2 \
  -DopenPMD_USE_PYTHON=OFF \
  -DBUILD_TESTING=OFF      \
  -DBUILD_EXAMPLES=OFF     \
  -DBUILD_CLI_TOOLS=OFF
