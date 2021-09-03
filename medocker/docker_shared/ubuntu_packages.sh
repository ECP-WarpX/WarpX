#!/bin/bash

set -e
set -x

export DEBIAN_FRONTEND=noninteractive

apt-get update -y

apt-get install -y \
    automake            \
    build-essential     \
    ca-certificates     \
    ccache              \
    clang               \
    cmake               \
    gcc                 \
    gfortran            \
    git                 \
    gnupg               \
    llvm-9-dev          \
    libblosc-dev        \
    libffi-dev          \
    libfreetype6-dev    \
    libhdf5-dev         \
    liblapack-dev       \
    liblapacke-dev      \
    libssl-dev          \
    libyaml-dev         \
    linux-libc-dev      \
    linux-tools-aws     \
    linux-tools-common  \
    locales             \
    lsof                \
    libblas-dev         \
    libboost-math-dev   \
    libfftw3-dev        \
    liblapack-dev       \
    make                \
    ninja-build         \
    pkg-config          \
    python3             \
    python3-pip         \
    python3-setuptools  \
    python3-venv        \
    wget                \
    sudo                \
    software-properties-common \
    ;

# Some python packages will fail to build with the default encoding
# due to non-English characters in their source
sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen
dpkg-reconfigure locales
update-locale LANG=en_US.UTF-8
echo 'export LANG=en_US.UTF-8' >> /etc/profile.d/10-set-locale.sh

echo 'export LLVM_CONFIG=/usr/bin/llvm-config-9' >> /etc/profile.d/11-configure-llvm.sh

if [ ${WARPXCOMPUTE}=CUDA ] ; then
    # install cuda toolkit
    wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
    sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
    sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/7fa2af80.pub
    sudo add-apt-repository "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/ /"
    sudo apt-get -y install cuda

    # so cmake knows where to find cuda compiler
    echo 'export CUDACXX=/usr/local/cuda-11.4/bin/nvcc' >> /etc/profile.d/12-configure-cuda.sh
fi
