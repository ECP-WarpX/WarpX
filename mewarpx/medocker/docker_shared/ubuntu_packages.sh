#!/bin/bash

set -e
set -x

export DEBIAN_FRONTEND=noninteractive

# Some python packages will fail to build with the default encoding
# due to non-English characters in their source
echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections
echo 'locales locales/default_environment_locale en_US.UTF-8' | debconf-set-selections
echo 'locales locales/locales_to_be_generated en_US.UTF-8' | debconf-set-selections
echo 'keyboard-configuration keyboard-configuration/layout English (US)' | debconf-set-selections

apt-get update -y

apt-get install -y -o Dpkg::Options::="--force-confdef" \
    automake            \
    build-essential     \
    ca-certificates     \
    ccache              \
    clang               \
    cmake               \
    curl                \
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

echo 'export LLVM_CONFIG=/usr/bin/llvm-config-9' >> /etc/profile.d/11-configure-llvm.sh

if [ ${WARPXCOMPUTE} == CUDA ] ; then
    # install cuda toolkit
    wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
    mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
    apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/3bf863cc.pub
    add-apt-repository "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/ /"
    apt-get update
    apt-get -y install cuda

    # so cmake knows where to find cuda compiler
    export PATH=/usr/local/cuda/bin${PATH:+:${PATH}}
    echo 'export CUDACXX=/usr/local/cuda/bin/nvcc' >> /etc/profile.d/12-configure-cuda.sh
    export LD_LIBRARY_PATH=/usr/local/cuda/lib64\
                         ${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

fi
