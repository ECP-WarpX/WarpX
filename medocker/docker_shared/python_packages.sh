#!/bin/bash

set -e
set -x

# We run `pip3 install --no-build-isolation foo` to ask it to build in
# the environment we already have, rather than in a freshly
# constructed one in which any declared build dependencies get
# installed afresh.

INSTALL_CMD="python3 -m pip install --no-build-isolation "

# Upgrade pip to get --no-build-isolation support

# Upgrade keyring and install keyring.alt to suppress annoying
# internal warnings
# https://github.com/Nike-Inc/gimme-aws-creds/issues/158
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade setuptools keyring keyrings.alt

# Build dependency of cmake
${INSTALL_CMD} scikit-build

# Build dependency of several downstream things
${INSTALL_CMD} cmake

${INSTALL_CMD} llvmlite

# The builds of scipy and h5py on aarch64 depends on having these already installed
${INSTALL_CMD} cython
${INSTALL_CMD} numpy
${INSTALL_CMD} pkgconfig

${INSTALL_CMD} numba

# Prevent building h5py from trying to pull in the oldest possible dependencies
H5PY_SETUP_REQUIRES=0 ${INSTALL_CMD} h5py

${INSTALL_CMD} \
     scipy \
     Forthon \
     matplotlib \
     pykern \
     yt

${INSTALL_CMD} scikit-image

# This is necessary for S3 copying
${INSTALL_CMD} awscli

${INSTALL_CMD} pandas

${INSTALL_CMD} mpi4py
