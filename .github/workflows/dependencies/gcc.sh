#!/usr/bin/env bash
#
# Copyright 2021 The WarpX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

sudo apt-get -qqq update
sudo apt-get install -y \
    build-essential     \
    ca-certificates     \
    ccache              \
    cmake               \
    gnupg               \
    ninja-build         \
    pkg-config          \
    wget
