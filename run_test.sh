#!/bin/bash

# Copyright 2018-2020 Axel Huebl, David Grote, Edoardo Zoni
# Luca Fedeli, Maxence Thevenet, Remi Lehe
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script runs some of WarpX's standard regression tests, but
# without comparing the output to previously run simulations.
# This checks that:
# - The code compiles and runs without error
# - For some of the tests, a Python script checks that the results are
# physically correct.

# The tests can be influenced by environment variables:
# Use `export WARPX_CI_DIM=3` or `export WARPX_CI_DIM=2` in order to
# select only the tests that correspond to this dimension
# Use `export WARPX_TEST_ARCH=CPU` or `export WARPX_TEST_ARCH=GPU` in order
# to run the tests on CPU or GPU respectively.

set -eu -o pipefail

# Parse command line arguments: if test names are given as command line arguments,
# store them in variable tests_arg and define new command line argument to call
# regtest.py with option --tests (works also for single test)
tests_arg=$*
tests_run=${tests_arg:+--tests=${tests_arg}}

# Remove contents and link to a previous test directory (intentionally two arguments)
rm -rf test_dir/* test_dir
# Create a temporary test directory
tmp_dir=$(mktemp --help >/dev/null 2>&1 && mktemp -d -t ci-XXXXXXXXXX || mktemp -d "${TMPDIR:-/tmp}"/ci-XXXXXXXXXX)
if [ $? -ne 0 ]; then
    echo "Cannot create temporary directory"
    exit 1
fi

# Copy WarpX into current test directory
mkdir ${tmp_dir}/warpx
cp -r ./* ${tmp_dir}/warpx

# Link the test directory
ln -s ${tmp_dir} test_dir

# Switch to the test directory
cd test_dir
echo "cd $PWD"

# Clone PICSAR, AMReX and warpx-data
git clone https://github.com/AMReX-Codes/amrex.git
cd amrex && git checkout 60578291de339964af938b81f2d0f7db5c6bf78c && cd -
# Use QED brach for QED tests
git clone https://github.com/ECP-WarpX/picsar.git
cd picsar && git checkout c16b642e3dcf860480dd1dd21cefa3874f395773 && cd -
# warpx-data contains various required data sets
git clone --depth 1 https://github.com/ECP-WarpX/warpx-data.git

# Clone the AMReX regression test utility
git clone https://github.com/ECP-WarpX/regression_testing.git

# Prepare regression tests
mkdir -p rt-WarpX/WarpX-benchmarks
cd warpx/Regression
echo "cd $PWD"
python prepare_file_travis.py
cp travis-tests.ini ../../rt-WarpX
cp -r Checksum ../../regression_testing/

# Run tests
cd ../../regression_testing/
echo "cd $PWD"
# run only tests specified in variable tests_arg (single test or multiple tests)
if [[ ! -z "${tests_arg}" ]]; then
  python regtest.py ../rt-WarpX/travis-tests.ini --no_update all "${tests_run}"
# run all tests (variables tests_arg and tests_run are empty)
else
  python regtest.py ../rt-WarpX/travis-tests.ini --no_update all
fi
