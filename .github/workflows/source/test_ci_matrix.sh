#!/usr/bin/env bash

set -eu -o pipefail

cp .github/workflows/source/ci_matrix.py Regression/
cd Regression/

# Put the name of all CI tests into a text file
python prepare_file_ci.py
grep "\[" ci-tests.ini > ci_all_tests.txt


export WARPX_CI_PSATD=TRUE

# Concatenate the names of all elements in CI matrix into another test file
WARPX_CI_REGULAR_CARTESIAN_2D=TRUE python prepare_file_ci.py
grep "\[" ci-tests.ini >  ci_matrix_elements.txt
WARPX_CI_REGULAR_CARTESIAN_3D=TRUE python prepare_file_ci.py
grep "\[" ci-tests.ini >>  ci_matrix_elements.txt
WARPX_CI_PYTHON_MAIN=TRUE       python prepare_file_ci.py
grep "\[" ci-tests.ini >> ci_matrix_elements.txt
WARPX_CI_SINGLE_PRECISION=TRUE  python prepare_file_ci.py
grep "\[" ci-tests.ini >> ci_matrix_elements.txt
WARPX_CI_RZ_OR_NOMPI=TRUE      python prepare_file_ci.py
grep "\[" ci-tests.ini >> ci_matrix_elements.txt
WARPX_CI_QED=TRUE               python prepare_file_ci.py
grep "\[" ci-tests.ini >> ci_matrix_elements.txt
WARPX_CI_EB=TRUE               python prepare_file_ci.py
grep "\[" ci-tests.ini >> ci_matrix_elements.txt

# Check that the resulting lists are equal
{
    python ci_matrix.py &&
        rm ci_all_tests.txt ci_matrix_elements.txt ci_matrix.py &&
        echo "test passed" &&
        exit 0
} || {
        rm ci_all_tests.txt ci_matrix_elements.txt ci_matrix.py &&
        echo "tests failed" &&
        exit 1
}
